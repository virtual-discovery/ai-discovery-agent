"""Reference solution for fresh-srtio3-soft-modes; runs as the Harbor oracle and in check_task.

It requests new pw.x / ph.x calculations from the lab. The functions at module level (request builders and
estimators) are also used by the construction scripts (benchmarks/construction/fresh-srtio3-soft-modes/) to plan
the reference calculations and to apply the same estimators to them.

Workflow: equation of state -> a0; frozen-phonon force constants of the z-polarized zone-centre modes on two k
meshes (the soft mode converges slowly with k) with the acoustic sum rule; electric-field response on two meshes
-> eps_inf, Z*; LO frequencies from the nonanalytic term; energy scan along the soft-mode eigendisplacement ->
polar well depth; octahedral-rotation scan in the R-point doubled cell with the undistorted reference in the SAME
cell and mesh -> tilt angle, depth and the R-mode harmonic frequency; biaxially strained cell with the
out-of-plane lattice parameter relaxed -> in-plane polar frequency.
"""
import json
import math
import os
import subprocess
import time
from pathlib import Path

import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
RY_CM = 109737.31568160          # 1 Ry in cm^-1 (frequency units: sqrt(Ry / (Ry mass units bohr^2)))
AMU_RY = 911.44424310865645      # 1 amu in Rydberg mass units
BOHR_A = 0.529177210903
RY_EV = 13.605693122994
MASS = {'Sr': 87.62, 'Ti': 47.867, 'O': 15.999}
ECUT = 60
EOS_A = (3.90, 3.92, 3.94, 3.96)
U_FP = 0.01                      # frozen-phonon displacement (Angstrom)
TILT_DEG = (2, 3, 4, 5)
STRAIN = 0.01                    # biaxial in-plane strain relative to a0
C_SCAN = (3.89, 3.91, 3.93)
POLAR_TI_AMPS = (0.01, 0.02, 0.03, 0.04)   # Ti displacement along the soft-mode eigendisplacement (Angstrom)
COLS_Z = {'Ti': [1], 'Oz': [4], 'Oxy': [2, 3]}      # displaced atoms (perovskite_5 order Sr, Ti, O_x, O_y, O_z)
COLS_X = {'Ti': [1], 'Ox': [2], 'Oy': [3], 'Oz': [4]}


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def compute(request, log):
    record = lab('evaluate', request)
    while record['status'] == 'running':
        time.sleep(5)
        record = lab('job/%d' % record['id'])
    if record['status'] != 'complete':
        raise RuntimeError('calculation %s failed: %s' % (request, record.get('error')))
    log.append(record['id'])
    print(json.dumps({'id': record['id'], 'charge': record['charge'], 'cumulative': record['cumulative_charge'],
                      'wall_s': record['result']['wall_seconds']}), flush=True)
    return record['result']


# ---------------------------------------------------------------- requests
def zeros(n):
    return [[0.0, 0.0, 0.0] for _ in range(n)]


def request(kind, cell, a, c, disp, k, ecut=ECUT):
    return {'kind': kind, 'cell': cell, 'lattice_A': round(a, 5), 'c_A': round(c, 5),
            'displacements_A': [[round(v, 5) + 0.0 for v in row] for row in disp], 'ecutwfc_Ry': ecut, 'kmesh': k}


def displaced(atoms, comp, u, n=5):
    d = zeros(n)
    for i in atoms:
        d[i][comp] = u
    return d


def tilt_pattern(a, theta_deg):
    """a0a0c- rotation about z in the rotated_10 cell: O_x and O_y of formula unit 1 (rotation sense +1) and of
    formula unit 2 (sense -1) move perpendicular to their Ti-O-Ti line by u = (a/2) tan(theta)."""
    u = a / 2 * math.tan(math.radians(theta_deg))
    d = zeros(10)
    for f, s in ((0, 1), (1, -1)):
        d[5 * f + 2] = [0.0, -s * u, 0.0]
        d[5 * f + 3] = [s * u, 0.0, 0.0]
    return d


def eos_requests():
    return [request('scf', 'perovskite_5', a, a, zeros(5), 6) for a in EOS_A]


def fp_z_requests(a, k, u=U_FP):
    return {name: request('scf', 'perovskite_5', a, a, displaced(atoms, 2, u), k) for name, atoms in COLS_Z.items()}


def dielectric_request(a, k):
    return request('dielectric', 'perovskite_5', a, a, zeros(5), k)


def cubic_request(a, k):
    return request('scf', 'perovskite_5', a, a, zeros(5), k)


def polar_requests(a, pattern, k):
    return [request('scf', 'perovskite_5', a, a, (np.array(pattern) * t).tolist(), k) for t in POLAR_TI_AMPS]


def tilt_requests(a, k, angles):
    return {th: request('scf', 'rotated_10', a, a, tilt_pattern(a, th), k) for th in angles}


def strained_a(a0):
    return round(a0 * (1 + STRAIN), 4)


def c_scan_requests(a0):
    return [request('scf', 'perovskite_5', strained_a(a0), c, zeros(5), 6) for c in C_SCAN]


def fp_x_requests(a0, c, k=6, u=U_FP):
    return {name: request('scf', 'perovskite_5', strained_a(a0), c, displaced(atoms, 0, u), k)
            for name, atoms in COLS_X.items()}


# ---------------------------------------------------------------- estimators (pure functions of lab results)
def a0_from_eos(results):
    a = np.array([r['settings']['lattice_A'] for r in results])
    e = np.array([r['total_energy_Ry'] for r in results])
    c = np.polyfit(a, e, 3)
    roots = [r.real for r in np.roots(np.polyder(c)) if abs(r.imag) < 1e-9 and a.min() < r.real < a.max()
             and np.polyval(np.polyder(c, 2), r.real) > 0]
    return float(roots[0])


def column(res, u_A, comp, atoms_moved):
    """-dF/du (Ry/bohr^2) for every atom's `comp` force when the listed atoms move together by u."""
    return -np.array(res['forces_Ry_per_bohr'])[:, comp] / (u_A / BOHR_A)


def z_block(fp, u=U_FP):
    """Force-constant matrix of the z-polar coordinates (Sr, Ti, O_z, s) where s moves O_x and O_y together.
    Columns from the displaced cells; the Sr column from the acoustic sum rule; symmetrized."""
    def gen(res):
        f = column(res, u, 2, None)
        return np.array([f[0], f[1], f[4], f[2] + f[3]])
    K = np.zeros((4, 4))
    K[:, 1], K[:, 2], K[:, 3] = gen(fp['Ti']), gen(fp['Oz']), gen(fp['Oxy'])
    K[:, 0] = -(K[:, 1] + K[:, 2] + K[:, 3])
    return 0.5 * (K + K.T)


Z_MASSES = np.array([MASS['Sr'], MASS['Ti'], MASS['O'], 2 * MASS['O']]) * AMU_RY


def neutral_charges(diel):
    z = np.array([diel['born_charges_raw'][lab] for lab in ('Sr', 'Ti', 'O_x', 'O_y', 'O_z')])
    return z - z.mean(axis=0)


def modes(K, masses, nac=None):
    """Frequencies (cm^-1, negative = imaginary) and mass-weighted eigenvectors; nac = (Zq vector, eps_qq,
    volume A^3) adds the macroscopic-field term 4 pi e^2 (q.Z_i)(q.Z_j)/(Omega q.eps.q) (e^2 = 2 in Ry units)."""
    K = np.array(K, float)
    if nac is not None:
        zq, eps, vol = nac
        zq = np.array(zq, float)
        K = K + 4 * np.pi * 2.0 / (vol / BOHR_A ** 3) * np.outer(zq, zq) / eps
    w2, v = np.linalg.eigh(K / np.sqrt(np.outer(masses, masses)))
    return np.sign(w2) * np.sqrt(np.abs(w2)) * RY_CM, v


def z_nac(diel):
    z = neutral_charges(diel)
    zq = [z[0][2][2], z[1][2][2], z[4][2][2], z[2][2][2] + z[3][2][2]]
    return zq, diel['dielectric_tensor_eps_inf'][2][2], diel['cell_volume_A3']


def soft_pattern(K):
    """Displacement pattern (5 atoms x 3, Angstrom units per unit Ti amplitude) of the lowest z-polar mode,
    with Sr fixed (a rigid translation does not change the energy) and Ti displaced by +1."""
    _, v = modes(K, Z_MASSES)
    e = v[:, 0] / np.sqrt(Z_MASSES)
    e = e - e[0]
    e = e / e[1]
    d = zeros(5)
    d[1][2], d[4][2], d[2][2], d[3][2] = e[1], e[2], e[3], e[3]
    return [[round(x, 4) for x in row] for row in d]


def well(amps, energies_meV):
    """Fit E = A x^2 + B x^4 (+ C x^6 with more than three points) through the origin; return the depth
    (meV, positive = energy gain), the minimizing amplitude and A."""
    x, y = np.array(amps, float), np.array(energies_meV, float)
    powers = (2, 4, 6) if len(x) > 3 else (2, 4)
    M = np.array([[xi ** p for p in powers] for xi in x])
    coef = np.linalg.lstsq(M, y, rcond=None)[0]
    grid = np.linspace(0, 1.3 * x.max(), 20001)
    e = sum(c * grid ** p for c, p in zip(coef, powers))
    i = int(np.argmin(e))
    return float(-e[i]), float(grid[i]), float(coef[0])


def tilt_analysis(ref, tilts, a):
    """Tilt angle (deg), depth (meV per formula unit) and the harmonic R-mode frequency (cm^-1, negative for an
    instability) from rotated_10 energies relative to the undistorted rotated_10 cell (same mesh)."""
    angles = sorted(tilts)
    u = [a / 2 * math.tan(math.radians(t)) for t in angles]
    de = [(tilts[t]['total_energy_Ry'] - ref['total_energy_Ry']) * RY_EV * 1000 / 2 for t in angles]
    depth, umin, A = well(u, de)
    # per formula unit two O atoms move by u: E = A u^2 = (1/2)(2 m_O) w^2 u^2  ->  w^2 = A / m_O
    A_ry = A / 1000 / RY_EV * BOHR_A ** 2
    w2 = A_ry / (MASS['O'] * AMU_RY)
    return math.degrees(math.atan(2 * umin / a)), depth, float(np.sign(w2) * math.sqrt(abs(w2)) * RY_CM)


def x_block(fp, u=U_FP):
    K = np.zeros((5, 5))
    for j, name in ((1, 'Ti'), (2, 'Ox'), (3, 'Oy'), (4, 'Oz')):
        K[:, j] = column(fp[name], u, 0, None)
    K[:, 0] = -K[:, 1:].sum(axis=1)
    return 0.5 * (K + K.T)


X_MASSES = np.array([MASS['Sr'], MASS['Ti'], MASS['O'], MASS['O'], MASS['O']]) * AMU_RY


def c_relaxed(results):
    c = np.array([r['settings']['c_A'] for r in results])
    e = np.array([r['total_energy_Ry'] for r in results])
    p = np.polyfit(c, e, 2)
    return float(-p[1] / (2 * p[0]))


# ---------------------------------------------------------------- workflow
def main():
    ev = []
    eos = [compute(r, ev) for r in eos_requests()]
    a0_fit = a0_from_eos(eos)
    a0 = round(a0_fit, 3)
    # zone-centre z-polar force constants on 6x6x6 and 8x8x8 meshes
    fp6 = {n: compute(r, ev) for n, r in fp_z_requests(a0, 6).items()}
    fp8 = {n: compute(r, ev) for n, r in fp_z_requests(a0, 8).items()}
    K6, K8 = z_block(fp6), z_block(fp8)
    w6, w8 = modes(K6, Z_MASSES)[0], modes(K8, Z_MASSES)[0]
    # electric-field response on 6x6x6 and 8x8x8 meshes
    d6, d8 = compute(dielectric_request(a0, 6), ev), compute(dielectric_request(a0, 8), ev)
    z8 = neutral_charges(d8)
    eps_inf = d8['dielectric_tensor_eps_inf'][0][0]
    lo8 = modes(K8, Z_MASSES, z_nac(d8))[0]
    to_soft, lo_max = float(w8[0]), float(lo8[-1])
    # polar double well along the soft-mode eigendisplacement (same cell and mesh for the reference)
    pattern = soft_pattern(K8)
    cub6 = compute(cubic_request(a0, 6), ev)
    pol = [compute(r, ev) for r in polar_requests(a0, pattern, 6)]
    pol_de = [(r['total_energy_Ry'] - cub6['total_energy_Ry']) * RY_EV * 1000 for r in pol]
    polar_depth, polar_amp, _ = well(POLAR_TI_AMPS, pol_de)
    # antiferrodistortive a0a0c- tilts in the R-point cell, referenced to the undistorted cell on the same mesh
    ref4 = compute(tilt_requests(a0, 4, [0])[0], ev)
    tilts = {th: compute(r, ev) for th, r in tilt_requests(a0, 4, TILT_DEG).items()}
    tilt_deg, tilt_depth, w_R = tilt_analysis(ref4, tilts, a0)
    # biaxial tension: relax c at a = a0 (1 + 1 %), then in-plane (x) polar force constants
    cs = [compute(r, ev) for r in c_scan_requests(a0)]
    c_rel = round(c_relaxed(cs), 3)
    fpx = {n: compute(r, ev) for n, r in fp_x_requests(a0, c_rel).items()}
    wx = modes(x_block(fpx), X_MASSES)[0]
    submission = {
        'a0_A': a0_fit, 'gamma_soft_mode_cm-1': to_soft, 'gamma_LO_highest_cm-1': lo_max,
        'eps_inf': float(eps_inf), 'born_charge_Ti': float(z8[1][0][0]),
        'polar_well_depth_meV_per_fu': polar_depth,
        'R_mode_cm-1': w_R, 'tilt_angle_deg': tilt_deg, 'tilt_well_depth_meV_per_fu': tilt_depth,
        'larger_energy_gain': 'polar' if polar_depth > tilt_depth else 'tilt',
        'strained_c_A': c_relaxed(cs), 'strained_inplane_polar_mode_cm-1': float(wx[0]),
        'evidence_ids': ev}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    zc = neutral_charges(d6)
    (APP/'report.md').write_text(f"""# SrTiO3 soft modes in PBE (reference solution)

**Equilibrium.** Four SCF energies (6x6x6, 60 Ry) fitted with a cubic polynomial give a0 = {a0_fit:.4f} A (PBE
overestimates the experimental 3.905 A by ~0.8%). Everything below is evaluated at a0 = {a0:.3f} A with 60 Ry
(SSSP recommends 75 Ry for the O PAW file; energy differences and forces at 60 Ry agree with 75 Ry here).

**Zone-centre force constants.** Finite displacements of 0.01 A (Ti, the apical O_z and the O_x/O_y pair together
along z; the Sr column from the acoustic sum rule, i.e. translational invariance) give the mass-weighted 4x4 matrix
of the z-polar coordinates (T1u block). Lowest mode: {w6[0]:.1f} cm^-1 (6x6x6) and {w8[0]:.1f} cm^-1 (8x8x8):
the soft mode converges slowly with k (4x4x4 gives ~50 cm^-1 more imaginary), so I use 8x8x8. Displacing Ti alone
(a Slater-mode guess) misses the O-cage counter-motion and gives a stable 'frequency'. Eigendisplacement (Ti = 1):
{pattern[1][2]:.2f} Ti, {pattern[4][2]:.2f} O_z, {pattern[2][2]:.2f} O_x/O_y, Sr fixed.

**Dielectric response and LO.** eps_inf = {d6['dielectric_tensor_eps_inf'][0][0]:.3f} (6x6x6),
{eps_inf:.3f} (8x8x8); Z*(Ti) = {zc[1][0][0]:.3f}, {z8[1][0][0]:.3f} after charge neutrality (4x4x4 overestimates both
strongly). The LO modes for q || z add 4 pi e^2 (q.Z)(q.Z)/(Omega q.eps_inf.q) to the full polar block (all T1u
modes mix); the highest LO is {lo_max:.1f} cm^-1; the soft TO has no LO partner of its own.

**Polar double well.** Energies along the soft eigendisplacement (Ti 0.01-0.04 A, same 5-atom cell and mesh as the
reference) fitted with E = A x^2 + B x^4: depth {polar_depth:.2f} meV per formula unit at Ti displacement
{polar_amp:.3f} A (fixed cubic cell).

**Octahedral tilts.** In the R-point doubled cell (2 f.u.), the a0a0c- rotation lowers the energy relative to the
undistorted doubled cell on the same mesh: depth {tilt_depth:.2f} meV/f.u. at {tilt_deg:.2f} deg; the quadratic
coefficient gives the harmonic R-mode frequency {w_R:.1f} cm^-1. (Construction check: 6x6x6 changes the 4 deg
energy by 0.1 meV/f.u.) Comparing the 10-atom energies with twice a 5-atom energy on a different mesh would
shift the depth by ~1 meV/f.u. Energies are per formula unit (the cell holds two).

**Competition.** At the PBE lattice the {'polar' if polar_depth > tilt_depth else 'tilt'} instability gives the larger energy gain.
Both are shallow (a few meV/f.u.); the two order parameters couple biquadratically, and the tilts suppress the polar
instability, so the gains are not additive. Zero-point motion (quantum paraelectricity) suppresses a ferroelectric
order with such a shallow well in real SrTiO3, while the tilt (105 K transition) survives.

**Strain.** At a = a0(1 + 1%) the energy minimum in c gives c = {c_relaxed(cs):.4f} A (Poisson contraction); the in-plane
polar mode is {wx[0]:.1f} cm^-1 (6x6x6), more unstable than in the cubic cell. Keeping c = a0 would overestimate the
volume and change this frequency.

**Limits.** Static-lattice PBE at one pseudopotential set; harmonic frequencies; fixed cubic cells for the wells (no
strain coupling); the polar instability is a PBE-volume artefact to a large degree: at the experimental lattice
parameter the soft mode is much less unstable.
""")
    lab('submit', submission)
    for _ in range(120):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
