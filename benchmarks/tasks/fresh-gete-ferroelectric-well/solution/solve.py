"""Reference solution for fresh-gete-ferroelectric-well; runs as the Harbor oracle and in check_task.

Route (no literature numbers are used; the only starting guess is a rhombohedral distortion of the protocol's own
rocksalt cell, refined by the computed gradients):
1. rocksalt pressures at three lattice parameters (k16) -> a0 where P = 0;
2. R3m energy surface from seven k12 SCFs (centre and +- steps in a, alpha, x) (energies + full gradients from forces and stress) -> quadratic model ->
   stationary point; one k20 SCF there gives the converged energy and a Newton correction of the structure;
3. rocksalt at a0 on k20/k24 meshes (the cubic energy converges slowly: small gap at L);
4. clamped cubic cell at a0: polar displacements x = 0.480, 0.485 (k20) with forces -> E(d) = E0 + A d^2 + B d^4 fit;
5. soft mode under compression: x = 0.495 at a_cub = 5.80, 5.75, 5.70 A (k20, MV 0.005; the compressed rocksalt
   phase is a semimetal) -> harmonic force constant from the force, zero crossing vs pressure;
6. spin-orbit: PSlibrary scalar and fully relativistic pairs at the two protocol structures (k12).
"""
import json
import math
import os
import subprocess
import time
from pathlib import Path

import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
BOHR_A = 0.529177210903
RY_MEV = 13605.693122994
GPA_PER_RY_BOHR3 = 14710.507848260711
S2 = math.sqrt(2)


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def calc(request):
    """Submit one calculation and poll until it is complete; returns (job id, result)."""
    rec = lab('evaluate', request)
    while rec['status'] == 'running':
        time.sleep(float(os.environ.get("POLL_SECONDS", "10")))
        rec = lab('job/%d' % rec['id'])
    if rec['status'] != 'complete':
        raise RuntimeError('job %d failed: %s' % (rec['id'], rec.get('error')))
    return rec['id'], rec['result']


def req(a_rh, alpha=60.0, x=0.5, kmesh=16, kshift=0, ecut=50, occ='fixed', degauss=None, pp='sssp'):
    return {'kind': 'scf', 'pseudopotential_set': pp, 'a_rh_A': round(a_rh, 6), 'alpha_deg': round(alpha, 5),
            'x_Te': round(x, 6), 'ecutwfc_Ry': ecut, 'kmesh': kmesh, 'kshift': kshift, 'occupations': occ,
            'degauss_Ry': degauss}


# ---------------------------------------------------------------- geometry and gradients
def cell_vectors(a, alpha_deg):
    c = math.cos(math.radians(alpha_deg))
    q = (math.sqrt(1 + 2 * c) + math.sqrt(1 - c)) / 3
    p = q - math.sqrt(1 - c)
    return a * np.array([[p, q, q], [q, p, q], [q, q, p]])


def volume(a, alpha_deg):
    c = math.cos(math.radians(alpha_deg))
    return a ** 3 * math.sqrt(1 - 3 * c * c + 2 * c ** 3)


def gradient(result):
    """(dE/da [Ry/A], dE/dalpha [Ry/deg], dE/dx [Ry]) from the pw.x stress (dE/deps = -Omega sigma) and forces."""
    s = result['settings']
    a, al = s['a_rh_A'], s['alpha_deg']
    dEde = -volume(a, al) / BOHR_A ** 3 * np.array(result['stress_GPa']) / GPA_PER_RY_BOHR3
    A0inv = np.linalg.inv(cell_vectors(a, al).T)
    g = []
    for da, dal in ((1e-5, 0), (0, 1e-4)):
        dF = (cell_vectors(a + da, al + dal).T - cell_vectors(a - da, al - dal).T) @ A0inv / (2 * (da + dal))
        g.append(float(np.sum(dEde * 0.5 * (dF + dF.T))))
    axis = cell_vectors(a, al).sum(axis=0) / BOHR_A
    g.append(float(-np.array(result['forces_Ry_per_bohr'][1]) @ axis))
    return np.array(g)


SCALE = np.array([0.05, 0.5, 0.005])        # natural steps of a (A), alpha (deg), x


def quadratic_model(results):
    """Least-squares quadratic E(a, alpha, x) through energies and gradients; returns (centre, coef, b, H) in scaled
    coordinates u = (p - centre) / SCALE, energies in meV."""
    pts = np.array([[r['settings']['a_rh_A'], r['settings']['alpha_deg'], r['settings']['x_Te']] for r in results])
    centre = pts.mean(axis=0)
    pairs = [(i, j) for i in range(3) for j in range(i, 3)]
    rows, rhs = [], []
    for r, p in zip(results, pts):
        u = (p - centre) / SCALE
        rows.append([1.0, *u, *[u[i] * u[j] for i, j in pairs]])
        rhs.append(r['total_energy_Ry'] * RY_MEV)
        g = gradient(r) * SCALE * RY_MEV
        for k in range(3):
            rows.append([0.0, *[1.0 if m == k else 0.0 for m in range(3)],
                         *[(u[j] if i == k else 0.0) + (u[i] if j == k else 0.0) for i, j in pairs]])
            rhs.append(g[k])
    coef = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)[0]
    H = np.zeros((3, 3))
    for c, (i, j) in zip(coef[4:], pairs):
        H[i, j] += c
        H[j, i] += c
    return centre, coef[0], coef[1:4], H


def model_minimum(results):
    centre, e0, b, H = quadratic_model(results)
    u = np.linalg.solve(H, -b)
    return centre + u * SCALE, H


def newton_refine(result, H):
    """Structure and energy at the stationary point from one SCF (energy + gradient) and a model Hessian (scaled)."""
    p = np.array([result['settings'][k] for k in ('a_rh_A', 'alpha_deg', 'x_Te')])
    g = gradient(result) * SCALE * RY_MEV
    u = np.linalg.solve(H, -g)
    return p + u * SCALE, result['total_energy_Ry'] * RY_MEV + 0.5 * g @ u


def zero_pressure_a(results):
    """a_cub where P = 0, from a quadratic fit of P(a_cub)."""
    a = np.array([r['settings']['a_rh_A'] * S2 for r in results])
    p = np.array([r['pressure_GPa'] for r in results])
    c = np.polyfit(a, p, 2)
    roots = [r.real for r in np.roots(c) if abs(r.imag) < 1e-12 and a.min() - 0.1 < r.real < a.max() + 0.1]
    return float(roots[0])


def clamped_depth(cubic, polar):
    """E(d) = E0 + A d^2 + B d^4 (d = 0.5 - x) fitted to the cubic energy and the polar energies and forces."""
    rows, rhs = [[1.0, 0.0, 0.0]], [cubic['total_energy_Ry'] * RY_MEV]
    for r in polar:
        d = 0.5 - r['settings']['x_Te']
        rows.append([1.0, d * d, d ** 4])
        rhs.append(r['total_energy_Ry'] * RY_MEV)
        rows.append([0.0, 2 * d, 4 * d ** 3])           # dE/dd = -dE/dx
        rhs.append(-gradient(r)[2] * RY_MEV)
    e0, A, B = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)[0]
    if A >= 0:
        return 0.0, 0.5
    return float(A * A / (4 * B)), float(0.5 - math.sqrt(-A / (2 * B)))


def soft_mode_pressure(polar):
    """Pressure where the harmonic force constant of x at the rocksalt site changes sign (linear interpolation of the
    force constant k = (dE/dx)/(x - 0.5), from forces at a small displacement, against the pressure)."""
    pts = sorted((r['pressure_GPa'], gradient(r)[2] / (r['settings']['x_Te'] - 0.5)) for r in polar)
    for (p1, k1), (p2, k2) in zip(pts, pts[1:]):
        if k1 < 0 <= k2:
            return float(p1 - k1 * (p2 - p1) / (k2 - k1)), pts
    raise RuntimeError('no sign change of the soft-mode force constant: %s' % pts)


R3M_START = (4.36, 58.3, 0.470)
R3M_STEP = (0.03, 0.4, 0.005)
R3M_POINTS = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
EOS_A = (5.98, 6.00, 6.02)
SOFT_A = (5.80, 5.75, 5.70)


def r3m_grid_requests():
    return [req(R3M_START[0] + u * R3M_STEP[0], R3M_START[1] + v * R3M_STEP[1], R3M_START[2] + w * R3M_STEP[2], kmesh=12)
            for u, v, w in R3M_POINTS]


def rounded_minimum(p):
    return round(float(p[0]), 3), round(float(p[1]), 2), round(float(p[2]), 4)


def main():
    evidence = []

    def buy(request):
        job, result = calc(request)
        evidence.append(job)
        return result
    # 1. rocksalt equation of state (pressures)
    eos = [buy(req(a / S2, kmesh=16)) for a in EOS_A]
    a0 = round(zero_pressure_a(eos), 3)
    # 2. R3m surface at k12, stationary point, converged energy there
    grid = [buy(r) for r in r3m_grid_requests()]
    pmin, H = model_minimum(grid)
    ra, ral, rx = rounded_minimum(pmin)
    r3m = buy(req(ra, ral, rx, kmesh=20))
    p_ref, e_r3m = newton_refine(r3m, H)
    # 3. rocksalt at a0 on dense meshes
    cub20 = buy(req(a0 / S2, kmesh=20))
    cub24 = buy(req(a0 / S2, kmesh=24))
    depth = cub24['total_energy_Ry'] * RY_MEV - e_r3m
    # 4. clamped cubic cell at a0 (k20 like the R3m point; the k20-k24 cubic difference is applied as a constant)
    polar = [buy(req(a0 / S2, x=x, kmesh=20)) for x in (0.480, 0.485)]
    clamped, x_clamped = clamped_depth(cub20, polar)
    # 5. soft mode under compression (semimetal: cold smearing on a dense mesh)
    soft = [buy(req(a / S2, x=0.495, kmesh=20, occ='smearing', degauss=0.005)) for a in SOFT_A]
    p_c, k_table = soft_mode_pressure(soft)
    # 6. spin-orbit at fixed protocol structures, same PSlibrary family with and without SOC
    soc = {}
    for pp in ('pslibrary', 'pslibrary_fullrel'):
        c = buy(req(a0 / S2, kmesh=12, pp=pp))
        r = buy(req(ra, ral, rx, kmesh=12, pp=pp))
        soc[pp] = (c['total_energy_Ry'] - r['total_energy_Ry']) * RY_MEV
    soc_change = soc['pslibrary_fullrel'] - soc['pslibrary']
    submission = {'cubic_a0_A': a0, 'r3m_a_rh_A': float(p_ref[0]), 'r3m_alpha_deg': float(p_ref[1]),
                  'r3m_x_Te': float(p_ref[2]), 'well_depth_meV_per_fu': float(depth),
                  'clamped_well_depth_meV_per_fu': clamped, 'soft_mode_critical_pressure_GPa': p_c,
                  'soc_well_depth_change_meV_per_fu': float(soc_change), 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    fraction = 1 - clamped / depth
    (APP/'report.md').write_text(f"""# GeTe ferroelectric double well (reference solution)

Protocol: PBE, SSSP pseudopotentials (GBRV Ge with 3d in valence, PSlibrary Te), 50 Ry, scalar-relativistic.
Energies per formula unit (the two-atom cell is one GeTe unit).

- Rocksalt: P(a) from k16 SCFs gives a0 = {a0:.3f} A. Its PBE gap at L is only ~0.35 eV, so its energy converges
  slowly with k (k8 and k12 overestimate the well depth); the cubic energy is taken from k24, the R3m energy from
  k20 (k-converged to ~0.2 meV).
- R3m: quadratic model of E(a, alpha, x) through seven k12 SCFs using energies and full gradients (stress and
  forces), refined by a Newton step from a k20 SCF at the model minimum: a = {p_ref[0]:.4f} A, alpha = {p_ref[1]:.3f}
  deg, x = {p_ref[2]:.4f} (each phase at its own zero pressure, i.e. equal enthalpy reference at P = 0).
- Well depth E(rocksalt, a0) - E(R3m, relaxed) = {depth:.2f} meV/f.u.
- Clamped cubic cell (alpha = 60 deg, a0; only the Ge/Te displacement along [111]): depth {clamped:.2f} meV/f.u. at
  x = {x_clamped:.4f}, from E0 + A d^2 + B d^4 through energies and forces. The rhombohedral shear and the volume
  relaxation (strain coupling ~ eta P^2, symmetry-allowed because the polar T1u mode squared contains the shear) supply
  {100 * fraction:.0f} % of the relaxed depth: the instability is a proper polar one (the clamped cell is already
  unstable) that is strongly deepened by strain.
- Under hydrostatic compression the rocksalt gap at L closes (semimetal from ~6-7 GPa), so the soft-mode force
  constant is computed with cold smearing on k20; it changes sign at P = {p_c:.2f} GPa (force constants vs pressure:
  {k_table}). This is the harmonic instability pressure of rocksalt; whether the transition is continuous needs the
  quartic coefficient including strain relaxation, which these data only bound.
- Spin-orbit (PSlibrary scalar vs fully relativistic pair at the protocol structures): the depth changes by
  {soc_change:+.2f} meV/f.u. Comparing the SOC depth with the SSSP depth would mix the pseudopotential difference (Ge 3d
  frozen vs in valence) into the 'SOC effect'.
- Limits: static lattice, PBE overestimates the volume, no SOC in the protocol, the depth is not a transition
  temperature (the paired GAP study's ~700 K transition involves anharmonic, order-disorder-like fluctuations).
""")
    lab('submit', submission)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
