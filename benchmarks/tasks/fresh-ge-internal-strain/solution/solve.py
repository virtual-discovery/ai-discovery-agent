"""Reference solution for fresh-ge-internal-strain; runs as the Harbor oracle and in check_task.

It requests new pw.x calculations from the lab (about 20 of them, ~35 minutes). The estimator functions at
module level are also applied by the construction script to the reference calculations.
"""
import json
import math
import os
import subprocess
import time
from pathlib import Path

APP = Path(os.environ.get('APP_DIR', '/app'))
BOHR_A = 0.529177210903
RY_BOHR3_GPA = 14710.507848260711          # 1 Ry/bohr^3 in GPa
AMU_KG, C_CM_S = 1.66053906660e-27, 2.99792458e10
EV_A2_SI = 16.021766340                    # 1 eV/A^2 in N/m
RY_BOHR2_EV_A2 = 13.605693122994 / BOHR_A ** 2
MASS_AMU = 72.63
GBRV = 'ge_pbe_v1.4.uspp.F.UPF'
ZERO = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
P_TARGET = 6.0                              # GPa
E4 = 0.01                                   # engineering shear strain
EXX = 0.01                                  # normal strain for C11, C12
U_PHONON = 0.01                             # A, relative displacement along x for the optical mode
U_SCAN = -0.012                             # A, trial relative displacement in the sheared cell


# ---------------------------------------------------------------- estimators (pure functions of lab results)
def shear(e4):
    return [[0, 0, 0], [0, 0, e4 / 2], [0, e4 / 2, 0]]


def normal(exx):
    return [[exx, 0, 0], [0, 0, 0], [0, 0, 0]]


def cauchy(res):
    """Cauchy stress (GPa, tensile positive) from the lab's pw.x-convention stress."""
    return [[-v for v in row] for row in res['stress_GPa']]


def omega_bohr3(a):
    return a ** 3 / 4 / BOHR_A ** 3


def root(points, target):
    """x where y(x) = target, from the quadratic (or line) through the points nearest the target."""
    pts = sorted(points, key=lambda p: abs(p[1] - target))[:3]
    pts.sort()
    if len(pts) == 2:
        (x0, y0), (x1, y1) = pts
        return x0 + (target - y0) * (x1 - x0) / (y1 - y0)
    (x0, y0), (x1, y1), (x2, y2) = pts
    c = ((y2 - y0) / (x2 - x0) - (y1 - y0) / (x1 - x0)) / (x2 - x1)
    b = (y1 - y0) / (x1 - x0) - c * (x0 + x1)
    a0 = y0 - b * x0 - c * x0 ** 2 - target
    if abs(c) < 1e-12:
        return -a0 / b
    xs = [(-b + s * math.sqrt(max(b * b - 4 * a0 * c, 0.0))) / (2 * c) for s in (1, -1)]
    return min(xs, key=lambda x: abs(x - pts[1][0]))


def lattice_at_pressure(runs, target):
    """runs: list of cubic scf results; returns the lattice parameter where pressure_GPa = target."""
    return root([(r['settings']['lattice_A'], r['pressure_GPa']) for r in runs], target)


def gamma_levels(res):
    """(E(Gamma2'), E(Gamma25')) from the band energies at the Gamma point of a Gamma-containing mesh."""
    point = res['mesh_point_nearest_gamma']
    assert max(abs(k) for k in point['k_cartesian_inverse_A']) < 1e-8, 'needs a mesh containing Gamma'
    levels = sorted(e for e in point['band_energies_eV'] if e > -10.0)[1:]      # drop the 3d shell and Gamma1
    groups = []
    for e in levels:
        if groups and abs(e - groups[-1][-1]) < 2e-3:
            groups[-1].append(e)
        else:
            groups.append([e])
    triplets = [g for g in groups if len(g) == 3]
    singlets = [g for g in groups if len(g) == 1 and g[0] < triplets[1][0]]
    return singlets[0][0], sum(triplets[0]) / 3


def gamma_gap(res):
    s, t = gamma_levels(res)
    return s - t


def interpolate(points, x):
    """Linear interpolation (or extrapolation) through the two points nearest x."""
    (x0, y0), (x1, y1) = sorted(sorted(points, key=lambda p: abs(p[0] - x))[:2])
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def internal_strain(e4_res, ph_res):
    """Clamped shear and zone-centre force constant at one lattice parameter.

    Lambda: force on atom 2 per unit engineering shear (Ry/bohr), from the clamped sheared cell;
    Phi: force constant of the relative coordinate (Ry/bohr^2), from the displacement along x.
    """
    a = e4_res['settings']['lattice_A']
    e4 = 2 * e4_res['settings']['strain'][1][2]
    u = ph_res['settings']['displacement_A'][0] / BOHR_A
    lam = -e4_res['forces_Ry_per_bohr'][1][0] / e4
    phi = -ph_res['forces_Ry_per_bohr'][1][0] / u
    c44_0 = cauchy(e4_res)[1][2] / e4
    soft = lam ** 2 / (omega_bohr3(a) * phi) * RY_BOHR3_GPA
    return {'Lambda_Ry_per_bohr': lam, 'Phi_Ry_per_bohr2': phi, 'C44_clamped': c44_0, 'C44_relaxed': c44_0 - soft,
            'zeta': 4 * lam / (a / BOHR_A * phi)}


def relaxed_from_scan(e4_res, scan_res):
    """Zero-force point of a displacement scan in the sheared cell (linear interpolation of force and stress)."""
    a = e4_res['settings']['lattice_A']
    e4 = 2 * e4_res['settings']['strain'][1][2]
    u1 = scan_res['settings']['displacement_A'][0]
    f0, f1 = e4_res['forces_Ry_per_bohr'][1][0], scan_res['forces_Ry_per_bohr'][1][0]
    s0, s1 = cauchy(e4_res)[1][2], cauchy(scan_res)[1][2]
    ustar = -f0 * u1 / (f1 - f0)
    return {'C44_relaxed': (s0 + (s1 - s0) * ustar / u1) / e4, 'zeta': -4 * ustar / (a * e4)}


def frequency_cm(phi):
    k = phi * RY_BOHR2_EV_A2 * EV_A2_SI
    return math.sqrt(k / (MASS_AMU / 2 * AMU_KG)) / (2 * math.pi * C_CM_S)


def normal_constants(plus, minus):
    e = plus['settings']['strain'][0][0]
    sp, sm = cauchy(plus), cauchy(minus)
    return (sp[0][0] - sm[0][0]) / (2 * e), (sp[1][1] - sm[1][1]) / (2 * e)


def gruneisen(ph_minus, ph_plus):
    am, ap = ph_minus['settings']['lattice_A'], ph_plus['settings']['lattice_A']
    wm = frequency_cm(-ph_minus['forces_Ry_per_bohr'][1][0] / (ph_minus['settings']['displacement_A'][0] / BOHR_A))
    wp = frequency_cm(-ph_plus['forces_Ry_per_bohr'][1][0] / (ph_plus['settings']['displacement_A'][0] / BOHR_A))
    return -(math.log(wp) - math.log(wm)) / (3 * math.log(ap / am))


# ---------------------------------------------------------------- lab access
def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


EVIDENCE = []


def calc(a, strain=ZERO, disp=(0, 0, 0), ecut=45, k=12, shift=0, kind='scf', thr=None):
    request = {'kind': kind, 'pseudopotential': GBRV, 'lattice_A': round(a, 4), 'strain': strain,
               'displacement_A': list(disp), 'ecutwfc_Ry': ecut, 'kmesh': k, 'kshift': shift,
               'occupations': 'fixed', 'degauss_Ry': None, 'forc_conv_thr_Ry_per_bohr': thr}
    job = lab('evaluate', request)
    while job['status'] == 'running':
        time.sleep(5)
        job = lab('job/%d' % job['id'])
    if job['status'] != 'complete':
        raise RuntimeError('calculation failed: %s' % job.get('error'))
    EVIDENCE.append(job['id'])
    return job['result']


def main():
    # Zero pressure: pressures at 60 Ry (the stress converges more slowly with the cutoff than the energy),
    # on a Gamma-containing mesh so that the band ordering at Gamma is visible.
    eos0 = [calc(a, ecut=60) for a in (5.75, 5.77)]
    a0 = lattice_at_pressure(eos0, 0.0)
    gap0 = interpolate([(r['settings']['lattice_A'], gamma_gap(r)) for r in eos0], a0)
    # 6 GPa
    eos6 = [calc(a, ecut=60) for a in (5.59, 5.61)]
    a6 = round(lattice_at_pressure(eos6, P_TARGET), 4)
    gap6 = interpolate([(r['settings']['lattice_A'], gamma_gap(r)) for r in eos6], a6)
    # Internal strain at 6 GPa on a Gamma-centred 12x12x12 mesh (shifted meshes miss the small direct gap at
    # Gamma, where the internal-strain response is concentrated; Gamma-centred k8/k12/k16 agree to <0.5 %).
    e4_12 = calc(a6, strain=shear(E4), k=12)
    ph_12 = calc(a6, disp=(U_PHONON, 0, 0), k=12)
    scan = calc(a6, strain=shear(E4), disp=(U_SCAN, 0, 0), k=12)
    I = internal_strain(e4_12, ph_12)
    R = relaxed_from_scan(e4_12, scan)
    # Maxwell check: Omega d(sigma_yz)/du from the displaced cell must equal Lambda from the sheared cell
    maxwell = cauchy(ph_12)[1][2] / RY_BOHR3_GPA * omega_bohr3(a6) / (U_PHONON / BOHR_A)
    exp_, exm = calc(a6, strain=normal(EXX)), calc(a6, strain=normal(-EXX))
    c11, c12 = normal_constants(exp_, exm)
    # +-0.5 % in a: on the expanded side the direct gap at Gamma closes quickly (0.31 eV at a6, ~0.03 eV at
    # +1 %) and the optical force constant softens anomalously, so larger steps are not local derivatives.
    phm, php = calc(a6 * 0.995, disp=(U_PHONON, 0, 0)), calc(a6 * 1.005, disp=(U_PHONON, 0, 0))
    gamma = gruneisen(phm, php)
    omega = frequency_cm(I['Phi_Ry_per_bohr2'])
    submission = {'a0_A': a0, 'gamma_gap_at_a0_eV': gap0, 'a_6GPa_A': a6, 'C11_GPa': c11, 'C12_GPa': c12,
                  'C44_clamped_ion_GPa': I['C44_clamped'], 'C44_relaxed_ion_GPa': I['C44_relaxed'],
                  'kleinman_zeta': I['zeta'], 'raman_frequency_cm-1': omega, 'raman_gruneisen': gamma,
                  'evidence_ids': list(EVIDENCE)}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# Ge internal strain under 6 GPa (reference solution)

**Protocol.** PBE, GBRV ultrasoft Ge with 3d in valence (the 3d shell is part of the converged core-valence
partition; the frozen-3d file shifts pressures by ~0.6 GPa and the optical force constant by ~4 %), fixed
occupations. Pressures were taken at 60 Ry: the stress of this pseudopotential carries a Pulay error of
~0.15 GPa at 45 Ry (and >1 GPa at 35 Ry), while forces and shear stresses are converged at 45 Ry.

**Zero pressure.** a0 = {a0:.4f} A from P(a) = 0. On a Gamma-containing mesh the s-like Gamma2' level lies
{-gap0:.3f} eV *below* the p-like Gamma25' triplet: PBE Ge without spin-orbit coupling is a band-inverted
semimetal at its own lattice parameter (a shifted mesh hides this and reports a small positive gap; a
Gamma-centred mesh with fixed occupations fills two thirds of the degenerate triplet). In that state the
internal-strain response has Fermi-surface contributions and does not converge with the lab's meshes, so the
elastic and Raman analysis is done at 6 GPa, where the direct gap at Gamma is open ({gap6:+.3f} eV).

**6 GPa.** a(6 GPa) = {a6:.4f} A. Clamped shear e4 = {E4} (pure shear, engineering strain, Cauchy stress):
C44^0 = {I['C44_clamped']:.2f} GPa; internal-strain force Lambda = {I['Lambda_Ry_per_bohr']:.4f} Ry/bohr per unit e4;
optical force constant of the relative coordinate Phi = {I['Phi_Ry_per_bohr2']:.5f} Ry/bohr^2, from a
displacement along x (no cubic anharmonic term along a cube axis). Maxwell check: Omega d(sigma_yz)/du_x =
{maxwell:.4f} Ry/bohr vs Lambda = {I['Lambda_Ry_per_bohr']:.4f}. Shifted meshes miss the small direct gap at Gamma and give
Lambda several per cent too small even at 12x12x12 (the Maxwell relation then fails), so Gamma-centred meshes
were used.
Closure: C44 = C44^0 - Lambda^2/(Omega Phi) = {I['C44_relaxed']:.2f} GPa and zeta = 4 Lambda/(a Phi) = {I['zeta']:.4f}
(u_x = -zeta a e4/4). The displacement scan in the sheared cell gives C44 = {R['C44_relaxed']:.2f} GPa and zeta = {R['zeta']:.4f}.
Internal relaxation removes {100 * (1 - I['C44_relaxed'] / I['C44_clamped']):.0f} % of the clamped shear stiffness: a large,
bond-bending (Keating-like) relaxation, not an artefact.
C11 = {c11:.1f}, C12 = {c12:.1f} GPa from +-{EXX} normal strains (stress-strain coefficients at 6 GPa; an
energy-strain fit must include the -P dV work term).
Raman (zone-centre optical) frequency with the reduced mass M/2: {omega:.1f} cm^-1; gamma = {gamma:.3f} from
+-0.5 % in a with identical settings. Steps of 1 % are not local: on the expanded side the direct gap at Gamma
almost closes (about 0.03 eV at +1 % in a) and the optical mode softens anomalously, so the Grueneisen parameter
of this protocol is only defined locally and changes rapidly on decompression.

**Limits.** Static lattice, no spin-orbit coupling (which would reduce the inversion at Gamma), harmonic
small-strain limit, one functional and pseudopotential; the zero-pressure internal-strain parameters of this
protocol are not supported by these calculations. Experiment (a = 5.658 A, C44 = 67 GPa, omega = 301 cm^-1 at
ambient pressure) is not directly comparable to the 6 GPa values.
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
