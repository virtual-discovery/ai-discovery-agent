"""Reference solution for fresh-h-chain-peierls-vs-afm; runs as the Harbor oracle and in check_task.

It requests new pw.x calculations from the lab, so it takes about an hour. The estimator functions
at module level are also used by the construction script on the reference calculations.
"""
import json
import math
import os
import subprocess
import time
from pathlib import Path

APP = Path(os.environ.get('APP_DIR', '/app'))
RY_MEV = 13605.693122994
STANDARD = dict(k_points=48, smearing='mv', degauss_Ry=0.005, box_A=10, ecutwfc_Ry=80)
CONVERGED = dict(STANDARD, degauss_Ry=0.002)       # with k_points set per request
AFM = [0.5, -0.5]
MAGNETIC_SITE_MOMENT = 0.002                         # muB; site moments are printed to 1e-4


# ---------- estimators (pure functions of lab results) ----------
def energy(res):
    """meV per atom."""
    return res['total_energy_Ry_per_cell'] / res['atoms_per_cell'] * RY_MEV


def site_moment(res):
    sites = res.get('site_moments_muB') or [0.0]
    return sum(abs(m) for m in sites) / len(sites)


def is_magnetic(res):
    return site_moment(res) > MAGNETIC_SITE_MOMENT


def parabola(points):
    """Coefficients (a, b, c) of a + b x + c x^2 through three (x, y) points."""
    (x0, y0), (x1, y1), (x2, y2) = points
    c = ((y2 - y0) / (x2 - x0) - (y1 - y0) / (x1 - x0)) / (x2 - x1)
    b = (y1 - y0) / (x1 - x0) - c * (x0 + x1)
    return y0 - b * x0 - c * x0 ** 2, b, c


def equilibrium(scan, R):
    """scan: list of (delta, result) for nonmagnetic 2-atom chains bracketing the minimum.

    delta_eq is the root of dE/d(delta) = -R F1 / 2 (per atom) from a parabola through the three points
    nearest the minimum; the minimum energy is the energy parabola evaluated there.
    """
    pts = sorted(scan, key=lambda p: energy(p[1]))[:3]
    pts.sort(key=lambda p: p[0])
    slope = [(d, -R * res['forces_eV_per_A_along_chain'][0] / 2) for d, res in pts]
    a, b, c = parabola(slope)
    roots = [r for r in ([-a / b] if abs(c) < 1e-12 else
                         [(-b + s * math.sqrt(max(b * b - 4 * a * c, 0))) / (2 * c) for s in (1, -1)])
             if pts[0][0] - 0.02 <= r <= pts[-1][0] + 0.02]
    d_eq = min(roots, key=lambda r: abs(r - pts[1][0]))
    ea, eb, ec = parabola([(d, energy(res)) for d, res in pts])
    return d_eq, ea + eb * d_eq + ec * d_eq ** 2


def critical_delta(points):
    """points: list of (delta, result) from AFM starts at one R. Returns the delta where the moment vanishes.

    Mean-field (gap-matching) form near the boundary: m^2 = A - B delta^2, fitted to the two magnetic
    points with the largest delta and kept inside the bracket set by the nearest nonmagnetic point.
    """
    mag = sorted((d, site_moment(r)) for d, r in points if is_magnetic(r))
    non = sorted(d for d, r in points if not is_magnetic(r))
    if not mag or mag[0][0] > 0:
        return 0.0
    upper = min([d for d in non if d > mag[-1][0]], default=None)
    if len(mag) < 2:
        return mag[-1][0] if upper is None else 0.5 * (mag[-1][0] + upper)
    (d1, m1), (d2, m2) = mag[-2], mag[-1]
    B = (m1 ** 2 - m2 ** 2) / (d2 ** 2 - d1 ** 2)
    dc = math.sqrt((m2 ** 2 + B * d2 ** 2) / B) if B > 0 else float('inf')
    if upper is not None and not d2 <= dc <= upper:
        dc = 0.5 * (d2 + upper)
    return dc


# ---------- lab access ----------
def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


EVIDENCE = []


def scf(R, delta, mag, **settings):
    request = {'kind': 'scf', 'bond_lengths_A': [round(R * (1 - delta), 6), round(R * (1 + delta), 6)],
               'starting_magnetization': mag, **dict(STANDARD, **settings)}
    job = lab('evaluate', request)
    while job['status'] == 'running':
        time.sleep(5)
        job = lab('job/%d' % job['id'])
    if job['status'] != 'complete':
        raise RuntimeError('calculation failed: %s' % job.get('error'))
    EVIDENCE.append(job['id'])
    return job['result']


def boundary(R, uniform_afm, first):
    """Bracket the critical dimerization with AFM starts on a 0.01 grid, then interpolate m^2 (3-5 calculations)."""
    step = 0.01
    on_grid = lambda x: round(round(x / step) * step, 4)
    pts = [(0.0, uniform_afm), (first, scf(R, first, AFM))]
    d = on_grid(critical_delta(pts + [(1.0, {'site_moments_muB': [0.0]})]))
    for _ in range(4):
        if d <= first or any(abs(d - p[0]) < 1e-9 for p in pts):
            break
        pts.append((d, scf(R, d, AFM)))
        mag = [p[0] for p in pts if is_magnetic(p[1])]
        non = [p[0] for p in pts if not is_magnetic(p[1])]
        if non and min(non) - max(mag) <= step + 1e-9:
            break
        d = on_grid(d + step) if is_magnetic(pts[-1][1]) else on_grid(d - step)
    return critical_delta(pts), pts


def main():
    out, notes = {}, {}
    # R = 1.25 A: weak-coupling regime. Uniform NM vs AFM, the dimerized minimum, and the converged SDW gap.
    nm125, afm125 = scf(1.25, 0, None), scf(1.25, 0, AFM)
    scan125 = [(d, scf(1.25, d, None)) for d in (0.39, 0.40, 0.41)]
    gap96 = scf(1.25, 0, AFM, **dict(CONVERGED, k_points=96))
    gap192 = scf(1.25, 0, AFM, **dict(CONVERGED, k_points=192))
    # R = 1.75 A: strong magnetism in the uniform chain; boundary of the magnetic window.
    nm175, afm175 = scf(1.75, 0, None), scf(1.75, 0, AFM)
    scan175 = [(d, scf(1.75, d, None)) for d in (0.56, 0.57, 0.58)]
    dc175, pts175 = boundary(1.75, afm175, 0.10)
    # R = 2.0 A
    nm200, afm200 = scf(2.0, 0, None), scf(2.0, 0, AFM)
    dc200, pts200 = boundary(2.0, afm200, 0.15)
    dim200_nm, dim200_afm = scf(2.0, 0.625, None), scf(2.0, 0.625, AFM)
    # R = 1.0 A
    afm100, dim100_afm = scf(1.0, 0, AFM), scf(1.0, 0.25, AFM)

    d125, e125 = equilibrium(scan125, 1.25)
    d175, e175 = equilibrium(scan175, 1.75)
    uniform125 = min(energy(nm125), energy(afm125))
    uniform175 = min(energy(nm175), energy(afm175))

    def label(uniform_states, dimer):
        best_uniform = min(uniform_states, key=energy)
        if energy(dimer) < energy(best_uniform):
            return 'dimerized_magnetic' if is_magnetic(dimer) else 'dimerized_nonmagnetic'
        return 'uniform_antiferromagnetic' if is_magnetic(best_uniform) else 'uniform_nonmagnetic'

    submission = {
        'delta_eq': {'1.25': d125, '1.75': d175},
        'dimerization_energy_meV_per_atom': {'1.25': uniform125 - e125, '1.75': uniform175 - e175},
        'uniform_magnetic_energy_gain_meV_per_atom': {'1.75': energy(nm175) - energy(afm175),
                                                      '2.0': energy(nm200) - energy(afm200)},
        'critical_dimerization': {'1.75': dc175, '2.0': dc200},
        'uniform_chain_band_gap_eV': {'1.25': gap192['band_gap_eV']},
        'ground_state': {'1.0': label([afm100], dim100_afm),
                         '2.0': label([nm200, afm200], min((dim200_nm, dim200_afm), key=energy))},
        'evidence_ids': list(EVIDENCE)}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    moments175 = ', '.join(f'{d:.3f}: {site_moment(r):.4f}' for d, r in sorted(pts175))
    moments200 = ', '.join(f'{d:.3f}: {site_moment(r):.4f}' for d, r in sorted(pts200))
    (APP/'report.md').write_text(f"""# Hydrogen chain in PBE: dimerization versus magnetism (reference solution)

**Setup.** 2-atom cells, bonds R(1-delta) and R(1+delta); energies per atom (total energy per cell / 2).
Standard settings: 48 k-points, Marzari-Vanderbilt 0.005 Ry, 10 A box, 80 Ry. A ferromagnetic start
(equal moments) is not enough: in this chain it relaxes to the nonmagnetic state, so magnetic solutions
were sought with antiparallel starting moments (+0.5, -0.5), the 2k_F = pi/R spin-density wave.

**Uniform chain.** The antiferromagnetic (SDW) solution exists at every R studied. At R = 1.25 A it is
weak: site moment {site_moment(gap192):.4f} muB, gain below 1 meV/atom. Its gap is set by nesting, so it needs
dense k-points and small smearing: {afm125['band_gap_eV']:.3f} eV at 48 k / 0.005 Ry, {gap96['band_gap_eV']:.3f} eV at 96 k /
0.002 Ry and {gap192['band_gap_eV']:.3f} eV at 192 k / 0.002 Ry (reported). Larger smearing (0.02 Ry MV, or
Fermi-Dirac 0.005 Ry, an electronic temperature above the mean-field ordering temperature) removes it, which
is why it is easy to mistake for an artefact; it strengthens, rather than disappears, as the smearing is
reduced. At 1.75 and 2.0 A the gain is large ({energy(nm175) - energy(afm175):.1f} and
{energy(nm200) - energy(afm200):.1f} meV/atom) and the gapped state is insensitive to the settings.

**Dimerization.** At fixed mean spacing the chain dimerizes all the way to H2-like molecules
(delta_eq = {d125:.4f} at 1.25 A and {d175:.4f} at 1.75 A, short bond about 0.75 A), gaining
{uniform125 - e125:.0f} and {uniform175 - e175:.0f} meV/atom relative to the lowest uniform chain (the AFM chain;
at 1.75 A using the nonmagnetic chain instead would overstate the gain by {energy(nm175) - energy(afm175):.0f} meV/atom).
The minimum was located from the zero of the force along the dimerization coordinate.

**Coexistence.** Magnetism survives in the dimerized chain only up to delta_c = {dc175:.3f} (R = 1.75 A) and
{dc200:.3f} (R = 2.0 A). Site moments (delta: muB): R = 1.75: {moments175}; R = 2.0: {moments200}.
The moment vanishes continuously, roughly as m^2 = A - B delta^2: the SDW survives while the Peierls gap of the
dimerized chain stays below the SDW gap of the uniform chain. The window grows with R but never reaches the
equilibrium dimerization, so the global minimum at every R from 1.0 to 2.0 A is the nonmagnetic dimerized
chain ({submission['ground_state']}). Antiferromagnetism is a property of the constrained (uniform or weakly
dimerized) chain only; dimerization and magnetism coexist only away from equilibrium.

**Interpretation and limits.** A broken-symmetry PBE solution is a mean-field way of representing static
correlation (growing U/t with R); it is not long-range order, which a 1D chain cannot have (Mermin-Wagner,
quantum fluctuations; the half-filled Hubbard chain is a Mott insulator with power-law spin correlations).
The results hold for a fixed collinear 1D geometry with clamped protons: no zero-point motion (large for H),
no noncollinear or longer-period orders beyond the 2-atom cell, PBE self-interaction and delocalization
errors; correlated benchmarks of the hydrogen chain (Motta et al. 2017, 2020) are needed for quantitative
statements about the uniform chain.
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
