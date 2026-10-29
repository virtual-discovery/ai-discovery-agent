"""Reference solution for sngeo2-coherent-miscibility; runs as the Harbor oracle and in check_task."""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, fsolve

APP = Path(os.environ.get('APP_DIR', '/app'))
CAT = json.loads((APP/'catalog.json').read_text())
KB = 8.617333e-5
X7 = [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]
EVIDENCE = []


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(1)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(str(r))
    EVIDENCE.append(r['id'])
    return r['result']


def vol(lattice):
    return abs(np.linalg.det(np.array(lattice)))


# ---------------- free-energy model, eV per rutile primitive cell (2 cations) ----------------
def subregular(x, y):
    x = np.asarray(x, float)
    return np.linalg.lstsq(np.c_[x * (1 - x), x * (1 - x) * (1 - 2 * x)], np.asarray(y, float), rcond=None)[0]


def g(x, a, b, T):
    return a * x * (1 - x) + b * x * (1 - x) * (1 - 2 * x) + 2 * KB * T * (x * np.log(x) + (1 - x) * np.log(1 - x))


def dg(x, a, b, T):
    return a * (1 - 2 * x) + b * (1 - 6 * x + 6 * x * x) + 2 * KB * T * np.log(x / (1 - x))


def spinodal_T(x, a, b):
    return (2 * a - b * (12 * x - 6)) * x * (1 - x) / (2 * KB)


def spinodal_Tc(a, b):
    x = np.linspace(1e-5, 1 - 1e-5, 200001)
    T = spinodal_T(x, a, b)
    return float(max(T.max(), 0.0))


def spinodal_at(a, b, T):
    x = np.linspace(1e-6, 1 - 1e-6, 200001)
    inside = x[spinodal_T(x, a, b) > T]
    return [float(inside.min()), float(inside.max())] if len(inside) else None


def binodal_at(a, b, T):
    sp = spinodal_at(a, b, T)
    if sp is None:
        return None
    def eqs(z):
        x1, x2 = 1 / (1 + np.exp(-z))
        return [dg(x1, a, b, T) - dg(x2, a, b, T), g(x2, a, b, T) - g(x1, a, b, T) - dg(x1, a, b, T) * (x2 - x1)]
    best = None
    for f in (0.2, 0.5, 0.8, 0.95):
        x0 = np.array([sp[0] * (1 - f) + 1e-9, sp[1] + (1 - sp[1]) * f - 1e-9])
        z, info, ok, msg = fsolve(eqs, np.log(x0 / (1 - x0)), full_output=True, xtol=1e-13)
        x1, x2 = sorted(1 / (1 + np.exp(-z)))
        if ok == 1 and x2 - x1 > sp[1] - sp[0] - 1e-9 and max(abs(np.array(eqs(z)))) < 1e-9:
            best = [float(x1), float(x2)]
            break
    return best


def binodal_Tc(a, b):
    lo, hi = 100.0, 1.05 * spinodal_Tc(a, b) + 1
    for _ in range(40):
        mid = (lo + hi) / 2
        if spinodal_at(a, b, mid) is not None and binodal_at(a, b, mid) is not None:
            lo = mid
        else:
            hi = mid
    return lo


# ---------------- coherency strain from the released scans ----------------
def scan_min(points, member, plane):
    """minimum over the stepped lattice parameter: parabola through the 5 points around the lowest energy"""
    def stepped(A):
        A = np.array(A)
        if plane == '001':
            return A[0, 2]
        if plane == '100':
            return np.linalg.norm(A[2])
        return np.linalg.norm(A[1] + A[2])
    pts = sorted((stepped(p['lattice_A']), p['energy_eV']) for p in points)
    s = np.array(pts); i = int(np.argmin(s[:, 1]))
    if i in (0, len(s) - 1):
        return float(s[i, 1])
    lo = max(0, min(i - 2, len(s) - 5)); w = s[lo:lo + 5]
    p = np.polyfit(w[:, 0], w[:, 1], 2)
    return float(np.polyval(p, -p[1] / (2 * p[0])))


def strain_energy(fn, plane, x):
    """Delta e_cs(x) per cell: end members held at the alloy's in-plane lattice (Vegard), relaxed out of plane"""
    E = CAT['end_members'][fn]
    a = x * E['GeO2']['a_A'] + (1 - x) * E['SnO2']['a_A']
    c = x * E['GeO2']['c_A'] + (1 - x) * E['SnO2']['c_A']
    ds = CAT['strain_scans']['datasets'][fn][plane]
    grids = list(ds['fixed'].values())
    if plane == '001':
        sid = '%02d' % int(np.argmin(np.abs(np.array(grids[0]) - a)))
    else:
        first = a if plane == '100' else np.sqrt(2) * a
        sid = '%02d-%02d' % (int(np.argmin(np.abs(np.array(grids[0]) - first))), int(np.argmin(np.abs(np.array(grids[1]) - c))))
    e = 0.0
    for member, w in (('GeO2', x), ('SnO2', 1 - x)):
        res = buy({'kind': 'strain_scan', 'functional': fn, 'id': f'{plane}-{member}-{sid}'})
        e += w * (scan_min(res['points'], member, plane) - E[member]['energy_eV_per_6atom_cell'])
    return e


def main():
    log = []
    # 1. Normalisation of the released SQS energies: compare with total energies run in one consistent setting.
    hse = {m: buy({'kind': 'hse_setting_run', 'name': f'{m}_PBE_values'}) for m in ('GeO2', 'SnO2')}
    alloy = buy({'kind': 'hse_setting_run', 'name': 'Alloy_PBE_Values/Sn_50_1'})
    nG, nS = alloy['species'].count('Ge'), alloy['species'].count('Sn')
    ef_cation = (alloy['energy_sigma0_eV'] - nG * hse['GeO2']['energy_sigma0_eV'] / 2 - nS * hse['SnO2']['energy_sigma0_eV'] / 2) / (nG + nS)
    # 2. PBE SQS
    pbe = []
    for x in X7:
        pbe.append(buy({'kind': 'sqs', 'functional': 'PBE', 'x_Ge': x, 'index': 0}))
    for i in (1, 2):
        pbe.append(buy({'kind': 'sqs', 'functional': 'PBE', 'x_Ge': 0.5, 'index': i}))
    released05 = [r['energy_as_released'] for r in pbe if r['x_Ge'] == 0.5 and r['index'] == 0][0]
    ratio = released05 / ef_cation          # ~2: the released value is per primitive cell (2 cations)
    per_cell = abs(ratio - 2) < 0.2
    assert per_cell, ratio
    # 3. r2SCAN SQS, screened against the Vegard volume line of the end members
    E = CAT['end_members']
    def vegard_dev(r, fn):
        v0 = {m: E[fn][m]['a_A'] ** 2 * E[fn][m]['c_A'] / 2 for m in E[fn]}
        return vol(r['lattice_A']) / 24 - (r['x_Ge'] * v0['GeO2'] + (1 - r['x_Ge']) * v0['SnO2'])
    r2 = []
    rejected = []
    for x in X7:
        for i in (0, 1, 2):
            r = buy({'kind': 'sqs', 'functional': 'r2SCAN', 'x_Ge': x, 'index': i})
            if abs(vegard_dev(r, 'r2SCAN')) < 0.3:
                r2.append(r)
                break
            rejected.append((x, i, round(vegard_dev(r, 'r2SCAN'), 2), r['energy_as_released']))
    fits = {}
    for fn, rows in (('PBE', pbe), ('r2SCAN', r2)):
        fits[fn] = subregular([r['x_Ge'] for r in rows], [r['energy_as_released'] for r in rows])
    a, b = fits['PBE']
    dh05 = float(np.mean([r['energy_as_released'] for r in pbe if r['x_Ge'] == 0.5])) / 2 * 1000
    # 4. coherency strain: all three released planes at x = 0.5 (PBE), (001) composition dependence for both functionals
    strain05 = {pl: strain_energy('PBE', pl, 0.5) / 2 * 1000 for pl in ('001', '100', '101')}
    habit = min(strain05, key=strain05.get)
    coh = {}
    for fn in ('PBE', 'r2SCAN'):
        xs = [0.25, 0.5, 0.75]
        ec = [strain05[habit] * 2 / 1000 if (fn == 'PBE' and x == 0.5) else strain_energy(fn, habit, x) for x in xs]
        ac, bc = subregular(xs, ec)
        coh[fn] = (fits[fn][0] - ac, fits[fn][1] - bc)
    tcs = {fn: spinodal_Tc(*coh[fn]) for fn in coh}
    bi, sp = binodal_at(*coh['PBE'], 700.0), spinodal_at(*coh['PBE'], 700.0)
    tinc = binodal_Tc(a, b)
    submission = {'mixing_enthalpy_x050_meV_per_cation': dh05,
                  'coherency_strain_energy_x050_meV_per_cation': strain05,
                  'habit_plane': habit,
                  'incoherent_Tc_K': tinc,
                  'coherent_spinodal_Tc_K': tcs,
                  'coherent_boundaries_700K_PBE': {'binodal': bi, 'spinodal': sp},
                  'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# Rutile Sn1-xGexO2: incoherent vs coherent miscibility (reference solution)

Units: x = Ge fraction; energies per cation unless stated; PBE throughout except the r2SCAN coherent Tc.

* Normalisation. The archive's SQS file holds total energies for x = 0 and 1 but small positive numbers for 0 < x < 1.
  PBE relaxations of the same SQS in one consistent setting (hybrid-functional settings) give a formation energy of
  {ef_cation*1000:.1f} meV/cation at x = 0.5; the released value is {released05:.4f} eV, a ratio of {ratio:.2f}, so the released
  numbers are formation energies per rutile primitive cell (2 cations), as in the paper's Eq. 1.
* Disordered alloy (PBE SQS, subregular fit per cell): alpha = {a:.3f} eV, beta = {b:.3f} eV; mixing enthalpy at x = 0.5
  {dh05:.1f} meV/cation, positive at every composition. The 80 small-cell orderings all have positive formation energies:
  no ordered ground state; the alloy is clustering-type (miscibility gap), not ordering-type.
* r2SCAN SQS screened against the end-member Vegard volume line: rejected {rejected} (volume about 0.7 A^3/cation too
  small, energies about 20 meV/cation above the PBE-scaled trend: relaxations that did not reach the equilibrium volume).
  Same decorations as the PBE SQS, which sit on the Vegard line.
* Incoherent miscibility gap (SQS enthalpy + ideal entropy, mean field): Tc = {tinc:.0f} K. Short-range order lowers it
  (cluster-expansion Monte Carlo gives roughly 10 % lower, and the value depends on the cluster truncation), but the
  gap stays far above growth temperatures: equilibrium (incoherent) solubility at 800-1000 K is only a few per cent.
* Coherency strain (Cahn approximation): end members held at the homogeneous alloy's in-plane lattice (Vegard, both
  in-plane parameters) and relaxed along the stepped out-of-plane parameter (parabolic minimum of each scan), relative
  to the relaxed end members. At x = 0.5: {', '.join(f'{k}: {v:.1f}' for k, v in strain05.items())} meV/cation. Habit
  plane = {habit} (the large c-axis mismatch relaxes out of plane). The archive's "101" dataset strains a2 - a3 and c,
  i.e. a (110)-type plane in standard rutile setting.
* Coherent spinodal (h_dis - e_cs - T s_ideal): Tc = {tcs['PBE']:.0f} K (PBE), {tcs['r2SCAN']:.0f} K (r2SCAN, screened SQS;
  keeping the unrelaxed SQS would give about 900 K). At 700 K (PBE): coherent binodal {bi}, spinodal {sp}.
* Explanation: the incoherent gap (Tc ~ 2500-2700 K) implies films grown at 800-1000 K are not equilibrium solid
  solutions; coherency strain lowers the gap to below ~900 K, so coherent (epitaxial) films can stay homogeneous and are
  metastable. A fixed substrate strain is linear in x and does not change the curvature, so it cannot explain it.
  Not supported by these data: vibrational entropy, interface energies, dislocations (semi-coherent interfaces),
  kinetics and amorphisation; the coherent result is an upper-bound-type estimate from end-member elasticity.
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
