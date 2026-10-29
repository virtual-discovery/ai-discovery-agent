"""Reference solution for cr2o3-magnetoelectric-sign; runs as the Harbor oracle and in check_task.

Route (no literature values): force-constant eigenmodes (units checked against the all-H phonons), Born charges
from +0.02 A displacements (displacements read from the OUTCAR positions, linearity checked at -0.02 A), induced
moments from Elk (+-0.01, +-0.02 A) and VASP (+-0.025 A), domain of each dataset from its converged Cr moments,
then alpha = mu0/Omega sum_m (dm/du_m)(p_m/k_m).
"""
import json
import os
import re
import subprocess
import time
from pathlib import Path

import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
MODES = ['Eu_1_x', 'Eu_2_x', 'Eu_3_x', 'Eu_4_x']
ATOMS = ['Cr1', 'Cr2', 'Cr3', 'Cr4', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6']
MU0, MUB = 1.25663706212e-6, 9.2740100783e-24
AMU, EV = 1.66053906660e-27, 1.602176634e-19
EVIDENCE = []


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(req):
    r = lab('evaluate', req)
    while r['status'] == 'running':
        time.sleep(1)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(str(r))
    EVIDENCE.append(r['id'])
    return r['result']


def vec(line, label):
    m = re.search(label + r'=\(\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s*\)', line)
    return np.array([float(x) for x in m.groups()]) if m else None


def dipole(run):
    """Physical dipole (e*A). VASP prints p in 'electrons Angst': p_ion = -ZVAL*r, so the charge unit is -|e|."""
    p = {}
    for l in run['OUTCAR_dipole_lines']:
        for lab_ in ('p\\[ion\\]', 'p\\[elc\\]'):
            v = vec(l, lab_)
            if v is not None:
                p[lab_] = v
    return -(p['p\\[ion\\]'] + p['p\\[elc\\]'])


def cart(poscar):
    return np.array(poscar['frac']) @ np.array(poscar['cell_A'])


def nearest_pair_domain(frac, cell, mz):
    pos = np.array(frac)[:4] @ cell
    shifts = [np.array([a, b, c]) @ cell for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]
    best = None
    for i in range(4):
        for j in range(i + 1, 4):
            for s in shifts:
                d = pos[j] + s - pos[i]
                if best is None or np.linalg.norm(d) < np.linalg.norm(best[2]) - 1e-9:
                    best = (i, j, d)
    i, j, d = best
    assert abs(d[0]) < 0.1 and abs(d[1]) < 0.1 and abs(d[2]) > 2.5   # face-sharing pair stacked along the trigonal axis
    si, sj = np.sign(mz[i] * d[2]), np.sign(mz[j] * -d[2])
    assert si == sj
    return 'out_pointing' if si < 0 else 'in_pointing', float(np.linalg.norm(d))


def elk_parse(res):
    blk = res['INFO.OUT_final_moments_block']
    tot = [float(x) for x in re.search(r'total moment\s+:\s+(\S+)\s+(\S+)\s+(\S+)', blk).groups()]
    atoms = re.findall(r'atom\s+\d+\s+:\s+(\S+)\s+(\S+)\s+(\S+)', blk)
    cr_mz = [float(a[2]) for a in atoms[:4]]
    rows = []
    for l in res['elk.in'][res['elk.in'].index('atoms'):].split('\n'):
        try:
            v = [float(x) for x in l.split()[:6]]
        except ValueError:
            continue
        if len(v) == 6:
            rows.append(v[:3])
    return np.array(tot), cr_mz, np.array(rows[:10])


def main():
    # ---- force-constant eigenmodes; the hdf5 values carry no units: check them against the all-H phonons
    fc = buy({'kind': 'force_constant_eigenmodes'})
    vals, vecs = fc['Cr_ph_eigenval_pol.hdf5'], fc['Cr_ph_eigenvec_pol.hdf5']
    ph = buy({'kind': 'phonon_gamma', 'masses': 'H'})
    mH = ph['masses_amu'][0]
    k_from_H = [mH * AMU * (2 * np.pi * f * 1e12) ** 2 / EV * 1e-20 for f in ph['frequencies_THz'] if f > 1]
    for m in MODES:
        assert min(abs(k / vals[m] - 1) for k in k_from_H) < 1e-4, m        # values are k in eV/A^2
    E = {m: np.array(vecs[m]) for m in MODES}
    K = {m: vals[m] for m in MODES}
    # ---- Born charges Z[kappa, a (polarization), b (displacement)]
    ref = buy({'kind': 'born_reference'})
    pos0 = np.array(ref['OUTCAR_ion_positions_cartesian_A'])
    cell = np.array(ref['POSCAR']['cell_A'])
    omega = abs(np.linalg.det(cell))
    d0 = dipole(ref)

    def born_rows(amplitude, atom):
        res = buy({'kind': 'born_displacement', 'amplitude': amplitude, 'atom': atom})
        k = ATOMS.index(atom)
        out = np.zeros((3, 3))
        for b, c in enumerate('xyz'):
            run = res['runs'][c]
            u_out = np.array(run['OUTCAR_ion_positions_cartesian_A']) - pos0
            u_pos = cart(run['POSCAR']) - cart(ref['POSCAR'])
            assert np.abs(u_out - u_pos).max() < 1e-4, 'POSCAR and OUTCAR positions disagree'
            out[:, b] = (dipole(run) - d0) / u_out[k, b]
        return out

    Z = np.array([born_rows('0_02', a) for a in ATOMS])
    # The reference run carries a dipole offset (z: -0.0074 e*A) absent from the displaced runs: one-sided
    # differences inherit offset/u in every atom's z row. The acoustic sum rule removes this common error;
    # the -0.02 A runs confirm linearity of the in-plane rows, where the offset is negligible.
    Zsum = Z.sum(axis=0)
    assert abs(Zsum[2]).max() > 1.0 and abs(Zsum[:2]).max() < 0.3
    for a in ('Cr1', 'O1'):
        assert np.abs(born_rows('-0_02', a)[:2] - Z[ATOMS.index(a)][:2]).max() < 0.06
    Z = Z - Z.mean(axis=0)                                                    # acoustic sum rule
    p = {m: np.einsum('kab,kb->a', Z, E[m]) for m in MODES}
    # ---- induced moments
    elk, vasp, doms = {}, {}, {}
    for m in MODES:
        pts = []
        for t in ('0_01', '-0_01', '0_02', '-0_02'):
            res = buy({'kind': 'elk_moment', 'mode': m, 'amplitude': t})
            tot, crz, frac = elk_parse(res)
            amp = float(t.replace('_', '.'))
            pts.append((amp, tot[0]))
            doms.setdefault('Elk', set()).add(nearest_pair_domain(frac, cell, crz)[0])
        a_, y_ = np.array(pts).T
        slope = np.polyfit(a_, y_, 1)[0]
        assert np.abs(y_ - slope * a_).max() < 0.05 * abs(slope) * 0.02 + 2e-5   # linear, symmetric
        elk[m] = slope
        pts = []
        for t in ('0_025', '-0_025'):
            res = buy({'kind': 'vasp_moment', 'mode': m, 'amplitude': t})
            line = res['OUTCAR_final_total_magnetization_line']
            mx = float(line.split('magnetization')[1].split()[0])
            pts.append((float(t.replace('_', '.')), mx))
            crz = [row[-1] for row in res['OUTCAR_last_ion_table_magnetization_z'][:4]]
            doms.setdefault('VASP', set()).add(nearest_pair_domain(res['POSCAR']['frac'], cell, crz)[0])
        vasp[m] = (pts[0][1] - pts[1][1]) / (pts[0][0] - pts[1][0])
    assert all(len(v) == 1 for v in doms.values())
    dom = {k: v.pop() for k, v in doms.items()}
    sign = {k: (1.0 if d == 'out_pointing' else -1.0) for k, d in dom.items()}
    d_out = {'Elk': {m: sign['Elk'] * elk[m] for m in MODES}, 'VASP': {m: sign['VASP'] * vasp[m] for m in MODES}}
    pref = MU0 * MUB * 1e-10 / (omega * 1e-30) * 1e12        # ps/m per (muB/A)(e)/(eV/A^2)
    res_alpha = {}
    for code in ('Elk', 'VASP'):
        contrib = {m: pref * p[m][0] * d_out[code][m] / K[m] for m in MODES}
        res_alpha[code] = (sum(contrib.values()), contrib)
    a_elk, c_elk = res_alpha['Elk']
    a_vasp = res_alpha['VASP'][0]
    assert abs(a_vasp / a_elk - 1) < 0.05, (a_elk, a_vasp)
    frac_c = {m: c_elk[m] / a_elk for m in MODES}
    # annealing: alpha_par(0K) has the sign of alpha_perp(0K); alpha_par changes sign below ~100 K, so near T_N
    # (where the domain is set) alpha_par has the opposite sign; F = -alpha E.H favours alpha_par > 0 for E || H.
    alpha_par_hot_out = -np.sign(a_elk)
    annealed = 'out_pointing' if alpha_par_hot_out > 0 else 'in_pointing'
    submission = {'alpha_perp_ps_per_m': {'out_pointing': float(a_elk), 'in_pointing': float(-a_elk)},
                  'mode_fractions': {m: float(v) for m, v in frac_c.items()},
                  'dominant_mode': max(frac_c, key=frac_c.get),
                  'dM_du_out_pointing_muB_per_A': {m: float(d_out['Elk'][m]) for m in ('Eu_3_x', 'Eu_4_x')},
                  'born_charge_Cr1': {'xx': float(Z[0, 0, 0]), 'xy': float(Z[0, 0, 1]), 'zz': float(Z[0, 2, 2])},
                  'dataset_domain': dom,
                  'domain_selected_by_parallel_field_annealing': annealed,
                  'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# Cr2O3: lattice-mediated spin magnetoelectric response (reference solution)

alpha_perp = mu0 dM_x/dE_x = (mu0/Omega) sum_m (dm_x/du_m) (p_m,x / k_m), with the eigenvectors e_m of the
force-constant matrix (not the mass-weighted phonons), k_m their eigenvalues (the hdf5 values; checked to equal
m_H (2 pi f)^2 of the all-H phonopy run, i.e. eV/A^2), p_m = sum_kappa Z*_kappa e_m,kappa (Z*_ab = (Omega/e)
dP_a/du_b) and dm/du_m the slope of the net cell spin moment (interstitial included). Omega = {omega:.3f} A^3
(rhombohedral cell).

Born charges: +0.02 A displacements against Disp_0_0, the displacement read from each OUTCAR position block and
checked against its POSCAR (the -0_005 folders do not pass this check: their OUTCARs belong to other runs), linear
at -0.02 A. The Disp_0_0 reference has an electronic dipole offset along z (-0.0074 e*A) that the displaced runs
lack, so one-sided differences violate the acoustic sum rule by ~3.7 e in the zz row; imposing the sum rule (or
central differences) removes it. VASP prints dipoles in electron units (p_ion = -ZVAL*r), so P = -(p_ion +
p_elc). Z*_Cr1 = xx {Z[0,0,0]:.3f}, xy {Z[0,0,1]:+.3f}, zz {Z[0,2,2]:.3f}.

Domains: the face-sharing Cr pair (2.62 A along the trigonal axis) defines in/out-pointing. From the converged
moments: Elk = {dom['Elk']} (its input fields bfcmt point the other way: Elk's field couples to spin), VASP =
{dom['VASP']}. The two codes therefore give opposite raw signs; after mapping to the out-pointing domain
alpha_perp(Elk) = {a_elk:.3f} ps/m and alpha_perp(VASP) = {a_vasp:.3f} ps/m. Mode fractions {', '.join(f'{m}: {v:.3f}' for m, v in frac_c.items())}:
the stiffer Eu_4 and Eu_3 modes carry the response; the softest force-constant mode contributes almost nothing.
The Elk Eu_3 +0.03 A run stopped early with 'Convergence targets achieved' at a moment off the linear trend and
should not be used. Sphere-projected moments (muffin-tin or VASP per-ion sums) miss ~10% of the net moment.

Temperature: at 0 K alpha_par and alpha_perp share the sign; alpha_par changes sign on warming, so near T_N the
out-pointing domain has alpha_par > 0 and annealing in parallel E and H (F = -alpha E.H) selects it:
{annealed}. Not supported by these data: clamped-ion, orbital and finite-temperature contributions, alpha_par.
""")
    lab('submit', submission)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(2)


if __name__ == '__main__':
    main()
