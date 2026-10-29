"""Reference solution for mg2x-midir-forbidden-transition; runs as the Harbor oracle and in check_task."""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
KB = 8.617333262e-5
HARTREE = 27.211386245988
BOHR = 0.529177210903
RY_MEV = 13605.693122994
COMPOUNDS = ['Mg2Si', 'Mg2Ge', 'Mg2Sn']
LEVELS = ['LDA', 'QSGW']


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


EVIDENCE = []


def buy(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(1)
        r = lab(f"job/{r['id']}")
    if r['status'] != 'complete':
        raise RuntimeError(r)
    EVIDENCE.append(r['id'])
    return r['result']


def mass(q, e, a, qfit=0.020):
    """Parabolic mass (m_e) from E(q) near X, q in units of 2 pi / a."""
    k2pi = 2 * np.pi / (a / BOHR)
    q, e = np.asarray(q), np.asarray(e)
    keep = q <= qfit + 1e-9
    return 1 / (2 * np.polyfit((q[keep] * k2pi) ** 2, (e[keep] - e[0]) / HARTREE, 1)[0])


def tote_rows(result):
    rows = []
    for line in result['rows_at_Gamma_and_X']:
        f = [float(x.replace('D', 'E')) for x in line.split()]
        rows.append({'q': f[:3], 'band': int(f[3]), 'e_in': f[5], 'c7': f[6], 'c8': f[7], 'Z': f[8]})
    return rows


def at_x(rows):
    return {r['band']: r for r in rows if np.allclose(r['q'], (0, 1, 0), atol=1e-6)}


def box_ratio(box, a, n_cm3, band):
    """<q_perp^2>(300 K) / <q_perp^2>(77 K) at density n, Fermi-Dirac, over the octant box (trapezoid weights)."""
    K = np.array(box['points'])
    e = np.array(box['E_eV'][f'band_{band}'])
    e = e - e[np.argmin(np.abs(K).sum(1))]
    h, qm = box['step_2pi_over_a'], box['q_max_2pi_over_a']
    upper = (qm['q_par'], qm['q_y'], qm['q_z'])
    w = np.ones(len(K)) * h ** 3
    for j in range(3):                       # half weight on the symmetry planes and on the outer faces
        w *= np.where(np.abs(K[:, j]) < 1e-9, 0.5, 1.0) * np.where(np.abs(K[:, j] - upper[j]) < 1e-9, 0.5, 1.0)
    G = 2 * np.pi / (a * 1e-10)
    norm = 2 * 3 * 8 * G ** 3 / (2 * np.pi) ** 3   # spin x valleys x octants, per m^3

    def mean_q2(T):
        kT = KB * T
        occ = lambda mu: w / (1 + np.exp(np.clip((e - mu) / kT, -700, 700)))
        lo, hi = -2.0, 1.0
        for _ in range(200):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if norm * occ(mid).sum() < n_cm3 * 1e6 else (lo, mid)
        f = occ((lo + hi) / 2)
        return (f * (K[:, 1] ** 2 + K[:, 2] ** 2)).sum() / f.sum()
    return float(mean_q2(300.0) / mean_q2(77.0))


def main():
    a = CATALOG['lattice_constants_A']
    nv = CATALOG['valence_bands_ecalj']
    # 1. X-point levels and the identity of the two states from their longitudinal masses (Mg2Si anchor: X1 lower).
    split, ident, notes = {c: {} for c in COMPOUNDS}, {}, []
    ml = {}
    for c in COMPOUNDS:
        for lev in LEVELS:
            s = buy({'kind': 'x_scan', 'compound': c, 'level': lev, 'direction': 'X-Gamma'})
            lo, hi = nv[c] + 1, nv[c] + 2
            q = s['q_from_X_2pi_over_a']
            m = {b: mass(q, s['E_eV'][f'band_{b}'], a[c]) for b in (lo, hi)}
            ex = {b: s['E_eV'][f'band_{b}'][0] for b in (lo, hi)}
            ml[(c, lev)] = m
            x1 = min(m, key=m.get)            # X1: the lighter longitudinal mass (the lower state of Mg2Si)
            x3 = hi if x1 == lo else lo
            ident[(c, lev)] = (x1, x3)
            split[c][lev] = ex[x3] - ex[x1]
            notes.append(f'{c} {lev}: lower band {lo} m_l = {m[lo]:.3f}, upper band {hi} m_l = {m[hi]:.3f} -> X1 = band {x1}; '
                         f'X3 - X1 = {1e3 * split[c][lev]:+.1f} meV')
    assert ident[('Mg2Si', 'LDA')][0] == nv['Mg2Si'] + 1 and ident[('Mg2Si', 'QSGW')][0] == nv['Mg2Si'] + 1
    # 2. One-shot G0W0 (TOTE): the column that equals e_in + Z (c8 - e_in) is the linearized solution.
    for c in COMPOUNDS:
        rows = tote_rows(buy({'kind': 'tote', 'run': f'g0w0_lda/{c}'}))
        for r in rows:
            assert abs((r['c7'] - r['e_in']) - r['Z'] * (r['c8'] - r['e_in'])) < 1e-6
        X = at_x(rows)
        x1, x3 = ident[(c, 'LDA')]              # TOTE bands are numbered as the LDA states it started from
        assert abs((X[x3]['e_in'] - X[x1]['e_in']) - split[c]['LDA']) < 2e-3
        split[c]['G0W0'] = X[x3]['c7'] - X[x1]['c7']
    # 3. Ladder correction: QSGW and QSGW + ladder in the same code and basis.
    qq = buy({'kind': 'questaal', 'run': 'prod/qsgw'})['eigenvalues_Ry']['X']
    ql = buy({'kind': 'questaal', 'run': 'prod/qsgw_ladder'})['eigenvalues_Ry']['X']
    vertex = ((ql[5] - ql[4]) - (qq[5] - qq[4])) * RY_MEV
    # 4. Selection rule from the PBE matrix elements of Mg2Si.
    pts = [[0, 0, 0], [0.02, 0, 0], [0.04, 0, 0], [0, 0.02, 0], [0, 0.04, 0], [0, 0, 0.02], [0, 0.02, 0.02],
           [0, 0.04, 0.04], [0.04, 0.02, 0], [0.08, 0.02, 0]]
    pm = {tuple(r[0]): r[3] for r in buy({'kind': 'momentum_matrix', 'points': pts})['rows']}
    assert pm[(0, 0, 0)] < 1e-12 and pm[(0.04, 0, 0)] < 1e-12            # forbidden at X and along the valley axis
    assert abs(pm[(0, 0.04, 0)] / pm[(0, 0.02, 0)] - 4) < 0.3              # ~ q_perp^2
    assert abs(pm[(0.08, 0.02, 0)] / pm[(0, 0.02, 0)] - 1) < 0.2           # no q_par dependence at leading order
    beta2 = pm[(0, 0.02, 0)] / (0.02 * 2 * np.pi / (a['Mg2Si'] / BOHR)) ** 2
    # 5. Infrared/Raman ratio on the QSGW boxes (lower conduction state; Mg2Sn: band nv+1 is X3).
    si_box = buy({'kind': 'qsgw_box', 'compound': 'Mg2Si', 'box': 'coarse'})
    sn_box = buy({'kind': 'qsgw_box', 'compound': 'Mg2Sn', 'box': 'coarse'})
    ratio = {'Mg2Si_1e16': box_ratio(si_box, a['Mg2Si'], 1e16, nv['Mg2Si'] + 1),
             'Mg2Si_3e18': box_ratio(si_box, a['Mg2Si'], 3e18, nv['Mg2Si'] + 1),
             'Mg2Sn_1e16': box_ratio(sn_box, a['Mg2Sn'], 1e16, nv['Mg2Sn'] + 1)}
    # 6. Electron-phonon (AHC, 24^3, non-adiabatic, eta = 10 meV): change of X3 - X1.
    def eph(c, T):
        req = {'kind': 'electron_phonon', 'compound': c, 'method': 'AHC', 'q_grid': 24, 'temperature_K': T}
        if c == 'Mg2Ge':
            req['ge_3d'] = 'core'
        s = buy(req)['renormalization_meV']['eta0.010']
        return s['X3'][0] - s['X1'][0]
    zpr_si, th_si, th_ge = eph('Mg2Si', 0), eph('Mg2Si', 300), eph('Mg2Ge', 300)
    # 7. QSGW at the 100 K lattice: first energy column of the last iteration (checked on the production run).
    prod = at_x(tote_rows(buy({'kind': 'tote', 'run': 'qsgw_final_iteration/Mg2Si'})))
    lo, hi = nv['Mg2Si'] + 1, nv['Mg2Si'] + 2
    assert abs((prod[hi]['e_in'] - prod[lo]['e_in']) - split['Mg2Si']['QSGW']) < 1e-3
    low = at_x(tote_rows(buy({'kind': 'tote', 'run': 'qsgw_low_temperature_lattice/Mg2Si'})))
    low_split = low[hi]['e_in'] - low[lo]['e_in']
    pred = {'Mg2Si_T0': low_split + 1e-3 * zpr_si, 'Mg2Si_300K': split['Mg2Si']['QSGW'] + 1e-3 * th_si,
            'Mg2Ge_300K': split['Mg2Ge']['QSGW'] + 1e-3 * th_ge}
    submission = {'x3_minus_x1_eV': {c: {lev: float(split[c][lev]) for lev in ('LDA', 'G0W0', 'QSGW')} for c in COMPOUNDS},
                  'vertex_correction_Mg2Si_meV': float(vertex),
                  'ir_raman_ratio_300K_over_77K': ratio,
                  'predicted_splitting_eV': {k: float(v) for k, v in pred.items()},
                  'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    s = split
    (APP/'report.md').write_text(f"""# Mid-infrared band of n-type Mg2X (reference solution)

**Assignment.** The band is the vertical transition between the two lowest conduction states at X, in which the
electrons sit (X1 -> X3 in Mg2Si and Mg2Ge, X3 -> X1 in Mg2Sn). Both states are even under inversion, so the
transition is dipole-forbidden at X: the PBE matrix element vanishes at X and along the valley axis, and grows as
q_perp^2 (|p|^2 = beta^2 q_perp^2 with beta^2 = {beta2:.1f} a.u. for Mg2Si). The measured bands are peaks near
0.40 eV (Mg2Si) and 0.56 eV (Mg2Ge), below or near the interband edges, not Drude tails; their position matches
the QSGW/G0W0 splitting and not the LDA one.

**Labels.** X1 and X3 were told apart by their dispersion: X1 has the lighter longitudinal mass (about 0.54) and
X3 the heavier (about 1.03). {'; '.join(notes)}. In Mg2Sn the order is inverted at every level, and the LDA
places X3 below the valence-band maximum (semimetal).

**Levels.** X3 - X1 (eV): Mg2Si LDA {s['Mg2Si']['LDA']:.3f}, G0W0 {s['Mg2Si']['G0W0']:.3f}, QSGW {s['Mg2Si']['QSGW']:.3f};
Mg2Ge {s['Mg2Ge']['LDA']:.3f}, {s['Mg2Ge']['G0W0']:.3f}, {s['Mg2Ge']['QSGW']:.3f}; Mg2Sn {s['Mg2Sn']['LDA']:.3f},
{s['Mg2Sn']['G0W0']:.3f}, {s['Mg2Sn']['QSGW']:.3f}. A scissor shifts both conduction states and cannot repair the LDA
splitting. The G0W0 value is the TOTE column that equals e_in + Z (eQP(Z=1) - e_in), whatever the README label
says; in the QSGW TOTE files the first energy column is the QSGW level. The ladder correction to the Mg2Si
splitting is {vertex:.1f} meV (Questaal QSGW vs QSGW+ladder, same basis), while it lowers the gap by about 90 meV;
QSGW80 instead lowers the splitting by about 60 meV, so the 0.8 scaling is not a correction for this splitting.
Comparing Questaal and ecalj directly would mix in a basis difference of about 28 meV.

**Temperature dependence of the strength.** At fixed n the IR/Raman ratio is <q_perp^2>. Over the QSGW boxes,
with Fermi-Dirac occupation (2 spins x 3 valleys x 8 octants), Q = {ratio['Mg2Si_1e16']:.2f} (Mg2Si, 1e16),
{ratio['Mg2Si_3e18']:.2f} (Mg2Si, 3e18) and {ratio['Mg2Sn_1e16']:.2f} (Mg2Sn, 1e16), against 300/77 = 3.90 for a
parabolic valley: the lower state flattens along the [011] diagonal of the transverse plane where it couples to
the upper one, which an X-W line or parabolic masses miss; degeneracy at 3e18 lowers the ratio. Beyond leading
order the matrix element saturates (mixing), which these numbers do not include.

**Predicted splittings.** QSGW plus the AHC renormalization of X3 - X1 ({zpr_si:.1f} meV at 0 K, {th_si:.1f} meV at
300 K for Mg2Si; {th_ge:.1f} meV at 300 K for Mg2Ge), with the QSGW splitting at the 100 K lattice
({low_split:.3f} eV) for T -> 0: Mg2Si {pred['Mg2Si_T0']:.3f} eV (T -> 0) and {pred['Mg2Si_300K']:.3f} eV (300 K),
Mg2Ge {pred['Mg2Ge_300K']:.3f} eV (300 K). Thermal expansion and electron-phonon coupling nearly cancel between
0 and 300 K.

**Not supported.** Electron-phonon and lattice effects for Mg2Sn; absolute strengths at the QSGW level; a
convergence of the ladder window; the measured temperature dependence of the band area, which mixes the
carrier density (freeze-out) and donor-bound electrons; spin-orbit effects.
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
