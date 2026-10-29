"""Reference solution for ti-codoping-exchange-tc; runs as the Harbor oracle and in check_task."""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
Z = CATALOG['dopants']
KB = 0.08617333262
EVIDENCE = []
ROWS, MOM = {}, {}


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    for attempt in range(5):
        out = subprocess.run(cmd, capture_output=True, text=True)
        try:
            result = json.loads(out.stdout)
            break
        except json.JSONDecodeError:
            time.sleep(3)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def evaluate(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(5)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(r)
    EVIDENCE.append(r['id'])
    return r['result']


def key(a, b):
    return f'{a}-{b}' if (Z[a], a) <= (Z[b], b) else f'{b}-{a}'


def buy(a, b, windows=('0-6', '7-16')):
    k = key(a, b)
    for w in windows:
        res = evaluate({'kind': 'pair', 'species': [a, b], 'shells': w})
        ROWS.setdefault(k, []).extend(res['rows'])
        for r in res['rows']:
            # single-impurity moments m01/m02; the dimer column m2 is impurity 2 only in shell 0
            MOM[k.split('-')[0]], MOM[k.split('-')[1]] = abs(r['m01']), abs(r['m02'])


def jlist(a, b):
    J = [None] * 17
    for r in ROWS[key(a, b)]:
        J[r['i_shell']] = r['J']
    return J


def mean_field(comp):
    names = list(comp)
    N = {r['i_shell']: r['num_neighbors'] for rows in ROWS.values() for r in rows}
    M = np.array([[comp[b] * sum(j * N[k] for k, j in enumerate(jlist(a, b))) for b in names] for a in names])
    return float(np.max(np.linalg.eigvals(M).real) / 3 / KB)


ETA = 0.0375   # 3D Heisenberg: <m^2> ~ L^-(1+eta) at T_C (weak disorder is irrelevant for Heisenberg order)


def mc(comp, size, samples, temps, seed=0):
    names = sorted(comp, key=lambda s: Z[s])
    ex = {key(a, b): [float(x) for x in jlist(a, b)] for a in names for b in names if (Z[a], a) <= (Z[b], b)}
    res = evaluate({'kind': 'monte_carlo', 'spins': comp, 'exchange_meV': ex, 'size': size, 'samples': samples,
                    'temperatures_K': temps, 'seed': seed})
    m2 = np.array(res['order_parameters']['total']['m2_per_sample'])
    m4 = np.array(res['order_parameters']['total']['m4_per_sample'])
    return np.array(res['temperatures_K']), m4.mean(0) / m2.mean(0) ** 2, m2.mean(0)


def crossing(T, a, b, tmin=2.0):
    """Highest temperature (T descending) where a - b turns from <= 0 to > 0 on cooling."""
    d = a - b
    for i in range(len(T) - 1):
        if d[i] <= 0 < d[i + 1] and T[i + 1] >= tmin:
            return float(T[i] + (T[i + 1] - T[i]) * d[i] / (d[i] - d[i + 1]))
    return None


def tc(comp):
    # coarse scan of the smaller cell to locate the transition region
    T, g, _ = mc(comp, 3, 16, [40, 30, 24, 19, 15, 11, 7, 3])
    star = None
    for i in range(len(T) - 1):
        if g[i] >= 1.2 > g[i + 1]:
            star = T[i] + (T[i + 1] - T[i]) * (g[i] - 1.2) / (g[i] - g[i + 1])
    if star is None or star < 3:
        grid = [8.0, 6.0, 5.0, 4.0, 3.0, 2.5, 2.0, 1.5]
    else:
        grid = sorted({round(float(star) * f, 2) for f in (1.35, 1.22, 1.12, 1.04, 0.96, 0.88, 0.78, 0.68)}, reverse=True)
    T, g3, m3 = mc(comp, 3, 32, grid)
    T, g4, m4 = mc(comp, 4, 32, grid)
    # finite-size scaling of the net magnetization: <m^2> L^(1+eta) is size independent at T_C
    x3, x4 = m3 * 3 ** (1 + ETA), m4 * 4 ** (1 + ETA)
    x = crossing(T, x4, x3)
    xb = crossing(T, g3, g4)
    # accept order only if the larger cell actually orders below the crossing (Binder ratio -> 1); a frustrated
    # network whose Binder ratio stays well above 1 at the lowest temperature has no long-range order there
    ordered = x is not None and g4[-1] < 1.12
    return (x if ordered else 0.0), {'g4_lowest_T': float(g4[-1]),'grid': [float(t) for t in T], 'binder_L3': g3.round(4).tolist(),
                                       'binder_L4': g4.round(4).tolist(), 'scaled_m2_L3': x3.round(4).tolist(),
                                       'scaled_m2_L4': x4.round(4).tolist(), 'binder_crossing': xb}


def main():
    pairs = [('Cr', 'Cr'), ('V', 'V'), ('V', 'Cr'), ('Nb', 'Nb'), ('Cr', 'Nb'), ('Mo', 'Mo'), ('Cr', 'Mo'),
             ('Mn', 'Mn'), ('Cr', 'Mn')]
    for a, b in pairs:
        buy(a, b)
    buy('Co', 'Co', windows=('0-6',))      # the archive has no short-range Co-Co dimers
    n_coco = len(ROWS.get(key('Co', 'Co'), []))
    assert n_coco == 0
    for a, b in pairs:
        assert all(j is not None for j in jlist(a, b)), (a, b)
    assert all(MOM[s] > 1e-3 for s in ('V', 'Cr', 'Nb', 'Mo', 'Mn'))
    comps = dict(CATALOG['compositions'])
    comps['Cr5Nb5'] = {'Cr': 0.05, 'Nb': 0.05}
    cid_of = {X: next(k for k, c in comps.items() if set(c) == {'Cr', X} and c['Cr'] == c[X] == 0.05)
              for X in CATALOG['codopant_candidates_with_Cr']}
    result, notes, mf = {}, {}, {}
    for cid in ('V5Cr5', 'Cr5Mo5', 'Cr5Mn5', 'Cr5Nb5', 'V2Cr2'):
        mf[cid] = mean_field(comps[cid])
        result[cid], notes[cid] = tc(comps[cid])
    result['Cr5Co5'] = None
    best = max(CATALOG['codopant_candidates_with_Cr'], key=lambda X: result[cid_of[X]])
    submission = {'curie_temperature_K': {k: result[k] for k in CATALOG['compositions']}, 'best_codopant_for_Cr': best,
                  'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    lines = '\n'.join(f'- {cid}: T_C = {result[cid]:.2f} K (crossing of <m^2> L^(1+eta), L = 3/4; 0 = none above 2 K; '
                      f'Binder-ratio crossing {notes[cid]["binder_crossing"]}); two-sublattice mean field {mf[cid]:.1f} K; '
                      f'grid {notes[cid]["grid"]} K, scaled <m^2> L3 {notes[cid]["scaled_m2_L3"]}, L4 {notes[cid]["scaled_m2_L4"]}'
                      for cid in notes)
    (APP/'report.md').write_text(f"""# Co-doped Bi2Te3: Curie temperatures from the archived exchange table (reference solution)

Model: classical unit spins on random Bi sites of both Bi layers, H = -(1/2) sum_ij J_ij e_i.e_j, with the
archived J of each of the 17 shells (4.38-12.17 A, both windows bought) for every species pair and nothing beyond.
All five dopants used carry moments (single-impurity m01/m02 = 1.2-4.4 muB). The archived dimer moment m2 equals
impurity 2's moment only in shell 0; beyond it reads ~0.01 muB for every pair (a host-site readout), which the
same species' m1 in its own rows contradicts, so it is not evidence of moment collapse.

T_C from finite-size scaling of the disorder-averaged net magnetization: a coarse L = 3 scan locates the
transition (Binder ratio g = [<m^4>]/[<m^2>]^2 near 1.2), then L = 3 and L = 4 with 32 random placements each.
T_C is the crossing of [<m^2>] L^(1+eta) (eta = 0.0375, 3D Heisenberg; dilution does not change the universality
class). Binder-ratio crossings agree but are noisier because the curves meet at a shallow angle:
{lines}

Mean field overestimates T_C three- to eightfold and fails qualitatively: at 5% + 5% the strong pairs within
6.5 A are rare, order percolates through weak long-range couplings, and antiferromagnetic shells frustrate it.
Nb has the largest mean-field T_C with Cr, but the antiferromagnetic Nb-Nb shells (6.10 A and 7.51 A) frustrate
the network: Cr5Nb5 shows at most weak, sample-dependent order (Binder ratio stays near 1.2), far below Cr5Mo5,
whose only strong antiferromagnetic shell (Mo-Mo, 6.49 A) has three neighbours. The best co-dopant for Cr is {best}. At 2% + 2% V/Cr nothing orders above 2 K although mean
field gives ~33 K. Co-Co dimers within 7.8 A are absent from the archive and the 9.99 A Co-Co run lost its moment,
so T_C of Cr5Co5 is not determined (null); filling the gap with zeros would be a guess.
Not supported: anisotropy, DMI and quantum effects (classical isotropic model), Curie temperatures of real films
(clustering, interstitials, carriers), and anything that needs couplings beyond 12.2 A.
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
