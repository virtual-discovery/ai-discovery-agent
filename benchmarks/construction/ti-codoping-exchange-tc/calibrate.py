"""Reference targets, ablation values, robustness values and the calibration table of ti-codoping-exchange-tc.

    python3 -B calibrate.py        (after the reference.py runs listed in README.md; numpy only)
Writes reference.json, ablations.json and calibration.json in this folder and the "calibration" block of the
task's design.json. make_task.py then stages the private targets and solution/ablation_values.json.
"""
import json
from pathlib import Path
import numpy as np
import model as M
import reference as R

HERE = Path(__file__).resolve().parent
DESIGN = HERE.parents[1]/'tasks/ti-codoping-exchange-tc/design.json'
COMPS = {'V5Cr5': {'V': 0.05, 'Cr': 0.05}, 'Cr5Nb5': {'Cr': 0.05, 'Nb': 0.05}, 'Cr5Mo5': {'Cr': 0.05, 'Mo': 0.05},
         'Cr5Mn5': {'Cr': 0.05, 'Mn': 0.05}, 'V2Cr2': {'V': 0.02, 'Cr': 0.02}, 'Cr5Co5': {'Cr': 0.05, 'Co': 0.05}}
ORDERED = ['V5Cr5', 'Cr5Mo5', 'Cr5Mn5']
CANDIDATES = ['V', 'Nb', 'Mo', 'Mn']
CID = {'V': 'V5Cr5', 'Nb': 'Cr5Nb5', 'Mo': 'Cr5Mo5', 'Mn': 'Cr5Mn5'}
TOL = {'relative': 0.2, 'no_order_max_K': 3.0}


def r(x, n=3):
    return None if x is None else float(round(x, n))


def scaled(tag):
    return R.scaled_crossings(tag)


def binder_x(tag):
    return R.crossings(tag)


def largest_scaled(tag):
    xs = scaled(tag)
    vals = [v['T'] for v in xs.values() if v['T'] is not None]
    return vals[-1] if vals else None


def chi_peak(tag, L=3):
    """Wrong path: T_C from the peak of the connected susceptibility of one supercell (no size comparison).
    chi = N (<m^2> - <|m|>^2) / T per sample, disorder averaged; parabolic refinement around the maximum."""
    d = R.load(tag)[L]
    arr = np.array([s['data'] for s in d['samples']])
    n = np.array([s['n_spins'] for s in d['samples']])[:, None]
    T = np.array(d['temps'])
    chi = (n * (arr[:, :, 2] - arr[:, :, 10] ** 2) / T[None, :]).mean(0)
    i = int(np.argmax(chi))
    if 0 < i < len(T) - 1:
        x, y = T[i - 1:i + 2], chi[i - 1:i + 2]
        a, b, _ = np.polyfit(x, y, 2)
        if a < 0:
            return float(-b / (2 * a))
    return float(T[i])


def single_placement(tag):
    """Wrong path: one random placement per size (the first sample of the L = 3 and L = 4 reference sets)."""
    runs = R.load(tag)
    T = np.array(runs[3]['temps'])
    x = R.crossing(T, R.scaled_m2(runs[4])[0], R.scaled_m2(runs[3])[0])
    return 0.0 if x is None else x


def no_order_evidence(tag):
    """For each size pair: the lowest temperature probed, and whether the scaled <m^2> still falls with size there."""
    runs = R.load(tag)
    sizes = sorted(runs)
    T = np.array(runs[sizes[0]]['temps'])
    out = {}
    for a, b in zip(sizes[:-1], sizes[1:]):
        xa, xb = R.scaled_m2(runs[a]).mean(0), R.scaled_m2(runs[b]).mean(0)
        out[f'{a}/{b}'] = {'falls_with_size_at_all_T': bool(np.all(xb < xa)),
                           'ratio_at_lowest_T': float(xb[-1] / xa[-1]), 'lowest_T': float(T[-1]),
                           'binder_L_large_at_lowest_T': float(R.binder(runs[b])[1][-1])}
    return out


def main():
    J0 = lambda a, b: [0.0 if j is None else j for j in M.pair_J(a, b)]
    targets, details, value = {}, {}, {}
    for cid in ORDERED:
        xs = scaled(cid)
        value[cid] = float(np.mean([xs['4/5']['T'], xs['5/6']['T']]))
        targets[cid] = {'type': 'value', 'K': r(value[cid])}
        details[cid] = {'scaled_crossings': xs, 'binder_crossings': binder_x(cid)}
    ev = no_order_evidence('V2Cr2')
    assert all(v['falls_with_size_at_all_T'] for v in ev.values()), ev
    assert all(v['T'] is None for v in scaled('V2Cr2').values())
    targets['V2Cr2'] = {'type': 'at_most', 'max_K': TOL['no_order_max_K']}
    details['V2Cr2'] = {'no_order': ev}
    targets['Cr5Co5'] = {'type': 'null'}
    # Cr5Nb5 is not graded as a number: its order is weak and strongly sample dependent (Binder ratio ~1.2 down to
    # 1-2 K); independent reference sets give no crossing (Cr5Nb5, L = 3/4/5) or a scaled crossing near 10.6 K
    # (Cr5Nb5_slow). Every estimate lies far below Cr5Mo5, so the best co-dopant is unique.
    nb = {'Cr5Nb5': {'scaled': scaled('Cr5Nb5'), 'binder': binder_x('Cr5Nb5')},
          'Cr5Nb5_slow': {'scaled': scaled('Cr5Nb5_slow'), 'binder': binder_x('Cr5Nb5_slow')}}
    def last_pair(xs):
        k = sorted(xs, key=lambda p: int(p.split('/')[1]))[-1]
        return xs[k]
    nb_last = [last_pair(d['scaled']) for d in nb.values()]
    nb_upper = max([(v['T'] + (v['err'] or 0)) for v in nb_last if v['T'] is not None] or [0.0])
    cand = {'V': value['V5Cr5'], 'Mo': value['Cr5Mo5'], 'Mn': value['Cr5Mn5'], 'Nb': nb_upper}
    ranked = sorted(CANDIDATES, key=lambda X: -cand[X])
    best = ranked[0]
    assert best == 'Mo' and value['Cr5Mo5'] * (1 - TOL['relative']) > max(cand[X] for X in ranked[1:]), cand
    details['Cr5Nb5'] = {k: {kk: {p: r(q['T']) for p, q in vv.items()} for kk, vv in d.items()} for k, d in nb.items()}

    # Wrong paths
    N = [s[2] for s in M.shells()]

    def pair_tc(a, b):
        return sum(j * n for j, n in zip(J0(a, b), N)) / 3 / M.KB
    mf = {cid: M.mean_field_tc(c, J0)[0] for cid, c in COMPS.items()}
    rule = {}
    for cid, c in COMPS.items():
        (a, ca), (b, cb) = list(c.items())
        rule[cid] = (pair_tc(a, a) + pair_tc(b, b) + pair_tc(a, b)) / 3 * ca
    t_cr5 = largest_scaled('Cr5')
    t_coco = largest_scaled('Cr5Co5_zero')
    near = {cid: largest_scaled(f'{cid}_near') for cid in ('V5Cr5', 'Cr5Mo5')}
    single = {cid: single_placement(cid) for cid in ORDERED}
    chi = {cid: chi_peak(f'{cid}_coarse') for cid in ORDERED}
    graded = ORDERED + ['V2Cr2', 'Cr5Co5']
    abl = {
        'mean_field_two_sublattice': {'curie_temperature_K': {k: r(mf[k]) for k in graded},
                                      'best_codopant_for_Cr': max(CANDIDATES, key=lambda X: mf[CID[X]])},
        'authors_averaged_mean_field': {'curie_temperature_K': {k: r(rule[k]) for k in graded},
                                        'best_codopant_for_Cr': max(CANDIDATES, key=lambda X: rule[CID[X]])},
        'single_placement_per_size': {'curie_temperature_K': {k: r(v) for k, v in single.items()}},
        'codopant_by_mean_field': {'best_codopant_for_Cr': max(CANDIDATES, key=lambda X: mf[CID[X]])},
        'linear_in_concentration': {'curie_temperature_K': {'V2Cr2': r(0.4 * value['V5Cr5'])}},
        'zero_filled_coco': {'curie_temperature_K': {'Cr5Co5': r(t_coco) if t_coco is not None else 0.0}},
        'near_shells_only': {'curie_temperature_K': {k: (r(v) if v is not None else 0.0) for k, v in near.items()}},
        'dimer_moment_collapse': {'curie_temperature_K': {k: (r(t_cr5) if t_cr5 is not None else 0.0)
                                                          for k in ('Cr5Mn5', 'Cr5Mo5')}},
    }
    rob = {}
    for name, fn in [('scaled_m2_crossing_L3_L4', lambda t: scaled(t)['3/4']['T']),
                     ('scaled_m2_crossing_L5_L6', lambda t: scaled(t)['5/6']['T']),
                     ('binder_crossing_L3_L4', lambda t: binder_x(t)['3/4']['T']),
                     ('binder_crossing_L5_L6', lambda t: binder_x(t)['5/6']['T'])]:
        rob[name] = {'curie_temperature_K': {cid: r(fn(cid)) for cid in ORDERED}}
    # V2Cr2: Binder curves at the paramagnetic value 5/3 cross by noise near 2 K for L = 4/5
    rob['V2Cr2_binder_noise_crossing'] = {'curie_temperature_K': {'V2Cr2': r(binder_x('V2Cr2')['4/5']['T'])}}

    def passes(cid, v):
        t = targets[cid]
        if t['type'] == 'null':
            return v is None
        if t['type'] == 'at_most':
            return v is not None and v <= t['max_K']
        return v is not None and abs(v / t['K'] - 1) <= TOL['relative']
    for name, entry in rob.items():
        for cid, v in entry['curie_temperature_K'].items():
            assert passes(cid, v), (name, cid, v)
    for name, entry in abl.items():
        fails = [cid for cid, v in entry.get('curie_temperature_K', {}).items() if not passes(cid, v)]
        fails += ['best'] if entry.get('best_codopant_for_Cr', best) != best else []
        assert fails, name
    calib = {}
    for cid in ORDERED:
        calib[f'curie_temperature_K[{cid}]'] = {
            'reference_mean_scaled_crossing_L4L5_L5L6': r(value[cid]),
            'reference_details': {'scaled': {k: r(v['T']) for k, v in details[cid]['scaled_crossings'].items()},
                                  'binder': {k: r(v['T']) for k, v in details[cid]['binder_crossings'].items()},
                                  'jackknife_err_scaled': {k: r(v['err']) for k, v in details[cid]['scaled_crossings'].items()}},
            'robustness': {n: e['curie_temperature_K'][cid] for n, e in rob.items() if cid in e['curie_temperature_K']},
            'ablations': {n: e['curie_temperature_K'][cid] for n, e in abl.items() if cid in e.get('curie_temperature_K', {})},
            'tolerance': f"within {int(TOL['relative'] * 100)}% ({value[cid] * (1 - TOL['relative']):.2f}-{value[cid] * (1 + TOL['relative']):.2f} K)"}
    calib['curie_temperature_K[V2Cr2]'] = {
        'reference': 'no order above 0.8 K (scaled <m^2> falls with size at every T for L = 3/4/5; Binder ratio ~1.6)',
        'evidence': details['V2Cr2']['no_order'],
        'ablations': {n: e['curie_temperature_K']['V2Cr2'] for n, e in abl.items() if 'V2Cr2' in e.get('curie_temperature_K', {})},
        'robustness': {'V2Cr2_binder_noise_crossing': rob['V2Cr2_binder_noise_crossing']['curie_temperature_K']['V2Cr2']},
        'tolerance': f"<= {TOL['no_order_max_K']} K (0 = no order above 2 K; 1 K allowance for noise-level Binder crossings)"}
    calib['curie_temperature_K[Cr5Co5]'] = {'reference': None, 'ablations': {
        'zero_filled_coco': abl['zero_filled_coco']['curie_temperature_K']['Cr5Co5'],
        'mean_field_two_sublattice': r(mf['Cr5Co5'])}, 'tolerance': 'must be null'}
    calib['not_an_ablation_single_cell_susceptibility_peak_L3'] = {k: r(v) for k, v in chi.items()}
    calib['ablation_notes'] = {
        'near_shells_only': 'shells 0-6 only: no order in the probed ranges (scaled <m^2> falls with size at every T, '
                            'Binder ratio ~1.65 down to 6 K for V5Cr5 and 10 K for Cr5Mo5), recorded as 0',
        'zero_filled_coco': 'missing Co-Co set to 0: no order above 2 K (scaled <m^2> falls with size); any number fails the null check',
        'dimer_moment_collapse': 'Mn/Mo dropped from the spin model: Cr 5% alone, reference T_C (L = 4/5 scaled crossing)',
        'single_placement_per_size': 'first random placement of the L = 3 and L = 4 reference sets only; 0 = no crossing'}
    calib['best_codopant_for_Cr'] = {'reference': best, 'monte_carlo_T_C_K': {X: r(cand[X]) for X in CANDIDATES},
                                     'Nb_note': 'Nb value is the largest crossing + jackknife error over two independent '
                                                'reference sets (upper bound); Cr5Nb5 is not graded as a number',
                                     'Cr5Nb5_reference_sets': details['Cr5Nb5'],
                                     'ranking_by_mean_field': sorted(CANDIDATES, key=lambda X: -mf[CID[X]]),
                                     'mean_field_K': {X: r(mf[CID[X]]) for X in CANDIDATES},
                                     'tolerance': 'exact'}
    (HERE/'reference.json').write_text(json.dumps({
        'targets': targets, 'tolerances': TOL, 'best_codopant': best,
        'details': {'estimator': 'mean of the <m^2> L^(1+eta) crossings of L = 4/5 and 5/6 (eta = 0.0375), reference runs of reference.py',
                    'Cr5Nb5': details['Cr5Nb5']},
        'hidden_numbers': [value[c] for c in ORDERED]}, indent=1) + '\n')
    (HERE/'ablations.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1) + '\n')
    (HERE/'calibration.json').write_text(json.dumps(calib, indent=1) + '\n')
    if (HERE/'oracle_run.json').exists():
        orc = json.loads((HERE/'oracle_run.json').read_text())
        for cid in ORDERED + ['V2Cr2', 'Cr5Co5']:
            calib[f'curie_temperature_K[{cid}]']['reference_oracle'] = orc['curie_temperature_K'][cid]
        calib['best_codopant_for_Cr']['reference_oracle'] = orc['best_codopant_for_Cr']
    design = json.loads(DESIGN.read_text())
    design['calibration'] = calib
    DESIGN.write_text(json.dumps(design, indent=1) + '\n')
    print(json.dumps({'targets': targets, 'best': best, 'cand': cand, 'ablations': abl, 'robustness': rob}, indent=1, default=float))


if __name__ == '__main__':
    main()
