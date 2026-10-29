"""Generate the staged, lab, private, ablation and robustness data for fluorite-superionic-finite-size.

    python3 -B fetch.py && python3 -B extract.py && python3 -B make_task.py

Reads the pool written by extract.py (~/.cache/ai-theorist/fluorite-superionic-finite-size/pool.npz).
Reference values use the full-resolution released files (analysis.py protocol). Ablation and robustness
values are computed the way an agent would, from the payloads the lab serves (served(), agent_*()).
"""
import json
import functools
import math
import os
from pathlib import Path
import numpy as np
import analysis as A

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/fluorite-superionic-finite-size'
CACHE = Path(os.path.expanduser('~/.cache/ai-theorist/fluorite-superionic-finite-size'))
META = {r['key']: r for r in json.loads((CACHE/'pool_meta.json').read_text())}
POOL = np.load(CACHE/'pool.npz')
FLUORITES = ('CaF2', 'PbF2', 'UO2')
BUDGET = 140
PREDICT = [('CaF2', 4, 1250), ('PbF2', 7, 600), ('PbF2', 7, 800), ('UO2', 4, 2450), ('UO2', 4, 3150),
           ('AgI', 6, 600), ('AgI', 6, 900)]
HIDDEN_KEYS = {f'{m}-{L}-{T}' for m, L, T in PREDICT}
THRESHOLDS = {'Tc_tol_K': 120.0, 'Ea_tol_eV': {'CaF2': {'ordered': 0.25, 'superionic': 0.25},
                                                'PbF2': {'ordered': 0.15, 'superionic': 0.10}},
              'logD_tol_dex': 0.15}
SPECIES = {'CaF2': ('Ca', 'F'), 'PbF2': ('Pb', 'F'), 'UO2': ('U', 'O'), 'AgI': ('I', 'Ag')}


def hidden(r):
    return (r['material'] == 'PbF2' and r['L'] >= 6) or (r['material'] == 'AgI' and r['L'] == 6) \
        or r['key'] in HIDDEN_KEYS


def cost(natoms):
    return max(1, math.ceil(natoms / 300))


def pid(m, L, T):
    return f'{m}-{L}x{L}x{L}-{T}K'


def msd_rows(n):
    return list(range(0, 50, 2)) + list(range(50, 500, 10)) + list(range(500, n, 50))


@functools.lru_cache(maxsize=None)
def served(key):
    """The lab's response for one run (what the agent sees)."""
    r = META[key]
    msd, log = POOL[key + '/msd'], POOL[key + '/log']
    rows = [[k] + [float(msd[k, c]) for c in (0, 2, 3, 4, 6)] for k in msd_rows(len(msd))]
    blocks = np.array_split(log, 20)
    sem = lambda c: float(np.std([b[:, c].mean() for b in blocks], ddof=1) / math.sqrt(len(blocks)))
    t1, t2 = SPECIES[r['material']]
    return {
        'material': r['material'], 'supercell': f"{r['L']}x{r['L']}x{r['L']}", 'atoms': r['natoms'],
        'temperature_K': r['T'], 'input_deck': r['deck'],
        'msd_table': {
            'origin': 'rows of the released msd_CM_B4.dat (ANALISI, 4 blocks), subset of lags',
            'columns': ['lag_frames', f'msd_type1_{t1}_A2', f'msd_type2_{t2}_A2', f'var_msd_type2_{t2}_A4',
                        f'msd_cm_type1_{t1}_A2', f'msd_cm_type2_{t2}_A2'],
            'notes': 'lag in saved trajectory frames; MSD of each atom type and of the centre of mass of each '
                     'type, in the frame of zero total momentum; var = ANALISI block-averaging variance estimate of msd_type2',
            'rows': rows},
        'thermo_summary': {
            'origin': 'released LAMMPS log of the production run (all samples)',
            'production_ps': round(float(log[-1, 0] - log[0, 0]), 3), 'samples': int(len(log)),
            'mean_temperature_K': round(float(log[:, 3].mean()), 3),
            'mean_total_energy_eV': round(float(log[:, 1].mean()), 5), 'sem_total_energy_eV': round(sem(1), 6),
            'var_total_energy_eV2': float('%.6g' % log[:, 1].var()),
            'mean_potential_energy_eV': round(float(log[:, 2].mean()), 5), 'sem_potential_energy_eV': round(sem(2), 6),
            'mean_pressure_bar': round(float(log[:, 4].mean()), 2), 'block_count_for_sem': 20},
        'source': 'retrieved author-released MD result (no new calculation)'}


# ---------------------------------------------------------------- reference (full-resolution data)
def ref_D(key, window=A.D_WINDOW):
    r = META[key]
    msd = POOL[key + '/msd']
    t = np.arange(len(msd)) * r['frame_ps']
    return A.diffusion(t, A.lab_msd(r['material'], msd), window)


def ref_plateau(key):
    r = META[key]
    msd = POOL[key + '/msd']
    t = np.arange(len(msd)) * r['frame_ps']
    return float(A.lab_msd(r['material'], msd)[(t >= 1) & (t <= 3)].mean())


def runs(m, L):
    return sorted([r for r in META.values() if r['material'] == m and r['L'] == L], key=lambda r: r['T'])


def ref_tc(m, L=4):
    rr = runs(m, L)
    Tm = np.array([POOL[r['key'] + '/log'][:, 3].mean() for r in rr])
    E = np.array([POOL[r['key'] + '/log'][:, 1].mean() / (r['natoms'] / A.FU_ATOMS[m]) for r in rr])
    return A.peak_temperature(Tm, A.heat_capacity(Tm, E))


def diffusive(D, plateau):
    """MSD gain over 100 ps at least twice the vibrational plateau (criterion used for the ordered window)."""
    return D > 0 and 6 * D * 100.0 >= 2 * plateau


def ref_ea(m, tc):
    rr = runs(m, 4)
    D = {r['T']: ref_D(r['key']) for r in rr}
    pl = {r['T']: ref_plateau(r['key']) for r in rr}
    lo = [T for T in D if T <= tc - 100 and diffusive(D[T], pl[T])]
    hi = [T for T in D if T >= tc + 100]
    return {'ordered': A.arrhenius(lo, [D[T] for T in lo])[0], 'superionic': A.arrhenius(hi, [D[T] for T in hi])[0],
            'ordered_T': lo, 'superionic_T': hi}


# ---------------------------------------------------------------- agent-side analysis (served payloads)
def frame_ps(deck):
    import re
    dt = float(re.search(r'variable\s+dt\s+equal\s+([\d.]+)', deck).group(1))
    every = int(re.search(r'variable\s+traj_print_interval\s+equal\s+(\d+)', deck).group(1))
    return dt * every


def agent_D(p, window=(20.0, 100.0), frame=None, frame_of_ref='lab'):
    rows = np.array(p['msd_table']['rows'])
    t = rows[:, 0] * (frame if frame else frame_ps(p['input_deck']))
    m1, m2 = A.MASS[p['material']]
    a = A.MOBILE_RATIO[p['material']] * m2 / m1
    msd = rows[:, 2] + a * (2 + a) * rows[:, 5] if frame_of_ref == 'lab' else rows[:, 2]
    sel = (t >= window[0]) & (t <= window[1])
    return float(np.polyfit(t[sel], msd[sel], 1)[0] / 6)


def agent_plateau(p):
    rows = np.array(p['msd_table']['rows'])
    t = rows[:, 0] * frame_ps(p['input_deck'])
    return float(rows[(t >= 1) & (t <= 3), 2].mean())


def agent_tc(ps, fluct=False):
    ps = sorted(ps, key=lambda p: p['temperature_K'])
    Tm = np.array([p['thermo_summary']['mean_temperature_K'] for p in ps])
    nfu = ps[0]['atoms'] / A.FU_ATOMS[ps[0]['material']]
    if fluct:
        cv = np.array([p['thermo_summary']['var_total_energy_eV2'] / (A.KB * T) ** 2 / nfu for p, T in zip(ps, Tm)])
    else:
        cv = A.heat_capacity(Tm, [p['thermo_summary']['mean_total_energy_eV'] / nfu for p in ps])
    return A.peak_temperature(Tm, cv)


def agent_ea(ps, tc, **kw):
    D = {p['temperature_K']: agent_D(p, **kw) for p in ps}
    pl = {p['temperature_K']: agent_plateau(p) for p in ps}
    lo = [T for T in D if T <= tc - 100 and diffusive(D[T], pl[T])]
    hi = [T for T in D if T >= tc + 100]
    out = {'ordered': A.arrhenius(lo, [D[T] for T in lo])[0]}
    if len(hi) >= 2:
        out['superionic'] = A.arrhenius(hi, [D[T] for T in hi])[0]
    return out


def interp_logD(T, pts):
    """ln D linear in 1/T between the two nearest bracketing temperatures."""
    below = max(t for t in pts if t < T)
    above = min(t for t in pts if t > T)
    x, x0, x1 = 1 / T, 1 / below, 1 / above
    y = np.log(pts[below]) + (np.log(pts[above]) - np.log(pts[below])) * (x - x0) / (x1 - x0)
    return float(np.exp(y))


def P(m, L, T):
    return served(f'{m}-{L}-{T}')


ORACLE_PLAN = {('CaF2', 4): [1100, 1200, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1750],
               ('PbF2', 4): [500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1100, 1250],
               ('PbF2', 5): [600, 800],
               ('UO2', 4): [2400, 2500, 2650, 2700, 2750, 2800, 2850, 2900, 3100, 3200],
               ('AgI', 2): [600, 900], ('AgI', 5): [600, 900],
               ('PbF2', 2): [600], ('PbF2', 3): [600], ('CaF2', 2): [1300], ('UO2', 3): [2750]}


def main():
    available = [r for r in META.values() if not hidden(r)]
    catalog_runs = [{'material': r['material'], 'supercell': r['L'], 'atoms': r['natoms'], 'temperature_K': r['T'],
                     'cost': cost(r['natoms'])} for r in sorted(available, key=lambda r: (r['material'], r['L'], r['T']))]
    predict = [{'id': pid(m, L, T), 'material': m, 'supercell': L, 'atoms': META[f'{m}-{L}-{T}']['natoms'],
                'temperature_K': T} for m, L, T in PREDICT]
    catalog = {
        'description': 'Released NVT molecular-dynamics runs (LAMMPS) of fluorite-structure CaF2, PbF2 and UO2 and of '
                       'alpha-AgI, F. Grasselli (2022). Supercells are n x n x n conventional cells at a fixed lattice '
                       'parameter; one 0.8 ns production run per (material, supercell, temperature).',
        'materials': {'CaF2': 'fluorite, rigid-ion Born-Mayer-Huggins (Fumi-Tosi) potential, 12 atoms per cell',
                      'PbF2': 'fluorite, rigid-ion Born-Mayer-Huggins potential, 12 atoms per cell',
                      'UO2': 'fluorite, partial-charge + EAM (many-body) potential, 12 atoms per cell',
                      'AgI': 'alpha-AgI (bcc iodide sublattice), Coulomb + Morse potential, 4 atoms per cell'},
        'request': {'kind': 'md_run', 'material': '<material>', 'supercell': '<n>', 'temperature_K': '<T>'},
        'budget': BUDGET, 'runs_available': catalog_runs, 'predict': predict}
    lab = {f"{r['material']}-{r['L']}-{r['T']}": served(r['key']) for r in available}

    # --- reference values (full-resolution released data)
    tc_ref = {m: ref_tc(m) for m in FLUORITES}
    ea_ref = {m: ref_ea(m, tc_ref[m]) for m in ('CaF2', 'PbF2')}
    d_ref = {pid(m, L, T): ref_D(f'{m}-{L}-{T}') for m, L, T in PREDICT}
    hidden_json = {
        'transition_temperature_K': tc_ref,
        'activation_energy_eV': {m: {k: ea_ref[m][k] for k in ('ordered', 'superionic')} for m in ea_ref},
        'diffusion_targets_A2_per_ps': d_ref, 'thresholds': THRESHOLDS,
        'reference_notes': {'Ea_windows_K': {m: {'ordered': ea_ref[m]['ordered_T'], 'superionic': ea_ref[m]['superionic_T']}
                                             for m in ea_ref},
                            'Tc_by_size_K': {m: {L: ref_tc(m, L) for L in (2, 3, 4)} for m in FLUORITES}}}

    # --- oracle-equivalent analysis on served data
    buy = {k: [P(k[0], k[1], T) for T in Ts] for k, Ts in ORACLE_PLAN.items()}
    o_tc = {m: agent_tc([p for p in buy[(m, 4)]]) for m in FLUORITES}
    o_ea = {m: agent_ea(buy[(m, 4)], o_tc[m]) for m in ('CaF2', 'PbF2')}
    D4 = {m: {p['temperature_K']: agent_D(p) for p in buy[(m, 4)]} for m in FLUORITES}
    o_pred = {pid('CaF2', 4, 1250): interp_logD(1250, D4['CaF2']),
              pid('UO2', 4, 2450): interp_logD(2450, D4['UO2']), pid('UO2', 4, 3150): interp_logD(3150, D4['UO2']),
              pid('PbF2', 7, 600): agent_D(P('PbF2', 5, 600)), pid('PbF2', 7, 800): agent_D(P('PbF2', 5, 800)),
              pid('AgI', 6, 600): agent_D(P('AgI', 5, 600)), pid('AgI', 6, 900): agent_D(P('AgI', 5, 900))}

    def preds(**changes):
        out = dict(o_pred)
        out.update(changes)
        return out

    # --- ablations (each a realistic wrong path, from served data)
    abl = {}
    abl['uniform_frame_time'] = preds(**{k: agent_D_scaled for k, agent_D_scaled in
                                         [(pid('UO2', 4, 2450), interp_logD(2450, {p['temperature_K']: agent_D(p, frame=0.04) for p in buy[('UO2', 4)]})),
                                          (pid('UO2', 4, 3150), interp_logD(3150, {p['temperature_K']: agent_D(p, frame=0.04) for p in buy[('UO2', 4)]}))]})
    two = {m: [P(m, 2, r['T']) for r in runs(m, 2)] for m in FLUORITES}
    D2 = {m: {p['temperature_K']: agent_D(p) for p in two[m]} for m in FLUORITES}
    abl['small_cell_Tc'] = {m: agent_tc(two[m]) for m in FLUORITES}
    abl['small_cell_Ea'] = {}
    for m in ('CaF2', 'PbF2'):
        lo = [T for T in D4[m] if T <= o_tc[m] - 100 and diffusive(D4[m][T], agent_plateau(P(m, 4, T)))]
        hi = [T for T in D4[m] if T >= o_tc[m] + 100]
        abl['small_cell_Ea'][m] = {'ordered': A.arrhenius(lo, [D2[m][T] for T in lo])[0],
                                   'superionic': A.arrhenius(hi, [D2[m][T] for T in hi])[0]}
    abl['small_cell_predictions'] = preds(**{pid('CaF2', 4, 1250): D2['CaF2'][1250], pid('UO2', 4, 2450): D2['UO2'][2450],
                                             pid('UO2', 4, 3150): D2['UO2'][3150], pid('PbF2', 7, 600): D2['PbF2'][600],
                                             pid('PbF2', 7, 800): D2['PbF2'][800], pid('AgI', 6, 600): agent_D(P('AgI', 2, 600)),
                                             pid('AgI', 6, 900): agent_D(P('AgI', 2, 900))})
    abl['768_cell_assumed_converged'] = preds(**{pid('PbF2', 7, 600): D4['PbF2'][600], pid('PbF2', 7, 800): D4['PbF2'][800]})

    def yeh_hummer(T, L_target=7):
        Ls = [2, 3, 4]
        Ds = [agent_D(P('PbF2', L, T)) for L in Ls]
        slope, intercept = np.polyfit(1 / np.array(Ls, float), Ds, 1)
        return float(intercept + slope / L_target)
    abl['yeh_hummer_1_over_L'] = preds(**{pid('PbF2', 7, 600): yeh_hummer(600), pid('PbF2', 7, 800): yeh_hummer(800)})
    abl['single_arrhenius'] = {}
    for m in ('CaF2', 'PbF2'):
        Ts = [T for T in D4[m] if diffusive(D4[m][T], agent_plateau(P(m, 4, T)))]
        e = A.arrhenius(Ts, [D4[m][T] for T in Ts])[0]
        abl['single_arrhenius'][m] = {'ordered': e, 'superionic': e}
    three = [P('UO2', 3, r['T']) for r in runs('UO2', 3)]
    abl['uo2_324_cell_Tc'] = agent_tc(three)
    abl['uo2_324_cell_predictions'] = preds(**{pid('UO2', 4, 2450): agent_D(P('UO2', 3, 2450)),
                                               pid('UO2', 4, 3150): agent_D(P('UO2', 3, 3150))})
    ratio = float(np.median([D4['CaF2'][T] / D2['CaF2'][T] for T in D4['CaF2'] if T >= o_tc['CaF2'] - 200]))
    abl['agi_fluorite_correction'] = preds(**{pid('AgI', 6, T): agent_D(P('AgI', 2, T)) * ratio for T in (600, 900)})
    abl['agi_fluorite_correction_ratio'] = ratio

    # --- robustness (defensible alternatives, from served data)
    rob = {}
    w = (10.0, 80.0)
    D4w = {m: {p['temperature_K']: agent_D(p, window=w) for p in buy[(m, 4)]} for m in FLUORITES}
    rob['window_10_80'] = {
        'activation_energy_eV': {m: agent_ea(buy[(m, 4)], o_tc[m], window=w) for m in ('CaF2', 'PbF2')},
        'predictions': {pid('CaF2', 4, 1250): interp_logD(1250, D4w['CaF2']), pid('UO2', 4, 2450): interp_logD(2450, D4w['UO2']),
                        pid('UO2', 4, 3150): interp_logD(3150, D4w['UO2']), pid('PbF2', 7, 600): agent_D(P('PbF2', 5, 600), window=w),
                        pid('PbF2', 7, 800): agent_D(P('PbF2', 5, 800), window=w), pid('AgI', 6, 600): agent_D(P('AgI', 5, 600), window=w),
                        pid('AgI', 6, 900): agent_D(P('AgI', 5, 900), window=w)}}
    f_tc = {m: agent_tc(buy[(m, 4)], fluct=True) for m in FLUORITES}
    rob['tc_from_fluctuations'] = {'transition_temperature_K': f_tc,
                                   'activation_energy_eV': {m: agent_ea(buy[(m, 4)], f_tc[m]) for m in ('CaF2', 'PbF2')}}
    d4, d5 = D4['PbF2'], {T: agent_D(P('PbF2', 5, T)) for T in (600, 800)}
    n4, n5, n7 = 768, 1500, 4116
    rob['pbf2_one_over_N'] = {pid('PbF2', 7, T): d5[T] + (d5[T] - d4[T]) * (1 / n7 - 1 / n5) / (1 / n5 - 1 / n4) for T in (600, 800)}
    rob['arrhenius_fit_predictions'] = {}
    for key, m, T, lo, hi in [(pid('CaF2', 4, 1250), 'CaF2', 1250, 1100, 1350), (pid('UO2', 4, 2450), 'UO2', 2450, 2400, 2650),
                              (pid('UO2', 4, 3150), 'UO2', 3150, 2900, 3200)]:
        Ts = [t for t in D4[m] if lo <= t <= hi]
        e, d0 = A.arrhenius(Ts, [D4[m][t] for t in Ts])
        rob['arrhenius_fit_predictions'][key] = d0 * math.exp(-e / (A.KB * T))
    rob['self_frame'] = {pid('AgI', 6, T): agent_D(P('AgI', 2, T), frame_of_ref='raw') for T in (600, 900)}

    # --- write files
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps(lab, separators=(',', ':')) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden_json, indent=1) + '\n')
    values = list(tc_ref.values()) + [ea_ref[m][k] for m in ea_ref for k in ('ordered', 'superionic')] + list(d_ref.values())
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1) + '\n')

    # --- calibration printout
    dex = lambda a, b: math.log10(a / b) if a > 0 else float('-inf')
    print(json.dumps({'Tc_ref': tc_ref, 'Tc_oracle': o_tc, 'Tc_fluct': f_tc, 'Tc_2x2x2': abl['small_cell_Tc'],
                      'Tc_UO2_3x3x3': abl['uo2_324_cell_Tc'],
                      'Ea_ref': {m: {k: ea_ref[m][k] for k in ('ordered', 'superionic', 'ordered_T', 'superionic_T')} for m in ea_ref},
                      'Ea_oracle': o_ea, 'Ea_fluctTc': rob['tc_from_fluctuations']['activation_energy_eV'],
                      'Ea_window_10_80': rob['window_10_80']['activation_energy_eV'],
                      'Ea_2x2x2': abl['small_cell_Ea'], 'Ea_single': abl['single_arrhenius'],
                      'D_ref': d_ref,
                      'dex_oracle': {k: dex(o_pred[k], v) for k, v in d_ref.items()},
                      'dex_window_10_80': {k: dex(rob['window_10_80']['predictions'][k], v) for k, v in d_ref.items()},
                      'dex_arrhenius_fit': {k: dex(v, d_ref[k]) for k, v in rob['arrhenius_fit_predictions'].items()},
                      'dex_one_over_N': {k: dex(v, d_ref[k]) for k, v in rob['pbf2_one_over_N'].items()},
                      'dex_self_frame_AgI': {k: dex(v, d_ref[k]) for k, v in rob['self_frame'].items()},
                      'dex_ablations': {name: {k: dex(v[k], d_ref[k]) for k in d_ref}
                                        for name, v in abl.items() if isinstance(v, dict) and set(d_ref) <= set(v)},
                      'agi_ratio': ratio,
                      'catalog_cost_total': sum(c['cost'] for c in catalog_runs),
                      'oracle_cost': sum(cost(META[f'{m}-{L}-{T}']['natoms']) for (m, L), Ts in ORACLE_PLAN.items() for T in Ts),
                      'oracle_evals': sum(len(Ts) for Ts in ORACLE_PLAN.values()),
                      'data_json_MB': (TASK/'environment/lab/data.json').stat().st_size / 2 ** 20}, indent=1, default=float))


if __name__ == '__main__':
    main()
