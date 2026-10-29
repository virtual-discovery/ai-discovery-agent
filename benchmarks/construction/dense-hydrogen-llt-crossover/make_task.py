"""Write the staged, lab, private and ablation data of dense-hydrogen-llt-crossover from the pool.

    docker run ... python3 -B make_task.py     (HCACHE = cache folder with pool_*.{npz,json}; TASK = task folder)

Prints the calibration numbers.
"""
import json
import os
import shutil
from pathlib import Path
import numpy as np
import analysis as A

HERE = Path(__file__).resolve().parent
TASK = Path(os.environ['TASK']) if 'TASK' in os.environ else HERE.parents[1]/'tasks/dense-hydrogen-llt-crossover'
CACHE = A.CACHE
BUDGET = 135
WIDOM_T = ['1100', '1200', '1400']
THRESHOLDS = {'widom_rs': 0.005, 'pimd_rs': 0.007, 'nqe_shift_GPa': 8.0}


def pimd_sets():
    sets = {}
    for key, m in A.META['pimd'].items():
        s = sets.setdefault(f"{m['T_label']}/{m['rs']}", {'replicas': {}})
        s['replicas'][str(m['replica'])] = m['data_key']
    return sets


def main():
    lab = TASK/'environment/lab'
    (lab/'data').mkdir(exist_ok=True)
    shutil.copy(CACHE/'pool_thermo.npz', lab/'data/thermo.npz')
    shutil.copy(CACHE/'pool_pimd.npz', lab/'data/pimd.npz')
    struct = json.loads((CACHE/'pool_structure.json').read_text())
    (lab/'data/structure.json').write_text(json.dumps(struct, separators=(',', ':')))
    runs = {k: {f: m[f] for f in ('N', 'T', 'rs', 'first_step', 'stride_steps', 'E_ref_eV', 'volume_A3', 'input_deck', 'jobs')}
            for k, m in A.META['runs'].items()}
    sets = pimd_sets()
    (lab/'data/meta.json').write_text(json.dumps({'runs': runs, 'pimd_sets': sets}, separators=(',', ':')))

    cost_t = {128: 1, 256: 1, 512: 2, 2048: 2}
    cost_s = {128: 2, 256: 2, 512: 3, 2048: 4}
    isotherms = []
    for n in (128, 256, 512, 2048):
        for t in sorted({m['T'] for m in A.META['runs'].values() if m['N'] == n}):
            rr = A.rs_list(n, t)
            isotherms.append({'N': n, 'T_K': t, 'thermo_cost': cost_t[n],
                              'runs': [{'rs': r, 'length_ps': A.META['runs'][f'{n}/{t}/{r}']['length_ps'],
                                        'structure_available': f'{n}/{t}/{r}' in struct} for r in rr]})
    pimd_cat = {}
    for key, s in sets.items():
        lab_t, rs = key.split('/')
        pimd_cat.setdefault(lab_t, []).append({'rs': rs, 'replicas': len(s['replicas'])})
    for v in pimd_cat.values():
        v.sort(key=lambda d: float(d['rs']))
    catalog = {
        'system': 'Dense fluid hydrogen; MACE machine-learning potential trained on PBE (Tenti, Jaeckl, Nakano, Rupp, Casula 2025).',
        'classical_md': 'LAMMPS NVT, classical nuclei, cubic cells of N atoms, timestep 0.2 fs, CSVR thermostat (tau 0.1 ps), '
                        'metal units. rs is the Wigner-Seitz radius in bohr: (4 pi / 3) (rs a0)^3 = V / N.',
        'pimd': 'Path-integral (PIOUD) runs with the same potential, N = 512, 12 beads, timestep 0.15 fs; several independent '
                'replicas per density started from decorrelated classical configurations. Archived under a temperature '
                'label (folder name); the archived columns are listed in the response.',
        'requests': {
            'thermo': {'request': {'kind': 'thermo', 'N': 2048, 'T_K': 1000, 'rs': '1.43'}, 'cost': 'thermo_cost of the isotherm',
                       'returns': 'thermodynamic time series (every 500 steps) of the run, its LAMMPS input deck and job history'},
            'structure': {'request': {'kind': 'structure', 'N': 2048, 'T_K': 900, 'rs': '1.425'},
                          'cost': {str(k): v for k, v in cost_s.items()},
                          'returns': 'observables from about 60 frames sampled from the released trajectory in six bursts of '
                                     'consecutive 20-fs dumps: per frame time, max_k S(k)/N (|k| <= 5.5 1/A, k = 2 pi n/L) and '
                                     'its |k|, fraction of atoms with a neighbour within 1.05 A; stable molecular fraction per '
                                     '80-fs window (pairs within 1.05 A in five consecutive dumps); g(r) per burst; MSD vs lag'},
            'pimd': {'request': {'kind': 'pimd', 'T_label': '1000K', 'rs': '1.48'}, 'cost': 1,
                     'returns': 'thermodynamic series (every 20 moves) of every replica at this density'}},
        'budget': BUDGET,
        'isotherms': isotherms,
        'pimd_available': pimd_cat}
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    for p in (lab/'public.json', TASK/'tests/public.json'):
        d = json.loads(p.read_text())
        d['solver_budget'] = BUDGET
        if p.parent.name == 'lab':
            d['requests'] = {'thermo': {'kind': 'thermo', 'N': '<int>', 'T_K': '<int>', 'rs': '<str>', 'cost': '1 (N<=256) or 2'},
                             'structure': {'kind': 'structure', 'N': '<int>', 'T_K': '<int>', 'rs': '<str>', 'cost': '2-4'},
                             'pimd': {'kind': 'pimd', 'T_label': '<str>', 'rs': '<str>', 'cost': 1}}
            d['catalog'] = '/app/catalog.json'
        p.write_text(json.dumps(d, indent=1) + '\n')

    # ---- references (from the served data) ----
    widom = {t: A.widom_rs(int(t)) for t in WIDOM_T}
    rs_cl, p_cl = A.classical_crossover(1000)
    rs_q, p_q = A.pimd_crossover('1000K')
    shift = p_q - p_cl
    targets = {'classification': {'900': 'first_order', '950': 'crossover', '1000': 'crossover', '1100': 'crossover'},
               'low_density_phase_900K': 'crystal', 'widom_rs_N2048': widom, 'pimd_crossover_rs_1000K': rs_q,
               'nqe_shift_GPa': {'1000': shift, '900': None}}
    # evidence asserted at construction
    loops = {}
    for t in (900, 950, 1000, 1100):
        for n in (128, 256, 512, 2048):
            rr = A.rs_list(n, t)
            x = np.array([float(r) for r in rr])
            p = np.array([A.pressure(n, t, r) for r in rr])
            loops[f'{n}/{t}'] = float((np.diff(p) / np.diff(x)).max())
    assert loops['2048/900'] > 0 and loops['2048/950'] < 0 and loops['2048/1000'] < 0 and loops['2048/1100'] < 0
    assert loops['128/1000'] > 0 and loops['256/1000'] > 0 and loops['256/950'] > 0
    sk = lambda k: float(np.median([r[1] for r in struct[k]['frames'] if r[0] > 30]))
    assert sk('2048/900/1.425') > 0.08 and sk('2048/1000/1.445') < 0.02 and sk('2048/1100/1.46') < 0.02
    dup = [k for k, m in A.META['pimd'].items() if m['T_label'] == '900K']
    assert all(A.META['pimd'][k]['data_key'].startswith('1000K') for k in dup)
    # ---- ablations / robustness ----
    abl = {'widom_whole_run': {t: A.widom_rs(int(t), t0=0.0) for t in WIDOM_T},
           'classification_small_cells': {'900': 'first_order', '950': 'first_order', '1000': 'first_order', '1100': 'crossover'},
           'classification_absolute_structure_factor': {'900': 'first_order', '950': 'crossover', '1000': 'crossover', '1100': 'first_order'},
           'phase_liquid_from_diffusion': 'liquid',
           'widom_N128': {t: A.widom_rs(int(t), n=128) for t in WIDOM_T},
           'nqe_ignored': {'rs': rs_cl, 'shift': {'1000': 0.0, '900': 0.0}},
           'nqe_fixed_density': {'1000': float(np.interp(rs_cl, *A.pimd_eos('1000K')) - p_cl), '900': None},
           'nqe_900_from_mislabelled_files': {'1000': shift, '900': A.pimd_crossover('900K')[1] - A.classical_crossover(900)[1]}}
    r_tanh, p_tanh = A.pimd_tanh_crossover('1000K')

    def widom_pe(t):
        rr = A.rs_list(2048, t)
        c = []
        for r in rr:
            time, _, pe, _, _ = A.series(2048, t, r)
            c.append(np.var(pe[time >= 10]) / (2048 * A.KB ** 2 * t * t) + 1.5)
        return A.peak([float(r) for r in rr], c)
    rob = {'widom_t0_5ps': {t: A.widom_rs(int(t), t0=5.0) for t in WIDOM_T},
           'widom_t0_20ps': {t: A.widom_rs(int(t), t0=20.0) for t in WIDOM_T},
           'widom_three_point_fit': {t: A.widom_rs(int(t), npts=3) for t in WIDOM_T},
           'widom_potential_energy_cv': {t: widom_pe(int(t)) for t in WIDOM_T},
           'widom_grid_maximum': {t: float(max(A.rs_list(2048, int(t)), key=lambda r: A.cv(2048, int(t), r))) for t in WIDOM_T},
           'pimd_tanh_fit': {'rs': r_tanh, 'shift': {'1000': p_tanh - p_cl, '900': None}},
           'classical_reference_N512': {'1000': p_q - A.classical_crossover(1000, 512)[1], '900': None},
           'pimd_no_transient_removal': {'rs': A.pimd_crossover('1000K', 0.0)[0],
                                         'shift': {'1000': A.pimd_crossover('1000K', 0.0)[1] - p_cl, '900': None}}}
    hidden = {'targets': targets, 'thresholds': THRESHOLDS,
              'reference': {'classical_crossover_1000K': [rs_cl, p_cl], 'pimd_crossover_1000K': [rs_q, p_q],
                            'pressure_slope_max_by_N_T': loops}}
    (lab/'private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = list(widom.values()) + [rs_q, shift, p_q, p_cl]   # rs_cl (1.4319) is not graded and coincides with grid labels
    (lab/'private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1) + '\n')
    print(json.dumps({'targets': targets, 'ablations': abl, 'robustness': rob, 'loops': loops}, indent=1, default=float))


if __name__ == '__main__':
    main()
