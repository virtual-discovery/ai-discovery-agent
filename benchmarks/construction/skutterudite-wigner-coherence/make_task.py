"""Generate the staged, lab, private and ablation data for skutterudite-wigner-coherence from pool/.

    python3 -B make_task.py [POOL]      (numpy only; run in the runtime image, see README)

Writes environment/agent/catalog.json, environment/lab/data/*, environment/lab/private/*,
solution/ablation_values.json and the budget in both public.json files; prints the calibration table.
"""
import json
import sys
from pathlib import Path
import numpy as np
import model as M

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/skutterudite-wigner-coherence'
POOL = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/skutterudite-wigner-coherence/pool'
BUDGET = 32
COSTS = {'mode_data': 1, 'velocity_operator': 4, 'archived_file': 1, 'harmonic': 1}
HIDDEN_RUNS = ['IrSb3/hiphive', 'YbIr4Sb12/hiphive']
THRESHOLDS = {'kappa_P_rel': 0.05, 'kappa_C_rel': 0.06, 'T_half_K': 25.0, 'irsb3_kappa_P_50K_rel': 0.04,
              'irsb3_kappa_TOT_300K_rel': 0.05, 'rattler_rel': 0.20}
FORCE_CONSTANTS = {
    'CoSb3/tc_calc_no_isotope': ('QE DFPT, 3x3x3 q-grid (dim-fc2 3 3 3)', 'finite differences (thirdorder.py displacements, folder 4_d3)'),
    'CoSb3/tc_calc_with_isotope': ('QE DFPT, 3x3x3 q-grid (dim-fc2 3 3 3)', 'finite differences (thirdorder.py displacements, folder 4_d3)'),
    'CoSb3/v_operator_mean_free_path_isotope': ('QE DFPT, 3x3x3 q-grid (dim-fc2 3 3 3)', 'not read (linewidths read from the previous run with --read-gamma)'),
    'YbCo4Sb12/hiphive': ('QE DFPT, imported (folder import_phQE_to_fc2), dim-fc2 2 2 2', 'hiPhive (30 structures, cutoffs 7.9/6.0 A)'),
    'FeSb3/hiphive': ('QE DFPT, imported (folder import_phQE_to_fc2), dim-fc2 2 2 2', 'hiPhive (30 structures, cutoffs 7.9/6.0 A)'),
    'YbFe4Sb12-NM/hiphive_fc2hip_fc3hip': ('hiPhive, dim-fc2 2 2 2', 'hiPhive'),
    'YbFe4Sb12-NM/thirdorder_no_isotope': ('QE DFPT, 3x3x3 q-grid (dim-fc2 3 3 3)', 'finite differences (ShengBTE thirdorder.py, 6th neighbours)'),
    'YbFe4Sb12-NM/thirdorder_isotope': ('QE DFPT, 3x3x3 q-grid (dim-fc2 3 3 3)', 'finite differences (ShengBTE thirdorder.py, 6th neighbours)'),
    'IrSb3/hiphive': ('QE DFPT, imported (folder import_phQE_to_fc2), dim-fc2 3 3 3', 'hiPhive (30 structures, cutoffs 7.9/6.0 A)'),
    'YbIr4Sb12/hiphive': ('QE DFPT, imported (folder import_phQE_to_fc2), dim-fc2 2 2 2', 'hiPhive (30 structures, cutoffs 7.9/6.0 A)'),
}


RATTLER = {}
KB, AMU = 1.380649e-23, 1.66053907e-27


def rattler_estimates(harmonic, meta, Y):
    """Yb-projected DOS maximum of YbIr4Sb12 on the archived axis, and on the corrected axis.

    The PDOS runs read eV/A^2 force constants with phonopy's QE calculator (factor 108.970772 instead of the
    15.633302 of the phono3py runs), so the archived frequency axis is inflated by their ratio."""
    h = harmonic['YbIr4Sb12']
    axis = np.array(h['pdos_frequency_axis']); yb = np.array(h['pdos_by_element']['Yb'])
    ratio = h['pdos_run']['phonopy_yaml']['frequency_unit_conversion_factor'] / meta['YbIr4Sb12/hiphive']['cell']['frequency_unit_conversion_factor']
    peak = axis[np.argmax(yb)]
    msd = h['msd_A2']['300']['Yb'] * 1e-20
    mass = [m for s, m in zip(meta['YbIr4Sb12/hiphive']['cell']['symbols'], meta['YbIr4Sb12/hiphive']['cell']['masses_amu']) if s == 'Yb'][0]
    einstein = np.sqrt(3 * KB * 300 / (mass * AMU * msd)) / (2 * np.pi) / 1e12
    gamma_optical = float(sorted(f for f in Y['frequency'][0] if f > 0.05)[0])
    return {'archived_axis_peak': float(peak), 'axis_ratio': float(ratio), 'pdos_peak': float(peak / ratio),
            'einstein_from_msd_300K': float(einstein), 'lowest_gamma_optical': gamma_optical}


def tag(run):
    return run.replace('/', '__')


def load(run):
    return dict(np.load(POOL/(tag(run) + '.npz')))


def series(d, vol, iso_mode='zero', temps=None, **kw):
    gi = d.get('gamma_isotope')
    if gi is not None and iso_mode is not None and iso_mode != 'archived':
        gi = M.repaired_isotope(gi, iso_mode)[0]
    G = d['G'] if 'G' in d else d['v2']
    T = d['temperature']
    out = {}
    for it, t in enumerate(T):
        if temps is not None and int(t) not in temps:
            continue
        out[int(t)] = M.kappa(d['frequency'], d['gamma'][it], gi, d['heat_capacity'][it], G, d['weight'], vol, **kw)
    return out


def targets(Y, Yv, I, Iv, y_kw=None, i_kw=None, thalf_mode='ratio', rattler=None):
    y_kw = y_kw or {}; i_kw = i_kw or {}
    sy = series(Y, Yv, **y_kw)
    si = series(I, Iv, temps={50, 300}, **i_kw)
    T = sorted(sy)
    kp = np.array([sy[t][0] for t in T]); kc = np.array([sy[t][1] for t in T])
    th = M.t_half(T, kp, kc, thalf_mode)
    if not np.isfinite(th):
        th = M.t_half(T, kp, kc, 'extrapolate')
    return {'YbIr4Sb12': {'kappa_P_W_mK': {'300': sy[300][0], '700': sy[700][0]},
                          'kappa_C_W_mK': {'100': sy[100][1], '300': sy[300][1], '700': sy[700][1]},
                          'T_half_K': th, 'rattler_frequency_THz': RATTLER['pdos_peak'] if rattler is None else rattler},
            'IrSb3': {'kappa_P_50K_W_mK': si[50][0], 'kappa_TOT_300K_W_mK': si[300][0] + si[300][1]}}


def flat(t):
    y, i = t['YbIr4Sb12'], t['IrSb3']
    return [y['kappa_P_W_mK']['300'], y['kappa_P_W_mK']['700'], y['kappa_C_W_mK']['100'], y['kappa_C_W_mK']['300'],
            y['kappa_C_W_mK']['700'], y['T_half_K'], y['rattler_frequency_THz'], i['kappa_P_50K_W_mK'], i['kappa_TOT_300K_W_mK']]


def errors(t, ref):
    a, b = flat(t), flat(ref)
    names = ['kP300', 'kP700', 'kC100', 'kC300', 'kC700', 'T_half', 'rattler', 'Ir_kP50', 'Ir_kTOT300']
    return {n: (round(x - y, 1) if n == 'T_half' else round(x / y - 1, 4)) for n, x, y in zip(names, a, b)}


def r7(a, sig=7):
    a = np.asarray(a, float)
    return np.vectorize(lambda x: float('%.*g' % (sig, x)))(a) if a.size else a


def main():
    meta = json.loads((POOL/'runs.json').read_text())
    harmonic = json.loads((POOL/'harmonic.json').read_text())
    runs = list(meta)
    data_dir = TASK/'environment/lab/data'
    data_dir.mkdir(parents=True, exist_ok=True)
    # 1. formula check: the reference evaluator reproduces every archived 2.1.0 table with the archived inputs
    for run in runs:
        d = load(run)
        if 'G' in d:
            vol = M.volume(meta[run]['cell']['lattice_A'])
            s = series(d, vol, iso_mode='archived')
            for it, t in enumerate(d['temperature']):
                kp, kc = s[int(t)]
                assert abs(kp / d['kappa_P_RTA'][it] - 1) < 1e-6 and abs(kc / d['kappa_C'][it] - 1) < 1e-6, (run, t)
        else:   # 1.19 runs: populations from archived group velocities equal the log K_P
            vol = M.volume(meta[run]['cell']['lattice_A'])
            logname = [f for f, v in meta[run]['files'].items() if v['kind'] == 'phono3py log'][0]
            kp = series(d, vol, iso_mode=None, temps={300}, iso=False)[300][0]
            logged = meta[run]['files'][logname]['tables']['K_P']['300']
            if run.endswith('with_isotope') or run.endswith('thirdorder_isotope'):
                # archive duplicate: the 'isotope' folder holds a byte-identical copy of the no-isotope kappa file
                sibling = run.replace('with_isotope', 'no_isotope').replace('thirdorder_isotope', 'thirdorder_no_isotope')
                assert meta[run]['kappa_md5_group'] == meta[sibling]['kappa_md5_group']
                iso_log = [v for f, v in meta[run]['files'].items() if v['kind'] == 'phono3py log' and
                           any('Mass variance' in h for h in v['header'])][0]['tables']['K_P']['300']
                assert kp - iso_log > 0.005, (run, kp, iso_log)
            else:
                assert abs(kp - logged) < 6e-3, run
    # 2. lab data (per-mode arrays; no archived kappa arrays) and lab metadata
    lab_meta = {'runs': {}, 'harmonic': {}}
    for run in runs:
        d = load(run)
        out = {k: d[k] for k in ('frequency', 'gamma', 'heat_capacity', 'qpoint', 'weight', 'temperature', 'mesh')}
        if 'gamma_isotope' in d:
            out['gamma_isotope'] = d['gamma_isotope']
        if 'G' in d:
            out['G'] = d['G']
            out['v2'] = np.einsum('qss->qs', d['G'])
        else:
            out['v2'] = d['v2']
        np.savez_compressed(data_dir/(tag(run) + '.npz'), **out)
        m = meta[run]
        entry = {'compound': m['compound'], 'velocity_operator': 'G' in d, 'withheld_outputs': run in HIDDEN_RUNS,
                 'files': {} if run in HIDDEN_RUNS else {f: v for f, v in m['files'].items() if v['kind'] != 'launch script'}}
        lab_meta['runs'][run] = entry
    for comp, h in harmonic.items():
        lab_meta['harmonic'][comp] = h
    (data_dir/'meta.json').write_text(json.dumps(lab_meta) + '\n')
    # 3. catalog
    cat_runs = []
    for run in runs:
        m = meta[run]
        fc2, fc3 = FORCE_CONSTANTS[run]
        launch = [c for f, v in m['files'].items() if v['kind'] == 'launch script' for c in v['commands']]
        cat_runs.append({
            'run': run, 'compound': m['compound'], 'archive_folder': m['archive_folder'],
            'kappa_file': m['kappa_file'], 'phono3py_version': m['cell']['version'],
            'frequency_unit_conversion_factor': m['cell']['frequency_unit_conversion_factor'], 'launch_commands': launch,
            'fc2': fc2, 'fc3': fc3, 'q_mesh': m['mesh'], 'n_irreducible_q': m['n_irreducible_q'], 'n_bands': m['n_bands'],
            'temperatures_K': m['temperatures_K'],
            'primitive_lattice_A': m['cell']['lattice_A'],
            'atoms': [{'symbol': s, 'mass_amu': w} for s, w in zip(m['cell']['symbols'], m['cell']['masses_amu'])],
            'velocity_operator_archived': m['has_velocity_operator'],
            'archived_files': [] if run in HIDDEN_RUNS else sorted(f for f, v in m['files'].items() if v['kind'] != 'launch script'),
            **({'note': 'Archived outputs (logs and thermal-conductivity tables) of this run are withheld; per-mode data and the velocity operator are served.'}
               if run in HIDDEN_RUNS else {})})
    catalog = {
        'source': 'E. Di Lucente, M. Simoncelli, N. Marzari, Phys. Rev. Research 5, 033125 (2023); Materials Cloud Archive 5jsfj-cdg70 (CC BY 4.0)',
        'archive_readme': (POOL.parent/'raw/README.txt').read_text() if (POOL.parent/'raw/README.txt').exists() else '',
        'compounds': ['CoSb3', 'YbCo4Sb12', 'FeSb3', 'YbFe4Sb12-NM', 'YbFe4Sb12-AFM', 'IrSb3', 'YbIr4Sb12'],
        'runs': cat_runs,
        'harmonic_available': {c: sorted(k for k in h if k != 'atoms') for c, h in harmonic.items()},
        'requests': {
            'mode_data': {'format': {'kind': 'mode_data', 'run': '<run>', 'temperature_K': '<one of temperatures_K>'}, 'cost': COSTS['mode_data'],
                          'returns': 'archived per-mode arrays on the irreducible q-points of the run mesh, indexed [q][band]: qpoints '
                                     '(reciprocal-lattice coordinates), weights, frequency_THz, gamma_THz (phono3py gamma at that '
                                     'temperature), gamma_isotope_THz (phono3py gamma_isotope; null if the run file has none), '
                                     'heat_capacity_eV_per_K (mode heat capacity), v2_A2THz2 = (1/3) sum_a |v_a|^2 of the band '
                                     'velocity (1.19 runs: archived group velocities; 2.1.0 runs: diagonal of the velocity operator).'},
            'velocity_operator': {'format': {'kind': 'velocity_operator', 'run': '<run with velocity_operator_archived>'}, 'cost': COSTS['velocity_operator'],
                                  'returns': 'G_A2THz2[q][s][s2] = (1/3) sum_a |V_a(q)_{s,s2}|^2 of the archived phono3py velocity '
                                             'operator at each irreducible q (diagonal = v2_A2THz2), with qpoints, weights and frequency_THz.'},
            'archived_file': {'format': {'kind': 'archived_file', 'run': '<run>', 'file': '<one of archived_files>'}, 'cost': COSTS['archived_file'],
                              'returns': 'phono3py logs: header lines and the printed K_P, K_C, K_T tables (xx component, W/m/K); '
                                         'other files: their text.'},
            'harmonic': {'format': {'kind': 'harmonic', 'compound': '<compound>'}, 'cost': COSTS['harmonic'],
                         'returns': 'the archived phonopy projected DOS (projected_dos.dat: its first column as written, and the '
                                    'per-atom columns summed by element) with the command and phonopy.yaml settings of that run; '
                                    'mean-square displacements (A^2, trace per atom, mean over the atoms of each element) at 300 and '
                                    '800 K with the command of that run; for YbFe4Sb12 the DFPT Gamma-point frequencies (cm^-1) of the '
                                    'nonmagnetic and antiferromagnetic (DFT+U) cells.'}},
        'budget': BUDGET,
        'units': {'frequency': 'THz (ordinary frequency, not angular)', 'velocity': 'A*THz = 100 m/s',
                  'heat_capacity': 'eV/K per mode (phono3py)', 'temperature': 'K'}}
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    # 4. reference values and wrong paths
    Y, I = load('YbIr4Sb12/hiphive'), load('IrSb3/hiphive')
    Yv, Iv = M.volume(meta['YbIr4Sb12/hiphive']['cell']['lattice_A']), M.volume(meta['IrSb3/hiphive']['cell']['lattice_A'])
    bad_y = M.repaired_isotope(Y['gamma_isotope'])[1]
    bad_i = M.repaired_isotope(I['gamma_isotope'])[1]
    assert bad_i.sum() == 0 and set(np.where(bad_y)[0]) == {1}
    RATTLER.update(rattler_estimates(harmonic, meta, Y))
    ref = targets(Y, Yv, I, Iv)
    corrupt_q = sorted(set(np.where(bad_y)[0].tolist()))
    abl = {
        'archived_isotope_rates_as_is': targets(Y, Yv, I, Iv, y_kw={'iso_mode': 'archived'}, i_kw={'iso_mode': 'archived'}),
        'isotope_scattering_omitted': targets(Y, Yv, I, Iv, y_kw={'iso': False}, i_kw={'iso': False}),
        'lifetime_convention_from_tc_data': targets(Y, Yv, I, Iv, y_kw={'tau_factor': 2.0}, i_kw={'tau_factor': 2.0}),
        'degenerate_pairs_counted_as_coherences': targets(Y, Yv, I, Iv, y_kw={'coh_degenerate': True}, i_kw={'coh_degenerate': True}),
        'gaussian_band_overlap': targets(Y, Yv, I, Iv, y_kw={'shape': 'gauss'}, i_kw={'shape': 'gauss'}),
        'half_width_lorentzian': targets(Y, Yv, I, Iv, y_kw={'width_factor': 1.0}, i_kw={'width_factor': 1.0}),
        'corrupted_qpoint_discarded': targets(Y, Yv, I, Iv, y_kw={'iso_mode': 'archived', 'drop_q': corrupt_q}),
        'pdos_axis_taken_as_THz': targets(Y, Yv, I, Iv, rattler=RATTLER['archived_axis_peak']),
    }
    rob = {
        'isotope_dropped_for_filled_compound_only': targets(Y, Yv, I, Iv, y_kw={'iso': False}),
        'isotope_repaired_with_band_median': targets(Y, Yv, I, Iv, y_kw={'iso_mode': 'band_median'}),
        'diagonal_populations': targets(Y, Yv, I, Iv, y_kw={'pop_block': False}, i_kw={'pop_block': False}),
        'degeneracy_threshold_0p1_cm1': targets(Y, Yv, I, Iv, y_kw={'deg_thr': 3.0e-3}, i_kw={'deg_thr': 3.0e-3}),
        'T_half_from_separate_interpolation': targets(Y, Yv, I, Iv, thalf_mode='separate'),
        'rattler_from_yb_msd_einstein': targets(Y, Yv, I, Iv, rattler=RATTLER['einstein_from_msd_300K']),
        'rattler_from_lowest_gamma_optical_mode': targets(Y, Yv, I, Iv, rattler=RATTLER['lowest_gamma_optical']),
    }
    hidden = {'reference': ref, 'thresholds': THRESHOLDS, 'rattler_estimates': RATTLER,
              'notes': {'isotope_repair': 'archived gamma_isotope > 1 THz (all at the q-point nearest Gamma, index %s) set to zero' % corrupt_q,
                        'archived_kappa_P_YbIr4Sb12_300K': float(Y.get('kappa_P_RTA', [np.nan] * 6)[5]) if 'kappa_P_RTA' in Y else None}}
    (TASK/'environment/lab/private').mkdir(parents=True, exist_ok=True)
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(flat(ref)) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'reference': ref, 'ablations': abl, 'robustness': rob}, indent=1) + '\n')
    for name in ('environment/lab/public.json', 'tests/public.json'):
        p = TASK/name
        pub = json.loads(p.read_text())
        pub['solver_budget'] = BUDGET
        if name.startswith('environment'):
            pub['requests'] = {k: {'cost': v} for k, v in COSTS.items()}
            pub['catalog'] = '/app/catalog.json'
        p.write_text(json.dumps(pub, indent=1) + '\n')
    print(json.dumps({'reference': ref, 'ablation_errors': {k: errors(v, ref) for k, v in abl.items()},
                      'robustness_errors': {k: errors(v, ref) for k, v in rob.items()}}, indent=1))


if __name__ == '__main__':
    main()
