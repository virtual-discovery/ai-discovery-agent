"""Write design.json (with the calibration table) from the private targets and the ablation values.

    python3 -B write_design.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/dense-hydrogen-llt-crossover'
H = json.loads((TASK/'environment/lab/private/hidden.json').read_text())
V = json.loads((TASK/'solution/ablation_values.json').read_text())
A, R, T, TH = V['ablations'], V['robustness'], H['targets'], H['thresholds']
r4 = lambda x: None if x is None else round(x, 4)
r2 = lambda x: None if x is None else round(x, 2)
ORACLE = json.loads((HERE/'oracle_values.json').read_text()) if (HERE/'oracle_values.json').exists() else {}

calibration = {}
for t in ('1100', '1200', '1400'):
    calibration[f'widom_rs_{t}K'] = {
        'reference': r4(T['widom_rs_N2048'][t]), 'oracle': r4(ORACLE.get('widom_rs_N2048', {}).get(t)),
        'robustness': {k: r4(R[k][t]) for k in ('widom_t0_5ps', 'widom_t0_20ps', 'widom_three_point_fit',
                                                 'widom_potential_energy_cv', 'widom_grid_maximum')},
        'nearest_ablation': {'whole_run_variance': r4(A['widom_whole_run'][t])},
        'tolerance': f"+-{TH['widom_rs']}"}
calibration['pimd_crossover_rs_1000K'] = {
    'reference': r4(T['pimd_crossover_rs_1000K']), 'oracle': r4(ORACLE.get('pimd_crossover_rs_1000K')),
    'robustness': {'tanh_fit': r4(R['pimd_tanh_fit']['rs']), 'no_transient_removal': r4(R['pimd_no_transient_removal']['rs'])},
    'nearest_ablation': {'nqe_ignored_classical_rs': r4(A['nqe_ignored']['rs'])}, 'tolerance': f"+-{TH['pimd_rs']}"}
calibration['nqe_shift_1000K_GPa'] = {
    'reference': r2(T['nqe_shift_GPa']['1000']), 'oracle': r2(ORACLE.get('nqe_crossover_pressure_shift_GPa', {}).get('1000')),
    'robustness': {'tanh_fit': r2(R['pimd_tanh_fit']['shift']['1000']), 'classical_N512': r2(R['classical_reference_N512']['1000']),
                   'no_transient_removal': r2(R['pimd_no_transient_removal']['shift']['1000'])},
    'nearest_ablation': {'fixed_density_difference': r2(A['nqe_fixed_density']['1000']), 'ignored': 0.0},
    'tolerance': f"+-{TH['nqe_shift_GPa']}"}
calibration['nqe_shift_900K_GPa'] = {'reference': None, 'nearest_ablation': {'mislabelled_files': r2(A['nqe_900_from_mislabelled_files']['900'])},
                                     'tolerance': 'must be null'}
calibration['classification'] = {'reference': T['classification'], 'nearest_ablation': {
    'small_cells': A['classification_small_cells'], 'absolute_S_k': A['classification_absolute_structure_factor']},
    'evidence_max_dP_drs_GPa_per_rs': {k: round(v) for k, v in H['reference']['pressure_slope_max_by_N_T'].items()},
    'tolerance': 'exact'}
calibration['low_density_phase_900K'] = {'reference': T['low_density_phase_900K'], 'nearest_ablation': 'liquid', 'tolerance': 'exact'}

design = {
 'task_id': 'dense-hydrogen-llt-crossover', 'catalog_id': 'dense-hydrogen-llt-crossover', 'start': 'released',
 'source_group': 'tenti_2025_hydrogen_llt', 'title': 'Dense hydrogen: liquid-liquid transition or crossover',
 'decision_points': [
  {'id': 'dp1', 'kind': 'evidence in the data (transient)',
   'fork': 'Whether to discard the start-up transient before computing c_v from energy fluctuations',
   'wrong_choice': 'Variance over whole runs (the 2048-atom runs at 1100-1800 K are only 20-60 ps long and start from a different state)',
   'consequence': 'c_v inflated 2-3x on the short runs; Widom-line rs too low by 0.010-0.020 (all three graded temperatures fail)',
   'ablation': 'widom_from_whole_runs'},
  {'id': 'dp2', 'kind': 'finite size / discriminating purchase',
   'fork': 'Read the van der Waals-like loops and large c_v peaks of N = 128/256 as first order, or test them at N = 512/2048',
   'wrong_choice': 'Classify from the cheap small cells',
   'consequence': 'First order at 950 and 1000 K (critical point near 1050 K) instead of crossover above ~950 K',
   'ablation': 'classification_from_small_cells'},
  {'id': 'dp3', 'kind': 'confounder (order parameter)',
   'fork': 'Liquid or solid on the molecular side at 900 K: diffusion (MSD) vs long-range order (max_k S(k)/N scaling)',
   'wrong_choice': 'Call it liquid because atoms diffuse (MSD grows by tens of A^2 over 40 ps)',
   'consequence': 'The 900 K first-order transition is misread as a liquid-liquid transition with a critical point',
   'ablation': 'molecular_side_liquid_because_it_diffuses'},
  {'id': 'dp4', 'kind': 'regime judgement (finite size of the order parameter)',
   'fork': 'Judge crystallinity from the absolute value of max_k S(k)/N or from its N dependence',
   'wrong_choice': 'Threshold S/N without scaling: N = 128 liquids give 0.05-0.3',
   'consequence': '1100 K molecular liquid classified as first order (solid-like)',
   'ablation': 'absolute_structure_factor_threshold'},
  {'id': 'dp5', 'kind': 'model comparison (NQE)',
   'fork': 'How to quantify nuclear quantum effects: compare crossover locations, compare pressures at fixed density, or ignore',
   'wrong_choice': 'Pressure difference at the classical crossover density, or no shift',
   'consequence': 'Shift of about -9 GPa or 0 instead of about -41 GPa; PIMD crossover rs taken as classical 1.432 instead of 1.488',
   'ablation': 'nqe_shift_at_fixed_density'},
  {'id': 'dp6', 'kind': 'evidence in metadata (mislabelled archive)',
   'fork': 'Trust the PIMD 900K folder label, or check the temperature column and file identity',
   'wrong_choice': 'Use the 900K-labelled files as 900 K path-integral data',
   'consequence': 'A spurious 900 K NQE shift (about -55 GPa); the files are byte-identical to 1000K files (Temp = 1000 K)',
   'ablation': 'nqe_900K_from_mislabelled_pimd_files'},
  {'id': 'dp7', 'kind': 'budget allocation',
   'fork': 'Where to spend 135 units: the 2048-atom isotherms cost 2 per density; the catalogue is worth about 1100 units',
   'wrong_choice': 'Coarse rs grids or small cells only',
   'consequence': 'Widom line from N = 128 (1100 K off by ~0.005) and small-cell loops read as first order',
   'ablation': 'small_cells_only_budget'}],
 'graded_fields': ['classification at 900/950/1000/1100 K (exact)', 'low_density_phase_900K (exact)',
                   'widom_rs_N2048 at 1100/1200/1400 K (+-0.005 rs)', 'pimd_crossover_rs_1000K (+-0.007 rs)',
                   'nqe_crossover_pressure_shift_GPa: 1000 K within 8 GPa; 900 K must be null'],
 'length_estimate': {'oracle_lab_evaluations': ORACLE.get('evaluations'), 'oracle_charge': ORACLE.get('charge'),
                     'budget': 135, 'catalog_total_cost': ORACLE.get('catalog_total'), 'tool_calls': '50-80',
                     'gpt6astra_cost_usd': '3-5', 'pilot': 'pending (coordinator)'},
 'runtime_dependencies': ['numpy (runtime image); no extra lab dependencies'],
 'data_provenance': {'archive_url': 'https://zenodo.org/records/14790037', 'doi': '10.5281/zenodo.14790037',
                     'paper_doi': '10.1103/pbrk-3zgd', 'arxiv': '2502.02447', 'license': 'CC BY 4.0',
                     'files_used': ['N_{128,256,512,2048}_T*K.tar (23 tars): input.lammps, log.lammps, slurm-*.out, '
                                    'six byte ranges of traj/trajectory.tar per run at 900-1100 K', 'PIMD.tar (all thermodata)'],
                     'member_checksums': ['members_thermo.json', 'members_frames.json', 'members_pimd.json'],
                     'not_used': ['aimd_*K.tar', 'mace_mod_loss_models.tar']},
 'exceptions': {},
 'notes': ('Design adaptations: the graded Widom line is the c_v maximum (the paper\'s estimator) rather than a compressibility '
           'maximum; the scaling of the peak is graded through the per-temperature classification (loops at N = 2048) rather '
           'than as a separate class; AIMD is not used (single-density 2-ps QE runs, ~3 GB each, cannot test an EOS kink); '
           'the PIMD 900K folder duplicates 1000K, so the 900 K NQE item is graded as "undeterminable" (null). '
           'dp7 (budget) fails mainly through its classification; its Widom values from N = 128 alone are within 0.006 rs.'),
 'calibration': calibration}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print(json.dumps(calibration, indent=1))
