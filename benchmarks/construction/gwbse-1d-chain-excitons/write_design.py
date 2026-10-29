"""Write design.json for gwbse-1d-chain-excitons from the private targets and the ablation values.

    cd benchmarks/construction/gwbse-1d-chain-excitons && python3 -B write_design.py
Run after make_task.py; oracle_run.json holds the oracle's submitted values and purchase counts (from the gate).
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/gwbse-1d-chain-excitons'
H = json.loads((TASK/'environment/lab/private/hidden.json').read_text())
V = json.loads((TASK/'solution/ablation_values.json').read_text())
O = json.loads((HERE/'oracle_run.json').read_text()) if (HERE/'oracle_run.json').exists() else {}
T, TH, A, R = H['targets'], H['thresholds'], V['ablations'], V['robustness']
r3 = lambda x: round(float(x), 3)
osub = O.get('submission', {})
DSOC_BI = H['reference_notes']['soc_gap_change_eV']['Bi2Te3']


def oracle(path):
    node = osub
    for key in path:
        node = node.get(key, {}) if isinstance(node, dict) else {}
    return r3(node) if isinstance(node, (int, float)) else None


calibration = {
    'scgw_direct_gap_eV Te3 with_soc': {
        'reference': r3(T['scgw_direct_gap_eV']['Te3']['with_soc']), 'oracle': oracle(['scgw_direct_gap_eV', 'Te3', 'with_soc']),
        'robustness': {'fundamental_gap_instead_of_direct (convention deviation)': r3(R['te3_fundamental_gap'])},
        'near_miss': {'first IP peak of the spectrum (bright onset above the weak band-edge transition)': r3(V['near_miss']['te3_gap_from_first_ip_peak'])},
        'nearest_ablation': {'soc_added_again': r3(A['soc_added_again_gaps']['Te3']['with_soc'])},
        'tolerance': f"+-{TH['gap_abs_eV']}"},
    'scgw_direct_gap_eV Te3 without_soc': {
        'reference': r3(T['scgw_direct_gap_eV']['Te3']['without_soc']), 'oracle': oracle(['scgw_direct_gap_eV', 'Te3', 'without_soc']),
        'robustness': {'k_resolved_corrections_on_soc_free_bands': r3(R['te3_without_soc_k_resolved']),
                       'soc_change_at_the_qp_minimum_k': r3(R['te3_without_soc_fixed_k'])},
        'nearest_ablation': {'soc_added_again': r3(A['soc_added_again_gaps']['Te3']['without_soc'])},
        'tolerance': f"+-{TH['gap_abs_eV']}"},
    'scgw_direct_gap_eV Bi2Te3 with_soc': {
        'reference': r3(T['scgw_direct_gap_eV']['Bi2Te3']['with_soc']), 'oracle': oracle(['scgw_direct_gap_eV', 'Bi2Te3', 'with_soc']),
        'robustness': {'ip_onset_1pct': r3(R['bi2te3_gap_ip_onset']), 'first_ip_peak': r3(R['bi2te3_gap_ip_first_peak'])},
        'nearest_ablation': {'soc_added_again': r3(A['soc_added_again_gaps']['Bi2Te3']['with_soc']),
                             'g0w0_gap': r3(A['g0w0_gap_for_bi2te3']),
                             'increment_transferred': r3(A['increment_transfer_gap'])},
        'tolerance': f"+-{TH['gap_abs_eV']}"},
    'scgw_direct_gap_eV Bi2Te3 without_soc': {
        'reference': r3(T['scgw_direct_gap_eV']['Bi2Te3']['without_soc']), 'oracle': oracle(['scgw_direct_gap_eV', 'Bi2Te3', 'without_soc']),
        'robustness': {'ip_onset_1pct': r3(R['bi2te3_gap_ip_onset'] - DSOC_BI), 'first_ip_peak': r3(R['bi2te3_gap_ip_first_peak'] - DSOC_BI),
                       'soc_change_at_the_qp_minimum_k': r3(R['bi2te3_without_soc_fixed_k'])},
        'nearest_ablation': {'soc_added_again': r3(A['soc_added_again_gaps']['Bi2Te3']['without_soc']),
                             'g0w0_gap': r3(A['g0w0_gap_for_bi2te3'] - DSOC_BI)},
        'tolerance': f"+-{TH['gap_abs_eV']}"},
    'scgw_optical_gap_eV': {
        'reference': {k: r3(v) for k, v in T['scgw_optical_gap_eV'].items()},
        'oracle': {k: oracle(['scgw_optical_gap_eV', k]) for k in T['scgw_optical_gap_eV']},
        'robustness': 'peak position on the 5 meV Haydock grid; any bright threshold between 3% and 30% gives the same peaks',
        'nearest_ablation': {'first_spectral_feature As2S3': r3(A['first_feature_optical_As2S3']),
                             'BSE_G0W0 run read for Te3': 2.716},
        'tolerance': f"+-{TH['optical_gap_abs_eV']}"},
    'screening_update_shift_eV': {
        'reference': {k: r3(v) for k, v in T['screening_update_shift_eV'].items()},
        'oracle': {k: oracle(['screening_update_shift_eV', k]) for k in T['screening_update_shift_eV']},
        'robustness': {'S3 with the IP-peak shift instead of the QP gap shift': r3(R['s3_dw_ip_shift'])},
        'nearest_ablation': {'run_labels_trusted Te3': r3(A['labels_trusted_dw']['Te3']),
                             'run_labels_trusted As2S3': r3(A['labels_trusted_dw']['As2S3']),
                             'S3 runs assumed to share energies': r3(A['s3_same_energies_dw']),
                             'first_spectral_feature As2S3': r3(A['first_feature_dw']['As2S3'])},
        'tolerance': f"+-{TH['screening_shift_abs_eV']}"},
    'bi2te3_binding_energy_eV': {
        'reference': r3(T['bi2te3_binding_energy_eV']), 'oracle': oracle(['bi2te3_binding_energy_eV']),
        'robustness': {'no_W_correction': r3(R['eb_no_w_correction']), 'additive_mean_shift': r3(R['eb_additive_mean_dw']),
                       'mean_ratio': r3(R['eb_mean_ratio']), 'ratio_scaled_by_self_consistency': r3(R['eb_ratio_scaled_by_self_consistency'])},
        'nearest_ablation': {'labels_trusted (G0W0 gap - BSE_G0W0 peak)': r3(A['labels_trusted_eb']),
                             'wire_model_calibrated_on_S3_Te3': r3(A['wire_model_calibrated_eb']),
                             'wire_model': r3(A['wire_model_eb']), 'mean_released': r3(A['released_mean_eb'])},
        'tolerance': f"window {TH['bi2te3_binding_window_eV']} (regime check)"},
}

design = {
    'task_id': 'gwbse-1d-chain-excitons', 'catalog_id': 'gwbse-1d-chain-excitons', 'start': 'released',
    'source_group': 'grillo_2026_1d_gwbse',
    'title': 'Exfoliable 1D chains: GW-BSE excitons, self-consistency and SOC',
    'decision_points': [
        {'id': 'dp1', 'fork': 'Which QP energies each archived BSE run used (BSE_G0W0 label versus the KfnQPdb header and the IP spectrum)',
         'wrong_choice': 'Trust the labels: treat the BSE_G0W0 runs of Te3 and As2S3 as G0W0-energy runs',
         'consequence': f"Delta_W Te3 {r3(A['labels_trusted_dw']['Te3'])} and As2S3 {r3(A['labels_trusted_dw']['As2S3'])} eV instead of "
                        f"{r3(T['screening_update_shift_eV']['Te3'])} and {r3(T['screening_update_shift_eV']['As2S3'])}; "
                        f"Bi2Te3 E_b {r3(A['labels_trusted_eb'])} eV from the G0W0 gap",
         'ablation': 'run_labels_trusted'},
        {'id': 'dp2', 'fork': 'Whether every chain\'s two runs share QP energies',
         'wrong_choice': 'Generalize the Te3 finding to S3 and take the raw difference of its two runs',
         'consequence': f"Delta_W S3 {r3(A['s3_same_energies_dw'])} eV (wrong sign) instead of {r3(T['screening_update_shift_eV']['S3'])}",
         'ablation': 's3_runs_assumed_same_energies'},
        {'id': 'dp3', 'fork': 'Whether the archived GW already includes SOC (Eo versus the SOC and SOC-free band structures)',
         'wrong_choice': 'Assume GW without SOC, as is common, and add the DFT SOC change to the archived gaps',
         'consequence': f"Te3 gaps low by {r3(-H['reference_notes']['soc_gap_change_eV']['Te3'])} eV and Bi2Te3 by "
                        f"{r3(-H['reference_notes']['soc_gap_change_eV']['Bi2Te3'])} eV, with and without SOC",
         'ablation': 'soc_added_again'},
        {'id': 'dp4', 'fork': 'How to obtain the withheld Bi2Te3 scGW gap: transfer the self-consistency increment or use Bi2Te3\'s own BSE_G0W0 run',
         'wrong_choice': 'Add the mean G0W0-to-scGW increment of S3, Te3 and As2S3 to the Bi2Te3 G0W0 gap',
         'consequence': f"Bi2Te3 gap {r3(A['increment_transfer_gap'])} eV instead of {r3(T['scgw_direct_gap_eV']['Bi2Te3']['with_soc'])} "
                        f"(its own increment is {r3(H['reference_notes']['self_consistency_gap_shift_eV']['Bi2Te3'])} eV)",
         'ablation': 'self_consistency_increment_transferred'},
        {'id': 'dp5', 'fork': 'Bright exciton versus weak or dark features below it',
         'wrong_choice': 'Take the first visible feature of the spectrum',
         'consequence': f"As2S3 optical gap {r3(A['first_feature_optical_As2S3'])} eV instead of "
                        f"{r3(T['scgw_optical_gap_eV']['As2S3'])}; As2S3 Delta_W {r3(A['first_feature_dw']['As2S3'])}",
         'ablation': 'first_spectral_feature'},
        {'id': 'dp6', 'fork': 'Transfer model for the Bi2Te3 binding energy',
         'wrong_choice': '1D screened-wire model with the DFPT supercell polarizability (raw or calibrated on S3 and Te3), or the mean of the released chains',
         'consequence': f"E_b {r3(A['wire_model_eb'])}, {r3(A['wire_model_calibrated_eb'])} or {r3(A['released_mean_eb'])} eV "
                        f"against {r3(T['bi2te3_binding_energy_eV'])}",
         'ablation': 'wire_model_transfer'},
    ],
    'graded_fields': [
        f"scgw_direct_gap_eV (Te3, Bi2Te3; with and without SOC) within +-{TH['gap_abs_eV']} eV",
        f"scgw_optical_gap_eV (S3, Te3, As2S3) within +-{TH['optical_gap_abs_eV']} eV",
        f"screening_update_shift_eV (S3, Te3, As2S3) within +-{TH['screening_shift_abs_eV']} eV",
        f"bi2te3_binding_energy_eV in {TH['bi2te3_binding_window_eV']} eV (regime window)",
    ],
    'length_estimate': {
        'oracle_lab_evaluations': O.get('evaluations'), 'oracle_charge': O.get('charge_used'), 'budget': 45,
        'catalog_units': 208, 'tool_calls': '40-70', 'gpt6astra_cost_usd': '2-4',
        'why': 'Parsing Yambo headers and spectra, cross-checking QP energies between runs and against the SOC and SOC-free bands, '
               'peak finding with the brightness rule, and a transfer estimate for Bi2Te3 under a budget that covers about a fifth of the catalog.'},
    'runtime_dependencies': ['numpy, scipy (runtime image); no extra layer'],
    'data_provenance': {'archive_url': 'https://archive.materialscloud.org/records/42t0v-jam38',
                        'doi': '10.24435/materialscloud:xq-dy', 'paper_doi': '10.1021/acsnano.5c14061',
                        'files_used': ['1D_ACS_Nano_Exciton.zip', 'README.txt'],
                        'md5': {'1D_ACS_Nano_Exciton.zip': 'a3cdbff7bbca1465fa0d2a8ec75ccec4', 'README.txt': '781fc14a3db2c37c787b60dac65d45cf'},
                        'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': ('Adapted from the candidate design (see the construction README). The archived GW and BSE already include SOC, so the '
              'design\'s "add DFT SOC to GW" step became the double-counting trap dp3. Cross-chain model transfer of the Bi2Te3 gap and '
              'binding energy is not uniquely recoverable (the 1D wire model overestimates every chain, and Bi2Te3\'s self-consistency '
              'increment is 0.11 eV against 0.75-0.81 eV for the others). The Bi2Te3 gap is therefore graded tightly through evidence in '
              'Bi2Te3\'s own BSE_G0W0 run, and its binding energy through a regime window. Phonon data are not staged.'),
    'calibration': calibration,
}
if O:
    design['length_estimate']['oracle_errors_eV'] = O.get('errors')
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print(json.dumps(calibration, indent=1))
