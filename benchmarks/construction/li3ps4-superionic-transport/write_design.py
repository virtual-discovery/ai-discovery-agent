"""Write design.json for li3ps4-superionic-transport from the private targets and ablation values.

    cd benchmarks/construction/li3ps4-superionic-transport && python3 -B write_design.py
Run after make_task.py (and after an oracle run, whose graded errors are copied from oracle_errors.json if present).
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/li3ps4-superionic-transport'
H = json.loads((TASK/'environment/lab/private/hidden.json').read_text())
V = json.loads((TASK/'solution/ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
T, TH = H['targets'], H['thresholds']
ORACLE = json.loads((HERE/'oracle_errors.json').read_text()) if (HERE/'oracle_errors.json').exists() else {}


def rd(x, n=4):
    return float('%.*g' % (n, x)) if isinstance(x, float) else x


def dex(d):
    return {k: rd(math.log10(v / T['sigma_S_per_m'][k]), 3) for k, v in d.items()}


design = {
    'task_id': 'li3ps4-superionic-transport',
    'catalog_id': 'li3ps4-superionic-transport',
    'start': 'released',
    'source_group': 'gigli_2024_li3ps4_transport',
    'title': 'Li3PS4: superionic transition, ionic conductivity and the paddle-wheel question',
    'decision_points': [
        {'id': 'dp1', 'fork': 'Which conductivity answers the question: Nernst-Einstein (tracer) or Green-Kubo (charge current)',
         'wrong_choice': 'Use the cheaper, smoother Nernst-Einstein series as sigma (H_R = 1)',
         'consequence': 'Predictions low by the Haven ratio (-0.26 to -0.29 dex); superionic E_a 0.28 eV instead of ~0.18 eV '
                        '(H_R rises with T); H_R = 1',
         'ablation': 'nernst_einstein_as_conductivity', 'difficulty_type': 'convention, with a cascade into E_a'},
        {'id': 'dp2', 'fork': 'Which Green-Kubo points below the transition are resolved',
         'wrong_choice': 'Fit every positive Green-Kubo value below the transition',
         'consequence': 'Below ~600 K the gamma-started Green-Kubo values sit at a ~0.03 S/m noise floor (apparent H_R > 1); '
                        'including them gives E_a ~0.16 eV instead of ~1.6 eV',
         'ablation': 'green_kubo_noise_floor_in_low_T_fit', 'difficulty_type': 'regime judgement discovered in the data'},
        {'id': 'dp3', 'fork': 'Which heat-capacity anomaly is the superionic transition',
         'wrong_choice': 'Take the largest Cp peak',
         'consequence': 'For every reference the largest Cp peak is melting of the PS4 sublattice (coincides with P diffusion '
                        'and free rotation): Tc 875/925 K instead of 737.5/787.5 K, and the withheld runs are then predicted '
                        'from the low-temperature branch',
         'ablation': 'largest_heat_capacity_peak_as_transition', 'difficulty_type': 'confounder that mimics the target; discriminating purchase (D_P, angle step)'},
        {'id': 'dp4', 'fork': 'Whether the transition depends on the DFT reference',
         'wrong_choice': "Transfer one reference's transition (PBEsol) or the experimental alpha transition (746 K) to all",
         'consequence': 'Tc off by 125-175 K (PBEsol transfer) or 41.5 K for r2SCAN (experimental value); budget must be split over references',
         'ablation': 'pbesol_transition_for_all_references', 'difficulty_type': 'discovered in the data under a binding budget'},
        {'id': 'dp5', 'fork': 'On which side of the jump each withheld run lies',
         'wrong_choice': 'Interpolate log sigma in 1/T between the neighbouring served temperatures',
         'consequence': 'Each withheld run is the first temperature above its transition; interpolation across the jump is '
                        'low by 0.27-0.52 dex. Resolving it requires buying a discriminating observable at that temperature '
                        '(angle fluctuation, Cp of the interval, order parameters or the 750 K MSD)',
         'ablation': 'interpolate_across_the_jump', 'difficulty_type': 'discriminating experiment; cascade from the transition location'},
        {'id': 'dp6', 'fork': 'Whether runs from different starting phases reach the same superionic state',
         'wrong_choice': 'Use alpha-started runs as the superionic reference for the gamma-started cells',
         'consequence': 'Alpha-started runs stay in a different, more conductive state: predictions high by 0.25-0.31 dex and '
                        'H_R 0.37 instead of 0.53',
         'ablation': 'alpha_started_runs_as_superionic_reference', 'difficulty_type': 'hysteresis discovered in the data'},
        {'id': 'dp7', 'fork': 'Frame for attributing conduction to Li-anion coupling',
         'wrong_choice': 'Read the Li-S term in the barycentric frame as dynamical coupling',
         'consequence': 'Li-S fraction 0.27 (equal to the kinematic recoil of the matrix, fixed by the masses) instead of ~0.05; '
                        'a spurious argument for the paddle-wheel picture',
         'ablation': 'barycentric_li_s_term_as_coupling', 'difficulty_type': 'confounder that mimics the target effect'},
        {'id': 'dp8', 'fork': 'Which PS4 motion to compare with Li hops',
         'wrong_choice': 'Compare with the PS4 libration rate implied by the angle fluctuation',
         'consequence': 'Rate ratio ~10^-0.3 (hops no faster than PS4 motion) instead of ~10^2.7; paddle-wheel wrongly supported',
         'ablation': 'librations_taken_for_reorientations', 'difficulty_type': 'discriminating observable between competing explanations'},
        {'id': 'dp9', 'fork': 'One Arrhenius law or separate phases',
         'wrong_choice': 'Fit one Arrhenius law through the transition for the gamma-started runs',
         'consequence': 'Single E_a ~0.65 eV for both regimes (true ~1.6 and ~0.18 eV); transition and superionic barrier missed',
         'ablation': 'single_arrhenius_through_the_transition', 'difficulty_type': 'model class'},
    ],
    'graded_fields': [
        f"transition_temperature_K (PBE0, r2SCAN) within +/-{TH['transition_K']:g} K",
        f"sigma_predictions_S_per_m: each of 3 withheld runs within +/-{TH['sigma_log10']:g} dex",
        f"haven_ratio_PBE0_gamma_750K within +/-{TH['haven_ratio']:g}",
        f"activation_energy_eV_PBE0_gamma: below_transition +/-{TH['ea_below_eV']:g} eV, superionic_to_1000K +/-{TH['ea_above_eV']:g} eV",
        f"li_s_cross_fraction within +/-{TH['li_s_fraction']:g}",
        f"log10_hop_to_flip_rate_ratio_PBE0_800K within +/-{TH['log10_rate_ratio']:g} dex",
        'paddle_wheel_supported must be false',
    ],
    'length_estimate': {'oracle_lab_evaluations': ORACLE.get('evaluations', 40), 'oracle_charge': ORACLE.get('charge', 50),
                        'budget': 60, 'tool_calls': '50-80', 'gpt6astra_cost_usd': '3-5', 'pilot': 'pending (coordinator)'},
    'runtime_dependencies': ['numpy (runtime image); the lab is a lookup of archived values'],
    'data_provenance': {'archive_url': 'https://archive.materialscloud.org/records/pkn2d-h8190', 'doi': '10.24435/materialscloud:qy-gv',
                        'paper_doi': '10.1021/acs.chemmater.3c02726', 'files_used': ['MaterialsCloudArchive-ChemMat2024.zip (47 members, see SOURCE.md)'],
                        'md5': {'MaterialsCloudArchive-ChemMat2024.zip': 'e838e8716fb42200c34f06119337e909'}, 'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': ('Adapted from the design: the archive holds no trajectories, only figure-level analyses, so the lab serves them '
              'per run and the withheld runs are the first temperature above each transition. PBE0 is graded with Green-Kubo '
              'activation energies (the paper tabulates Nernst-Einstein ones, 1.43/0.269 eV). The hop/flip ratio is graded at 800 K '
              'in the superionic phase (600 K, as in the design, lies below the PBE0 transition, where D_omega is librational). '
              'Beta-started runs above 625 K are not served because the archived r2SCAN and PBE0 beta series coincide there.'),
    'calibration': {
        'transition_temperature_K': {'reference': T['transition_temperature_K'],
                                     'robustness': {'upper_bracket': R['tc_upper_bracket'], 'lower_bracket': R['tc_lower_bracket']},
                                     'oracle': ORACLE.get('transition_K'),
                                     'nearest_ablation': {'experimental_746K (r2SCAN)': 746.0, 'melting_peak': A['tc_melting_peak'],
                                                          'pbesol_transfer': A['tc_pbesol_for_all']},
                                     'tolerance': f"+/-{TH['transition_K']:g} K"},
        'sigma_predictions_log10_error': {'reference_S_per_m': {k: rd(v) for k, v in T['sigma_S_per_m'].items()},
                                          'robustness': {'tracer_branch_over_HR': dex(R['sigma_tracer_branch_over_haven']),
                                                         'gk_superionic_fit_to_1000K': dex(R['sigma_superionic_gk_fit']),
                                                         'msd_route_PBE0_750': rd(math.log10(R['sigma_PBE0_750_msd_route'] / T['sigma_S_per_m']['PBE0/gamma/750']), 3)},
                                          'oracle': ORACLE.get('sigma_log10'),
                                          'nearest_ablation': {'alpha_branch': dex(A['sigma_alpha_branch']), 'tracer_as_sigma': dex(A['sigma_tracer']),
                                                               'interpolate_across_jump': dex(A['sigma_interpolated_across_jump']),
                                                               'low_branch': dex(A['sigma_low_branch'])},
                                          'tolerance': f"+/-{TH['sigma_log10']:g} dex each"},
        'haven_ratio_PBE0_gamma_750K': {'reference': rd(T['haven_ratio_PBE0_gamma_750K']),
                                        'robustness': {'first_three_superionic': rd(R['haven_first_three']), 'all_superionic_to_1100K': rd(R['haven_all_superionic_to_1100K'])},
                                        'oracle_error': ORACLE.get('haven_ratio'),
                                        'nearest_ablation': {'alpha_branch': rd(A['haven_alpha_branch']), 'tracer_equals_sigma': 1.0},
                                        'tolerance': f"+/-{TH['haven_ratio']:g}"},
        'activation_energy_below_transition_eV': {'reference': rd(T['activation_energy_eV']['below_transition']),
                                                  'robustness': {'unweighted': rd(R['ea_below_unweighted']), 'tracer': rd(R['ea_below_tracer'])},
                                                  'oracle_error': ORACLE.get('ea_below_eV'),
                                                  'nearest_ablation': {'single_arrhenius': rd(A['ea_single_arrhenius']), 'noise_floor': rd(A['ea_noise_floor']['below_transition'])},
                                                  'tolerance': f"+/-{TH['ea_below_eV']:g} eV"},
        'activation_energy_superionic_eV': {'reference': rd(T['activation_energy_eV']['superionic_to_1000K']),
                                            'robustness': {'window_to_875K': rd(R['ea_above_to_875K']), 'without_750K': rd(R['ea_above_without_750K'])},
                                            'oracle_error': ORACLE.get('ea_above_eV'),
                                            'nearest_ablation': {'tracer': rd(A['ea_tracer']['superionic_to_1000K']), 'single_arrhenius': rd(A['ea_single_arrhenius'])},
                                            'tolerance': f"+/-{TH['ea_above_eV']:g} eV"},
        'li_s_cross_fraction': {'reference': rd(T['li_s_cross_fraction']), 'robustness': {'zero_by_construction': 0.0},
                                'nearest_ablation': {'barycentric': rd(A['li_s_barycentric']),
                                                     'recoil_prediction': rd(H['reference_details']['barycentric_recoil_prediction'])},
                                'tolerance': f"+/-{TH['li_s_fraction']:g}"},
        'log10_hop_to_flip_rate_ratio_PBE0_800K': {'reference': rd(T['log10_hop_to_flip_rate_ratio_PBE0_800K']),
                                                   'robustness': {'domega_as_rate': rd(R['log_ratio_domega_as_rate']), 'msd_750K_hop_3A': rd(R['log_ratio_msd_750_hop3A'])},
                                                   'nearest_ablation': {'libration': rd(A['log_ratio_libration'])},
                                                   'tolerance': f"+/-{TH['log10_rate_ratio']:g} dex"},
    },
}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print('wrote design.json')
