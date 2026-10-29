"""Write design.json for bab2-hole-doping-two-gap from calibration.json (make_task.py) and oracle_run.json (record_oracle.py).

    cd benchmarks/construction/bab2-hole-doping-two-gap && python3 -B write_design.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/bab2-hole-doping-two-gap'
C = json.loads((HERE/'calibration.json').read_text())
ORACLE = json.loads((HERE/'oracle_run.json').read_text()) if (HERE/'oracle_run.json').exists() else {}
ref, rob, abl, th, info = C['references'], C['robustness'], C['ablations'], C['thresholds'], C['info']
r3 = lambda v, n=3: round(v, n)

DPS = [
    {'id': 'dp1', 'fork': 'Which doping each archive folder holds: README order versus each run\'s tot_charge',
     'wrong_choice': 'Map folders 0-6 onto the README list (0.0, 0.02, 0.04, 0.06, 0.065, 0.08, 0.1) in order',
     'consequence': f"x_c {r3(abl['xc_readme_order'], 4)} instead of {r3(ref['critical_hole_doping'], 4)}; the 'x = 0.065' and 'x = 0.10' Tc "
                    f"come from x = 0.10 and x = 0 ({r3(abl['iso_readme_order']['0.065'], 2)} and {r3(abl['iso_readme_order']['0.10'], 2)} K)",
     'ablation': 'readme_order_mapping'},
    {'id': 'dp2', 'fork': 'How to locate x_c when the bracketing runs (0.065, 0.08) are withheld: interpolate omega or omega^2',
     'wrong_choice': 'Linear interpolation of the signed soft-mode frequency between x = 0.06 and 0.10',
     'consequence': f"x_c {r3(abl['xc_omega_linear'], 4)} (outside +-{th['critical_doping_abs']}); homogeneous 1/9 soft mode "
                    f"{r3(abl['hom_omega_linear'], 1)} instead of {r3(ref['soft_mode_cm1']['homogeneous_one_ninth'], 1)} cm-1",
     'ablation': 'omega_linear_crossing'},
    {'id': 'dp3', 'fork': 'What the archived supercells contain: the README calls 331 homogeneously 1-hole doped, its output shows 144 electrons '
                          '(undoped); the primitive H point folds onto supercell A, not Gamma',
     'wrong_choice': 'Read the 331 folded-H mode as the homogeneous 1/9 h/f.u. answer, or read the Cs cell at Gamma',
     'consequence': f"homogeneous {r3(abl['hom_331_as_homogeneous'], 2)} cm-1 (unstable) instead of +{r3(ref['soft_mode_cm1']['homogeneous_one_ninth'], 1)}; "
                    f"Cs {r3(abl['cs_gamma_point'], 2)} instead of {r3(ref['soft_mode_cm1']['cs_substituted'], 2)} cm-1",
     'ablation': 'supercell_331_as_homogeneous'},
    {'id': 'dp4', 'fork': 'Which isotropic calculation answers: the archived mu* = 0.13 logs and their Allen-Dynes header estimate, or new '
                          'Eliashberg solutions at the stated mu* = 0.10',
     'wrong_choice': 'Extrapolate the archived mu* = 0.13 gaps (15-30 K) to zero, or use Allen-Dynes with f1 f2',
     'consequence': f"Tc {r3(abl['iso_archived_mu013']['0.065'], 2)}/{r3(abl['iso_archived_mu013']['0.10'], 2)} K (archived logs) or "
                    f"{r3(abl['iso_allen_dynes']['0.065'], 2)}/{r3(abl['iso_allen_dynes']['0.10'], 2)} K (Allen-Dynes) instead of "
                    f"{r3(ref['iso_tc_K']['0.065'], 2)}/{r3(ref['iso_tc_K']['0.10'], 2)} K",
     'ablation': 'archived_mu013_logs'},
    {'id': 'dp5', 'fork': 'Convergence near Tc: critical slowing down leaves finite, unconverged gaps above Tc (as in the archived EPW log at '
                          'x = 0.08, which ends in NaN at 30 K)',
     'wrong_choice': 'Scan T with the default iteration limit and take the first T where the reported gap drops below 0.01 meV',
     'consequence': f"Tc {abl['iso_unconverged_scan']['0.065']:g}/{abl['iso_unconverged_scan']['0.10']:g} K instead of "
                    f"{r3(ref['iso_tc_K']['0.065'], 2)}/{r3(ref['iso_tc_K']['0.10'], 2)} K",
     'ablation': 'unconverged_gap_scan'},
    {'id': 'dp6', 'fork': 'Anisotropic Tc when the archived gaps are still open at the last temperature (60 K)',
     'wrong_choice': 'Report the isotropic Tc as the superconducting Tc, or extrapolate Delta (not Delta^2) linearly',
     'consequence': f"{r3(abl['aniso_iso_reported'], 1)} K or {r3(abl['aniso_linear_gap'], 1)} K instead of {r3(ref['aniso_tc_K'], 1)} K",
     'ablation': 'linear_gap_extrapolation'},
    {'id': 'dp7', 'fork': 'Gap structure from the imaginary-axis distribution or from the Pade real-axis continuation',
     'wrong_choice': 'Take the cluster gaps from the real-axis (pade_aniso_gap0) file',
     'consequence': f"ratio {r3(abl['ratio_real_axis'], 2)} instead of {r3(ref['gap_ratio_25K'], 2)} (+-{int(th['gap_ratio_rel'] * 100)}%)",
     'ablation': 'real_axis_gap_ratio'},
]

CALIB = {
    'critical_hole_doping': {
        'reference': r3(ref['critical_hole_doping'], 5),
        'definition': 'zero of omega^2 at H, piecewise linear between the archived x = 0.06 and 0.065 runs (all archived DFPT)',
        'oracle': ORACLE.get('submission', {}).get('critical_hole_doping'),
        'robustness': {k: r3(v, 5) for k, v in rob.items() if k.startswith('xc_')},
        'nearest_ablation': {'omega_linear': r3(abl['xc_omega_linear'], 5), 'readme_order': r3(abl['xc_readme_order'], 5)},
        'not_ablations': {'paper_value_0.065': 0.065},
        'tolerance': f"+-{th['critical_doping_abs']}"},
    'iso_tc_K': {
        'reference': {k: r3(v, 2) for k, v in ref['iso_tc_K'].items()},
        'definition': 'leading eigenvalue of the linearized equations of the lab solver = 1 (mu* 0.10, 500 meV, a2f column 1)',
        'oracle': ORACLE.get('submission', {}).get('iso_tc_K'),
        'robustness': {'converged_delta2_extrapolation': {k: r3(v, 2) for k, v in rob['iso_converged_extrapolation'].items()},
                       'a2f_column_10': {k: r3(v, 2) for k, v in rob['iso_a2f_column_10'].items()}},
        'nearest_ablation': {'archived_mu013_logs': {k: r3(v, 2) for k, v in abl['iso_archived_mu013'].items()},
                             'unconverged_scan': abl['iso_unconverged_scan'],
                             'allen_dynes': {k: r3(v, 2) for k, v in abl['iso_allen_dynes'].items()},
                             'readme_order': {k: r3(v, 2) for k, v in abl['iso_readme_order'].items()}},
        'not_ablations': {'paper_Tc_iso_mu013': 32.0,
                          'linearized_mu013': {k: r3(v, 2) for k, v in info['iso_mu013_linearized'].items()}},
        'tolerance': f"+-{int(th['iso_tc_rel'] * 100)}% each"},
    'aniso_tc_K': {
        'reference': r3(ref['aniso_tc_K'], 2),
        'definition': 'Delta_max(imaginary axis)^2 extrapolated linearly to zero from the archived 55 and 60 K (x = 0.065)',
        'oracle': ORACLE.get('submission', {}).get('aniso_tc_K'),
        'robustness': {k: r3(v, 2) for k, v in rob.items() if k.startswith('aniso_')},
        'nearest_ablation': {'linear_gap_extrapolation': r3(abl['aniso_linear_gap'], 2), 'isotropic_reported': r3(abl['aniso_iso_reported'], 2)},
        'not_ablations': {'last_archived_temperature': 60.0, 'paper_value': 61.6},
        'tolerance': f"+-{int(th['aniso_tc_rel'] * 100)}%"},
    'gap_ratio_25K': {
        'reference': r3(ref['gap_ratio_25K'], 3),
        'definition': 'weight-averaged large-gap cluster / small-gap cluster, archived imaginary-axis distribution, x = 0.065, 25 K',
        'oracle': ORACLE.get('submission', {}).get('gap_ratio_25K'),
        'robustness': {'median': r3(rob['ratio_median'], 3), 'unweighted_bins': r3(rob['ratio_unweighted_bins'], 3)},
        'nearest_ablation': {'real_axis': r3(abl['ratio_real_axis'], 3)},
        'not_ablations': {'mode_of_each_cluster (violates the stated weighting)': r3(info['ratio_mode_convention_violation'], 3),
                          'paper_gaps_11.54_over_1.50': r3(info['paper_values']['gap_ratio'], 3)},
        'tolerance': f"+-{int(th['gap_ratio_rel'] * 100)}%"},
    'soft_mode_cm1.homogeneous_one_ninth': {
        'reference': r3(ref['soft_mode_cm1']['homogeneous_one_ninth'], 2),
        'definition': 'omega^2 at H extrapolated linearly from the archived x = 0.08 and 0.10 runs to x = 1/9',
        'oracle': ORACLE.get('submission', {}).get('soft_mode_cm1', {}).get('homogeneous_one_ninth'),
        'robustness': {k: r3(v, 2) for k, v in rob.items() if k.startswith('hom_')},
        'nearest_ablation': {'omega_linear': r3(abl['hom_omega_linear'], 2), '331_as_homogeneous': r3(abl['hom_331_as_homogeneous'], 2)},
        'not_ablations': {'x_0.10_value_unextrapolated': r3(info['soft_modes']['0.1'], 2)},
        'tolerance': f"+-{th['homogeneous_soft_mode_abs_cm1']} cm-1"},
    'soft_mode_cm1.cs_substituted': {
        'reference': r3(ref['soft_mode_cm1']['cs_substituted'], 2),
        'definition': 'lowest archived 331_cs frequency at the supercell A point (= folded primitive H)',
        'oracle': ORACLE.get('submission', {}).get('soft_mode_cm1', {}).get('cs_substituted'),
        'nearest_ablation': {'gamma_point': r3(abl['cs_gamma_point'], 2)},
        'tolerance': f"+-{th['cs_soft_mode_abs_cm1']} cm-1"},
}

design = {
    'task_id': 'bab2-hole-doping-two-gap', 'catalog_id': 'bab2-hole-doping-two-gap', 'start': 'released',
    'source_group': 'bercx_2025_bcs_landscape', 'title': 'Hole-doped BaB2: stability threshold, Eliashberg Tc and two gaps',
    'decision_points': DPS,
    'graded_fields': [f"critical_hole_doping within +-{th['critical_doping_abs']} h/f.u.",
                      f"iso_tc_K (x = 0.065, 0.10; mu* 0.10, 500 meV) within +-{int(th['iso_tc_rel'] * 100)}%",
                      f"aniso_tc_K (x = 0.065, archived anisotropic run) within +-{int(th['aniso_tc_rel'] * 100)}%",
                      f"gap_ratio_25K (x = 0.065) within +-{int(th['gap_ratio_rel'] * 100)}%",
                      f"soft_mode_cm1.homogeneous_one_ninth within +-{th['homogeneous_soft_mode_abs_cm1']} cm-1",
                      f"soft_mode_cm1.cs_substituted within +-{th['cs_soft_mode_abs_cm1']} cm-1"],
    'length_estimate': {'oracle_lab_evaluations': ORACLE.get('evaluations'), 'oracle_charge': ORACLE.get('charge_used'),
                        'budget': ORACLE.get('budget', 34), 'catalog_units': 'about 300 (DFPT 10, a2f 7, logs 7, 136 gap files x 2 units, supercells 4) plus any number of solver runs',
                        'tool_calls': '45-75', 'gpt6astra_cost_usd': '2.5-4.5',
                        'why': ('Folder-to-doping mapping from headers, DFPT soft-mode scan and interpolation, solver validation and '
                                'Tc scans with convergence control at two dopings, budgeted anisotropic gap files with an '
                                'extrapolation, and a supercell cross-check; each wrong fork moves a graded value outside its tolerance.'),
                        'oracle_errors': ORACLE.get('errors')},
    'runtime_dependencies': ['numpy (runtime image); no extra layer'],
    'data_provenance': {'archive_url': 'https://archive.materialscloud.org/records/0kqen-w8h90', 'doi': '10.24435/materialscloud:x2-ap',
                        'paper_doi': '10.1103/sb28-fjc9', 'files_used': ['bab2-analysis.zip', 'README.md', 'files_description.md'],
                        'md5': {'bab2-analysis.zip': 'a63a59c6c913b4cab7d5050eb0f1710c', 'README.md': '57c5c2ef290197ab4035051e52fcc5d1',
                                'files_description.md': '405192de6ce3f5f9b25d2f376ba92694'},
                        'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': ('Adapted from the candidate design (see the construction README). Hidden dopings keep their archived alpha2F and '
              'anisotropic data purchasable; only their DFPT is withheld, because Tc or gaps transferred to a withheld doping are '
              'not uniquely recoverable (Tc and gap ratios are non-monotonic in x near threshold). The isotropic target uses '
              'mu* = 0.10 at the archived 500 meV cutoff so that it requires new solutions, and is graded at x = 0.065 and 0.10. '
              'The memorised paper numbers x_c = 0.065 and Tc_aniso = 61.6 K fall inside their tolerances; Tc_iso = 32.0 K and the '
              'paper gap ratio 7.7 do not.'),
    'calibration': CALIB,
}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print(json.dumps({'decision_points': len(DPS), 'oracle': ORACLE.get('errors')}, indent=1))
