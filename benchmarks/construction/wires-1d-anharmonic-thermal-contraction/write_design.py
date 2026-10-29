"""Write design.json (decision points, calibration table, length estimate) from calibration.json and check.json."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/wires-1d-anharmonic-thermal-contraction'
cal = json.loads((HERE/'calibration.json').read_text())
check = json.loads((HERE/'check.json').read_text()) if (HERE/'check.json').exists() else {}
ref, rob, abl, info = cal['ref'], cal['rob'], cal['abl'], cal['info']
th = json.loads((TASK/'environment/lab/private/hidden.json').read_text())['thresholds']
r = lambda x, n=4: float(f'{x:.{n}g}')
pct = lambda a, b: f'{100 * (a / b - 1):+.1f}%'

DPS = [
    ('dp1', 'Rigid-body modes of a 1D wire: four zero modes at Gamma (3 translations + axial rotation); the archived Hessians carry a small imaginary rotational value at Gamma at every T',
     'Treat the imaginary Gamma value as a dynamical instability', 'TaSe3 called unstable from 0 K', 'gamma_rotation_mode_counted_as_instability'),
    ('dp2', 'Auxiliary SSCHA matrices vs free-energy Hessian', 'Judge stability from auxiliary frequencies (always positive) or compute C_V from Hessian frequencies',
     f'No instability found up to 1000 K; C_V {pct(abl["cv_cuc2_hessian_freq"], ref["cv_cuc2"])} (CuC2), {pct(abl["cv_tase3_hessian_freq"], ref["cv_tase3"])} (TaSe3)',
     'stability_from_auxiliary_frequencies'),
    ('dp3', 'Sign convention of the archived SSCHA stress (pressure-like, as pw.x); established from the kinetic SSCHA-minus-average difference or the pw.x P line',
     'Read the printed stress as tension-positive', 'alpha_L positive: expansion predicted instead of contraction', 'stress_read_as_tension_positive'),
    ('dp4', 'SSCHA free-energy stress vs ab initio ensemble-average stress (both printed in each log)', 'Differentiate the ab initio average stress',
     f'alpha_L {pct(abl["alpha_average_stress"], ref["alpha"])}', 'ab_initio_average_stress'),
    ('dp5', 'Acoustic sum rule in the long-wavelength stiffness: plain ph.x matrices vs the 1D-ASR set', 'Fourier-interpolate the plain DFPT matrices to small q',
     f'k {pct(abl["k_plain_asr"], ref["k"])} at q = 0.01 (2 pi/c) and divergent as q -> 0; alpha_L {pct(abl["alpha_plain_asr"], ref["alpha"])}', 'plain_asr_long_wavelength_stiffness'),
    ('dp6', 'Harmonic (infinitesimal) stiffness vs the single archived strained calculation (+4.7 %)', 'Use the finite-strain energy secant 2 dE / eps^2',
     f'k {pct(abl["k_finite_strain"], ref["k"])}; alpha_L {pct(abl["alpha_finite_strain"], ref["alpha"])}', 'finite_strain_secant_stiffness'),
    ('dp7', 'Normalization of a 1D stress: cell volume incl. vacuum consistent with the per-cell stiffness; k = d2E/deps2 (not E = k eps^2)',
     'Use E = k eps^2 (stiffness halved, as the paper appears to) or the supercell volume', 'alpha_L x2 or x6', 'stiffness_from_E_equals_k_eps2'),
    ('dp8', 'Sign of the frequency-derivative term in C_V = sum c (1 - T/w dw/dT) (the archived notebook and the paper write +)', 'Use the + sign',
     f'C_V {pct(abl["cv_cuc2_paper_sign"], ref["cv_cuc2"])} (CuC2), {pct(abl["cv_tase3_paper_sign"], ref["cv_tase3"])} (TaSe3): the deviation from harmonic has the wrong sign',
     'derivative_term_with_plus_sign'),
    ('dp9', 'C_V route: entropy derivative vs harmonic C_V with temperature-dependent (or DFPT) frequencies', 'Harmonic formula with SSCHA or DFPT frequencies',
     f'C_V {pct(abl["cv_cuc2_harmonic_sscha_freq"], ref["cv_cuc2"])} (CuC2), {pct(abl["cv_tase3_harmonic_sscha_freq"], ref["cv_tase3"])} (TaSe3)',
     'harmonic_cv_with_sscha_frequencies'),
    ('dp10', 'Withheld 650/700 K Hessians: the onset must be interpolated between released temperatures', 'Report the first released unstable temperature (800 K)',
     f'{abl["T_inst_first_released_unstable"] - ref["T_inst"]:+.0f} K', 'first_released_unstable_temperature'),
]

design = {
    'task_id': TASK.name, 'catalog_id': TASK.name, 'start': 'released', 'source_group': 'cignarella_2025_1d_wires',
    'title': '1D wires: thermal contraction, anharmonic heat capacity and instability from SSCHA data',
    'decision_points': [{'id': i, 'fork': f, 'wrong_choice': w, 'consequence': c, 'ablation': a} for i, f, w, c, a in DPS],
    'graded_fields': [
        f'TaSe3.axial_stiffness_eV within {th["k_rel"]:.0%} of {r(ref["k"])}',
        f'TaSe3.alpha_L_300K_1e-6_per_K within {th["alpha_rel"]:.0%} of {r(ref["alpha"])}',
        f'CuC2.cv_ratio_300K within {th["cv_cuc2_rel"]:.1%} of {r(ref["cv_cuc2"])}',
        f'TaSe3.cv_ratio_300K within {th["cv_tase3_rel"]:.1%} of {r(ref["cv_tase3"])}',
        f'TaSe3.instability_temperature_K within {th["T_inst_abs_K"]:.0f} K of {r(ref["T_inst"])}',
    ],
    'length_estimate': {
        'oracle_lab_evaluations': (check.get('summary', {}).get('oracle', {}) or {}).get('evaluations'),
        'budget': 36, 'catalog_units': 137, 'tool_calls': '45-75', 'gpt6astra_cost_usd': '2.5-4.5',
        'why': 'Three wires, two kinds of SSCHA matrices plus logs per temperature, DFPT with two sum-rule variants and one strained '
               'calculation; the agent parses QE dynamical matrices, Fourier-interpolates force constants, differentiates noisy '
               'stresses and mode frequencies, and must pick converged runs within a budget of about a quarter of the catalog.'},
    'runtime_dependencies': ['numpy (runtime image); lab serves archived text only'],
    'data_provenance': {'archive_url': 'https://archive.materialscloud.org/records/k45dq-ebf35', 'doi': '10.24435/materialscloud:fj-20',
                        'files_used': ['Extreme-anharmonicity-thermal-contration-1D-wires.zip', 'README.txt'],
                        'md5': {'Extreme-anharmonicity-thermal-contration-1D-wires.zip': '1de3fb44caf32b1ed403cbe0e7b687ce',
                                'README.txt': 'c866d88afdf77c93202d26d1c30e0b38'},
                        'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': ('Graded values are computed from archived ingredients with a validated numpy pipeline (calib.py), as in the '
              'skutterudite task; the archive holds no alpha or C_V outputs. The TaSe3 instability is bracketed by withheld '
              'Hessians (650/700 K). QHA alpha (design target) dropped: the QHA AiiDA archive holds one strained TaSe3 '
              'calculation only. AuSe2 is served but not graded (no Hessians; runs at T >= 150 K end at the Kong-Liu limit).'),
    'calibration': {
        'TaSe3.axial_stiffness_eV': {'reference': r(ref['k']), 'robustness': {'q_0.05': r(rob['k_q005']), 'T0_sscha_hessian': r(rob['k_T0_hessian'])},
                                     'nearest_ablation': {'finite_strain': r(abl['k_finite_strain']), 'plain_asr_q0.01': r(abl['k_plain_asr'])},
                                     'other': {'T300_sscha_hessian (isothermal, excluded by the stated definition)': r(info['k_T300_hessian'])},
                                     'tolerance': f'+-{th["k_rel"]:.0%}'},
        'TaSe3.alpha_L_300K_1e-6_per_K': {'reference': r(ref['alpha']),
                                          'robustness': {k: r(rob[k]) for k in ('alpha_wide_quadratic', 'alpha_global_cubic_T0_hessian_k', 'alpha_k_q005')},
                                          'nearest_ablation': {k: r(abl[k]) for k in ('alpha_finite_strain', 'alpha_plain_asr', 'alpha_average_stress', 'alpha_half_k', 'alpha_sign')},
                                          'other': {'stiffness from T300 Hessian': r(info['alpha_T300_hessian_k'])},
                                          'tolerance': f'+-{th["alpha_rel"]:.0%}'},
        'CuC2.cv_ratio_300K': {'reference': r(ref['cv_cuc2']),
                               'robustness': {k: r(rob[k]) for k in ('cv_cuc2_cubic', 'cv_cuc2_linear_200_500', 'cv_cuc2_all_converged_cubic')},
                               'nearest_ablation': {k: r(abl[k]) for k in ('cv_cuc2_harmonic_sscha_freq', 'cv_cuc2_dfpt_harmonic', 'cv_cuc2_per_atom_3N', 'cv_cuc2_paper_sign', 'cv_cuc2_hessian_freq', 'cv_cuc2_with_unconverged_400K')},
                               'other': {'unconverged 400 K included, quadratic (passes)': r(info['cv_cuc2_with_unconverged_400K_quadratic'])},
                               'tolerance': f'+-{th["cv_cuc2_rel"]:.1%}'},
        'TaSe3.cv_ratio_300K': {'reference': r(ref['cv_tase3']),
                                'robustness': {k: r(rob[k]) for k in ('cv_tase3_cubic', 'cv_tase3_linear_200_500', 'cv_tase3_all_converged_cubic')},
                                'nearest_ablation': {k: r(abl[k]) for k in ('cv_tase3_harmonic_sscha_freq', 'cv_tase3_dfpt_harmonic', 'cv_tase3_paper_sign', 'cv_tase3_per_atom_3N', 'cv_tase3_hessian_freq')},
                                'tolerance': f'+-{th["cv_tase3_rel"]:.1%}'},
        'TaSe3.instability_temperature_K': {'reference': r(ref['T_inst']),
                                            'reference_basis': {'100k_650_700_crossing': r(info['T_inst_100k']), '10k_650_700_crossing': r(info['T_inst_10k'])},
                                            'robustness': {k: r(rob[k]) for k in ('T_inst_released_600_800_10k', 'T_inst_released_600_800_100k', 'T_inst_quadratic_fit')},
                                            'nearest_ablation': {'first_released_unstable': 800, 'gamma_rotation': 0, 'auxiliary': 'none'},
                                            'tolerance': f'+-{th["T_inst_abs_K"]:.0f} K'},
    },
}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print('design.json written;', len(DPS), 'decision points')
