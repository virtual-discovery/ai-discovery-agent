"""Write the task's design.json from calibration.json (make_task.py) and oracle_run.json (record_oracle.py).

    cd benchmarks/construction/metal-resistivity-te-vs-anharmonicity && python3 -B write_design.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/metal-resistivity-te-vs-anharmonicity'
cal = json.loads((HERE/'calibration.json').read_text())
oracle = json.loads((HERE/'oracle_run.json').read_text()) if (HERE/'oracle_run.json').exists() else {}
t, A, R, th = cal['truth'], cal['ablations'], cal['robustness'], cal['thresholds']
r4 = lambda x: round(x, 4) if isinstance(x, float) else {k: r4(v) for k, v in x.items()}
pct = lambda v, ref: round(100 * (v / ref - 1), 2)

DPS = [
    ('dp1', 'Which harmonic reference isolates thermal expansion for Al: the archive holds 1-DFPT runs at 50^3, 60^3 '
            'and 70^3 fine meshes, the TE runs use 50^3 (visible only in each run\'s metadata)',
     'Use the finest (70^3, "converged") harmonic run', 'Al Delta_TE(481 K) %.2f%% instead of %.2f%% (+2 points of pure '
     'mesh error, outside +-1)' % (A['al_reference_mesh70'], t['al_te_change_percent_481K']), 'al_reference_finest_mesh70'),
    ('dp2', 'Combined TE + anharmonic resistivity from the combined (TDEP at the expanded lattice) physics, or by '
            'combining separately estimated TE and 0 K-lattice anharmonic changes',
     'Multiply (or add) the separate effects', 'Pb rho(510 K) %.1f (multiplicative, %+.1f%%) or %.1f (additive, %+.1f%%) '
     'instead of %.2f' % (A['te_anh_510_separable_multiplicative'], pct(A['te_anh_510_separable_multiplicative'], t['pb_te_anharm_resistivity_510K']),
                          A['te_anh_510_separable_additive'], pct(A['te_anh_510_separable_additive'], t['pb_te_anharm_resistivity_510K']),
                          t['pb_te_anharm_resistivity_510K']), 'separable_multiplicative'),
    ('dp3', 'Whether the anharmonic correction is a property of T alone: the 0 K-lattice TDEP runs are withheld, the '
            'expanded-lattice TE+anh/TE pairs are released',
     'Use the expanded-lattice anharmonic correction for the 0 K lattice', 'Pb Delta_anh %.1f / %.1f%% instead of '
     '%.1f / %.1f%% at 410 / 610 K' % (A['anharm_from_expanded_lattice']['410'], A['anharm_from_expanded_lattice']['610'],
                                     t['pb_anharmonic_change_percent_0K_lattice']['410'], t['pb_anharmonic_change_percent_0K_lattice']['610']),
     'anharmonic_correction_from_expanded_lattice'),
    ('dp4', 'Which phonon surrogate transfers EPW resistivity ratios: the lab prints the lambda-like <omega^-2> moment; '
            'a transport-weighted kernel (|q|^2 (x/sinh x)^2 / omega^2) must be built from the grid and validated on '
            'released same-lattice pairs',
     'Use the printed <omega^-2> moment ratio', 'Pb Delta_anh(410 K) %.1f%% (outside +-3) and rho_TE+anh(510 K) %.1f '
     '(%+.1f%%)' % (A['anharm_inverse_square_moment']['410'], A['te_anh_510_inverse_square_moment'],
                    pct(A['te_anh_510_inverse_square_moment'], t['pb_te_anharm_resistivity_510K'])), 'lambda_like_inverse_square_moment'),
    ('dp5', 'Transfer to Nb: the Pb/Al cancellation of TE and anharmonicity and the Pb-validated surrogate calibration '
            'do not hold for d-band Nb; its own 1000 K TE+anh run must calibrate the estimate',
     'Assume cancellation (TE+anh = harmonic), or reuse the Pb calibration, or hold TE+anh/TE fixed',
     'Nb rho(2000 K) %.1f (%+.1f%%), %.1f (%+.1f%%) or %.1f (%+.1f%%) instead of %.2f' % (
         A['nb_cancellation'], pct(A['nb_cancellation'], t['nb_te_anharm_resistivity_2000K']),
         A['nb_pb_calibrated_surrogate'], pct(A['nb_pb_calibrated_surrogate'], t['nb_te_anharm_resistivity_2000K']),
         A['nb_constant_anharmonic_fraction'], pct(A['nb_constant_anharmonic_fraction'], t['nb_te_anharm_resistivity_2000K']),
         t['nb_te_anharm_resistivity_2000K']), 'nb_cancellation_as_in_pb_al'),
    ('dp6', 'Pb thermal expansion near melting: superlinear growth, with an electronic/deformation-potential residual '
            'beyond phonon softening that grows with T',
     'Extrapolate rho_TE linearly, or scale by phonon softening alone', 'Pb rho_TE(610 K) %.1f (%+.1f%%) or %.1f (%+.1f%%) '
     'instead of %.2f' % (A['te_610_linear_extrapolation'], pct(A['te_610_linear_extrapolation'], t['pb_te_resistivity_610K']),
                          A['te_610_phonon_only'], pct(A['te_610_phonon_only'], t['pb_te_resistivity_610K']), t['pb_te_resistivity_610K']),
     'te_linear_extrapolation'),
    ('dp7', 'SERTA (printed first) vs the converged iterative BTE', 'Read the SERTA tensors',
     'Pb rho_TE+anh(510 K) %+.1f%%, Nb %+.1f%%' % (pct(A['serta']['ta510'], t['pb_te_anharm_resistivity_510K']),
                                                   pct(A['serta']['nb'], t['nb_te_anharm_resistivity_2000K'])), 'serta_instead_of_ibte'),
    ('dp8', 'The 10 K results are grid-limited (rho 1.57 vs 2.95 micro-Ohm cm for lattices 0.1% apart; Bloch-Gruneisen '
            'regime): exclude them from fits and surrogate calibration',
     'Fit TE+anh linearly including 10 K, or calibrate the surrogate including the 10 K pair',
     'Pb rho_TE+anh(510 K) %.1f (%+.1f%%) or %.1f (%+.1f%%)' % (
         A['te_anh_510_linear_fit_with_10K'], pct(A['te_anh_510_linear_fit_with_10K'], t['pb_te_anharm_resistivity_510K']),
         A['calibration_with_10K']['ta510'], pct(A['calibration_with_10K']['ta510'], t['pb_te_anharm_resistivity_510K'])),
     'linear_fit_including_10K'),
]

design = {
    'task_id': 'metal-resistivity-te-vs-anharmonicity', 'catalog_id': 'metal-resistivity-te-vs-anharmonicity',
    'start': 'released', 'source_group': 'wang_2026_metal_resistivity',
    'title': 'Pb, Nb, Al: thermal expansion versus anharmonicity in phonon-limited resistivity',
    'decision_points': [{'id': i, 'fork': f, 'wrong_choice': w, 'consequence': c, 'ablation': a} for i, f, w, c, a in DPS],
    'graded_fields': ['al_te_change_percent_481K within +-%.1f percentage points' % th['al_te_change_abs_pp'],
                      'pb_te_resistivity_610K_uohm_cm within +-%g%%' % (100 * th['pb_te_610_rel']),
                      'pb_anharmonic_change_percent_0K_lattice (410, 610 K) within +-%.1f percentage points each' % th['pb_anharm_abs_pp'],
                      'pb_te_anharm_resistivity_510K_uohm_cm within +-%g%%' % (100 * th['pb_te_anharm_510_rel']),
                      'nb_te_anharm_resistivity_2000K_uohm_cm within +-%g%%' % (100 * th['nb_te_anharm_2000_rel'])],
    'length_estimate': {'oracle_lab_evaluations': oracle.get('evaluations'), 'oracle_charge': oracle.get('charge_used'),
                        'budget': oracle.get('budget'),
                        'catalog_units': '86 (27 served transport runs x 2 + 32 force-constant sets x 1; phonon requests '
                                         'can be repeated at other grids)',
                        'tool_calls': '45-75', 'gpt6astra_cost_usd': '2.5-4.5',
                        'why': 'Transport lookups across three metals and four conditions, phonon grids for about 15 '
                               'force-constant sets, building and validating a transport-weighted surrogate, an '
                               'extrapolation near melting and a Nb-specific calibration.',
                        'oracle_errors': oracle.get('errors')},
    'runtime_dependencies': ['numpy (runtime image); no extra layer'],
    'data_provenance': {'archive_url': 'https://archive.materialscloud.org/records/zxf0b-rs406',
                        'doi': '10.24435/materialscloud:4h-r1', 'paper_doi': '10.1103/yk57-bt6t',
                        'files_used': ['data_for_paper.tar.gz (150 members, see members.md5)', 'README.md', 'Acknowledgements.md'],
                        'md5': {'data_for_paper.tar.gz': '83486659df1249b85d0025c2231d06b6',
                                'README.md': '089442c4cc8bfd798f062c178048897c',
                                'Acknowledgements.md': '968797ce577c1c08e4ab2b2195796e48'},
                        'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': 'Adapted from the candidate design (see the construction README): the Nb 1500 K TE+anh output is missing '
             'from the archive, so Nb is graded at 2000 K; phonon-moment surrogates miss EPW ratios by 5-17%, so '
             'targets are withheld archived values whose defensible estimators (transport-weighted surrogate validated '
             'on released pairs, ratio and polynomial extrapolations) agree within the tolerances.',
    'calibration': {
        'al_te_change_percent_481K': {'reference': r4(t['al_te_change_percent_481K']),
                                      'definition': 'TE run (50^3) / harmonic 50^3 run interpolated quadratically to 481 K',
                                      'robustness': {'linear_interpolation': r4(R['al_linear_interpolation'])},
                                      'nearest_ablation': {'mesh70_reference': r4(A['al_reference_mesh70']),
                                                           'mesh60_reference': r4(A['al_reference_mesh60'])},
                                      'not_ablations': {'serta': r4(A['serta']['al'])},
                                      'tolerance': '+-%.1f pp' % th['al_te_change_abs_pp']},
        'pb_te_resistivity_610K_uohm_cm': {'reference': r4(t['pb_te_resistivity_610K']),
                                           'oracle_route': 'q2 surrogate DFPT(a_T)/DFPT(a_0) x residual e(T) linear from 410/510 K',
                                           'oracle_model': r4(cal['oracle_model']['te610']),
                                           'robustness': {k: r4(R[k]) for k in ('te_610_ratio_quadratic', 'te_610_rho_cubic',
                                                                                 'te_610_rho_quadratic', 'te_610_residual_constant',
                                                                                 'te_610_uniform_weight_with_residual')} |
                                                         {'longitudinal_q2': r4(R['longitudinal_q2']['te610'])},
                                           'nearest_ablation': {'linear_extrapolation': r4(A['te_610_linear_extrapolation']),
                                                                'phonon_only': r4(A['te_610_phonon_only'])},
                                           'tolerance': '+-%g%%' % (100 * th['pb_te_610_rel'])},
        'pb_anharmonic_change_percent_0K_lattice': {
            'reference': r4(t['pb_anharmonic_change_percent_0K_lattice']),
            'oracle_route': 'q2 surrogate TDEP(a_0)/DFPT(a_0), calibrated on released 210/310 K TE+anh vs TE pairs',
            'oracle_model': r4(cal['oracle_model']['anh']),
            'robustness': {'longitudinal_q2': r4(R['longitudinal_q2']['anh']), 'calibration_110_310': r4(R['calibration_110_310']['anh']),
                           'uncalibrated_q2': r4(R['uncalibrated_q2_pb']['anh'])},
            'nearest_ablation': {'expanded_lattice_transfer': r4(A['anharm_from_expanded_lattice']),
                                 'inverse_square_moment': r4(A['anharm_inverse_square_moment'])},
            'not_ablations': {'calibration_including_10K (passes here, fails at 510 K)': r4(A['calibration_with_10K']['anh']),
                              'serta': r4(A['serta']['anh'])},
            'tolerance': '+-%.1f pp each' % th['pb_anharm_abs_pp']},
        'pb_te_anharm_resistivity_510K_uohm_cm': {
            'reference': r4(t['pb_te_anharm_resistivity_510K']),
            'oracle_route': 'TE(510 K) x q2 surrogate TDEP(510)/DFPT(a_510) / c',
            'oracle_model': r4(cal['oracle_model']['ta510']),
            'robustness': {'longitudinal_q2': r4(R['longitudinal_q2']['ta510']), 'calibration_110_310': r4(R['calibration_110_310']['ta510']),
                           'ratio_to_harmonic_extrapolated': r4(R['te_anh_510_ratio_to_harmonic_linear']),
                           'linear_fit_110_310': r4(R['te_anh_510_linear_fit_110_310'])},
            'nearest_ablation': {'calibration_including_10K': r4(A['calibration_with_10K']['ta510']), 'serta': r4(A['serta']['ta510']),
                                 'inverse_square_moment': r4(A['te_anh_510_inverse_square_moment']),
                                 'linear_fit_including_10K': r4(A['te_anh_510_linear_fit_with_10K']),
                                 'separable_multiplicative': r4(A['te_anh_510_separable_multiplicative']),
                                 'separable_additive': r4(A['te_anh_510_separable_additive'])},
            'tolerance': '+-%g%%' % (100 * th['pb_te_anharm_510_rel'])},
        'nb_te_anharm_resistivity_2000K_uohm_cm': {
            'reference': r4(t['nb_te_anharm_resistivity_2000K']),
            'oracle_route': 'TE(2000 K) x q2 surrogate TDEP/DFPT(a_2000) / c_Nb (Nb 1000 K pair)',
            'oracle_model': r4(cal['oracle_model']['nb']),
            'robustness': {'longitudinal_q2': r4(R['longitudinal_q2']['nb']), 'ratio_to_harmonic': r4(R['nb_ratio_to_harmonic']),
                           'proportional_T': r4(R['nb_proportional_T'])},
            'nearest_ablation': {'pb_calibrated_surrogate': r4(A['nb_pb_calibrated_surrogate']),
                                 'uncalibrated_surrogate': r4(A['nb_uncalibrated_surrogate']), 'serta': r4(A['serta']['nb']),
                                 'constant_anharmonic_fraction': r4(A['nb_constant_anharmonic_fraction']),
                                 'cancellation': r4(A['nb_cancellation'])},
            'tolerance': '+-%g%%' % (100 * th['nb_te_anharm_2000_rel'])}},
}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print('wrote design.json with', len(DPS), 'decision points')
