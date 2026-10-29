"""Plausible wrong paths for metal-resistivity-te-vs-anharmonicity; each must fail grading.

Values come from solution/ablation_values.json, computed from the released data by
benchmarks/construction/metal-resistivity-te-vs-anharmonicity/make_task.py (model.py holds the estimators).
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
AL, TE610, ANH, TA510, NB = ('al_te_change_percent_481K', 'pb_te_resistivity_610K_uohm_cm',
                             'pb_anharmonic_change_percent_0K_lattice', 'pb_te_anharm_resistivity_510K_uohm_cm',
                             'nb_te_anharm_resistivity_2000K_uohm_cm')


def setting(**changes):
    def mutate(payload, trace):
        for field, value in changes.items():
            payload[field] = value
        return payload, trace
    return mutate


S = A['serta']
ABLATIONS = [
    ('dp1', 'al_reference_finest_mesh70', setting(**{AL: A['al_reference_mesh70']})),
    ('dp1', 'al_reference_mesh60', setting(**{AL: A['al_reference_mesh60']})),
    ('dp2', 'separable_multiplicative', setting(**{TA510: A['te_anh_510_separable_multiplicative']})),
    ('dp2', 'separable_additive', setting(**{TA510: A['te_anh_510_separable_additive']})),
    ('dp3', 'anharmonic_correction_from_expanded_lattice', setting(**{ANH: A['anharm_from_expanded_lattice']})),
    ('dp4', 'lambda_like_inverse_square_moment', setting(**{ANH: A['anharm_inverse_square_moment'],
                                                            TA510: A['te_anh_510_inverse_square_moment']})),
    ('dp5', 'nb_cancellation_as_in_pb_al', setting(**{NB: A['nb_cancellation']})),
    ('dp5', 'nb_surrogate_calibrated_on_pb', setting(**{NB: A['nb_pb_calibrated_surrogate']})),
    ('dp5', 'nb_surrogate_uncalibrated', setting(**{NB: A['nb_uncalibrated_surrogate']})),
    ('dp5', 'nb_constant_anharmonic_fraction', setting(**{NB: A['nb_constant_anharmonic_fraction']})),
    ('dp6', 'te_linear_extrapolation', setting(**{TE610: A['te_610_linear_extrapolation']})),
    ('dp6', 'te_phonon_softening_only', setting(**{TE610: A['te_610_phonon_only']})),
    ('dp7', 'serta_instead_of_ibte', setting(**{AL: S['al'], TE610: S['te610'], ANH: S['anh'], TA510: S['ta510'],
                                                 NB: S['nb']})),
    ('dp8', 'linear_fit_including_10K', setting(**{TA510: A['te_anh_510_linear_fit_with_10K']})),
    ('dp8', 'surrogate_calibration_including_10K', setting(**{ANH: A['calibration_with_10K']['anh'],
                                                              TA510: A['calibration_with_10K']['ta510']})),
]

L = R['longitudinal_q2']
C = R['calibration_110_310']
# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('longitudinal_projected_q2_surrogate', setting(**{ANH: L['anh'], TA510: L['ta510'], NB: L['nb'], TE610: L['te610']})),
    ('calibration_on_110_210_310K', setting(**{ANH: C['anh'], TA510: C['ta510']})),
    ('uncalibrated_q2_surrogate_for_pb', setting(**{ANH: R['uncalibrated_q2_pb']['anh']})),
    ('te_610_ratio_quadratic_extrapolation', setting(**{TE610: R['te_610_ratio_quadratic']})),
    ('te_610_cubic_fit', setting(**{TE610: R['te_610_rho_cubic']})),
    ('te_610_quadratic_fit', setting(**{TE610: R['te_610_rho_quadratic']})),
    ('te_610_constant_residual', setting(**{TE610: R['te_610_residual_constant']})),
    ('te_610_uniform_weight_with_residual', setting(**{TE610: R['te_610_uniform_weight_with_residual']})),
    ('nb_ratio_to_harmonic_from_1000K', setting(**{NB: R['nb_ratio_to_harmonic']})),
    ('nb_proportional_to_T', setting(**{NB: R['nb_proportional_T']})),
    ('te_anh_510_ratio_to_harmonic_extrapolated', setting(**{TA510: R['te_anh_510_ratio_to_harmonic_linear']})),
    ('te_anh_510_linear_fit_110_310K', setting(**{TA510: R['te_anh_510_linear_fit_110_310']})),
    ('al_linear_interpolation', setting(**{AL: R['al_linear_interpolation']})),
]
