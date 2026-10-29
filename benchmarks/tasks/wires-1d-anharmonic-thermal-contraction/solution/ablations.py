"""Plausible wrong paths for wires-1d-anharmonic-thermal-contraction; each must fail grading.

Values come from solution/ablation_values.json, computed from the archived data by
benchmarks/construction/wires-1d-anharmonic-thermal-contraction/calib.py.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(**changes):
    """changes: field path 'TaSe3.alpha_L_300K_1e-6_per_K' -> value."""
    def mutate(payload, trace):
        for path, value in changes.items():
            wire, field = path.split('.', 1)
            payload[wire][field] = value
        return payload, trace
    return mutate


def paths(**kw):
    names = {'k': 'TaSe3.axial_stiffness_eV', 'alpha': 'TaSe3.alpha_L_300K_1e-6_per_K', 'cv_ta': 'TaSe3.cv_ratio_300K',
             'T': 'TaSe3.instability_temperature_K', 'cv_cu': 'CuC2.cv_ratio_300K'}
    return setting(**{names[k]: v for k, v in kw.items()})


ABLATIONS = [
    ('dp1', 'gamma_rotation_mode_counted_as_instability', paths(T=A['T_inst_rotation_mode'])),
    ('dp2', 'stability_from_auxiliary_frequencies', paths(T=None)),
    ('dp2', 'heat_capacity_from_hessian_frequencies', paths(cv_cu=A['cv_cuc2_hessian_freq'], cv_ta=A['cv_tase3_hessian_freq'])),
    ('dp3', 'stress_read_as_tension_positive', paths(alpha=A['alpha_sign'])),
    ('dp4', 'ab_initio_average_stress', paths(alpha=A['alpha_average_stress'])),
    ('dp5', 'plain_asr_long_wavelength_stiffness', paths(k=A['k_plain_asr'], alpha=A['alpha_plain_asr'])),
    ('dp6', 'finite_strain_secant_stiffness', paths(k=A['k_finite_strain'], alpha=A['alpha_finite_strain'])),
    ('dp7', 'stiffness_from_E_equals_k_eps2', paths(alpha=A['alpha_half_k'])),
    ('dp7', 'supercell_volume_with_cell_stiffness', paths(alpha=A['alpha_supercell_volume'])),
    ('dp8', 'derivative_term_with_plus_sign', paths(cv_cu=A['cv_cuc2_paper_sign'], cv_ta=A['cv_tase3_paper_sign'])),
    ('dp9', 'harmonic_cv_with_sscha_frequencies', paths(cv_cu=A['cv_cuc2_harmonic_sscha_freq'], cv_ta=A['cv_tase3_harmonic_sscha_freq'])),
    ('dp9', 'harmonic_cv_with_dfpt_frequencies', paths(cv_cu=A['cv_cuc2_dfpt_harmonic'], cv_ta=A['cv_tase3_dfpt_harmonic'])),
    ('dp10', 'first_released_unstable_temperature', paths(T=A['T_inst_first_released_unstable'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('stiffness_q_0.05', paths(k=R['k_q005'], alpha=R['alpha_k_q005'])),
    ('stiffness_from_T0_sscha_hessian_global_cubic_stress', paths(k=R['k_T0_hessian'], alpha=R['alpha_global_cubic_T0_hessian_k'])),
    ('stress_quadratic_100_500', paths(alpha=R['alpha_wide_quadratic'])),
    ('cv_cubic_fits', paths(cv_cu=R['cv_cuc2_cubic'], cv_ta=R['cv_tase3_cubic'])),
    ('cv_linear_200_500', paths(cv_cu=R['cv_cuc2_linear_200_500'], cv_ta=R['cv_tase3_linear_200_500'])),
    ('cv_all_converged_cubic', paths(cv_cu=R['cv_cuc2_all_converged_cubic'], cv_ta=R['cv_tase3_all_converged_cubic'])),
    ('instability_600_800_100k', paths(T=R['T_inst_released_600_800_100k'])),
    ('instability_quadratic_fit_500_1000', paths(T=R['T_inst_quadratic_fit'])),
]
