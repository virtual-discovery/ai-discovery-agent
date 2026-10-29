"""Plausible wrong paths for bab2-hole-doping-two-gap; each must fail grading.

Values come from solution/ablation_values.json, computed from the archived data by
benchmarks/construction/bab2-hole-doping-two-gap/make_task.py.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(**changes):
    def mutate(payload, trace):
        for field, value in changes.items():
            if isinstance(value, dict) and isinstance(payload.get(field), dict):
                payload[field] = {**payload[field], **value}
            else:
                payload[field] = value
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'readme_order_mapping', setting(critical_hole_doping=A['xc_readme_order'], iso_tc_K=A['iso_readme_order'])),
    ('dp2', 'omega_linear_crossing', setting(critical_hole_doping=A['xc_omega_linear'],
                                             soft_mode_cm1={'homogeneous_one_ninth': A['hom_omega_linear']})),
    ('dp3', 'supercell_331_as_homogeneous', setting(soft_mode_cm1={'homogeneous_one_ninth': A['hom_331_as_homogeneous']})),
    ('dp3', 'supercell_gamma_point', setting(soft_mode_cm1={'cs_substituted': A['cs_gamma_point']})),
    ('dp4', 'archived_mu013_logs', setting(iso_tc_K=A['iso_archived_mu013'])),
    ('dp4', 'allen_dynes_tc', setting(iso_tc_K=A['iso_allen_dynes'])),
    ('dp5', 'unconverged_gap_scan', setting(iso_tc_K=A['iso_unconverged_scan'])),
    ('dp6', 'isotropic_tc_as_anisotropic', setting(aniso_tc_K=A['aniso_iso_reported'])),
    ('dp6', 'linear_gap_extrapolation', setting(aniso_tc_K=A['aniso_linear_gap'])),
    ('dp7', 'real_axis_gap_ratio', setting(gap_ratio_25K=A['ratio_real_axis'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('xc_omega2_piecewise_released', setting(critical_hole_doping=R['xc_w2_linear_released_0.06_0.10'],
                                             soft_mode_cm1={'homogeneous_one_ninth': R['hom_w2_linear_released_0.06_0.10']})),
    ('xc_omega2_quadratic_released', setting(critical_hole_doping=R['xc_w2_quadratic_released'],
                                             soft_mode_cm1={'homogeneous_one_ninth': R['hom_w2_quadratic_released']})),
    ('xc_omega2_linear_regression_released', setting(critical_hole_doping=R['xc_w2_linear_regression_released'],
                                                     soft_mode_cm1={'homogeneous_one_ninth': R['hom_w2_linear_regression_released']})),
    ('iso_linearized_reference_column_10', setting(iso_tc_K=R['iso_a2f_column_10'])),
    ('iso_converged_extrapolation_grid', setting(iso_tc_K=R['iso_converged_extrapolation'])),
    ('aniso_large_cluster_squared', setting(aniso_tc_K=R['aniso_large_cluster_squared'])),
    ('aniso_sqrt_free_energy', setting(aniso_tc_K=R['aniso_sqrt_free_energy'])),
    ('aniso_three_point_fit', setting(aniso_tc_K=R['aniso_three_point_fit'])),
    ('aniso_real_axis_max', setting(aniso_tc_K=R['aniso_real_axis_max_squared'])),
    ('gap_ratio_median', setting(gap_ratio_25K=R['ratio_median'])),
    ('gap_ratio_unweighted_bins', setting(gap_ratio_25K=R['ratio_unweighted_bins'])),
]
