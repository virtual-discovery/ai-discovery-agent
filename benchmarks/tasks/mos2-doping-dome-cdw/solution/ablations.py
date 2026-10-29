"""Plausible wrong paths for mos2-doping-dome-cdw; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/mos2-doping-dome-cdw/make_task.py from the archived data (the hidden-doping
values of the wrong paths use the withheld spectra, i.e. what each wrong method would really give).
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(**changes):
    def mutate(payload, trace):
        for field, value in changes.items():
            payload[field] = value
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'tc_linear_interpolation_across_lambda_jump', setting(tc_1x1_K=A['tc_linear_interpolation'])),
    ('dp1', 'unscaled_average_of_neighbour_spectra', setting(tc_1x1_K=A['tc_unscaled_average_spectrum'])),
    ('dp2', 'nonadiabatic_shifts_only_NA', setting(tc_1x1_K=A['tc_NA_spectra'], tc_1x1_soc_at_critical_doping_K=A['soc_NA_spectra'])),
    ('dp2', 'adiabatic_spectra', setting(tc_1x1_K=A['tc_A_spectra'])),
    ('dp3', 'constant_dos_for_1x1', setting(tc_1x1_K=A['tc_constant_dos'], tc_1x1_soc_at_critical_doping_K=A['soc_constant_dos'])),
    ('dp4', 'allen_dynes_formula', setting(tc_1x1_K=A['tc_allen_dynes'], tc_1x1_soc_at_critical_doping_K=A['soc_allen_dynes'],
                                           tc_cdw_high_K=A['cdw_allen_dynes'])),
    ('dp4', 'mcmillan_formula', setting(tc_1x1_K=A['tc_mcmillan'], tc_1x1_soc_at_critical_doping_K=A['soc_mcmillan'],
                                        tc_cdw_high_K=A['cdw_mcmillan'])),
    ('dp5', 'boundary_at_soft_mode_instability', setting(critical_doping_e_per_cell=A['xc_soft_mode_instability'],
                                                         tc_1x1_soc_at_critical_doping_K=A['soc_at_soft_mode_instability'])),
    ('dp5', 'boundary_at_cdw_spinodal', setting(critical_doping_e_per_cell=A['xc_cdw_spinodal'],
                                                tc_1x1_soc_at_critical_doping_K=A['soc_at_cdw_spinodal'])),
    ('dp6', 'doping_per_2x2_cell', setting(critical_doping_e_per_cell=A['xc_per_2x2_cell'])),
    ('dp6', 'doping_as_areal_density', setting(critical_doping_e_per_cell=A['xc_areal_density_1e14'])),
    ('dp7', 'soc_omitted', setting(tc_1x1_soc_at_critical_doping_K=A['soc_omitted'])),
    ('dp8', 'cdw_dopings_from_caption', setting(tc_cdw_high_K=A['cdw_caption_dopings'])),
    ('dp8', 'cdw_labels_in_table_order', setting(tc_cdw_high_K=A['cdw_table_order'])),
    ('dp9', 'archived_tables_taken_as_eliashberg', setting(tc_1x1_K=A['tc_archived_table'], tc_cdw_high_K=A['cdw_archived_table'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('custom_spectrum_rigid_band_lower_dos', setting(tc_1x1_K=R['tc_custom_rigid_band'])),
    ('custom_spectrum_rigid_band_upper_dos', setting(tc_1x1_K=R['tc_custom_rigid_band_upper_dos'])),
    ('moments_allen_dynes_without_shape_correction', setting(tc_1x1_K=R['tc_moments_no_shape'])),
    ('lambda_weighted_mix_of_neighbour_spectra', setting(tc_1x1_K=R['tc_lambda_weighted_mix'])),
    ('quadratic_pes_crossing', setting(critical_doping_e_per_cell=R['xc_quadratic'],
                                       tc_1x1_soc_at_critical_doping_K=R['soc_at_xc_quadratic'])),
    ('cdw_well_minimum_crossing', setting(critical_doping_e_per_cell=R['xc_cdw_well_minimum'],
                                          tc_1x1_soc_at_critical_doping_K=R['soc_at_xc_well_minimum'])),
    ('soc_ratio_to_scalar_relativistic_curve', setting(tc_1x1_soc_at_critical_doping_K=R['soc_ratio_method'])),
]
