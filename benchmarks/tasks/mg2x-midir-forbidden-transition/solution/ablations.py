"""Plausible wrong paths for mg2x-midir-forbidden-transition; each must fail grading.

Values come from solution/ablation_values.json, computed from the archived data by
benchmarks/construction/mg2x-midir-forbidden-transition/make_task.py.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def merge(base, new):
    out = dict(base)
    for k, v in new.items():
        out[k] = merge(base.get(k, {}), v) if isinstance(v, dict) and isinstance(base.get(k), dict) else v
    return out


def setting(**changes):
    def mutate(payload, trace):
        for field, value in changes.items():
            payload[field] = merge(payload[field], value) if isinstance(value, dict) else value
        return payload, trace
    return mutate


ABLATIONS = [
    # dp1: X1/X3 labelled by energy order (Mg2Sn sign wrong at every level)
    ('dp1', 'energy_order_labels', setting(x3_minus_x1_eV=A['split']['energy_order_labels'])),
    # dp2: G0W0 read from the TOTE column described as eQP without Z (Z = 1), or from the input eigenvalues
    ('dp2', 'g0w0_z_equal_1_column', setting(x3_minus_x1_eV=A['split']['g0w0_z_equal_1_column'])),
    ('dp2', 'g0w0_input_eigenvalue_column', setting(x3_minus_x1_eV=A['split']['g0w0_input_eigenvalue_column'])),
    ('dp2', 'qsgw_tote_qp_column_at_100K_lattice', setting(predicted_splitting_eV=A['pred']['qsgw_column7_low_temperature'])),
    # dp3: ladder correction across codes, or the QSGW80 scaling taken as the vertex correction
    ('dp3', 'vertex_cross_code', setting(vertex_correction_Mg2Si_meV=A['vertex']['cross_code'])),
    ('dp3', 'vertex_from_qsgw80', setting(vertex_correction_Mg2Si_meV=A['vertex']['qsgw80_scaling'])),
    ('dp3', 'qsgw80_as_qsgw', setting(x3_minus_x1_eV=A['split']['qsgw80_as_qsgw'])),
    # dp4: matrix element taken as allowed (constant), as isotropic q^2, or as the full PBE matrix element
    ('dp4', 'constant_matrix_element', setting(ir_raman_ratio_300K_over_77K=A['ratio']['constant_matrix_element'])),
    ('dp4', 'isotropic_q2_matrix_element', setting(ir_raman_ratio_300K_over_77K=A['ratio']['isotropic_q2'])),
    ('dp4', 'full_pbe_matrix_element', setting(ir_raman_ratio_300K_over_77K=A['ratio']['full_pbe_matrix_element'])),
    # dp5: parabolic valley from the fitted masses, or the X-W line taken as isotropic
    ('dp5', 'parabolic_valley', setting(ir_raman_ratio_300K_over_77K=A['ratio']['parabolic'])),
    ('dp5', 'xw_line_isotropic', setting(ir_raman_ratio_300K_over_77K=A['ratio']['xw_line_isotropic'])),
    # dp6: degeneracy and counting on the box
    ('dp6', 'boltzmann_statistics', setting(ir_raman_ratio_300K_over_77K=A['ratio']['boltzmann'])),
    ('dp6', 'spin_not_counted', setting(ir_raman_ratio_300K_over_77K=A['ratio']['spin_not_counted'])),
    ('dp6', 'one_octant_only', setting(ir_raman_ratio_300K_over_77K=A['ratio']['one_octant_only'])),
    ('dp6', 'no_symmetry_plane_weights', setting(ir_raman_ratio_300K_over_77K=A['ratio']['no_symmetry_plane_weights'])),
    # dp7: Mg2Sn ratio transferred from Mg2Si, or the X1 state assumed occupied in Mg2Sn
    ('dp7', 'mg2si_ratio_for_mg2sn', setting(ir_raman_ratio_300K_over_77K=A['ratio']['mg2si_value_for_mg2sn'])),
    ('dp7', 'x1_occupied_in_mg2sn', setting(ir_raman_ratio_300K_over_77K=A['ratio']['x1_occupied_in_mg2sn'])),
    # dp8: lattice constant not matched to the temperature
    ('dp8', 'room_temperature_lattice_at_T0', setting(predicted_splitting_eV=A['pred']['room_temperature_lattice_at_T0'])),
    ('dp8', 'low_temperature_lattice_at_300K', setting(predicted_splitting_eV=A['pred']['low_temperature_lattice_at_300K'])),
    # dp9: electron-phonon renormalization ignored, zero-point only, or with the sign of X1 - X3
    ('dp9', 'no_electron_phonon', setting(predicted_splitting_eV=A['pred']['no_electron_phonon'])),
    ('dp9', 'zero_point_only_at_300K', setting(predicted_splitting_eV=A['pred']['zero_point_only_at_300K'])),
    ('dp9', 'electron_phonon_sign_reversed', setting(predicted_splitting_eV=A['pred']['electron_phonon_sign_reversed'])),
    # also a data-selection error: G0W0 from the PBE start
    ('dp2', 'g0w0_from_pbe_start', setting(x3_minus_x1_eV=A['split']['g0w0_from_pbe'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('ratio_nested_boxes_at_both_temperatures', setting(ir_raman_ratio_300K_over_77K=R['ratio']['nested_boxes_both_T'])),
    ('ratio_octant_unfolded_plain_sum', setting(ir_raman_ratio_300K_over_77K=R['ratio']['unfolded_plain_sum'])),
    ('ratio_coarse_box_at_both_temperatures', setting(ir_raman_ratio_300K_over_77K=R['ratio']['coarse_box_both_T'])),
    ('eph_special_displacement_4x4x4', setting(predicted_splitting_eV=R['pred']['zg_4x4x4'])),
    ('eph_ahc_adiabatic', setting(predicted_splitting_eV=R['pred']['ahc_adiabatic'])),
    ('eph_ahc_q16_or_ge3d_valence', setting(predicted_splitting_eV=R['pred']['ahc_q16_or_ge3d'])),
    ('vertex_ladder_window_nc12', setting(vertex_correction_Mg2Si_meV=R['vertex']['nc12'])),
    ('vertex_larger_basis', setting(vertex_correction_Mg2Si_meV=R['vertex']['basis2'])),
]
