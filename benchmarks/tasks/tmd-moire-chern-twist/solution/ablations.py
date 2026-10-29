"""Plausible wrong paths for tmd-moire-chern-twist; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/tmd-moire-chern-twist/make_task.py from the archived data and reference calculations.
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
    ('dp1', 'spin_down_sector', setting(dft_chern_K_up=A['dft_spin_down_sector'])),
    ('dp2', 'indicator_without_minus_sign', setting(dft_chern_K_up=A['dft_no_minus_sign'])),
    ('dp3', 'kappa_prime_as_conjugate', setting(dft_chern_K_up=A['dft_kappa_prime_conjugate'])),
    ('dp4', 'state_index_as_band_index', setting(dft_chern_K_up=A['dft_state_index_as_band'])),
    ('dp5', 'continuum_plot_script_basis', setting(continuum_chern_3p89=A['continuum_plot_script_basis'])),
    ('dp5', 'continuum_left_handed_basis', setting(continuum_chern_3p89=A['continuum_left_handed_basis'])),
    ('dp5', 'continuum_layers_swapped', setting(continuum_chern_3p89=A['continuum_layers_swapped'])),
    ('dp6', 'literature_local_stacking_set', setting(continuum_set_matching_dft=A['continuum_set_literature_prior'])),
    ('dp7', 'fixed_continuum_transferred_to_dft', setting(dft_chern_K_up=A['dft_from_fixed_continuum_2p14'])),
    ('dp7', 'continuum_assumed_to_track_dft', setting(continuum_P_chern_2p14=A['continuum_2p14_equals_dft'])),
    ('dp8', 'xm_mx_swapped', setting(hartree_xm_meV=A['hartree_xm_mx_swapped'],
                                     hartree_first_star_xm_meV=A['first_star_xm_mx_swapped'])),
    ('dp9', 'first_star_reconstruction_as_xm', setting(hartree_xm_meV=A['hartree_xm_first_star_only'])),
    ('dp9', 'full_value_as_first_star', setting(hartree_first_star_xm_meV=A['first_star_equals_full_xm'])),
    ('dp10', 'fourier_magnitudes_with_large_angle_sign', setting(hartree_first_star_xm_meV=A['first_star_from_magnitudes_3p15_sign'])),
    ('dp10', 'no_sign_flip_assumed', setting(hartree_xm_meV=A['hartree_xm_no_sign_flip'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('xm_nearest_grid_point', setting(hartree_xm_meV=R['hartree_xm_nearest_grid'])),
    ('xm_bilinear_interpolation', setting(hartree_xm_meV=R['hartree_xm_bilinear'])),
    ('xm_from_archived_line_cuts', setting(hartree_xm_meV=R['hartree_xm_archived_linecut'])),
    ('first_star_fft_with_duplicated_endpoint', setting(hartree_first_star_xm_meV=R['first_star_full_grid_fft'])),
]
