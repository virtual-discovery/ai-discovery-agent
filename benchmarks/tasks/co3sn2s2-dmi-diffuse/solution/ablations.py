"""Plausible wrong paths for co3sn2s2-dmi-diffuse; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/co3sn2s2-dmi-diffuse/make_task.py from the archived data and the archived Hamiltonian.
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
    ('dp1', 'narrow_scans_only', setting(narrow_003_profiles=A['narrow_scans_only'])),
    ('dp1', 'narrow_scans_only_bin_interpolation', setting(narrow_003_profiles=A['narrow_scans_only_bin_interpolation'])),
    ('dp1', 'narrow_scans_only_independent_fits', setting(narrow_003_profiles=A['narrow_scans_only_independent_fits'])),
    ('dp2', 'critical_scattering_mirrored_about_tc', setting(narrow_003_profiles=A['mirror_about_tc'])),
    ('dp3', 'wide_101_background_from_1p5K', setting(wide_101_profile_180K=A['wide_101_background_from_1p5K'])),
    ('dp3', 'wide_101_background_from_300K', setting(wide_101_profile_180K=A['wide_101_background_from_300K'])),
    ('dp4', 'wide_101_equal_to_003', setting(wide_101_profile_180K=A['wide_101_equal_to_003'])),
    ('dp4', 'wide_101_c_axis_polarisation_only', setting(wide_101_profile_180K=A['wide_101_c_axis_polarisation_only'])),
    ('dp5', 'wide_101_in_plane_afm_correlations', setting(wide_101_profile_180K=A['wide_101_in_plane_afm_correlations'])),
    ('dp6', 'dm_matrix_copied_to_all_bonds', setting(canting_angle_deg=A['tilt_copy_matrix'], spin_wave_gap_meV=A['gap_copy_matrix'])),
    ('dp7', 'out_of_plane_dm_only', setting(canting_angle_deg=A['tilt_dz_only'], spin_wave_gap_meV=A['gap_dz_only'])),
    ('dp7', 'in_plane_dm_only', setting(canting_angle_deg=A['tilt_in_plane_only'], spin_wave_gap_meV=A['gap_in_plane_only'])),
    ('dp8', 'dm_matrix_on_reversed_bond', setting(canting_angle_deg=A['tilt_reversed_bond'], spin_wave_gap_meV=A['gap_reversed_bond'])),
    ('dp9', 'spin_length_from_ordered_moment', setting(spin_wave_gap_meV=A['gap_spin_from_ordered_moment'])),
    ('dp10', 'bonds_counted_twice', setting(spin_wave_gap_meV=A['gap_bonds_counted_twice'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('narrow_fit_with_250K_and_300K_data', setting(narrow_003_profiles=R['narrow_more_data'])),
    ('narrow_per_annulus_interpolation', setting(narrow_003_profiles=R['narrow_bin_interpolation'])),
    ('narrow_independent_fits', setting(narrow_003_profiles=R['narrow_independent_fits'])),
    ('wide_101_background_only', setting(wide_101_profile_180K=R['wide_101_background_only'])),
    ('wide_101_xy_polarisation_ratio', setting(wide_101_profile_180K=R['wide_101_xy_ratio'])),
    ('wide_101_unsmoothed_background', setting(wide_101_profile_180K=R['wide_101_raw_background'])),
    ('wide_101_background_165K_only', setting(wide_101_profile_180K=R['wide_101_background_165K_only'])),
    ('tilt_rounded_and_conventional_cell_gap', setting(canting_angle_deg=R['tilt_rounded_ground_state'],
                                                       spin_wave_gap_meV=R['gap_from_conventional_cell'])),
]
