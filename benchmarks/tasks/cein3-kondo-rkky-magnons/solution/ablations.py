"""Plausible wrong paths for cein3-kondo-rkky-magnons; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/cein3-kondo-rkky-magnons/calibrate.py and make_task.py from the released data
(fits of each model family to all purchasable data, predictions made with the same binning).
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(values):
    def mutate(payload, trace):
        for field, value in values.items():
            payload[field] = value
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'nearest_neighbour_model_from_12meV_data_only', setting(A['J1_HEonly'])),
    ('dp2', 'short_range_J1_J3_model', setting(A['J123_short'])),
    ('dp2', 'short_range_six_shell_model', setting(A['J6shell_short'])),
    ('dp3', 'momentum_integration_window_ignored', setting(A['noperp_J1LR'])),
    ('dp4', 'nearest_neighbour_model_forced_through_all_data', setting(A['J1_forced_all'])),
    ('dp5', 'free_ion_J_5_2_as_effective_spin', setting(A['free_ion_spin'])),
    ('dp6', 'broad_spectra_read_as_incoherent_continuum', setting(A['incoherent_continuum'])),
    ('dp7', 'ten_free_shells_overfit', setting(A['J10shell_free'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('J1_J3_plus_long_range_term', setting(R['J123LR'])),
    ('integration_window_read_as_total_width_0p08', setting(R['J1LR_h004'])),
    ('J1_plus_long_range_term_all_visible_data', setting(R['J1LR'])),
    ('J1_and_W_from_bandwidth_2p75_meV', setting(R['J1_from_bandwidth_2p75'])),
    ('eta_from_cut_box_widths', setting(R['eta_from_cut_box_widths'])),
]
