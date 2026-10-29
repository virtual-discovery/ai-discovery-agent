"""Plausible wrong paths for hgte-gated-subband-densities; each must fail grading.

Values come from solution/ablation_values.json, computed from the released data by
benchmarks/construction/hgte-gated-subband-densities/make_task.py.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
F, A, R = V['fields'], V['ablations'], V['robustness']


def setting(values):
    def mutate(payload, trace):
        for key, value in values.items():
            section, name = F[key]
            (payload if section is None else payload[section])[name] = value
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'bulk_spin_degenerate', setting(A['bulk_spin_degenerate'])),
    ('dp2', 'hole_strongest_peak', setting(A['hole_strongest_peak'])),
    ('dp2', 'hole_peak_nearest_model', setting(A['hole_peak_nearest_model'])),
    ('dp3', 'top_is_larger_density', setting(A['top_is_larger_density'])),
    ('dp4', 'four_equivalent_pockets', setting(A['four_equivalent_pockets'])),
    ('dp4', 'single_hole_pocket', setting(A['single_hole_pocket'])),
    ('dp5', 'wedge_own_fermi_level', setting(A['wedge_own_fermi_level'])),
    ('dp6', 'band_minus1_split_in_two', setting(A['band_minus1_split_in_two'])),
    ('dp6', 'band_minus1_split_in_four', setting(A['band_minus1_split_in_four'])),
    ('dp6', 'wedge_not_normalised', setting(A['wedge_not_normalised'])),
    ('dp7', 'rigid_band_from_neutral_run', setting(A['rigid_band_from_neutral_run'])),
    ('dp8', 'holes_appear_at_neutrality', setting(A['holes_appear_at_neutrality'])),
    ('dp8', 'holes_appear_at_hall_sign_change', setting(A['holes_appear_at_hall_sign_change'])),
    ('dp9', 'steep_ridge_as_bulk_harmonic_as_top', setting(A['steep_ridge_as_bulk_harmonic_as_top'])),
    ('dp9', 'steep_ridge_as_bulk_sum_as_top', setting(A['steep_ridge_as_bulk_sum_as_top'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, setting(values)) for name, values in R.items()]
