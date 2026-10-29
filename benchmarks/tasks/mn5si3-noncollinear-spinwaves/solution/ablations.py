"""Plausible wrong paths for mn5si3-noncollinear-spinwaves; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/mn5si3-noncollinear-spinwaves/calibrate.py from purchasable data only.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())


def route(name):
    row = V[name]

    def mutate(payload, trace):
        payload['l_branch_positions_rlu'] = row['pred']
        payload['zone_centre_excitation_meV'] = round(row['optical'], 2)
        payload['mn1_carries_moment'] = row['mn1']
        payload['mn2_canting_deg'] = round(row['canting'], 2)
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'mn2_only_sites', route('dp1_mn2_only_sites')),
    ('dp2', 'bonds_double_counted', route('dp2_bonds_double_counted')),
    ('dp3', 'free_ion_spin_5_2', route('dp3_free_ion_spin_5_2')),
    ('dp4', 'lowest_branch_only', route('dp4_lowest_branch_only')),
    ('dp4', 'branch_weights_ignored', route('dp4_branch_weights_ignored')),
    ('dp5', 'collinear_ground_state_not_relaxed', route('dp5_collinear_ground_state_not_relaxed')),
    ('dp6', 'angle_to_mn1_axis', route('dp6_angle_to_mn1_axis')),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('collinear_reference_exchange_table', route('robust_collinear_reference_table')),
    ('long_range_17A_exchange_table', route('robust_long_range_r17_table')),
    ('momentum_resolution_0.04', route('robust_sigQ_0.04')),
    ('momentum_resolution_0.08', route('robust_sigQ_0.08')),
    ('energy_scale_0.9', route('robust_scale_0.9')),
    ('energy_scale_1.1', route('robust_scale_1.1')),
    ('resolution_from_all_cuts', route('robust_all_cuts_for_resolution')),
    ('mn1_spin_0.6', route('robust_mn1_spin_0.6')),
]
