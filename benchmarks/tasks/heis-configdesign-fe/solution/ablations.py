"""Plausible wrong paths for heis-configdesign-fe; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/heis-configdesign-fe/make_task.py from the released data.
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
    ('dp1', 'spiral_001_only', setting(spiral_predictions_meV_per_atom=A['spiral_001_only'], J1_meV=A['J1_001_only'])),
    ('dp2', 'flat_spiral_convention', setting(J1_meV=A['J1_flat_spiral'], spin_wave_stiffness_meV_A2=A['stiffness_flat_spiral'])),
    ('dp3', 'four_shells', setting(spiral_predictions_meV_per_atom=A['spiral_four_shells'])),
    ('dp4', 'spiral_model_for_random_states', setting(random_predictions_meV_per_atom=A['random_from_spiral_model'])),
    ('dp5', 'random_341_only', setting(random_predictions_meV_per_atom=A['random_341_only'])),
    ('dp5', 'random_221_three_shells', setting(random_predictions_meV_per_atom=A['random_221_three_shells'])),
    ('dp6', 'random_wrong_lattice_constant', setting(random_predictions_meV_per_atom=A['random_wrong_lattice'])),
    ('dp7', 'stiffness_without_2pi', setting(spin_wave_stiffness_meV_A2=A['stiffness_no_2pi'])),
    ('dp8', 'stiffness_from_truncated_model', setting(spin_wave_stiffness_meV_A2=A['stiffness_from_truncated_model'])),
    ('dp9', 'J1_from_random_states', setting(J1_meV=A['J1_from_random_states'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('ten_shell_spiral_model', setting(spiral_predictions_meV_per_atom=R['spiral_ten_shells'], J1_meV=R['J1_ten_shells'])),
    ('eight_shell_random_model_all_cells', setting(random_predictions_meV_per_atom=R['random_eight_shells'])),
    ('stiffness_three_smallest_q', setting(spin_wave_stiffness_meV_A2=R['stiffness_kmax3'])),
]
