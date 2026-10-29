"""Plausible wrong paths for ti-codoping-exchange-tc; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/ti-codoping-exchange-tc/reference.py and make_task.py from the archived exchange table
(mean-field formulas, or Monte Carlo of the wrong model with the lab's own code path).
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(entry):
    def mutate(payload, trace):
        payload['curie_temperature_K'].update(entry.get('curie_temperature_K', {}))
        if 'best_codopant_for_Cr' in entry:
            payload['best_codopant_for_Cr'] = entry['best_codopant_for_Cr']
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'mean_field_two_sublattice', setting(A['mean_field_two_sublattice'])),
    ('dp1', 'authors_averaged_mean_field', setting(A['authors_averaged_mean_field'])),
    ('dp2', 'single_placement_per_size', setting(A['single_placement_per_size'])),
    ('dp3', 'codopant_by_mean_field', setting(A['codopant_by_mean_field'])),
    ('dp4', 'linear_in_concentration', setting(A['linear_in_concentration'])),
    ('dp5', 'zero_filled_coco', setting(A['zero_filled_coco'])),
    ('dp6', 'near_shells_only', setting(A['near_shells_only'])),
    ('dp7', 'dimer_moment_collapse', setting(A['dimer_moment_collapse'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, setting(entry)) for name, entry in R.items()]
