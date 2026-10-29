"""Plausible wrong paths for tan-pressure-kappa-channels; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/tan-pressure-kappa-channels/reference.py with the lab's own code path from the
archived data (each is the full submission the wrong path produces).
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(values):
    def mutate(payload, trace):
        payload.update(values)
        return payload, trace
    return mutate


ABLATIONS = [(dp, name, setting(A[name])) for dp, name in V['ablation_decision_points']]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, setting(values)) for name, values in R.items()]
