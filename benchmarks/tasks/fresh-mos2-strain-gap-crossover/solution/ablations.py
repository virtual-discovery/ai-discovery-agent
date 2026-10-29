"""Plausible wrong paths for fresh-mos2-strain-gap-crossover; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/fresh-mos2-strain-gap-crossover/make_task.py from the reference pw.x calculations
(the same lab code path the agent uses) with the estimators of solve.py.
"""
import copy
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def merge(target, changes):
    for key, value in changes.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            merge(target[key], value)
        else:
            target[key] = copy.deepcopy(value)


def setting(changes):
    def mutate(payload, trace):
        merge(payload, changes)
        return payload, trace
    return mutate


ABLATIONS = [(entry['decision_point'], name, setting(entry['changes'])) for name, entry in A.items()]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, setting(changes)) for name, changes in R.items()]
