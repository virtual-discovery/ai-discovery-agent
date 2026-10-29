"""Plausible wrong paths for fresh-nio-exchange-dftu; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/fresh-nio-exchange-dftu/make_task.py from the reference pw.x calculations
(the same lab code path the agent uses), with the estimators of solution/solve.py.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(changes):
    def mutate(payload, trace):
        for field, value in changes.items():
            payload[field] = dict(payload[field], **value) if isinstance(value, dict) else value
        return payload, trace
    return mutate


ABLATIONS = [(dp, name, setting(A[name])) for dp, name in V['ablation_order']]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, setting(changes)) for name, changes in R.items()]
