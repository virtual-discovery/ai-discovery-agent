"""Plausible wrong paths for fresh-srtio3-soft-modes; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/fresh-srtio3-soft-modes/make_task.py from the reference pw.x/ph.x calculations (the same
lab code path the agent uses) with the estimators of solve.py; none is an arbitrary perturbation.
"""
import copy
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(changes):
    def mutate(payload, trace):
        for key, value in changes.items():
            payload[key] = copy.deepcopy(value)
        return payload, trace
    return mutate


# dp1 k4 force constants and polar well (cascade into the competition verdict); dp2 Ti-only (Slater) displacement;
# dp3 experimental lattice parameter; dp4 4x4x4 dielectric tensors; dp5 analytic or uncoupled LO;
# dp6 tilt energies referenced to a 5-atom cell on another mesh, or per 10-atom cell; dp7 no Poisson relaxation.
ABLATIONS = [(entry['decision_point'], name, setting(entry['changes'])) for name, entry in A.items()]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, setting(changes)) for name, changes in R.items()]
