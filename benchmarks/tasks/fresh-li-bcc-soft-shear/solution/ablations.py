"""Plausible wrong paths for fresh-li-bcc-soft-shear; each must fail grading.

Every value comes from solution/ablation_values.json, computed by
benchmarks/construction/fresh-li-bcc-soft-shear/make_task.py from real calculations made through the task's own
lab code (task_lab.run), not from arbitrary perturbations. Each entry is a complete submission as that path would
produce it (fields the path does not affect keep the reference analysis).
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
FIELDS = ['a0_A', 'bulk_modulus_GPa', 'Cprime_GPa', 'C44_GPa', 'zener_anisotropy', 'mechanically_stable_at_a0',
          'compressed_pressure_GPa', 'compressed_Bprime_GPa', 'compressed_B44_GPa']


def use(values, fields=FIELDS):
    """Replace the listed submission fields by the values a wrong (or alternative) path produced."""
    def mutate(payload, trace):
        for field in fields:
            payload[field] = values[field]
        return payload, trace
    return mutate


ABLATIONS = [(dp, name, use(A[name])) for dp, name in V['ablation_decision_points']]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, use(values)) for name, values in R.items()]
