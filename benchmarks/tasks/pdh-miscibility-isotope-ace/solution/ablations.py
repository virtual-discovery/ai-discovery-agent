"""Plausible wrong paths for pdh-miscibility-isotope-ace; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/pdh-miscibility-isotope-ace/make_task.py from the reference calculations.
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
    ('dp1', 'fixed_volume_configurational', setting(A['fixed_volume_configurational'])),
    ('dp2', 'no_coherency', setting(A['no_coherency'])),
    ('dp2', 'coherency_111', setting(A['coherency_111'])),
    ('dp3', 'small_cell_dilute', setting(A['small_cell_dilute'])),
    ('dp3', 'small_cell_fixed_volume_dilute', setting(A['small_cell_fixed_volume_dilute'])),
    ('dp4', 'harmonic_zpe', setting(A['harmonic_zpe'])),
    ('dp5', 'zpe_only_no_thermal', setting(A['zpe_only_no_thermal'])),
    ('dp6', 'static_site_preference', setting(A['static_site_preference'])),
    ('dp6', 'tetrahedral_site_as_dilute', setting(A['tetrahedral_site_as_dilute'])),
    ('dp7', 'elastic_only_wagner_horner', setting(A['elastic_only_wagner_horner'])),
    ('dp8', 'ordered_small_cells', setting(A['ordered_small_cells'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, setting(values)) for name, values in R.items()]
