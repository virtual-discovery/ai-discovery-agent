"""Plausible wrong paths for fresh-gete-ferroelectric-well; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/fresh-gete-ferroelectric-well/make_task.py from the reference pw.x calculations (the same
lab code path the agent uses) with the estimators of solve.py.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(changes):
    changes = {k: v for k, v in changes.items() if k != 'note'}

    def mutate(payload, trace):
        payload.update(changes)
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'clamped_reported_as_well_depth', setting(A['clamped_reported_as_well_depth'])),
    ('dp2', 'cubic_at_r3m_volume', setting(A['cubic_at_r3m_volume'])),
    ('dp2', 'x_relaxed_in_experimental_cell', setting(A['x_relaxed_in_experimental_cell'])),
    ('dp3', 'k12_mesh_both_phases', setting(A['k12_mesh_both_phases'])),
    ('dp4', 'one_sided_quadratic_model', setting(A['one_sided_quadratic_model'])),
    ('dp5', 'per_atom_energies', setting(A['per_atom_energies'])),
    ('dp6', 'fixed_occupations_soft_mode_k16', setting(A['fixed_occupations_soft_mode_k16'])),
    ('dp7', 'soft_mode_k12_mv0.01', setting(A['soft_mode_k12_mv0.01'])),
    ('dp8', 'soft_mode_from_dx0.02', setting(A['soft_mode_from_dx0.02'])),
    ('dp9', 'soc_compared_with_sssp', setting(A['soc_compared_with_sssp'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, setting(values)) for name, values in R.items()]
