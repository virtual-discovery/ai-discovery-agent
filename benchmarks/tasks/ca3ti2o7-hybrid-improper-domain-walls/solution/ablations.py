"""Plausible wrong paths for ca3ti2o7-hybrid-improper-domain-walls; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/ca3ti2o7-hybrid-improper-domain-walls/make_task.py from relaxations
run with the lab's own NEP code path (construction runs/*.json).
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


ABLATIONS = [
    ('dp1', 'swap_psi_phi_construction', setting(A['swap_psi_phi_construction'])),
    ('dp2', 'antiphase_boundary_as_psi_wall', setting(A['antiphase_boundary_as_psi_wall'])),
    ('dp3', 'polar_axis_swapped', setting(A['polar_axis_swapped'])),
    ('dp4', 'layer_wall_in_bilayer', setting(A['layer_wall_in_bilayer'])),
    ('dp5', 'charged_slope_from_short_cells', setting(A['charged_slope_from_short_cells'])),
    ('dp5', 'charged_wall_treated_as_converged', setting(A['charged_wall_treated_as_converged'])),
    ('dp5', 'slope_per_supercell_length', setting(A['slope_per_supercell_length'])),
    ('dp6', 'short_range_NEP', setting(A['short_range_NEP'])),
    ('dp7', 'one_wall_per_cell_normalization', setting(A['one_wall_per_cell_normalization'])),
    ('dp8', 'profile_from_polarization', setting(A['profile_from_polarization'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, setting(vals)) for name, vals in R.items()]
