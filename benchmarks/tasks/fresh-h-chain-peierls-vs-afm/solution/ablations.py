"""Plausible wrong paths for fresh-h-chain-peierls-vs-afm; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/fresh-h-chain-peierls-vs-afm/make_task.py from the reference pw.x calculations
(the same lab code path the agent uses).
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(changes):
    def mutate(payload, trace):
        for field, values in changes.items():
            payload[field] = dict(payload[field], **values)
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'ferromagnetic_start_only', setting(A['ferromagnetic_start_only'])),
    ('dp1', 'one_atom_cell_spin_polarized', setting(A['one_atom_cell_spin_polarized'])),
    ('dp2', 'large_smearing_mv002', setting(A['large_smearing_mv002'])),
    ('dp2', 'fermi_dirac_0005', setting(A['fermi_dirac_0005'])),
    ('dp3', 'coarse_k_sdw_gap', setting(A['coarse_k_sdw_gap'])),
    ('dp4', 'coarse_boundary_grid', setting(A['coarse_boundary_grid'])),
    ('dp4', 'linear_moment_extrapolation', setting(A['linear_moment_extrapolation'])),
    ('dp5', 'nonmagnetic_uniform_reference', setting(A['nonmagnetic_uniform_reference'])),
    ('dp6', 'afm_wins_at_large_spacing', setting(A['afm_wins_at_large_spacing'])),
    ('dp7', 'grid_minimum_dimerization', setting(A['grid_minimum_dimerization'])),
    ('dp8', 'energies_per_cell', setting(A['energies_per_cell'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, setting(changes)) for name, changes in R.items()]
