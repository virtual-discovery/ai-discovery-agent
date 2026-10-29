"""Plausible wrong paths for 2d-qsh-shc-plateau; each must fail grading.

Values come from solution/ablation_values.json, computed by benchmarks/construction/2d-qsh-shc-plateau
(reference.py with the lab's own evaluator on the lab's models, and make_task.py) from the archived data.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(changes):
    def mutate(payload, trace):
        for field, value in changes.items():
            payload[field] = {**payload[field], **value}
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'plateau_at_scf_fermi', setting(A['plateau_at_scf_fermi'])),
    ('dp1', 'plateau_at_win_fermi_energy', setting(A['plateau_at_win_fermi_energy'])),
    ('dp2', 'plateau_at_dft_mid_gap', setting(A['plateau_at_dft_mid_gap'])),
    ('dp3', 'mesh_120_unrefined', setting(A['mesh_120_unrefined'])),
    ('dp3', 'mesh_60_unrefined', setting(A['mesh_60_unrefined'])),
    ('dp4', 'lattice_position_convention', setting(A['lattice_position_convention'])),
    ('dp5', 'relaxed_archive_values', setting(A['relaxed_archive_values'])),
    ('dp5', 'paper_table_values', setting(A['paper_table_values'])),
    ('dp6', 'z2_from_plateau_magnitude', setting(A['z2_from_plateau_magnitude'])),
    ('dp7', 'z2_from_coarse_wcc', setting(A['z2_from_coarse_wcc'])),
    ('dp8', 'atomic_thickness_normalisation', setting(A['atomic_thickness_normalisation'])),
    ('dp8', 'e_over_4pi_units', setting(A['e_over_4pi_units'])),
    ('dp9', 'sigma_z_only_spin_direction', setting(A['sigma_z_only_spin_direction'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, setting(changes)) for name, changes in R.items()]
