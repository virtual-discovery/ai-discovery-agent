"""Plausible wrong paths for cu2wse4-orbital-edelstein; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/cu2wse4-orbital-edelstein/make_task.py from the lab's own cached evaluations.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(changes):
    def mutate(payload, trace):
        for field, value in changes.items():
            payload[field] = value
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'rotated_axes_setting', setting(A['rotated_axes_setting'])),
    ('dp1', 'chiral_d2_assumption', setting(A['chiral_d2_assumption'])),
    ('dp2', 'intra_atomic_orbital_moment', setting(A['intra_atomic_orbital_moment'])),
    ('dp2', 'berry_plus_intra_atomic_moment', setting(A['berry_plus_intra_atomic_moment'])),
    ('dp3', 'no_spin_factor_in_nosoc_alpha', setting(A['no_spin_factor_in_nosoc_alpha'])),
    ('dp4', 'no_spin_factor_in_nosoc_holes', setting(A['no_spin_factor_in_nosoc_holes'])),
    ('dp4', 'no_spin_factor_anywhere_nosoc', setting(A['no_spin_factor_anywhere_nosoc'])),
    ('dp5', 'common_header_energy_axis', setting(A['common_header_energy_axis'])),
    ('dp6', 'coarse_mesh_32', setting(A['coarse_mesh_32'])),
    ('dp7', 'tau_hbar_over_gamma', setting(A['tau_hbar_over_gamma'])),
    ('dp7', 'positive_charge_sign', setting(A['positive_charge_sign'])),
    ('dp8', 'conduction_ratio_from_valence', setting(A['conduction_ratio_from_valence'])),
    ('dp9', 'nlh_full_permutation_symmetry', setting(A['nlh_full_permutation_symmetry'])),
    ('dp10', 'chirality_from_chiral_class', setting(A['chirality_from_chiral_class'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, setting(changes)) for name, changes in R.items()]
