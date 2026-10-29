"""Plausible wrong paths for inmnti2-weyl-census; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/inmnti2-weyl-census/make_task.py from the archived data with the lab's own evaluator.
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
    ('dp1', 'standard_model_as_magnetic', setting(A['standard_model_as_magnetic'])),
    ('dp2', 'band_index_from_exclude_bands', setting(A['band_index_from_exclude_bands'])),
    ('dp2', 'same_band_index_both_models', setting(A['same_band_index_both_models'])),
    ('dp3', 'wannier_nodes_trusted', setting(A['wannier_nodes_trusted'])),
    ('dp3', 'grid_validated_nearest_pair', setting(A['grid_validated_nearest_pair'])),
    ('dp4', 'cubic_single_orbit', setting(A['cubic_single_orbit'])),
    ('dp4', 'unitary_operations_only', setting(A['unitary_operations_only'])),
    ('dp4', 'planar_families', setting(A['planar_families'])),
    ('dp4', 'magnetization_axis_split', setting(A['magnetization_axis_split'])),
    ('dp5', 'chirality_of_upper_band', setting(A['chirality_of_upper_band'])),
    ('dp5', 'cubic_symmetric_node', setting(A['cubic_symmetric_node'])),
    ('dp6', 'kubo_folder_as_broadening', setting(A['kubo_folder_as_broadening'])),
    ('dp6', 'kubo_mesh_58', setting(A['kubo_mesh_58'])),
    ('dp7', 'raw_hall_components', setting(A['raw_hall_components'])),
    ('dp8', 'nscf_path_fermi_level', setting(A['nscf_path_fermi_level'])),
    ('dp9', 'pbe_count_with_avoided_crossings', setting(A['pbe_count_with_avoided_crossings'])),
    ('dp9', 'pbe_count_coarse_grid', setting(A['pbe_count_coarse_grid'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, setting(changes)) for name, changes in R.items()]
