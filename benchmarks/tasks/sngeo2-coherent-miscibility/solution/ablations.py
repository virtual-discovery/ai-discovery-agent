"""Plausible wrong paths for sngeo2-coherent-miscibility; each must fail grading.

Each wrong path is the complete set of graded numbers that analysis would produce, computed from the
released data by benchmarks/construction/sngeo2-coherent-miscibility/make_task.py (model.analysis with
the corresponding choice changed) and stored in solution/ablation_values.json.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
GRADED = ['mixing_enthalpy_x050_meV_per_cation', 'coherency_strain_energy_x050_meV_per_cation', 'habit_plane',
          'incoherent_Tc_K', 'coherent_spinodal_Tc_K', 'coherent_boundaries_700K_PBE']


def path(values):
    def mutate(payload, trace):
        for field in GRADED:
            payload[field] = values[field]
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'sqs_values_read_per_cation', path(A['sqs_values_read_per_cation'])),
    ('dp1', 'strain_per_cell_mixed_with_per_cation', path(A['strain_per_cell_mixed_with_per_cation'])),
    ('dp2', 'r2scan_unrelaxed_sqs_kept', path(A['r2scan_unrelaxed_sqs_kept'])),
    ('dp2', 'r2scan_index0_sqs_unscreened', path(A['r2scan_index0_sqs_unscreened'])),
    ('dp3', 'inplane_c_not_constrained', path(A['inplane_c_not_constrained'])),
    ('dp4', 'no_out_of_plane_relaxation', path(A['no_out_of_plane_relaxation'])),
    ('dp5', 'habit_plane_100', path(A['habit_plane_100'])),
    ('dp6', 'incoherent_model_reported_as_coherent', path(A['incoherent_model_reported_as_coherent'])),
    ('dp7', 'binodal_spinodal_swapped', path(A['binodal_spinodal_swapped'])),
    ('dp8', 'functional_mixing', path(A['functional_mixing'])),
    ('dp9', 'overtruncated_cluster_expansion_Tc', path(A['overtruncated_cluster_expansion_Tc'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, path(values)) for name, values in R.items()]
