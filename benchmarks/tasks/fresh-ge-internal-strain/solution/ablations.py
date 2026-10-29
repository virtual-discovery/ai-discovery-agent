"""Plausible wrong paths for fresh-ge-internal-strain; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/fresh-ge-internal-strain/make_task.py from the reference pw.x calculations
(the same lab code path the agent uses) with the estimators of solve.py.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(changes):
    def mutate(payload, trace):
        payload.update(changes)
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'gap_from_shifted_mesh_band_edges', setting(A['gap_from_shifted_mesh_band_edges'])),
    ('dp1', 'gap_from_gamma_mesh_band_edges', setting(A['gap_from_gamma_mesh_band_edges'])),
    ('dp2', 'pressures_at_35_Ry', setting(A['pressures_at_35_Ry'])),
    ('dp2', 'experimental_lattice_parameter', setting(A['experimental_lattice_parameter'])),
    ('dp3', 'frozen_3d_core_pseudopotential', setting(A['frozen_3d_core_pseudopotential'])),
    ('dp4', 'shifted_k8_mesh', setting(A['shifted_k8_mesh'])),
    ('dp4', 'shifted_k8_relaxation', setting(A['shifted_k8_relaxation'])),
    ('dp5', 'wide_smearing_mv002', setting(A['wide_smearing_mv002'])),
    ('dp5', 'smearing_0.01_gamma_k12', setting(A['smearing_0.01_gamma_k12'])),
    ('dp6', 'clamped_reported_as_relaxed', setting(A['clamped_reported_as_relaxed'])),
    ('dp6', 'loose_relaxation_at_small_strain', setting(A['loose_relaxation_at_small_strain'])),
    ('dp7', 'one_sided_111_displacement_0.02', setting(A['one_sided_111_displacement_0.02'])),
    ('dp7', 'one_sided_111_displacement_0.05', setting(A['one_sided_111_displacement_0.05'])),
    ('dp8', 'energy_route_ignoring_prestress', setting(A['energy_route_ignoring_prestress'])),
    ('dp8', 'one_sided_normal_strain_0.02', setting(A['one_sided_normal_strain_0.02'])),
    ('dp10', 'gruneisen_from_1pct_steps', setting(A['gruneisen_from_1pct_steps'])),
    ('dp10', 'gruneisen_from_1pct_steps_k8', setting(A['gruneisen_from_1pct_steps_k8'])),
    ('dp9', 'atomic_mass_in_frequency', setting(A['atomic_mass_in_frequency'])),
    ('dp9', 'tensor_shear_strain_as_engineering', setting(A['tensor_shear_strain_as_engineering'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, setting(changes)) for name, changes in R.items()]
