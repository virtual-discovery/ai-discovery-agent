"""Plausible wrong paths for polar-phonons-1d-loto-collapse; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/polar-phonons-1d-loto-collapse/make_task.py from the archived data.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, A3, R, R3 = V['ablations'], V['ablations_3d'], V['robustness'], V['robustness_3d']


def tubes(name):
    def mutate(payload, trace):
        payload['BN55_LO_1D_cm-1'] = A[name]['BN55']
        payload['BN66_LO_1D_cm-1'] = A[name]['BN66']
        return payload, trace
    return mutate


def three_d(value):
    def mutate(payload, trace):
        payload['BN55_LO_3D_q0_cm-1'] = value
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'three_d_boundary_condition_dispersion', tubes('three_d_boundary_condition_dispersion')),
    ('dp1', 'analytic_gamma_frequency_as_3d_limit', three_d(A3['abl_analytic_gamma_frequency'])),
    ('dp2', 'linear_2d_collapse', tubes('linear_2d_collapse')),
    ('dp3', 'standard_fourier_interpolation', tubes('standard_fourier_interpolation')),
    ('dp4', 'atomic_shell_radius', tubes('atomic_shell_radius')),
    ('dp4', 'cell_dielectric_as_cylinder_dielectric', tubes('cell_dielectric_as_cylinder_dielectric')),
    ('dp5', 'asymptotic_kernel_in_interpolation', tubes('asymptotic_kernel_in_interpolation')),
    ('dp6', 'gamma_only_long_range_model', tubes('gamma_only_long_range_model')),
    ('dp7', 'q_fraction_of_2pi_over_L', tubes('q_fraction_of_2pi_over_L')),
    ('dp7', 'e2_hartree_in_long_range', tubes('e2_hartree_in_long_range')),
    ('dp8', 'cylinder_dielectric_for_3d_artefact', three_d(A3['abl_cylinder_dielectric'])),
    ('dp8', 'cylinder_dielectric_and_volume_for_3d_artefact', three_d(A3['abl_cylinder_dielectric_and_volume'])),
    ('dp8', 'perpendicular_dielectric_for_3d_artefact', three_d(A3['abl_perpendicular_dielectric'])),
    ('dp8', 'unscreened_3d_artefact', three_d(A3['abl_unscreened'])),
    ('dp9', 'family_transfer_from_44', tubes('family_transfer_from_44')),
]


def robust(name):
    def mutate(payload, trace):
        payload['BN55_LO_1D_cm-1'] = R[name]['BN55']
        payload['BN66_LO_1D_cm-1'] = R[name]['BN66']
        return payload, trace
    return mutate


# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('authors_radius_d_over_4', robust('authors_radius_d_over_4')),
    ('no_ewald_damping', robust('no_ewald_damping')),
    ('ewald_alpha_2', robust('ewald_alpha_2')),
    ('radius_1p5_atomic', robust('radius_1p5_atomic')),
    ('3d_limit_from_1d_code_gamma_run', three_d(R3['1D_code_gamma'])),
]
