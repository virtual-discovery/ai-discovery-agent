"""Plausible wrong paths for zrsix-butterfly-amr; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/zrsix-butterfly-amr/make_task.py through the lab's own evaluators.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
GRADED = ['neutral_mu_eV', 'carrier_density_cm3', 'butterfly_peak_angle_deg', 'mr_b_over_c', 'mr_exponent_b',
          'zrsite_mr_c_over_b', 'mr_ratio_se_over_te']


def setting(values):
    def mutate(payload, trace):
        for field, value in values.items():
            if field in GRADED:
                payload[field] = value
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'input_e_fermi', setting(A['input_e_fermi'])),
    ('dp2', 'zrsis_copied_input', setting(A['zrsis_copied_input'])),
    ('dp3', 'thirteen_filled_bands', setting(A['thirteen_filled_bands'])),
    ('dp3', 'carrier_density_one_spin', setting(A['carrier_density_one_spin'])),
    ('dp4', 'mesh_80', setting(A['mesh_80'])),
    ('dp4', 'mesh_40', setting(A['mesh_40'])),
    ('dp5', 'closed_orbit_mass_model', setting(A['closed_orbit_mass_model'])),
    ('dp6', 'samples_at_equal_btau', setting(A['samples_at_equal_btau'])),
    ('dp6', 'samples_at_equal_btau_swapped', setting(A['samples_at_equal_btau_swapped'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('mesh_120_potentials_and_densities', setting(R['mesh120_potentials'])),
    ('cubic_peak_finder', setting(R['cubic_peak'])),
    ('oracle_angle_subset', setting(R['oracle_angles'])),
]
