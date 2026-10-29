"""Plausible wrong paths for mnwo4-halfmoon-frustration; each must fail grading.

Values come from solution/ablation_values.json, computed from the archived data by
benchmarks/construction/mnwo4-halfmoon-frustration/calibrate.py and make_task.py.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(values, **extra):
    def mutate(payload, trace):
        for field, value in values.items():
            payload[field] = value
        for field, value in extra.items():
            payload[field] = value
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'dispersion_only_min_norm_split', setting(A['abl_dispersion_only_min_norm'])),
    ('dp2', 'spin_length_one_in_dynamics', setting(A['abl_spin_length_one'])),
    ('dp3', 'mean_field_rpa_nominal_T', setting(A['abl_rpa_nominal_T'])),
    ('dp3', 'mean_field_rpa_lowest_valid_T', setting(A['abl_rpa_lowest_valid_T'])),
    ('dp4', 'isotropic_model_no_anisotropy', setting(A['abl_isotropic'])),
    ('dp5', 'magnetic_form_factor_applied', setting(A['abl_form_factor'])),
    ('dp6', 'temperature_insensitive_spin_liquid', setting(A['abl_temperature_insensitive'],
                                                           half_moon_origin='spin_liquid_local_constraint')),
    ('dp7', 'quasi_1d_chain_model', setting(A['abl_chain_only'], half_moon_origin='quasi_1d_zigzag_chain')),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('archived_hamiltonian_scga', setting(R['reference_model'])),
    ('34_spectra_null_from_20K_map_only', setting(R['robust_wide_spectra_20K_only'])),
    ('null_from_20_30_60K_maps', setting(R['robust_three_maps'])),
    ('null_from_11K_zone_peak', setting(R['robust_zone_peak_null'])),
    ('no_literature_joint_de_best_of_4_seeds', setting(R['robust_nonrecall_best_of_seeds'])),
]
