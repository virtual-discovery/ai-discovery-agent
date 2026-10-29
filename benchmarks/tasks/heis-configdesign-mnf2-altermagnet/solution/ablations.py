"""Plausible wrong paths for heis-configdesign-mnf2-altermagnet; each must fail grading.

Values come from solution/ablation_values.json, computed from the released data by
benchmarks/construction/heis-configdesign-mnf2-altermagnet/make_task.py.
"""
import copy
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(**changes):
    def mutate(payload, trace):
        for field, value in changes.items():
            payload[field] = copy.deepcopy(value)
        return payload, trace
    return mutate


def model(entry, alt_mean=None, tn=None):
    """Replace the fitted-model outputs (predictions, J1-J3, and optionally T_N and the split-shell mean)."""
    def mutate(payload, trace):
        payload['spiral_predictions_meV_per_cell'] = dict(entry['predictions'])
        payload['J1_meV'], payload['J2_meV'], payload['J3_meV'] = entry['J']
        if tn is not None:
            payload['neel_temperature_K'] = tn
        if alt_mean is not None:
            payload['altermagnetic_shell'] = dict(payload['altermagnetic_shell'], sublattice_mean_J_meV=alt_mean)
        return payload, trace
    return mutate


def alt(entry):
    return setting(altermagnetic_shell={'distance_A': entry['distance_A'], 'sublattice_mean_J_meV': entry['mean_meV'],
                                        'sublattice_difference_J_meV': entry['difference_meV']})


ABLATIONS = [
    ('dp1', 'per_atom_label', model(A['per_atom_label'], A['per_atom_label']['split_mean'], A['per_atom_label']['TN'])),
    ('dp1', 'per_mn_label', model(A['per_mn_label'], A['per_mn_label']['split_mean'], A['per_mn_label']['TN'])),
    ('dp2', 'include_truncated_214', model(A['include_truncated_214'])),
    ('dp2', 'dense_214_with_truncated', model(A['dense_214_with_truncated'])),
    ('dp3', 'no_214_ten_shells', model(A['no_214_ten_shells'])),
    ('dp3', '001_only', model(A['001_only'])),
    ('dp4', 'coarse_commensurate_grid', model(A['coarse_commensurate_grid'])),
    ('dp5', 'six_shells', model(A['6_shells'], A['6_shells']['split_mean'])),
    ('dp5', 'three_shells', model(A['3_shells'], A['3_shells']['split_mean'])),
    ('dp6', 'neel_reference', model(A['neel_reference'])),
    ('dp6', 'lattice_phase_only', model(A['lattice_phase_only'])),
    ('dp7', 'TN_mean_field', setting(neel_temperature_K=A['TN_mean_field'])),
    ('dp8', 'alt_zero_difference', alt(A['alt_zero_difference'])),
    ('dp8', 'alt_orthorhombic_fit', alt(A['alt_orthorhombic_fit'])),
    ('dp8', 'alt_wrong_shell', alt(A['alt_wrong_shell'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('twelve_shell_model', model(R['12_shells'], R['12_shells']['split_mean'])),
    ('eight_shell_model', model(R['8_shells'], R['8_shells']['split_mean'])),
    ('all_clean_available_spirals', model(R['all_clean_available'], R['all_clean_available']['split_mean'])),
    ('TN_coarse_scan_L8_L12', setting(neel_temperature_K=R['TN_coarse_L8_L12'])),
]
