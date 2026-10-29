"""Plausible wrong paths for cr2o3-magnetoelectric-sign; each must fail grading.

Each entry replaces the graded fields with the values the wrong path produces from the archived data. The values
come from solution/ablation_values.json, computed by
benchmarks/construction/cr2o3-magnetoelectric-sign/make_task.py.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def use(values):
    def mutate(payload, trace):
        payload.update(json.loads(json.dumps(values)))
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'phonon_eigenvectors_as_modes', use(A['phonon_eigenvectors_as_modes'])),
    ('dp1', 'hdf5_values_as_THz', use(A['hdf5_values_as_THz'])),
    ('dp2', 'elk_domain_from_input_field', use(A['elk_domain_from_input_field'])),
    ('dp3', 'axial_second_pair_as_neighbours', use(A['axial_second_pair_as_neighbours'])),
    ('dp4', 'codes_assumed_same_domain_averaged', use(A['codes_assumed_same_domain_averaged'])),
    ('dp5', 'born_from_mislabeled_0005_pair', use(A['born_from_mislabeled_0005_pair'])),
    ('dp6', 'born_tensor_transposed', use(A['born_tensor_transposed'])),
    ('dp6', 'born_species_averaged', use(A['born_species_averaged'])),
    ('dp11', 'born_one_sided_no_asr', use(A['born_one_sided_no_asr'])),
    ('dp7', 'elk_eu3_from_003_pair', use(A['elk_eu3_from_003_pair'])),
    ('dp8', 'muffin_tin_moments', use(A['muffin_tin_moments'])),
    ('dp8', 'vasp_ion_projected_moments', use(A['vasp_ion_projected_moments'])),
    ('dp9', 'hexagonal_cell_volume', use(A['hexagonal_cell_volume'])),
    ('dp10', 'annealing_ignores_sign_change', use(A['annealing_ignores_sign_change'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('vasp_data_with_domain_flip', use(R['vasp_flipped'])),
    ('elk_window_and_born_005', use(R['elk_window_010_015_born_005'])),
    ('elk_vasp_mean', use(R['elk_vasp_mean'])),
    ('born_one_sided_002', use(R['born_one_sided_002'])),
]
