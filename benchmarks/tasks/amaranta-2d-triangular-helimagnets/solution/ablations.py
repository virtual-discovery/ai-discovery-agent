"""Plausible wrong paths for amaranta-2d-triangular-helimagnets; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/amaranta-2d-triangular-helimagnets/make_task.py from the teacher model
(the authors' archived AMaRaNTA parameters).
"""
import copy
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def replace(ground_states=None, predictions=None, axis=None):
    def mutate(payload, trace):
        payload = copy.deepcopy(payload)
        if ground_states:
            payload['ground_states'].update(copy.deepcopy(ground_states))
        if predictions:
            payload['energy_predictions_meV_per_site'].update(predictions)
        if axis is not None:
            payload['VI2_moment_direction'] = list(axis)
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'six_equivalent_nearest_neighbours', replace(A['ideal_triangular']['ground_states'], A['ideal_triangular']['predictions'])),
    ('dp2', 'third_neighbour_aliased_in_4_periodic_cell', replace(A['J3_doubled']['ground_states'], A['J3_doubled']['predictions'])),
    ('dp3', 'exchange_truncated_after_J2', replace(A['no_J3']['ground_states'], A['no_J3']['predictions'])),
    ('dp4', 'only_FM_120_stripe_compared', replace(A['commensurate_only']['ground_states'])),
    ('dp5', 'VI2_axis_from_cartesian_axes_or_SIA', replace(axis=A['VI2_axis_cartesian_or_SIA_only'])),
    ('dp6', 'VI2_tensor_copied_to_all_bonds', replace(axis=A['VI2_axis_tensor_copied'])),
    ('dp7', 'ordered_pairs_double_counted', replace(predictions=A['pairs_double_counted'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('isotropic_model_predictions', replace(predictions=R['isotropic_predictions'])),
    ('anisotropic_luttinger_tisza_q', replace(ground_states=R['anisotropic_LT_q'])),
    ('collinear_zz_mapping_q', replace(ground_states=R['collinear_zz_mapping_q'])),
    ('VI2_axis_from_six_noisy_probes', replace(axis=R['VI2_axis_six_probes'])),
]
