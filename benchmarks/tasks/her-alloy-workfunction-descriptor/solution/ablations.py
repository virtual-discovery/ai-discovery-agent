"""Plausible wrong paths for her-alloy-workfunction-descriptor; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/her-alloy-workfunction-descriptor/make_task.py from the released data
(same purchases as the reference solution unless stated).
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
FIELD = {'lib3': 'lib3_log10_abs_current', 'volcano_centre_eV': 'volcano_centre_eV',
         'coverage_shift_eV': 'coverage_shift_eV', 'alpha_wf': 'alpha_wf'}


def use(values):
    def mutate(payload, trace):
        for key, value in values.items():
            payload[FIELD[key]] = value
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'coverage_total_slope', use(A['coverage_total_slope'])),
    ('dp1', 'coverage_first_last', use(A['coverage_first_last'])),
    ('dp2', 'literature_centre', use(A['literature_centre'])),
    ('dp3', 'mean_energy_descriptor', use(A['mean_energy_descriptor'])),
    ('dp4', 'no_transport_cap', use(A['no_transport_cap'])),
    ('dp5', 'no_library_calibration', use(A['no_library_calibration'])),
    ('dp5', 'prefactor_only_calibration', use(A['prefactor_only_calibration'])),
    ('dp6', 'wf_only_coefficient', use(A['wf_only_coefficient'])),
    ('dp7', 'wf_with_dilute_binding', use(A['wf_with_dilute_binding'])),
    ('dp8', 'pooled_without_library_offsets', use(A['pooled_without_library_offsets'])),
    ('dp9', 'linear_combination_descriptor', use(A['linear_combination_descriptor'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, use(values)) for name, values in R.items()]
