"""Plausible wrong paths for gan-hole-piezomobility; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/gan-hole-piezomobility/make_task.py from the released data.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
PRED = 'mobility_predictions_cm2_per_Vs'


def setting(**changes):
    def mutate(payload, trace):
        for field, value in changes.items():
            payload[field] = value
        return payload, trace
    return mutate


def best(value, mu_key):
    return {'strain_voigt': value['strain_voigt'], mu_key: value[mu_key]}


ABLATIONS = [
    ('dp1', 'W4_same_sign_series', setting(**{PRED: A['W4_same_sign_series']})),
    ('dp1', 'shear_offdiag_sign_not_flipped', setting(**{PRED: A['shear_offdiag_sign_not_flipped']})),
    ('dp2', 'W3_swap_from_e1', setting(**{PRED: A['W3_swap_from_e1']})),
    ('dp3', 'linear_piezomobility_from_1pct', setting(**{PRED: A['linear_piezomobility_from_1pct']})),
    ('dp4', 'W1_precrossing_extrapolation', setting(**{PRED: A['W1_precrossing_extrapolation']})),
    ('dp5', 'best_mu_zz_from_serta', setting(best_mu_zz_2pct=best(A['best_mu_zz_from_serta'], 'mu_zz_cm2_per_Vs'))),
    ('dp6', 'best_inplane_diagonal_only',
     setting(best_inplane_2pct=best(A['best_inplane_diagonal_only'], 'mu_inplane_max_cm2_per_Vs'))),
    ('dp7', 'best_inplane_label_as_voigt',
     setting(best_inplane_2pct=best(A['best_inplane_label_as_voigt'], 'mu_inplane_max_cm2_per_Vs'))),
    ('dp8', 'flagged_runs_replaced_by_scaled_serta', setting(**{PRED: A['flagged_runs_replaced_by_scaled_serta']})),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('W3_linear_extrapolation_of_e2_series', setting(**{PRED: R['W3_linear_series']})),
    ('W3_anisotropy_ratio_extrapolation', setting(**{PRED: R['W3_anisotropy_ratio']})),
    ('W5_linear_extrapolation_of_compress_e6_series', setting(**{PRED: R['W5_linear_compress_series']})),
    ('W1_lower_bound_saturated_at_2pct', setting(**{PRED: R['W1_saturated_at_2pct']})),
    ('W1_upper_bound_linear_1_to_2pct', setting(**{PRED: R['W1_linear_1_to_2pct']})),
    ('W2_flat_from_expand_1pct', setting(**{PRED: R['W2_flat_from_1pct']})),
    ('best_inplane_from_expand_2pct_e6',
     setting(best_inplane_2pct=best(R['best_inplane_from_expand_e6'], 'mu_inplane_max_cm2_per_Vs'))),
]
