"""Plausible wrong paths for li3ps4-superionic-transport; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/li3ps4-superionic-transport/make_task.py from the archived data.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(**changes):
    def mutate(payload, trace):
        for field, value in changes.items():
            if isinstance(value, dict) and isinstance(payload.get(field), dict):
                payload[field] = {**payload[field], **value}
            else:
                payload[field] = value
        return payload, trace
    return mutate


def ea(below=None, above=None):
    out = {}
    if below is not None:
        out['below_transition'] = below
    if above is not None:
        out['superionic_to_1000K'] = above
    return out


ABLATIONS = [
    ('dp1', 'nernst_einstein_as_conductivity',
     setting(sigma_predictions_S_per_m=A['sigma_tracer'], haven_ratio_PBE0_gamma_750K=A['haven_one'],
             activation_energy_eV_PBE0_gamma=ea(A['ea_tracer']['below_transition'], A['ea_tracer']['superionic_to_1000K']))),
    ('dp2', 'green_kubo_noise_floor_in_low_T_fit',
     setting(activation_energy_eV_PBE0_gamma=ea(below=A['ea_noise_floor']['below_transition']))),
    ('dp3', 'largest_heat_capacity_peak_as_transition',
     setting(transition_temperature_K=A['tc_melting_peak'], sigma_predictions_S_per_m=A['sigma_low_branch'])),
    ('dp4', 'pbesol_transition_for_all_references', setting(transition_temperature_K=A['tc_pbesol_for_all'])),
    ('dp4', 'experimental_alpha_transition_for_all_references', setting(transition_temperature_K=A['tc_experimental_alpha'])),
    ('dp5', 'interpolate_across_the_jump', setting(sigma_predictions_S_per_m=A['sigma_interpolated_across_jump'])),
    ('dp6', 'alpha_started_runs_as_superionic_reference',
     setting(sigma_predictions_S_per_m=A['sigma_alpha_branch'], haven_ratio_PBE0_gamma_750K=A['haven_alpha_branch'])),
    ('dp7', 'barycentric_li_s_term_as_coupling', setting(li_s_cross_fraction=A['li_s_barycentric'])),
    ('dp8', 'librations_taken_for_reorientations',
     setting(log10_hop_to_flip_rate_ratio_PBE0_800K=A['log_ratio_libration'], paddle_wheel_supported=True)),
    ('dp9', 'single_arrhenius_through_the_transition',
     setting(activation_energy_eV_PBE0_gamma=ea(A['ea_single_arrhenius'], A['ea_single_arrhenius']))),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('tracer_branch_over_haven_ratio',
     setting(sigma_predictions_S_per_m=R['sigma_tracer_branch_over_haven'], haven_ratio_PBE0_gamma_750K=R['haven_first_three'],
             activation_energy_eV_PBE0_gamma=ea(below=R['ea_below_tracer']))),
    ('superionic_fit_to_1000K_haven_to_1100K_window_to_875K',
     setting(sigma_predictions_S_per_m=R['sigma_superionic_gk_fit'], haven_ratio_PBE0_gamma_750K=R['haven_all_superionic_to_1100K'],
             activation_energy_eV_PBE0_gamma=ea(R['ea_below_unweighted'], R['ea_above_to_875K']))),
    ('upper_bracket_domega_as_rate_zero_coupling',
     setting(transition_temperature_K=R['tc_upper_bracket'], log10_hop_to_flip_rate_ratio_PBE0_800K=R['log_ratio_domega_as_rate'],
             li_s_cross_fraction=0.0)),
    ('lower_bracket_msd_route',
     setting(transition_temperature_K=R['tc_lower_bracket'],
             sigma_predictions_S_per_m={'PBE0/gamma/750': R['sigma_PBE0_750_msd_route']},
             log10_hop_to_flip_rate_ratio_PBE0_800K=R['log_ratio_msd_750_hop3A'],
             activation_energy_eV_PBE0_gamma=ea(above=R['ea_above_without_750K']))),
]
