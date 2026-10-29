"""Plausible wrong paths for fluorite-superionic-finite-size; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/fluorite-superionic-finite-size/make_task.py from the served (released) data.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
D, TC, EA = 'diffusion_predictions_A2_per_ps', 'transition_temperature_K', 'activation_energy_eV'


def setting(**changes):
    """Replace whole fields, or merge into a dict field when the value is given as ('merge', dict)."""
    def mutate(payload, trace):
        for field, value in changes.items():
            field = {'D': D, 'TC': TC, 'EA': EA}[field]
            if isinstance(value, tuple):
                payload[field] = {**payload[field], **value[1]}
            else:
                payload[field] = value
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'uniform_frame_time_for_UO2', setting(D=A['uniform_frame_time'])),
    ('dp2', 'small_cell_transition_temperatures', setting(TC=A['small_cell_Tc'])),
    ('dp2', 'small_cell_activation_energies', setting(EA=A['small_cell_Ea'])),
    ('dp2', 'small_cell_predictions', setting(D=A['small_cell_predictions'])),
    ('dp3', '768_atom_cell_assumed_converged', setting(D=A['768_cell_assumed_converged'])),
    ('dp4', 'yeh_hummer_one_over_L_extrapolation', setting(D=A['yeh_hummer_1_over_L'])),
    ('dp5', 'single_arrhenius_across_Tc', setting(EA=A['single_arrhenius'])),
    ('dp6', 'uo2_324_atom_cell_transition', setting(TC=('merge', {'UO2': A['uo2_324_cell_Tc']}))),
    ('dp6', 'uo2_324_atom_cell_predictions', setting(D=A['uo2_324_cell_predictions'])),
    ('dp7', 'fluorite_size_correction_applied_to_AgI', setting(D=A['agi_fluorite_correction'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('msd_fit_window_10_80_ps', setting(EA=R['window_10_80']['activation_energy_eV'], D=R['window_10_80']['predictions'])),
    ('Tc_from_energy_fluctuations_with_cascaded_windows',
     setting(TC=R['tc_from_fluctuations']['transition_temperature_K'], EA=R['tc_from_fluctuations']['activation_energy_eV'])),
    ('PbF2_large_cell_by_one_over_N_extrapolation', setting(D=('merge', R['pbf2_one_over_N']))),
    ('Arrhenius_fit_instead_of_interpolation', setting(D=('merge', R['arrhenius_fit_predictions']))),
    ('AgI_from_32_atom_cell_in_momentum_frame', setting(D=('merge', R['self_frame']))),
]
