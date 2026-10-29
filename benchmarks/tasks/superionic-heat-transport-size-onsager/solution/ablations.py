"""Plausible wrong paths for superionic-heat-transport-size-onsager; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/superionic-heat-transport-size-onsager/make_task.py from the released (served-format) runs.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def path(entry):
    def mutate(payload, trace):
        payload['kappa_predictions_W_mK'] = dict(entry['pred'])
        if 'runs' in entry:
            payload['kappa_runs_W_mK'] = dict(entry['runs'])
        if 'bias' in entry:
            payload['thermostat_biases_kappa'] = entry['bias']
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'raw_heat_flux_single_component_gk', path(A['plain_flux'])),
    ('dp2', 'smallest_96_atom_cells', path(A['smallest_cells'])),
    ('dp3', 'spectrum_normalisation_factor_two', path(A['factor_two'])),
    ('dp4', 'phonon_gas_one_over_T_extrapolation', path(A['phonon_gas'])),
    ('dp5', 'uo2_324_atom_cell_assumed_converged', path(A['cell3_converged'])),
    ('dp6', 'thermostat_bias_assumed_nvt_discarded', path(A['thermostat_bias'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('cepstral_lowpass_3THz', path(R['lowpass_3THz'])),
    ('nve_runs_only', path(R['nve_only'])),
    ('single_neighbour_temperature', path(R['single_neighbour'])),
    ('linear_extrapolation_both_ensembles', path(R['linear_extrapolation_both_ensembles'])),
]
