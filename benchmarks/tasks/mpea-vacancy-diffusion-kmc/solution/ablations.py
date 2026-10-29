"""Plausible wrong paths for mpea-vacancy-diffusion-kmc; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/mpea-vacancy-diffusion-kmc/calibrate.py from the archived hops with the lab kMC engine.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
GRADED = ['barrier_predictions_eV', 'log10_eta_vacancy', 'fastest_tracer_B', 'log10_tracer_ratio_B',
          'activation_energy_C_eV', 'log10_eta_thermo_ideal_B', 'log10_eta_kra_simple_B']


def setting(values):
    def mutate(payload, trace):
        for field, value in values.items():
            payload[field] = value
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'thermo_ideal_as_real', setting(A['thermo_ideal_as_real'])),
    ('dp1', 'direct_forward_barrier_fit', setting(A['direct_forward_barrier_fit'])),
    ('dp2', 'species_mean_barriers', setting(A['species_mean_barriers'])),
    ('dp3', 'jump_count_tracers', setting(A['jump_count_tracers'])),
    ('dp4', 'rom_from_linear_model_pure_limit', setting(A['rom_from_linear_model_pure_limit'])),
    ('dp5', 'ea_from_rom_barrier', setting(A['ea_from_rom_barrier'])),
    ('dp5', 'ea_from_mean_accepted_barrier', setting(A['ea_from_mean_accepted_barrier'])),
    ('dp6', 'untested_attribution', setting(A['untested_attribution'])),
]

# Defensible alternative analyses (all graded fields replaced); rounding to 4 significant digits is checked automatically.
ROBUSTNESS = [(name, setting({k: v for k, v in values.items() if k in GRADED})) for name, values in R.items()]
