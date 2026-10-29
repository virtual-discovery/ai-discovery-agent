"""Plausible wrong paths for fresh-al-surface-workfunction; each must fail grading.

Every value comes from solution/ablation_values.json, computed by
benchmarks/construction/fresh-al-surface-workfunction/make_task.py from real construction runs made through the
lab's own QE code path (for example the default-band slabs, the k = 12 ladder, the primitive-cell bulk energies and
the total local potentials are actual calculations, not perturbations).
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
GRADED = ('surface_energy_J_per_m2', 'work_function_eV', 'work_function_ordering', 'relaxation_d12_percent_111',
          'thick_slab_predictions')


def use(values):
    def mutate(payload, trace):
        for field in GRADED:
            payload[field] = json.loads(json.dumps(values[field]))
        return payload, trace
    return mutate


def ordering(value):
    def mutate(payload, trace):
        payload['work_function_ordering'] = value
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'default_bands_ladder', use(A['default_bands_ladder'])),
    ('dp1', 'default_bands_slope_fit', use(A['default_bands_slope_fit'])),
    ('dp2', 'primitive_bulk_reference', use(A['primitive_bulk_reference'])),
    ('dp3', 'slope_fit_window_N7_12', use(A['slope_fit_window_N7_12'])),
    ('dp3', 'slope_fit_window_N8_12', use(A['slope_fit_window_N8_12'])),
    ('dp4', 'total_local_potential_vacuum_level', use(A['total_local_potential_vacuum_level'])),
    ('dp5', 'thin_slab_ladder_N3_6', use(A['thin_slab_ladder_N3_6'])),
    ('dp5', 'single_7_layer_slab', use(A['single_7_layer_slab'])),
    ('dp6', 'relaxation_from_7_layers', use(A['relaxation_from_7_layers'])),
    ('dp7', 'k12_ladder_for_sigma_phi', use(A['k12_ladder_for_sigma_phi'])),
    ('dp8', 'area_counted_once', use(A['area_counted_once'])),
    ('dp5', 'ordering_called_unresolved', ordering('unresolved')),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, use(values)) for name, values in R.items()]
