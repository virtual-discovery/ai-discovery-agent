"""Plausible wrong paths for ge-pressure-polymorph-pathways; each must fail grading.

Values come from solution/ablation_values.json, computed from the released data by
benchmarks/construction/ge-pressure-polymorph-pathways/make_task.py.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']
FIELDS = ('transition_pressures_GPa', 'zero_pressure_energy_differences_meV_per_atom', 'metastable_order_at_0GPa')


def use(values):
    def mutate(payload, trace):
        for f in FIELDS:
            payload[f] = json.loads(json.dumps(values[f]))
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'rattled_frames_included', use(A['rattled_frames_included'])),
    ('dp2', 'enthalpy_from_printed_stress', use(A['enthalpy_from_printed_stress'])),
    ('dp3', 'linear_dE_over_dV', use(A['linear_dE_over_dV'])),
    ('dp4', 'archived_structures_uncorrected', use(A['archived_structures_uncorrected'])),
    ('dp5', 'surrogate_only', use(A['surrogate_only'])),
    ('dp5', 'surrogate_for_metastable_phases', use(A['surrogate_for_metastable_phases'])),
    ('dp6', 'force_filter_discards_compressed_frames', use(A['force_filter_discards_compressed_frames'])),
    ('dp7', 'pstress_as_pressure', use(A['pstress_as_pressure'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, use(R[name])) for name in ('vinet_eos', 'oracle_records_only', 'oracle_records_without_highest_pressure_run',
                                                'correct_all_phases', 'stricter_hydrostatic_filter', 'relaxation_runs_only')]
