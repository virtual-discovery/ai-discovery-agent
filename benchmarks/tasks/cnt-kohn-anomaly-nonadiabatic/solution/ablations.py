"""Plausible wrong paths for cnt-kohn-anomaly-nonadiabatic; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/cnt-kohn-anomaly-nonadiabatic/make_task.py from the archived tables, the
archived coupling data (through the lab's own evaluator) and the archived threshold program.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def setting(fields):
    def mutate(payload, trace):
        for field, value in fields.items():
            payload[field] = value
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'no_dfpt_subtraction', setting(A['no_dfpt_subtraction'])),
    ('dp2', 'same_temperature_reference', setting(A['same_temperature_reference'])),
    ('dp2', 'smearing_read_as_eV', setting(A['smearing_read_as_eV'])),
    ('dp3', 'adiabatic_for_nonadiabatic', setting(A['adiabatic_for_nonadiabatic'])),
    ('dp3', 'nonadiabatic_for_adiabatic', setting(A['nonadiabatic_for_adiabatic'])),
    ('dp4', 'header_column_order', setting(A['header_column_order'])),
    ('dp5', 'both_processes_summed', setting(A['both_processes_summed'])),
    ('dp6', 'sorted_intraband_only', setting(A['sorted_intraband_only'])),
    ('dp7', 'coarse_mesh_1e5', setting(A['coarse_mesh_1e5'])),
    ('dp7', 'native_archived_grid', setting(A['native_archived_grid'])),
    ('dp8', 'fd_scan_reference_frequency', setting(A['fd_scan_reference_frequency'])),
    ('dp8', 'dfpt_scan_as_temperature_proxy', setting(A['dfpt_scan_as_temperature_proxy'])),
    ('dp9', 'inverted_mass_parameter', setting(A['inverted_mass_parameter'])),
    ('dp10', 'interpolate_released_tables', setting(A['interpolate_released_tables'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, setting(fields)) for name, fields in R.items()]
