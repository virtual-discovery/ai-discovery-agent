"""Plausible wrong paths for skutterudite-wigner-coherence; each must fail grading.

Values come from solution/ablation_values.json, computed by
benchmarks/construction/skutterudite-wigner-coherence/make_task.py from the archived per-mode data.
Each wrong path replaces all graded numbers by what that analysis produces.
"""
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def answers(values):
    def mutate(payload, trace):
        payload['YbIr4Sb12'] = json.loads(json.dumps(values['YbIr4Sb12']))
        payload['IrSb3'] = json.loads(json.dumps(values['IrSb3']))
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'archived_isotope_rates_as_is', answers(A['archived_isotope_rates_as_is'])),
    ('dp2', 'isotope_scattering_omitted', answers(A['isotope_scattering_omitted'])),
    ('dp3', 'lifetime_convention_from_tc_data', answers(A['lifetime_convention_from_tc_data'])),
    ('dp4', 'degenerate_pairs_counted_as_coherences', answers(A['degenerate_pairs_counted_as_coherences'])),
    ('dp5', 'gaussian_band_overlap', answers(A['gaussian_band_overlap'])),
    ('dp5', 'half_width_lorentzian', answers(A['half_width_lorentzian'])),
    ('dp6', 'corrupted_qpoint_discarded', answers(A['corrupted_qpoint_discarded'])),
    ('dp7', 'pdos_axis_taken_as_THz', answers(A['pdos_axis_taken_as_THz'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [(name, answers(values)) for name, values in R.items()]
