"""Plausible wrong paths for hydrogen-plasma-facing-materials; each must fail grading.

Values come from solution/ablation_values.json, computed from the archived runs by
benchmarks/construction/hydrogen-plasma-facing-materials/make_task.py.
ABLATIONS lists (decision_point_id, name, mutate); mutate(payload, trace) returns the payload the wrong path
would have produced (only the affected entries change).
"""
import copy
import json
from pathlib import Path

V = json.loads(Path(__file__).with_name('ablation_values.json').read_text())
A, R = V['ablations'], V['robustness']


def patch(changes):
    def mutate(payload, trace):
        payload = copy.deepcopy(payload)
        for field, entries in changes.items():
            for key, value in entries.items():
                if isinstance(value, dict) and isinstance(payload[field].get(key), dict) and field == 'surface_binding_energy_eV':
                    payload[field][key].update(value)
                else:
                    payload[field][key] = value
        return payload, trace
    return mutate


ABLATIONS = [
    ('dp1', 'atom_reference', patch(A['atom_reference'])),
    ('dp2', 'textbook_site', patch(A['textbook_site'])),
    ('dp3', 'archived_runs_as_sites', patch(A['archived_runs_as_sites'])),
    ('dp3', 'energy_binning_0p1eV', patch(A['energy_binning_0p1eV'])),
    ('dp4', 'w_atom_from_w_set', patch(A['w_atom_from_w_set'])),
    ('dp4', 'w_atom_from_wc_set', patch(A['w_atom_from_wc_set'])),
    ('dp5', 'other_termination', patch(A['other_termination'])),
    ('dp5', 'paper_protocol_mean', patch(A['paper_protocol_mean'])),
    ('dp5', 'site_mean_lowest_termination', patch(A['site_mean_lowest_termination'])),
    ('dp6', 'cohesive_bulk_from_hife_supercell', patch(A['cohesive_bulk_from_hife_supercell'])),
    ('dp7', 'starting_geometry', patch(A['starting_geometry'])),
    ('dp8', 'cohesive_energy_as_sbe', patch(A['cohesive_energy_as_sbe'])),
]

# Defensible variations that must still pass (rounding to 4 significant digits is checked automatically).
ROBUSTNESS = [
    ('duplicate_cluster_mean', patch(R['duplicate_cluster_mean'])),
    ('carbon_atom_at_60Ry', patch(R['carbon_atom_at_60Ry'])),
    ('distance_from_another_duplicate', patch(R['distance_from_another_duplicate'])),
    ('rydberg_13p6057', patch(R['rydberg_13p6057'])),
    ('tib2_ti_reported_as_undetermined', patch({'surface_binding_energy_eV': {'TiB2-h-191': {'Ti': None}}})),
]
