"""Lab logic for heis-configdesign-mnf2-altermagnet: retrieval of author-released DFT spin spirals and grading.

Purchases return the released calculation for a catalog spiral (VASP total energy of the 6-atom cell, the
archive's CSV energy, constrained moments, INCAR, k-mesh, OUTCAR warnings). The [111] and [113] spirals are
never served; their energies, and reference values derived from the released data, grade the submission.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())


def validate(request):
    if request.get('kind') != 'spiral':
        raise ValueError('kind must be "spiral"')
    if set(request) != {'kind', 'key'} or request['key'] not in DATA['spirals']:
        raise ValueError('Use {"kind": "spiral", "key": <key from spirals_available>}')
    return 1


def run(request, workdir):
    return {'key': request['key'], **DATA['spirals'][request['key']],
            'source': 'retrieved author-released DFT result (no new calculation)'}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    t, ref = PRIVATE['thresholds'], PRIVATE['reference']
    fields = ['spiral_predictions_meV_per_cell', 'J1_meV', 'J2_meV', 'J3_meV', 'neel_temperature_K',
              'altermagnetic_shell', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    sp, targets = payload['spiral_predictions_meV_per_cell'], PRIVATE['spiral_targets']
    if not isinstance(sp, dict) or set(sp) != set(targets) or not all(finite(v) for v in sp.values()):
        return False, {'error': 'spiral_predictions_meV_per_cell must give a finite energy for every spirals_to_predict id'}
    scalars = ['J1_meV', 'J2_meV', 'J3_meV', 'neel_temperature_K']
    if not all(finite(payload[k]) for k in scalars):
        return False, {'error': 'J1_meV, J2_meV, J3_meV and neel_temperature_K must be finite numbers'}
    alt = payload['altermagnetic_shell']
    if (not isinstance(alt, dict) or set(alt) != {'distance_A', 'sublattice_mean_J_meV', 'sublattice_difference_J_meV'}
            or not finite(alt['distance_A']) or not finite(alt['sublattice_mean_J_meV'])
            or not (alt['sublattice_difference_J_meV'] is None or finite(alt['sublattice_difference_J_meV']))):
        return False, {'error': 'altermagnetic_shell must hold distance_A, sublattice_mean_J_meV (numbers) and '
                                'sublattice_difference_J_meV (number or null)'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    rmse = math.sqrt(sum((sp[k] - v) ** 2 for k, v in targets.items()) / len(targets))
    ground = PRIVATE['ground_state_ids']
    ground_error = sum(sp[k] - targets[k] for k in ground) / len(ground)
    tn_rel = payload['neel_temperature_K'] / ref['TN_K'] - 1
    in_window = lambda x, w: w[0] <= x <= w[1]
    checks = {'evidence_executed': evidence,
              'spiral_predictions': rmse <= t['spiral_rmse_meV'],
              'ground_state_energy': abs(ground_error) <= t['ground_state_tol_meV'],
              'J1': in_window(payload['J1_meV'], t['J1_window_meV']),
              'J2': in_window(payload['J2_meV'], t['J2_window_meV']),
              'J3': in_window(payload['J3_meV'], t['J3_window_meV']),
              'neel_temperature': abs(tn_rel) <= t['TN_rel_tol'],
              'altermagnetic_shell_distance': abs(alt['distance_A'] - ref['split_shell_distance_A']) <= t['split_shell_distance_tol_A'],
              'altermagnetic_shell_mean': in_window(alt['sublattice_mean_J_meV'], t['split_shell_mean_window_meV']),
              'altermagnetic_difference_identifiable': alt['sublattice_difference_J_meV'] is None}
    report = {'checks': checks, 'spiral_rmse_meV_per_cell': rmse, 'ground_state_error_meV_per_cell': ground_error,
              'J_meV': [payload['J1_meV'], payload['J2_meV'], payload['J3_meV']],
              'neel_temperature_relative_error': tn_rel}
    return all(checks.values()), report
