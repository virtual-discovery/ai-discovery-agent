"""Lab logic for heis-configdesign-fe: retrieval of author-released DFT energies and grading.

Purchases return released energies (and moments) for catalog configurations. The [214] spirals
and the 3x3x4 random states are never served; their energies grade the submission.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())


def validate(request):
    kind = request.get('kind')
    if kind == 'spiral':
        if set(request) != {'kind', 'key'} or request['key'] not in DATA['spirals']:
            raise ValueError('Use {"kind": "spiral", "key": <key from spirals_available>}')
        return 1
    if kind == 'random_state':
        if set(request) != {'kind', 'label'} or request['label'] not in DATA['random_states']:
            raise ValueError('Use {"kind": "random_state", "label": <label from random_states_available>}')
        return 1
    raise ValueError('kind must be "spiral" or "random_state"')


def run(request, workdir):
    if request['kind'] == 'spiral':
        return {'key': request['key'], **DATA['spirals'][request['key']],
                'source': 'retrieved author-released DFT result (no new calculation)'}
    return {'label': request['label'], **DATA['random_states'][request['label']],
            'source': 'retrieved author-released DFT result (no new calculation)'}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    t, ref = PRIVATE['thresholds'], PRIVATE['reference']
    fields = ['spiral_predictions_meV_per_atom', 'random_predictions_meV_per_atom', 'J1_meV',
              'spin_wave_stiffness_meV_A2', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    sp, rp = payload['spiral_predictions_meV_per_atom'], payload['random_predictions_meV_per_atom']
    targets_s, targets_r = PRIVATE['spiral_targets'], PRIVATE['random_targets']
    if not isinstance(sp, dict) or set(sp) != set(targets_s) or not all(finite(v) for v in sp.values()):
        return False, {'error': 'spiral_predictions_meV_per_atom must give a finite energy for every spirals_to_predict id'}
    if not isinstance(rp, dict) or set(rp) != set(targets_r) or not all(finite(v) for v in rp.values()):
        return False, {'error': 'random_predictions_meV_per_atom must give a finite energy for every random_states_to_predict id'}
    if not finite(payload['J1_meV']) or not finite(payload['spin_wave_stiffness_meV_A2']):
        return False, {'error': 'J1_meV and spin_wave_stiffness_meV_A2 must be finite numbers'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    spiral_rmse = math.sqrt(sum((sp[k] - v) ** 2 for k, v in targets_s.items()) / len(targets_s))
    random_rmse = math.sqrt(sum((rp[k] - v) ** 2 for k, v in targets_r.items()) / len(targets_r))
    random_mean_error = sum(rp.values()) / len(rp) - ref['random_mean_meV']
    stiffness_rel = payload['spin_wave_stiffness_meV_A2'] / ref['stiffness_meV_A2'] - 1
    low, high = t['J1_window_meV']
    checks = {'evidence_executed': evidence,
              'spiral_predictions': spiral_rmse <= t['spiral_rmse_meV'],
              'random_state_mean_energy': abs(random_mean_error) <= t['random_mean_tol_meV'],
              'random_state_rmse': random_rmse <= t['random_rmse_meV'],
              'J1': low <= payload['J1_meV'] <= high,
              'stiffness': abs(stiffness_rel) <= t['stiffness_rel_tol']}
    report = {'checks': checks, 'spiral_rmse_meV': spiral_rmse, 'random_rmse_meV': random_rmse,
              'random_mean_error_meV': random_mean_error, 'J1_meV': payload['J1_meV'],
              'stiffness_relative_error': stiffness_rel}
    return all(checks.values()), report
