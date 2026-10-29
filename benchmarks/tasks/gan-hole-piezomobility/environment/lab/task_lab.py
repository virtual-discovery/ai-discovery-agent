"""Lab logic for gan-hole-piezomobility: retrieval of author-archived GW band energies and EPW hole mobilities.

Purchases return archived results for the 37 strain states (no new calculation). The mobilities of the states in
catalog.json `to_predict` are never served; their 300 K IBTE tensors, and the two design answers derived from the
2 % states, grade the submission.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
BANDS = json.loads((HERE/'bands.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
COSTS = DATA['costs']
STATES = set(BANDS['energies_eV'])
WITHHELD = set(DATA['withheld'])
SOURCE = 'retrieved author-archived result (no new calculation)'
KEYS = {'xx': (0, 0), 'yy': (1, 1), 'zz': (2, 2), 'yz': (1, 2), 'xz': (0, 2), 'xy': (0, 1)}


def check_state(request):
    if request.get('state') not in STATES:
        raise ValueError('state must be one of the labels listed under "states" in /app/catalog.json')


def validate(request):
    kind = request.get('kind')
    if kind == 'bands':
        if set(request) != {'kind', 'state'}:
            raise ValueError('Use {"kind": "bands", "state": <label>}')
        check_state(request)
        return COSTS['bands']
    if kind == 'mobility':
        if set(request) != {'kind', 'state', 'solver'}:
            raise ValueError('Use {"kind": "mobility", "state": <label>, "solver": "serta" or "ibte"}')
        check_state(request)
        if request['solver'] not in ('serta', 'ibte'):
            raise ValueError('solver must be "serta" or "ibte"')
        if request['state'] in WITHHELD:
            raise ValueError(f"the mobility of {request['state']} is withheld for grading (listed in to_predict); "
                             'its bands can be purchased')
        return COSTS['mobility_' + request['solver']]
    raise ValueError('kind must be "bands" or "mobility"')


def run(request, workdir):
    state = request['state']
    if request['kind'] == 'bands':
        return {'state': state, 'band_indices': list(range(BANDS['band_indices'][0], BANDS['band_indices'][1] + 1)),
                'vbm_eV': DATA['vbm_eV'][state], 'units': 'eV (absolute GW quasiparticle energies)',
                'k_points': 'rows follow kgrid_crystal in /app/catalog.json (crystal coordinates of this state)',
                'energies_eV': BANDS['energies_eV'][state], 'source': SOURCE}
    return {**DATA['mobility'][state][request['solver']], 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def strain_matches(vec, index, value, either_sign, tol):
    if not isinstance(vec, list) or len(vec) != 6 or not all(finite(v) for v in vec):
        return False
    others = all(abs(v) <= tol for i, v in enumerate(vec) if i != index)
    hit = abs(vec[index] - value) <= tol or (either_sign and abs(vec[index] + value) <= tol)
    return others and hit


def grade(payload, trace, validation):
    th, targets, graded = PRIVATE['thresholds'], PRIVATE['targets'], PRIVATE['graded']
    fields = ['mobility_predictions_cm2_per_Vs', 'best_mu_zz_2pct', 'best_inplane_2pct', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    preds = payload['mobility_predictions_cm2_per_Vs']
    if not isinstance(preds, dict) or set(preds) != set(targets):
        return False, {'error': 'mobility_predictions_cm2_per_Vs must give a tensor for every state in to_predict'}
    for lab, tensor in preds.items():
        if not isinstance(tensor, dict) or set(tensor) != set(KEYS) or not all(finite(v) for v in tensor.values()):
            return False, {'error': f'{lab}: give finite xx, yy, zz, yz, xz and xy (cm^2/Vs)'}
    bz, bi = payload['best_mu_zz_2pct'], payload['best_inplane_2pct']
    if not isinstance(bz, dict) or set(bz) != {'strain_voigt', 'mu_zz_cm2_per_Vs'} or not finite(bz['mu_zz_cm2_per_Vs']):
        return False, {'error': 'best_mu_zz_2pct must be {"strain_voigt": [6 numbers], "mu_zz_cm2_per_Vs": number}'}
    if not isinstance(bi, dict) or set(bi) != {'strain_voigt', 'mu_inplane_max_cm2_per_Vs'} or not finite(bi['mu_inplane_max_cm2_per_Vs']):
        return False, {'error': 'best_inplane_2pct must be {"strain_voigt": [6 numbers], "mu_inplane_max_cm2_per_Vs": number}'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    checks, errors = {'evidence_executed': evidence}, {}
    for lab, keys in graded.items():
        for k in keys:
            pred, true = preds[lab][k], targets[lab][k]
            window = th['window'].get(lab, {}).get(k)
            if window:
                checks[f'{lab}:{k}'] = window[0] <= pred <= window[1]
                errors[f'{lab}:{k}'] = pred
            else:
                rel = pred / true - 1
                checks[f'{lab}:{k}'] = abs(rel) <= (th['diag_rel'] if k in ('xx', 'yy', 'zz') else th['offdiag_rel'])
                errors[f'{lab}:{k}'] = rel
    ref_z, ref_i = PRIVATE['best_mu_zz_2pct'], PRIVATE['best_inplane_2pct']
    zi = max(range(6), key=lambda i: abs(ref_z['strain_voigt'][i]))
    checks['best_mu_zz_2pct_strain'] = strain_matches(bz['strain_voigt'], zi, ref_z['strain_voigt'][zi], False, th['strain_abs'])
    errors['best_mu_zz_2pct_mu'] = bz['mu_zz_cm2_per_Vs'] / ref_z['mu_zz'] - 1
    checks['best_mu_zz_2pct_mu'] = abs(errors['best_mu_zz_2pct_mu']) <= th['best_mu_zz_rel']
    checks['best_inplane_2pct_strain'] = strain_matches(bi['strain_voigt'], ref_i['component_index'], ref_i['abs_value'],
                                                        True, th['strain_abs'])
    errors['best_inplane_2pct_mu'] = bi['mu_inplane_max_cm2_per_Vs'] / ref_i['mu_ref'] - 1
    checks['best_inplane_2pct_mu'] = abs(errors['best_inplane_2pct_mu']) <= th['best_inplane_rel']
    return all(checks.values()), {'checks': checks, 'relative_errors_or_values': errors}
