"""Lab logic for mg2x-midir-forbidden-transition: retrieval of author-archived results and grading.

Purchases return archived outputs of Sugawara (2026), Zenodo record 23319151 (ecalj LDA/G0W0/QSGW bands, scans
and boxes around X; TOTE files; Questaal eigenvalues; BerkeleyGW quasiparticle energies; PBE momentum matrix
elements; electron-phonon renormalizations). Nothing is recomputed. The graded quantities are derived values
(signed X3 - X1 splittings, the ladder correction, infrared/Raman ratios and predicted splittings) held in
private/hidden.json and never served.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
COSTS = DATA['costs']
COMPOUNDS = ('Mg2Si', 'Mg2Ge', 'Mg2Sn')
SOURCE = 'retrieved author-archived result (no new calculation)'


def need(request, keys, usage):
    if set(request) != set(keys):
        raise ValueError('Use ' + usage)


def grid_key(point):
    if not (isinstance(point, list) and len(point) == 3 and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in point)):
        raise ValueError('each point must be [q_par, q_y, q_z]')
    key = '|'.join(f'{v:.2f}' for v in point)
    if any(abs(v * 50 - round(v * 50)) > 1e-6 for v in point) or key not in DATA['momentum_matrix']:
        raise ValueError(f'point {point} is not on the archived grid (step 0.02; 0 <= q_par <= 0.24; 0 <= q_y, q_z <= 0.16)')
    return key


def validate(request):
    kind = request.get('kind')
    if kind == 'bands':
        need(request, ('kind', 'compound', 'level'), '{"kind": "bands", "compound": "Mg2Si|Mg2Ge|Mg2Sn", "level": "LDA|QSGW"}')
        if f"{request['compound']}|{request['level']}" not in DATA['bands']:
            raise ValueError('compound must be Mg2Si, Mg2Ge or Mg2Sn and level LDA or QSGW')
        return COSTS['bands']
    if kind == 'x_scan':
        need(request, ('kind', 'compound', 'level', 'direction'),
             '{"kind": "x_scan", "compound": "Mg2Si|Mg2Ge|Mg2Sn", "level": "LDA|QSGW", "direction": "X-Gamma|X-W"}')
        if f"{request['compound']}|{request['level']}|{request['direction']}" not in DATA['x_scan']:
            raise ValueError('compound must be Mg2Si, Mg2Ge or Mg2Sn, level LDA or QSGW, direction X-Gamma or X-W')
        return COSTS['x_scan']
    if kind == 'qsgw_box':
        need(request, ('kind', 'compound', 'box'), '{"kind": "qsgw_box", "compound": "Mg2Si|Mg2Ge|Mg2Sn", "box": "coarse|fine_0.01|fine_0.005"}')
        if f"{request['compound']}|{request['box']}" not in DATA['qsgw_box']:
            raise ValueError('available boxes: coarse for every compound; fine_0.01 and fine_0.005 for Mg2Si')
        return COSTS['qsgw_box']
    if kind == 'momentum_matrix':
        need(request, ('kind', 'points'), '{"kind": "momentum_matrix", "points": [[q_par, q_y, q_z], ...]}')
        pts = request['points']
        if not isinstance(pts, list) or not 0 < len(pts) <= DATA['max_momentum_matrix_points']:
            raise ValueError(f"points must be a list of 1 to {DATA['max_momentum_matrix_points']} grid points")
        keys = [grid_key(p) for p in pts]
        if len(set(keys)) != len(keys):
            raise ValueError('points must be distinct')
        return COSTS['momentum_matrix']
    if kind in ('tote', 'questaal', 'pseudopotential_g0w0'):
        need(request, ('kind', 'run'), '{"kind": "%s", "run": "<run from the catalog>"}' % kind)
        if request['run'] not in DATA[kind]:
            raise ValueError(f'unknown {kind} run; see catalog items.{kind}')
        return COSTS[kind]
    if kind == 'electron_phonon':
        method = request.get('method')
        if method == 'AHC':
            c = request.get('compound')
            keys = ('kind', 'compound', 'method', 'q_grid', 'temperature_K') + (('ge_3d',) if c == 'Mg2Ge' else ())
            need(request, keys, '{"kind": "electron_phonon", "compound": "Mg2Si|Mg2Ge", "method": "AHC", "q_grid": <int>, '
                                '"temperature_K": 0|300} (Mg2Ge also "ge_3d": "core|valence")')
            key = f"{c}|AHC|{request['q_grid']}|{request['temperature_K']}|{request.get('ge_3d', 'core')}"
        elif method == 'ZG':
            need(request, ('kind', 'compound', 'method', 'supercell', 'temperature_K'),
                 '{"kind": "electron_phonon", "compound": "Mg2Si|Mg2Ge", "method": "ZG", "supercell": <int>, "temperature_K": 0|300}')
            key = f"{request['compound']}|ZG|{request['supercell']}|{request['temperature_K']}"
        else:
            raise ValueError('method must be "AHC" or "ZG"')
        if key not in DATA['electron_phonon']:
            raise ValueError('no such archived electron-phonon calculation; see catalog items.electron_phonon')
        return COSTS['electron_phonon']
    raise ValueError('kind must be bands, x_scan, qsgw_box, momentum_matrix, tote, questaal, pseudopotential_g0w0 or electron_phonon')


def run(request, workdir):
    kind = request['kind']
    if kind == 'bands':
        out = DATA['bands'][f"{request['compound']}|{request['level']}"]
    elif kind == 'x_scan':
        out = DATA['x_scan'][f"{request['compound']}|{request['level']}|{request['direction']}"]
    elif kind == 'qsgw_box':
        out = DATA['qsgw_box'][f"{request['compound']}|{request['box']}"]
    elif kind == 'momentum_matrix':
        rows = [[p] + DATA['momentum_matrix'][grid_key(p)] for p in request['points']]
        out = {'compound': 'Mg2Si', 'level': 'PBE (pseudopotential, Quantum ESPRESSO)',
               'columns': ['q_par, q_y, q_z (2 pi / a, from X = (1, 0, 0))', 'E_band13_eV', 'E_band14_eV', '|p_13,14|^2 (atomic units)'],
               'rows': rows}
    elif kind in ('tote', 'questaal', 'pseudopotential_g0w0'):
        out = DATA[kind][request['run']]
    else:
        if request['method'] == 'AHC':
            key = f"{request['compound']}|AHC|{request['q_grid']}|{request['temperature_K']}|{request.get('ge_3d', 'core')}"
        else:
            key = f"{request['compound']}|ZG|{request['supercell']}|{request['temperature_K']}"
        out = DATA['electron_phonon'][key]
    return {**out, 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    t, tol = PRIVATE['targets'], PRIVATE['tolerances']
    fields = ['x3_minus_x1_eV', 'vertex_correction_Mg2Si_meV', 'ir_raman_ratio_300K_over_77K', 'predicted_splitting_eV',
              'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    sp = payload['x3_minus_x1_eV']
    if (not isinstance(sp, dict) or set(sp) != set(COMPOUNDS)
            or not all(isinstance(v, dict) and set(v) == {'LDA', 'G0W0', 'QSGW'} and all(finite(x) for x in v.values())
                       for v in sp.values())):
        return False, {'error': 'x3_minus_x1_eV must give finite LDA, G0W0 and QSGW values for Mg2Si, Mg2Ge and Mg2Sn'}
    for name in ('ir_raman_ratio_300K_over_77K', 'predicted_splitting_eV'):
        f = payload[name]
        if not isinstance(f, dict) or set(f) != set(t[name]) or not all(finite(v) for v in f.values()):
            return False, {'error': f'{name} must give a finite value for each of {sorted(t[name])}'}
    if not finite(payload['vertex_correction_Mg2Si_meV']):
        return False, {'error': 'vertex_correction_Mg2Si_meV must be a finite number'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    for c in COMPOUNDS:
        for lev in ('LDA', 'G0W0', 'QSGW'):
            k = f'x3_minus_x1:{c}:{lev}'
            errors[k] = sp[c][lev] - t['x3_minus_x1_eV'][c][lev]
            checks[k] = abs(errors[k]) <= tol['x3_minus_x1_eV']
    errors['vertex_correction'] = payload['vertex_correction_Mg2Si_meV'] - t['vertex_correction_Mg2Si_meV']
    checks['vertex_correction'] = abs(errors['vertex_correction']) <= tol['vertex_meV']
    for k, ref in t['ir_raman_ratio_300K_over_77K'].items():
        errors[f'ir_raman_ratio:{k}'] = payload['ir_raman_ratio_300K_over_77K'][k] / ref - 1
        checks[f'ir_raman_ratio:{k}'] = abs(errors[f'ir_raman_ratio:{k}']) <= tol['ratio_rel']
    for k, ref in t['predicted_splitting_eV'].items():
        errors[f'predicted_splitting:{k}'] = payload['predicted_splitting_eV'][k] - ref
        checks[f'predicted_splitting:{k}'] = abs(errors[f'predicted_splitting:{k}']) <= tol['predicted_eV']
    return all(checks.values()), {'checks': checks, 'errors': errors}
