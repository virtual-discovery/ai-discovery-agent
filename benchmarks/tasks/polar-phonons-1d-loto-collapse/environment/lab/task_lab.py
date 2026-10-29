"""Lab logic for polar-phonons-1d-loto-collapse: retrieval of author-archived phonon data and grading.

Purchases return archived files of Materials Cloud record 08gsg-gpd20 (Rivano, Marzari, Sohier 2023):
ph.x dynamical-matrix files (dfpt), the archived q2r/matdyn interpolated dispersions (interpolated, parsed
frequencies and axial weights) and archived program outputs (archived_file). Nothing is recomputed. The explicit
small-q 1D runs of BN(5,5) and BN(6,6) below 0.15 |Gamma-Z| and the interpolated dispersions of both tubes are
never served; their LO frequencies (and the 3D-periodic q->0 LO of BN(5,5)) grade the submission.
"""
import gzip
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
INDEX = json.loads((HERE/'data/index.json').read_text())
PUBLIC = json.loads((HERE/'public.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
COSTS = {k: v['cost'] for k, v in PUBLIC['requests'].items()}
SRC = 'retrieved author-archived data (no new calculation)'


def blob(key):
    return gzip.decompress((HERE/'data/blobs'/INDEX[key]['blob']).read_bytes())


def validate(request):
    kind = request.get('kind')
    if kind == 'dfpt':
        if set(request) != {'kind', 'file'} or not isinstance(request['file'], str):
            raise ValueError('Use {"kind": "dfpt", "file": <id from dfpt_runs>}')
        if INDEX.get(request['file'], {}).get('kind') != 'dfpt':
            raise ValueError('No such dynamical-matrix file in the catalog (it is withheld or not archived)')
    elif kind == 'interpolated':
        if set(request) != {'kind', 'folder'} or not isinstance(request['folder'], str):
            raise ValueError('Use {"kind": "interpolated", "folder": <folder from interpolated_dispersions>}')
        if 'interpolated:' + request['folder'] not in INDEX:
            raise ValueError('No archived interpolated dispersion is served for this folder')
    elif kind == 'archived_file':
        if set(request) != {'kind', 'file'} or not isinstance(request['file'], str):
            raise ValueError('Use {"kind": "archived_file", "file": <path from archived_outputs>}')
        if 'file:' + request['file'] not in INDEX:
            raise ValueError('No such archived output in the catalog')
    else:
        raise ValueError('kind must be "dfpt", "interpolated" or "archived_file"')
    return COSTS[kind]


def run(request, workdir):
    kind = request['kind']
    if kind == 'dfpt':
        return {'file': request['file'], 'text': blob(request['file']).decode(), 'source': SRC}
    if kind == 'interpolated':
        return json.loads(blob('interpolated:' + request['folder']))
    return json.loads(blob('file:' + request['file']))


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    targets, t = PRIVATE['targets_1D_LO_cm-1'], PRIVATE['thresholds']
    fields = {'BN55_LO_1D_cm-1': 'BN55', 'BN66_LO_1D_cm-1': 'BN66'}
    missing = [f for f in list(fields) + ['BN55_LO_3D_q0_cm-1', 'evidence_ids'] if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    for field, tube in fields.items():
        sub = payload[field]
        if not isinstance(sub, dict) or set(sub) != set(targets[tube]) or not all(finite(v) for v in sub.values()):
            return False, {'error': f'{field} must give a finite frequency for each of the keys {sorted(targets[tube])}'}
    if not finite(payload['BN55_LO_3D_q0_cm-1']):
        return False, {'error': 'BN55_LO_3D_q0_cm-1 must be a finite number'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    checks, errors = {'evidence_executed': evidence}, {}
    for field, tube in fields.items():
        for q, ref in targets[tube].items():
            err = payload[field][q] - ref
            errors[f'{tube}@{q}'] = round(err, 2)
            checks[f'{tube}_LO_1D_q{q}'] = abs(err) <= t['tolerance_cm-1'][q]
    err3 = payload['BN55_LO_3D_q0_cm-1'] - PRIVATE['target_3D_LO_q0_cm-1']
    errors['BN55_3D_q0'] = round(err3, 2)
    checks['BN55_LO_3D_q0'] = abs(err3) <= t['tolerance_3d_cm-1']
    return all(checks.values()), {'checks': checks, 'errors_cm-1': errors}
