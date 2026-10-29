"""Lab logic for sngeo2-coherent-miscibility: retrieval of author-released DFT results and grading.

Purchases return archived results (SQS, coherency-strain scans, r2SCAN orderings, HSE-setting PBE runs).
Nothing graded is served: the targets are derived quantities recomputed from the full archive.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
COST = {'PBE': 1, 'r2SCAN': 2}
SOURCE = 'retrieved author-released DFT result (no new calculation)'


def _sqs_key(request):
    x = request.get('x_Ge')
    if not isinstance(x, (int, float)) or isinstance(x, bool):
        return None
    return f"{request.get('functional')}|{float(x):.3f}|{request.get('index')}"


def validate(request):
    kind = request.get('kind')
    if kind == 'sqs':
        if set(request) != {'kind', 'functional', 'x_Ge', 'index'} or _sqs_key(request) not in DATA['sqs']:
            raise ValueError('Use {"kind": "sqs", "functional": "PBE"|"r2SCAN", "x_Ge": <x>, "index": <i>} from sqs_available')
        return COST[request['functional']]
    if kind == 'strain_scan':
        if set(request) != {'kind', 'functional', 'id'} or f"{request.get('functional')}|{request.get('id')}" not in DATA['strain']:
            raise ValueError('Use {"kind": "strain_scan", "functional": "PBE"|"r2SCAN", "id": "<plane>-<member>-<i>[-<j>]"} '
                             '(see strain_scans in the catalog)')
        return COST[request['functional']]
    if kind == 'ordering_r2SCAN':
        if set(request) != {'kind', 'id'} or request.get('id') not in DATA['orderings_r2SCAN']:
            raise ValueError('Use {"kind": "ordering_r2SCAN", "id": "<id from orderings_r2SCAN_available>"}')
        return 2
    if kind == 'hse_setting_run':
        if set(request) != {'kind', 'name'} or request.get('name') not in DATA['hse']:
            raise ValueError('Use {"kind": "hse_setting_run", "name": "<name from hse_setting_runs_available>"}')
        return 1
    raise ValueError('kind must be "sqs", "strain_scan", "ordering_r2SCAN" or "hse_setting_run"')


def run(request, workdir):
    kind = request['kind']
    if kind == 'sqs':
        return {**DATA['sqs'][_sqs_key(request)], 'source': SOURCE}
    if kind == 'strain_scan':
        return {**DATA['strain'][f"{request['functional']}|{request['id']}"], 'energy_unit': 'eV per 6-atom cell', 'source': SOURCE}
    if kind == 'ordering_r2SCAN':
        return {'id': request['id'], **DATA['orderings_r2SCAN'][request['id']], 'source': SOURCE}
    return {**DATA['hse'][request['name']], 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def pair(v):
    return isinstance(v, list) and len(v) == 2 and all(finite(q) for q in v)


def grade(payload, trace, validation):
    ref, t = PRIVATE['reference'], PRIVATE['thresholds']
    fields = ['mixing_enthalpy_x050_meV_per_cation', 'coherency_strain_energy_x050_meV_per_cation', 'habit_plane',
              'incoherent_Tc_K', 'coherent_spinodal_Tc_K', 'coherent_boundaries_700K_PBE', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    strain, tcs, bnd = payload['coherency_strain_energy_x050_meV_per_cation'], payload['coherent_spinodal_Tc_K'], payload['coherent_boundaries_700K_PBE']
    if not isinstance(strain, dict) or set(strain) != {'001', '100', '101'} or not all(finite(v) for v in strain.values()):
        return False, {'error': 'coherency_strain_energy_x050_meV_per_cation needs finite values for "001", "100" and "101"'}
    if not isinstance(tcs, dict) or set(tcs) != {'PBE', 'r2SCAN'} or not all(finite(v) for v in tcs.values()):
        return False, {'error': 'coherent_spinodal_Tc_K needs finite values for "PBE" and "r2SCAN"'}
    if not finite(payload['mixing_enthalpy_x050_meV_per_cation']) or not finite(payload['incoherent_Tc_K']):
        return False, {'error': 'mixing_enthalpy_x050_meV_per_cation and incoherent_Tc_K must be finite numbers'}
    if not isinstance(bnd, dict) or set(bnd) != {'binodal', 'spinodal'}:
        return False, {'error': 'coherent_boundaries_700K_PBE needs "binodal" and "spinodal"'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    lo, hi = t['mixing_enthalpy_window']
    tlo, thi = t['incoherent_Tc_window_K']
    tol = t['boundary_abs_tol']
    def bounds_ok(name):
        v, r = bnd[name], ref['coherent_boundaries_700K_PBE'][name]
        return pair(v) and all(abs(a - b) <= tol for a, b in zip(sorted(v), r))
    checks = {'evidence_executed': evidence,
              'mixing_enthalpy': lo <= payload['mixing_enthalpy_x050_meV_per_cation'] <= hi,
              **{f'strain_energy_{pl}': abs(strain[pl] / ref['coherency_strain_energy_x050_meV_per_cation'][pl] - 1) <= t['strain_rel_tol']
                 for pl in ('001', '100', '101')},
              'habit_plane': payload['habit_plane'] == ref['habit_plane'],
              'incoherent_Tc': tlo <= payload['incoherent_Tc_K'] <= thi,
              **{f'coherent_Tc_{fn}': abs(tcs[fn] / ref['coherent_spinodal_Tc_K'][fn] - 1) <= t['coherent_Tc_rel_tol'][fn]
                 for fn in ('PBE', 'r2SCAN')},
              'coherent_binodal_700K': bounds_ok('binodal'),
              'coherent_spinodal_700K': bounds_ok('spinodal')}
    report = {'checks': checks,
              'strain_relative_error': {pl: strain[pl] / ref['coherency_strain_energy_x050_meV_per_cation'][pl] - 1 for pl in strain},
              'coherent_Tc_relative_error': {fn: tcs[fn] / ref['coherent_spinodal_Tc_K'][fn] - 1 for fn in tcs}}
    return all(checks.values()), report
