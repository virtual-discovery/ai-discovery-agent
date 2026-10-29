"""Lab logic for mn5si3-noncollinear-spinwaves: retrieval of author-archived Mn5Si3 data and grading.

Purchases return archived items of Biniskos, dos Santos et al. (Materials Cloud paxx3-a1k09 and 4tgt6-1p378):
KKR exchange tables, the authors' Mn1 couplings, 10 K dispersion tables and constant-energy cuts, energy scans,
the gap table and constant-energy Q scans. Everything measured along (0,3,Ql), the 5-20 meV energy scan at (2,1,0)
and the authors' relaxed ground state are never served; they grade the submission.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
SOURCE = 'retrieved author-archived data (no new measurement or calculation)'
KINDS = ('exchange_set', 'dispersion_table', 'constant_energy_cut', 'energy_scan', 'gap_table', 'q_scan')


def resolve(request):
    kind = request.get('kind')
    if kind == 'exchange_set':
        if set(request) != {'kind', 'name'} or request['name'] not in DATA['exchange']:
            raise ValueError('Use {"kind": "exchange_set", "name": <one of %s>}' % sorted(DATA['exchange']))
        return kind, request['name'], 1
    if kind == 'dispersion_table':
        if set(request) != {'kind', 'direction'} or request['direction'] not in DATA['dispersion']:
            raise ValueError('Use {"kind": "dispersion_table", "direction": "h" | "k"}; (0,3,Ql) data are not sold')
        return kind, request['direction'], 1
    if kind == 'constant_energy_cut':
        d, E = request.get('direction'), request.get('energy_meV')
        if set(request) != {'kind', 'direction', 'energy_meV'} or isinstance(E, bool) or not isinstance(E, (int, float)):
            raise ValueError('Use {"kind": "constant_energy_cut", "direction": "h" | "k", "energy_meV": <number>}')
        key = '%s_%g' % (d, E)
        if key not in DATA['cuts']:
            raise ValueError('No such cut for sale; see constant_energy_cuts in /app/catalog.json')
        return kind, key, 2
    if kind == 'energy_scan':
        if set(request) != {'kind', 'name'} or request['name'] not in DATA['escans']:
            raise ValueError('Use {"kind": "energy_scan", "name": <name from energy_scans in /app/catalog.json>}')
        return kind, request['name'], 1
    if kind == 'gap_table':
        if set(request) != {'kind'}:
            raise ValueError('Use {"kind": "gap_table"}')
        return kind, None, 1
    if kind == 'q_scan':
        if set(request) != {'kind', 'name'} or request['name'] not in DATA['qscans']:
            raise ValueError('Use {"kind": "q_scan", "name": <name from q_scans in /app/catalog.json>}')
        return kind, request['name'], 1
    raise ValueError('kind must be one of: ' + ', '.join(KINDS))


def validate(request):
    return resolve(request)[2]


def run(request, workdir):
    kind, key, _ = resolve(request)
    if kind == 'exchange_set':
        out = DATA['exchange'][key]
    elif kind == 'dispersion_table':
        out = DATA['dispersion'][key]
    elif kind == 'constant_energy_cut':
        out = DATA['cuts'][key]
    elif kind == 'energy_scan':
        out = {'name': key, **DATA['escans'][key]}
    elif kind == 'gap_table':
        out = DATA['gap_table']
    else:
        out = {'name': key, **DATA['qscans'][key]}
    return {**out, 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    P, th = PRIVATE, PRIVATE['thresholds']
    fields = ['l_branch_positions_rlu', 'zone_centre_excitation_meV', 'mn1_carries_moment', 'mn2_canting_deg', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    pred = payload['l_branch_positions_rlu']
    targets = P['l_targets']
    if not isinstance(pred, dict) or set(pred) != set(targets) or not all(
            isinstance(v, list) and len(v) == 2 and all(finite(x) for x in v) for v in pred.values()):
        return False, {'error': 'l_branch_positions_rlu must map each of %s to [Ql below 1, Ql above 1]' % sorted(targets, key=float)}
    if not isinstance(payload['mn1_carries_moment'], bool):
        return False, {'error': 'mn1_carries_moment must be true or false'}
    for f in ('zone_centre_excitation_meV', 'mn2_canting_deg'):
        if not finite(payload[f]):
            return False, {'error': f + ' must be a finite number'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    dev = []
    for E, t in targets.items():
        lo, hi = sorted(pred[E])
        dev += [lo - t['low']['q_rlu'], hi - t['high']['q_rlu']]
    rms = math.sqrt(sum(d * d for d in dev) / len(dev))
    lo_o, hi_o = th['optical_window_meV']
    lo_c, hi_c = th['canting_window_deg']
    checks = {'evidence_executed': evidence,
              'l_branch_positions': rms <= th['l_rms_rlu'],
              'zone_centre_excitation': lo_o <= payload['zone_centre_excitation_meV'] <= hi_o,
              'mn1_carries_moment': payload['mn1_carries_moment'] is P['mn1_carries_moment'],
              'mn2_canting': lo_c <= payload['mn2_canting_deg'] <= hi_c}
    return all(checks.values()), {'checks': checks, 'l_rms_rlu': rms}
