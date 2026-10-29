"""Lab logic for mnte-altermagnet-chiral-magnons: retrieval of archived polarized-neutron data and grading.

Purchases return the authors' archived HYSPEC slices (in four tiles), cuts and polarized diffraction, and the raw IN20
polarized diffraction files. The IN20 constant-energy inelastic scans (states A and B) are never served;
their peak positions and chiral signs grade the submission.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
COSTS = {'hyspec_tile': 1, 'hyspec_cut': 1, 'pnd': 1}
TILES = ('u_low_v_low', 'u_low_v_high', 'u_high_v_low', 'u_high_v_high')


def key_of(request):
    kind = request.get('kind')
    if kind == 'hyspec_tile':
        if set(request) != {'kind', 'quantity', 'tile'} or request['quantity'] not in ('symmetric', 'chiral') or request['tile'] not in TILES:
            raise ValueError('Use {"kind": "hyspec_tile", "quantity": "symmetric" | "chiral", "tile": one of %s}' % list(TILES))
        return 'hyspec_tile:%s:%s' % (request['quantity'], request['tile'])
    if kind == 'hyspec_cut':
        if set(request) != {'kind', 'cut'} or request['cut'] not in ('g', 'h', 'i', 'j'):
            raise ValueError('Use {"kind": "hyspec_cut", "cut": "g" | "h" | "i" | "j"}')
        return 'hyspec_cut:' + request['cut']
    if kind == 'pnd':
        if request.get('instrument') == 'HYSPEC' and set(request) == {'kind', 'instrument'}:
            return 'pnd:HYSPEC'
        if request.get('instrument') == 'IN20' and set(request) == {'kind', 'instrument', 'state'} and request['state'] in ('A', 'B'):
            return 'pnd:IN20:' + request['state']
        raise ValueError('Use {"kind": "pnd", "instrument": "HYSPEC"} or {"kind": "pnd", "instrument": "IN20", "state": "A" | "B"}')
    raise ValueError('kind must be "hyspec_tile", "hyspec_cut" or "pnd"')


def validate(request):
    key_of(request)
    return COSTS[request['kind']]


def run(request, workdir):
    return {'item': key_of(request), **DATA[key_of(request)],
            'source': 'retrieved author-archived measurement (no new experiment)'}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def q_score(q):
    """|g-wave form factor| of the 9.83 A altermagnetic bond pair at q, relative to its maximum (3 sqrt(3)/2)."""
    h, k, l = q
    s = lambda x: math.sin(2 * math.pi * x)
    g = s(2 * h + k) + s(k - h) - s(h + 2 * k)
    return abs(g * s(l)) / (1.5 * math.sqrt(3))


def grade(payload, trace, validation):
    t, ref = PRIVATE['thresholds'], PRIVATE['in20']
    fields = ['split_exchange_difference_meV', 'in20_peak_qh', 'in20_chirality', 'max_splitting_q_rlu', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    dj, pk, ch, q = (payload[f] for f in fields[:4])
    if not finite(dj):
        return False, {'error': 'split_exchange_difference_meV must be a finite number'}
    if not isinstance(pk, dict) or set(pk) != {'low', 'high'} or not all(finite(v) for v in pk.values()):
        return False, {'error': 'in20_peak_qh must be {"low": number, "high": number}'}
    if not (isinstance(ch, dict) and set(ch) == {'A', 'B'} and all(isinstance(ch[s], dict) and set(ch[s]) == {'low', 'high'}
                                                                    and all(v in ('+', '-') for v in ch[s].values()) for s in ch)):
        return False, {'error': 'in20_chirality must be {"A": {"low": "+"|"-", "high": "+"|"-"}, "B": {...}}'}
    if not (isinstance(q, list) and len(q) == 3 and all(finite(x) for x in q)):
        return False, {'error': 'max_splitting_q_rlu must be a list of three numbers'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    lo, hi = t['dJ_window_meV']
    mean = (pk['low'] + pk['high']) / 2
    sep = pk['high'] - pk['low']
    sep_rel = sep / ref['separation_qh'] - 1
    score = q_score(q)
    checks = {'evidence_executed': evidence,
              'split_exchange_difference': lo <= abs(dj) <= hi,
              'in20_peak_mean': abs(mean - ref['mean_qh']) <= t['in20_mean_tol_qh'],
              'in20_peak_separation': abs(sep_rel) <= t['in20_separation_rel_tol'],
              'in20_chirality_state_A': ch['A'] == ref['chirality']['A'],
              'in20_chirality_state_B': ch['B'] == ref['chirality']['B'],
              'max_splitting_wave_vector': score >= t['q_score_min']}
    report = {'checks': checks, 'peak_mean_error_qh': mean - ref['mean_qh'], 'peak_separation_relative_error': sep_rel,
              'q_score': score}
    return all(checks.values()), report
