"""Lab logic for ti-codoping-exchange-tc.

Two request kinds:
  pair         - retrieval of the authors' archived Jij_table.csv rows for one dopant pair and one shell window;
  monte_carlo  - classical Monte Carlo of an agent-specified dilute exchange model on the Bi sublattice of Bi2Te3
                 (dilute.py builds the random configurations, mc_dilute.c runs heat-bath + over-relaxation sweeps).
The reference ordering temperatures are never served; they grade the submission.
"""
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import dilute

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
SPECIES = DATA['species']
WINDOWS = DATA['windows']
MC = DATA['monte_carlo']


def pair_key(a, b):
    za, zb = SPECIES[a], SPECIES[b]
    return f'{a}-{b}' if (za, a) <= (zb, b) else f'{b}-{a}'


def check_mc(request):
    allowed = {'kind', 'spins', 'exchange_meV', 'size', 'samples', 'temperatures_K', 'seed'}
    if not set(request) <= allowed or not {'spins', 'exchange_meV', 'size', 'samples', 'temperatures_K'} <= set(request):
        raise ValueError('monte_carlo needs spins, exchange_meV, size, samples, temperatures_K (optional seed); see catalog')
    spins = request['spins']
    if not isinstance(spins, dict) or not 1 <= len(spins) <= 2 or not all(s in SPECIES for s in spins):
        raise ValueError('spins: {species: concentration} for one or two dopant species of the catalog')
    if not all(isinstance(c, (int, float)) and not isinstance(c, bool) and 0.005 <= c <= 0.11 for c in spins.values()):
        raise ValueError('each concentration (fraction of Bi sites) must lie in [0.005, 0.11]')
    if sum(spins.values()) > 0.11 + 1e-12:
        raise ValueError('total spin concentration must not exceed 0.11')
    ex = request['exchange_meV']
    need = {pair_key(a, b) for a in spins for b in spins}
    if not isinstance(ex, dict) or set(ex) != need:
        raise ValueError('exchange_meV must have exactly the keys %s (species ordered as in the archive)' % sorted(need))
    for k, v in ex.items():
        if not (isinstance(v, list) and len(v) == 17 and all(isinstance(x, (int, float)) and not isinstance(x, bool)
                                                            and math.isfinite(x) and abs(x) <= 100 for x in v)):
            raise ValueError(f'exchange_meV[{k}] must be a list of 17 finite numbers (meV, archive shells 0-16)')
    size, samples = request['size'], request['samples']
    if type(size) is not int or str(size) not in MC['cost_per_16_samples']:
        raise ValueError('size must be one of %s' % sorted(int(s) for s in MC['cost_per_16_samples']))
    if type(samples) is not int or samples < 1 or samples > MC['max_samples'][str(size)]:
        raise ValueError(f'samples must be an integer from 1 to {MC["max_samples"][str(size)]} for size {size}')
    temps = request['temperatures_K']
    if not (isinstance(temps, list) and 1 <= len(temps) <= MC['max_temperatures']
            and all(isinstance(t, (int, float)) and not isinstance(t, bool) and 0.5 <= t <= 500 for t in temps)
            and len(set(temps)) == len(temps)):
        raise ValueError(f'temperatures_K: 1-{MC["max_temperatures"]} distinct values between 0.5 and 500 K')
    seed = request.get('seed', 0)
    if type(seed) is not int or not 0 <= seed < 10**6:
        raise ValueError('seed must be an integer in [0, 1000000)')
    return MC['cost_per_16_samples'][str(size)] * math.ceil(samples / 16)


def validate(request):
    kind = request.get('kind')
    if kind == 'pair':
        if set(request) != {'kind', 'species', 'shells'}:
            raise ValueError('Use {"kind": "pair", "species": [A, B], "shells": "0-6" or "7-16"}')
        sp = request['species']
        if not (isinstance(sp, list) and len(sp) == 2 and all(isinstance(s, str) and s in SPECIES for s in sp)):
            raise ValueError('species: two dopant symbols from the catalog (may be equal)')
        if request['shells'] not in WINDOWS:
            raise ValueError('shells must be one of %s' % sorted(WINDOWS))
        return 1
    if kind == 'monte_carlo':
        return check_mc(request)
    raise ValueError('kind must be "pair" or "monte_carlo"')


def run(request, workdir):
    if request['kind'] == 'pair':
        key = pair_key(*request['species'])
        lo, hi = WINDOWS[request['shells']]
        rows = [r for r in DATA['rows'].get(key, []) if lo <= r['i_shell'] <= hi]
        return {'pair': key, 'shells': request['shells'], 'rows': rows, 'n_rows': len(rows),
                'source': 'rows of the archived Jij_table.csv, verbatim (no new calculation)'}
    spins = request['spins']
    names = sorted(spins, key=lambda s: SPECIES[s])
    conc = {n: float(spins[n]) for n in names}
    J = {(a, b): [float(x) for x in request['exchange_meV'][pair_key(a, b)]] for a in names for b in names}
    size, samples = request['size'], request['samples']
    seed = request.get('seed', 0)
    base = int(hashlib.sha256(f'{seed}:{size}'.encode()).hexdigest()[:8], 16)
    seeds = [base + i for i in range(samples)]
    res, temps = dilute.run_many(size, conc, J, seeds, request['temperatures_K'], MC['equilibration_sweeps'],
                                 MC['measurement_sweeps'], MC['overrelaxation'], workdir=str(workdir),
                                 nworkers=min(MC['threads'], dilute.workers()))
    data = np.array([r['data'] for r in res])               # (samples, T, 14)
    labels = {'total': 0, names[0]: 2}
    if len(names) == 2:
        labels[names[1]] = 3
    box = dilute.bi_sites(size)[1]
    out = {'size': size, 'n_bi_sites': int(336 * size ** 3), 'box_A': box.round(4).tolist(), 'samples': samples, 'seed': seed,
           'spins_per_sample': [r['n_spins'] for r in res], 'bonds_per_sample': [r['n_bonds'] for r in res],
           'temperatures_K': temps, 'energy_meV_per_spin': data[:, :, 1].mean(0).round(6).tolist(),
           'order_parameters': {}}
    for name, q in labels.items():
        m2, m4, m1 = data[:, :, 2 + 2 * q], data[:, :, 3 + 2 * q], data[:, :, 10 + q]
        out['order_parameters'][name] = {
            'm2_per_sample': m2.round(8).tolist(), 'm4_per_sample': m4.round(10).tolist(),
            'abs_m_per_sample': m1.round(8).tolist()}
    out['order_parameters']['total']['m2_kmin_per_sample'] = data[:, :, 14].round(8).tolist()
    if len(names) == 2:
        out['order_parameters'][f'{names[0]}-minus-{names[1]}'] = {
            'm2_per_sample': data[:, :, 4].round(8).tolist(), 'm4_per_sample': data[:, :, 5].round(10).tolist(),
            'abs_m_per_sample': data[:, :, 11].round(8).tolist()}
    out['definitions'] = DATA['mc_output_definitions']
    return out


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    targets, tol = PRIVATE['targets'], PRIVATE['tolerances']
    missing = [f for f in ('curie_temperature_K', 'best_codopant_for_Cr', 'evidence_ids') if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    tc = payload['curie_temperature_K']
    if not isinstance(tc, dict) or set(tc) != set(targets):
        return False, {'error': 'curie_temperature_K must give a value (number or null) for every composition id'}
    if not all(v is None or (finite(v) and v >= 0) for v in tc.values()):
        return False, {'error': 'curie_temperature_K values must be non-negative numbers or null'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    report = {}
    for cid, t in targets.items():
        v = tc[cid]
        if t['type'] == 'null':
            ok = v is None
        elif t['type'] == 'at_most':
            ok = v is not None and v <= t['max_K']
        else:
            ok = v is not None and abs(v / t['K'] - 1) <= tol['relative']
        checks[f'tc:{cid}'] = ok
        report[cid] = {'submitted': v, 'relative_error': (v / t['K'] - 1) if (t['type'] == 'value' and v is not None) else None}
    checks['best_codopant_for_Cr'] = payload['best_codopant_for_Cr'] == PRIVATE['best_codopant_for_Cr']
    report['checks'] = checks
    return all(checks.values()), report
