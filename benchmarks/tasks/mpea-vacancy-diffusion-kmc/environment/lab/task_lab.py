"""Lab logic for mpea-vacancy-diffusion-kmc.

Two request kinds:
  hops  - retrieval of author-archived DFT energies (two end states and the transition state) of up to 20 vacancy
          hops, with each hop's occupation of the 128-site bcc supercell (no new calculation);
  kmc   - kinetic Monte Carlo (kmc.c) of a barrier model the agent specifies (local-environment KRA expansion plus
          an optional pair Hamiltonian for the end states), on an L x L x L bcc supercell with one vacancy.
The 100 hops in `hops_to_predict` are never served; their barriers, and the reference diffusion quantities computed
by the builders from a model trained on every clean archived hop, grade the submission.
"""
import json
import math
import os
from pathlib import Path
import kmcservice

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
ELEMENTS = kmcservice.ELEMENTS
HOPS = DATA['hops']
LIM = DATA['kmc_limits']


def num(x, bound=50.0):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x) and abs(x) <= bound


def intval(x, lo, hi):
    return type(x) is int and lo <= x <= hi


def check_model(model):
    if not isinstance(model, dict) or set(model) - {'kra', 'end_state'} or 'kra' not in model:
        raise ValueError('model must be {"kra": {...}, "end_state": {...} or null}; see kmc_service in the catalog')
    kra = model['kra']
    if not isinstance(kra, dict) or set(kra) - {'orbits', 'K0', 'W', 'Q'} or 'K0' not in kra:
        raise ValueError('model.kra needs K0 and optional orbits, W, Q')
    n_orb = kra.get('orbits', 0)
    if not intval(n_orb, 0, len(DATA['orbits'])):
        raise ValueError(f'model.kra.orbits must be an integer from 0 to {len(DATA["orbits"])}')
    if not isinstance(kra['K0'], dict) or set(kra['K0']) != set(ELEMENTS) or not all(num(v) for v in kra['K0'].values()):
        raise ValueError('model.kra.K0 must give a finite value (eV) for each of %s' % ELEMENTS)
    if n_orb:
        W = kra.get('W')
        if not isinstance(W, dict) or set(W) != set(ELEMENTS) or not all(
                isinstance(W[e], list) and len(W[e]) == n_orb and all(
                    isinstance(r, list) and len(r) == 6 and all(num(v) for v in r) for r in W[e]) for e in ELEMENTS):
            raise ValueError('model.kra.W must map each hopping element to `orbits` rows of 6 finite numbers (eV)')
    elif kra.get('W'):
        raise ValueError('model.kra.W given but orbits is 0')
    if kra.get('Q') is not None:
        Q = kra['Q']
        if n_orb < 3:
            raise ValueError('model.kra.Q needs orbits >= 3 (it acts on the composition of orbits 1-3)')
        if not isinstance(Q, dict) or set(Q) != set(ELEMENTS) or not all(
                isinstance(Q[e], list) and len(Q[e]) == 6 and all(
                    isinstance(r, list) and len(r) == 6 and all(num(v) for v in r) for r in Q[e]) for e in ELEMENTS):
            raise ValueError('model.kra.Q must map each hopping element to a 6 x 6 matrix of finite numbers (eV)')
    es = model.get('end_state')
    if es is not None:
        if not isinstance(es, dict) or set(es) != {'shells', 'V'} or not intval(es['shells'], 1, len(DATA['shells'])):
            raise ValueError(f'model.end_state must be null or {{"shells": 1-{len(DATA["shells"])}, "V": [...]}}')
        V = es['V']
        if not (isinstance(V, list) and len(V) == es['shells'] and all(
                isinstance(m, list) and len(m) == 7 and all(isinstance(r, list) and len(r) == 7 and all(num(v) for v in r)
                                                            for r in m) for m in V)):
            raise ValueError('model.end_state.V must be `shells` matrices of 7 x 7 finite numbers (eV), order %s'
                             % kmcservice.SPECIES7)
        if any(abs(m[a][b] - m[b][a]) > 1e-9 for m in V for a in range(7) for b in range(7)):
            raise ValueError('model.end_state.V matrices must be symmetric')


def kmc_cost(req):
    n = 2 * req['cell'] ** 3
    work = req['trajectories'] * (req['steps'] + req.get('equilibration_sweeps', 0) * n / LIM['sweep_site_per_step'])
    return max(1, math.ceil(work / LIM['work_per_unit']))


def validate(request):
    kind = request.get('kind')
    if kind == 'hops':
        ids = request.get('ids')
        if set(request) != {'kind', 'ids'} or not isinstance(ids, list) or not 1 <= len(ids) <= DATA['hops_per_purchase'] \
                or len(set(ids)) != len(ids) or not all(isinstance(i, str) and i in HOPS for i in ids):
            raise ValueError('Use {"kind": "hops", "ids": [1-%d distinct ids from hops_available]}' % DATA['hops_per_purchase'])
        return 1
    if kind == 'kmc':
        allowed = {'kind', 'composition', 'temperature_K', 'cell', 'steps', 'trajectories', 'equilibration_sweeps',
                   'checkpoints', 'seed', 'model'}
        need = {'composition', 'temperature_K', 'cell', 'steps', 'trajectories', 'model'}
        if set(request) - allowed or not need <= set(request):
            raise ValueError('kmc needs composition, temperature_K, cell, steps, trajectories, model '
                             '(optional equilibration_sweeps, checkpoints, seed); see kmc_service in the catalog')
        comp = request['composition']
        if not isinstance(comp, dict) or not comp or set(comp) - set(ELEMENTS) or not all(
                num(v, 1.0) and v >= 0 for v in comp.values()) or sum(comp.values()) <= 0:
            raise ValueError('composition: {element: fraction} with elements from %s' % ELEMENTS)
        if not num(request['temperature_K'], 1e5) or not LIM['temperature_K'][0] <= request['temperature_K'] <= LIM['temperature_K'][1]:
            raise ValueError('temperature_K must lie in %s' % LIM['temperature_K'])
        if not intval(request['cell'], *LIM['cell']):
            raise ValueError('cell (conventional cubes per edge) must be an integer in %s' % LIM['cell'])
        if not intval(request['steps'], *LIM['steps']) or not intval(request['trajectories'], *LIM['trajectories']):
            raise ValueError('steps must be in %s and trajectories in %s' % (LIM['steps'], LIM['trajectories']))
        if not intval(request.get('equilibration_sweeps', 0), *LIM['equilibration_sweeps']):
            raise ValueError('equilibration_sweeps must be an integer in %s' % LIM['equilibration_sweeps'])
        if not intval(request.get('checkpoints', 0), *LIM['checkpoints']):
            raise ValueError('checkpoints must be an integer in %s' % LIM['checkpoints'])
        if not intval(request.get('seed', 0), 0, 10 ** 6):
            raise ValueError('seed must be an integer in [0, 1000000]')
        counts = kmcservice.counts_for(comp, request['cell'])
        if min(c for c, e in zip(counts, ELEMENTS) if comp.get(e, 0) > 0) < 1:
            raise ValueError('a listed element rounds to zero atoms in this cell')
        check_model(request['model'])
        cost = kmc_cost(request)
        if cost > LIM['max_cost']:
            raise ValueError(f'this run would cost {cost} units; the maximum per run is {LIM["max_cost"]}')
        return cost
    raise ValueError('kind must be "hops" or "kmc"')


def run(request, workdir):
    if request['kind'] == 'hops':
        return {'hops': {i: HOPS[i] for i in request['ids']},
                'source': 'author-archived DFT energies (Lee & Natarajan 2026), retrieved; no new calculation'}
    workers = max(1, min(4, os.cpu_count() or 1))
    traj, pairs = kmcservice.run(request, DATA['orbits'], DATA['shells'], workers)
    return {'composition_counts': dict(zip(ELEMENTS, kmcservice.counts_for(request['composition'], request['cell']))),
            'sites': 2 * request['cell'] ** 3, 'temperature_K': request['temperature_K'],
            'trajectories': traj, 'nn_pair_counts_after_equilibration': {'order': ELEMENTS, 'counts': pairs},
            'definitions': DATA['kmc_output_definitions']}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    tg, tol = PRIVATE['targets'], PRIVATE['tolerances']
    fields = ['barrier_predictions_eV', 'log10_eta_vacancy', 'fastest_tracer_B', 'log10_tracer_ratio_B', 'activation_energy_C_eV', 'log10_eta_thermo_ideal_B', 'log10_eta_kra_simple_B',
              'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    bp, eta = payload['barrier_predictions_eV'], payload['log10_eta_vacancy']
    hb = PRIVATE['hidden_barriers_eV']
    if not isinstance(bp, dict) or set(bp) != set(hb) or not all(finite(v) for v in bp.values()):
        return False, {'error': 'barrier_predictions_eV must give a finite barrier (eV) for every id in hops_to_predict'}
    if not isinstance(eta, dict) or set(eta) != set(tg['log10_eta_vacancy']) or not all(finite(v) for v in eta.values()):
        return False, {'error': 'log10_eta_vacancy must give a finite number for each composition id %s'
                       % sorted(tg['log10_eta_vacancy'])}
    scalars = ['log10_tracer_ratio_B', 'activation_energy_C_eV', 'log10_eta_thermo_ideal_B', 'log10_eta_kra_simple_B']
    if not all(finite(payload[k]) for k in scalars):
        return False, {'error': '%s must be finite numbers' % scalars}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    rmse = math.sqrt(sum((bp[k] - v) ** 2 for k, v in hb.items()) / len(hb))
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed,
              'barrier_predictions': rmse <= tol['barrier_rmse_eV']}
    report = {'barrier_rmse_eV': rmse}
    for c, ref in tg['log10_eta_vacancy'].items():
        checks[f'log10_eta_vacancy:{c}'] = abs(eta[c] - ref) <= tol['log10_eta']
        report[f'log10_eta_vacancy:{c}'] = eta[c] - ref
    checks['fastest_tracer_B'] = payload['fastest_tracer_B'] == tg['fastest_tracer_B']
    for k in ('log10_tracer_ratio_B', 'log10_eta_thermo_ideal_B', 'log10_eta_kra_simple_B'):
        checks[k] = abs(payload[k] - tg[k]) <= tol[k]
        report[k] = payload[k] - tg[k]
    checks['activation_energy_C_eV'] = abs(payload['activation_energy_C_eV'] - tg['activation_energy_C_eV']) <= tol['activation_energy_eV']
    report['activation_energy_C_eV'] = payload['activation_energy_C_eV'] - tg['activation_energy_C_eV']
    report['checks'] = checks
    return all(checks.values()), report
