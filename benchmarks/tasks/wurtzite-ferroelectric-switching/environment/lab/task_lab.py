"""Lab logic for wurtzite-ferroelectric-switching: retrieval of author-archived QE results and grading.

Purchases return released neb.x run summaries, final NEB image structures, Berry-phase (lberry) outputs,
the authors' Born-charge file and vc-relax results. No new calculation is run. Barriers, polarizations and
mechanism descriptors are never served; they grade the submission.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
COSTS = {'neb_run': 2, 'neb_image': 1, 'berry_phase': 1, 'born_charges': 1, 'vc_relax': 1}
SOURCE = 'retrieved author-archived Quantum ESPRESSO result (no new calculation)'


def _system(request):
    s = request.get('system')
    if s not in ('AlN', 'AlBN'):
        raise ValueError('system must be "AlN" or "AlBN"')
    return s


def validate(request):
    kind = request.get('kind')
    if kind not in COSTS:
        raise ValueError('kind must be one of %s' % sorted(COSTS))
    if kind == 'neb_run':
        if set(request) != {'kind', 'system', 'run'}:
            raise ValueError('Use {"kind": "neb_run", "system": ..., "run": ...}')
        if request['run'] not in DATA['neb_runs'][_system(request)]:
            raise ValueError('run must be one of %s' % sorted(DATA['neb_runs']['AlN']))
    elif kind == 'neb_image':
        if set(request) != {'kind', 'system', 'run', 'image'}:
            raise ValueError('Use {"kind": "neb_image", "system": ..., "run": ..., "image": <1..N>}')
        runs = DATA['neb_images'][_system(request)]
        if request['run'] not in runs:
            raise ValueError('run must be one of %s' % sorted(runs))
        n = len(runs[request['run']])
        if type(request['image']) is not int or not 1 <= request['image'] <= n:
            raise ValueError('image must be an integer from 1 to %d for this run' % n)
    elif kind == 'berry_phase':
        if set(request) != {'kind', 'system', 'folder'}:
            raise ValueError('Use {"kind": "berry_phase", "system": ..., "folder": "NN_image"}')
        if request['folder'] not in DATA['berry'][_system(request)]:
            raise ValueError('folder must be one of berry_phase_folders in the catalog')
    elif kind == 'born_charges':
        if set(request) != {'kind', 'system'} or request.get('system') != 'AlBN':
            raise ValueError('Use {"kind": "born_charges", "system": "AlBN"} (only AlBN is archived)')
    elif kind == 'vc_relax':
        if set(request) != {'kind', 'key'} or request['key'] not in DATA['vc_relax']:
            raise ValueError('key must be one of %s' % sorted(DATA['vc_relax']))
    return COSTS[kind]


def run(request, workdir):
    kind = request['kind']
    if kind == 'neb_run':
        out = DATA['neb_runs'][request['system']][request['run']]
        return {'system': request['system'], 'run': request['run'], **out, 'source': SOURCE}
    if kind == 'neb_image':
        s, r, i = request['system'], request['run'], request['image']
        img = DATA['neb_images'][s][r][i - 1]
        return {'system': s, 'run': r, 'image': i, 'cell_A': DATA['neb_runs'][s][r]['cell_A'], **img,
                'units': 'Cartesian positions in angstrom', 'source': SOURCE}
    if kind == 'berry_phase':
        out = DATA['berry'][request['system']][request['folder']]
        return {'system': request['system'], 'folder': request['folder'], **out, 'source': SOURCE}
    if kind == 'born_charges':
        return {'system': 'AlBN', 'file': 'Born effective charge (AlBN).csv', 'rows': DATA['born_charges']['AlBN'],
                'note': 'Atom order as in the AlBN NEB inputs. The archive does not state for which structure these '
                        'charges were computed.', 'source': SOURCE}
    out = DATA['vc_relax'][request['key']]
    return {'key': request['key'], **out, 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    ref, t = PRIVATE['reference'], PRIVATE['thresholds']
    fields = ['barrier_meV_per_fu', 'AlBN_collective_barrier_meV_per_fu', 'spontaneous_polarization_C_m2',
              'AlBN_saddle_polarization_change_C_m2', 'AlBN_intermediate', 'AlBN_switching_mode', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    b, ps, inter = payload['barrier_meV_per_fu'], payload['spontaneous_polarization_C_m2'], payload['AlBN_intermediate']
    pair_ok = lambda d: isinstance(d, dict) and set(d) == {'AlN', 'AlBN'} and all(finite(v) for v in d.values())
    if not pair_ok(b) or not pair_ok(ps):
        return False, {'error': 'barrier_meV_per_fu and spontaneous_polarization_C_m2 need finite "AlN" and "AlBN" values'}
    if not (isinstance(inter, dict) and set(inter) == {'energy_meV_per_fu', 'inverted_cations'}
            and finite(inter['energy_meV_per_fu']) and finite(inter['inverted_cations'])):
        return False, {'error': 'AlBN_intermediate needs finite energy_meV_per_fu and inverted_cations'}
    if not finite(payload['AlBN_collective_barrier_meV_per_fu']) or not finite(payload['AlBN_saddle_polarization_change_C_m2']):
        return False, {'error': 'AlBN_collective_barrier_meV_per_fu and AlBN_saddle_polarization_change_C_m2 must be finite numbers'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    rel = lambda x, r: abs(x / r - 1)
    errs = {'barrier_AlN_rel': rel(b['AlN'], ref['barrier_meV_per_fu']['AlN']),
            'barrier_AlBN_rel': rel(b['AlBN'], ref['barrier_meV_per_fu']['AlBN']),
            'collective_barrier_AlBN_rel': rel(payload['AlBN_collective_barrier_meV_per_fu'], ref['AlBN_collective_barrier_meV_per_fu']),
            'Ps_AlN_abs': abs(ps['AlN'] - ref['spontaneous_polarization_C_m2']['AlN']),
            'Ps_AlBN_abs': abs(ps['AlBN'] - ref['spontaneous_polarization_C_m2']['AlBN']),
            'saddle_dP_abs': abs(payload['AlBN_saddle_polarization_change_C_m2'] - ref['AlBN_saddle_polarization_change_C_m2']),
            'intermediate_energy_rel': rel(inter['energy_meV_per_fu'], ref['AlBN_intermediate']['energy_meV_per_fu'])}
    mode = payload['AlBN_switching_mode']
    checks = {'evidence_executed': evidence,
              'barrier_AlN': errs['barrier_AlN_rel'] <= t['barrier_rel_tol'],
              'barrier_AlBN': errs['barrier_AlBN_rel'] <= t['barrier_rel_tol'],
              'collective_barrier_AlBN': errs['collective_barrier_AlBN_rel'] <= t['barrier_rel_tol'],
              'Ps_AlN': errs['Ps_AlN_abs'] <= t['polarization_abs_tol_C_m2'],
              'Ps_AlBN': errs['Ps_AlBN_abs'] <= t['polarization_abs_tol_C_m2'],
              'saddle_polarization_change': errs['saddle_dP_abs'] <= t['saddle_dP_abs_tol_C_m2'],
              'intermediate_energy': errs['intermediate_energy_rel'] <= t['intermediate_energy_rel_tol'],
              'intermediate_inverted_cations': inter['inverted_cations'] == ref['AlBN_intermediate']['inverted_cations'],
              'switching_mode': isinstance(mode, str) and mode.strip().lower() == ref['AlBN_switching_mode']}
    return all(checks.values()), {'checks': checks, 'errors': errs}
