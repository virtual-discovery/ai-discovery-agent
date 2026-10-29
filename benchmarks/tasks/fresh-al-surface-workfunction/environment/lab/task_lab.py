"""Lab logic for fresh-al-surface-workfunction: fresh Quantum ESPRESSO 6.7 calculations for Al slabs and bulk Al.

validate() checks a structured request against the allowed ranges in public.json and returns its charge, which
grows with the cost of the calculation. run() writes the pw.x input (and the pp.x inputs for the planar-averaged
potentials), runs `mpirun -np <task cpus> pw.x` in the job folder (slab_qe.py) and returns the parsed physical
results. Nothing is looked up. grade() compares the terminal submission with reference values computed at
construction time through this same code path (private/hidden.json; see
benchmarks/construction/fresh-al-surface-workfunction/).
"""
import json
import math
from pathlib import Path

import slab_qe

HERE = Path(__file__).parent
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
PUBLIC = json.loads((HERE/'public.json').read_text())
LIM = PUBLIC['limits']
COST = PUBLIC['costs']
SLAB_KEYS = {'kind', 'face', 'layers', 'vacuum_A', 'kpoints', 'smearing_Ry', 'nbnd', 'relax'}
BULK_KEYS = {'kind', 'cell', 'face', 'kpoints', 'kz', 'smearing_Ry'}


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def integer(x):
    return type(x) is int


def check_smearing(request):
    sm = request.get('smearing_Ry', 0.02)
    if isinstance(sm, bool) or sm not in LIM['smearing_Ry']:
        raise ValueError('smearing_Ry (Marzari-Vanderbilt width) must be one of %s' % LIM['smearing_Ry'])


def validate(request):
    kind = request.get('kind')
    if kind == 'slab':
        extra = set(request) - SLAB_KEYS
        if extra:
            raise ValueError('unknown slab fields %s; allowed: %s' % (sorted(extra), sorted(SLAB_KEYS)))
        for key in ('face', 'layers', 'kpoints'):
            if key not in request:
                raise ValueError('missing field %s' % key)
        face, n, k = request['face'], request['layers'], request['kpoints']
        if face not in LIM['faces']:
            raise ValueError('face must be one of %s' % LIM['faces'])
        if not integer(k) or str(k) not in LIM['max_layers_by_kpoints']:
            raise ValueError('kpoints (n of the in-plane n x n x 1 Gamma-centred mesh) must be one of %s'
                             % sorted(map(int, LIM['max_layers_by_kpoints'])))
        nmax = LIM['max_layers_by_kpoints'][str(k)]
        if not integer(n) or not LIM['min_layers'] <= n <= nmax:
            raise ValueError('layers must be an integer from %d to %d at kpoints %d' % (LIM['min_layers'], nmax, k))
        vac = request.get('vacuum_A', 15.0)
        if not number(vac) or not LIM['vacuum_A'][0] <= vac <= LIM['vacuum_A'][1]:
            raise ValueError('vacuum_A must be a number in %s' % LIM['vacuum_A'])
        check_smearing(request)
        if 'nbnd' in request:
            nb = request['nbnd']
            if not integer(nb) or not 1.5 * n + 1 <= nb <= 4 * n + 8:
                raise ValueError('nbnd must be an integer from %d to %d for %d layers'
                                 % (math.ceil(1.5 * n + 1), 4 * n + 8, n))
        relax = request.get('relax', 'none')
        if relax not in LIM['relax']:
            raise ValueError('relax must be one of %s' % LIM['relax'])
        if relax != 'none' and (k not in LIM['relax_kpoints'] or not LIM['relax_layers'][0] <= n <= LIM['relax_layers'][1]):
            raise ValueError('relaxations need kpoints in %s and %d to %d layers'
                             % (LIM['relax_kpoints'], LIM['relax_layers'][0], LIM['relax_layers'][1]))
        units = max(1, int(math.ceil(n * COST['slab_kpoint_factor'][str(k)] / 4 - 1e-9)))
        return units * (COST['relax_factor'] if relax != 'none' else 1)
    if kind == 'bulk':
        extra = set(request) - BULK_KEYS
        if extra:
            raise ValueError('unknown bulk fields %s; allowed: %s' % (sorted(extra), sorted(BULK_KEYS)))
        cell, k = request.get('cell'), request.get('kpoints')
        lo, hi = LIM['bulk_kpoints']
        if cell == 'primitive':
            if 'face' in request or 'kz' in request:
                raise ValueError('a primitive bulk cell takes no face or kz')
            if not integer(k) or not lo <= k <= hi:
                raise ValueError('kpoints must be an integer from %d to %d (n x n x n Gamma-centred mesh)' % (lo, hi))
        elif cell == 'oriented':
            if request.get('face') not in LIM['faces']:
                raise ValueError('an oriented bulk cell needs face in %s' % LIM['faces'])
            if not integer(k) or str(k) not in LIM['max_layers_by_kpoints']:
                raise ValueError('kpoints (in-plane) must be one of %s' % sorted(map(int, LIM['max_layers_by_kpoints'])))
            if not integer(request.get('kz')) or not lo <= request['kz'] <= hi:
                raise ValueError('kz must be an integer from %d to %d' % (lo, hi))
        else:
            raise ValueError('cell must be "primitive" or "oriented"')
        check_smearing(request)
        return COST['bulk']
    raise ValueError('kind must be "slab" or "bulk"')


def run(request, workdir):
    req = dict(request)
    if req['kind'] == 'slab':
        req.setdefault('vacuum_A', 15.0)
        req.setdefault('relax', 'none')
    req.setdefault('smearing_Ry', 0.02)
    req['ecutwfc_Ry'] = 30
    result = slab_qe.calculate(req, workdir)
    result['source'] = 'new Quantum ESPRESSO 6.7 calculation (PBE, SSSP 1.3 precision Al PAW)'
    return result


FACES = ('111', '100')
THICK = ('111_N18', '100_N18')
ORDERINGS = ('phi111>phi100', 'phi100>phi111', 'unresolved')


def grade(payload, trace, validation):
    ref, tol = PRIVATE['reference'], PRIVATE['tolerances']
    fields = ['surface_energy_J_per_m2', 'work_function_eV', 'work_function_ordering', 'relaxation_d12_percent_111',
              'thick_slab_predictions', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    sig, phi, thick = payload['surface_energy_J_per_m2'], payload['work_function_eV'], payload['thick_slab_predictions']
    if not (isinstance(sig, dict) and isinstance(phi, dict) and all(number(sig.get(f)) and number(phi.get(f)) for f in FACES)):
        return False, {'error': 'surface_energy_J_per_m2 and work_function_eV need finite numbers for "111" and "100"'}
    if not (isinstance(thick, dict) and all(isinstance(thick.get(k), dict) and number(thick[k].get('free_energy_offset_eV'))
                                            and number(thick[k].get('work_function_eV')) for k in THICK)):
        return False, {'error': 'thick_slab_predictions needs free_energy_offset_eV and work_function_eV for %s' % list(THICK)}
    if payload['work_function_ordering'] not in ORDERINGS or not number(payload['relaxation_d12_percent_111']):
        return False, {'error': 'work_function_ordering must be one of %s and relaxation_d12_percent_111 a number'
                                % list(ORDERINGS)}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    err = {'sigma_' + f: sig[f] - ref['sigma'][f] for f in FACES}
    err.update({'phi_' + f: phi[f] - ref['phi'][f] for f in FACES})
    err['d12_111'] = payload['relaxation_d12_percent_111'] - ref['d12_111']
    for k in THICK:
        err['E_' + k] = thick[k]['free_energy_offset_eV'] - ref['thick'][k]['free_energy_offset_eV']
        err['phi_' + k] = thick[k]['work_function_eV'] - ref['thick'][k]['work_function_eV']
    checks = {'evidence_executed': evidence}
    for f in FACES:
        checks['surface_energy_' + f] = abs(err['sigma_' + f]) <= tol['sigma_J_per_m2']
        checks['work_function_' + f] = abs(err['phi_' + f]) <= tol['phi_eV']
    checks['work_function_ordering'] = payload['work_function_ordering'] == ref['ordering']
    checks['relaxation_d12_111'] = abs(err['d12_111']) <= tol['d12_percent']
    for k in THICK:
        checks['thick_slab_energy_' + k] = abs(err['E_' + k]) <= tol['thick_energy_eV'][k.split('_')[0]]
        checks['thick_slab_work_function_' + k] = abs(err['phi_' + k]) <= tol['thick_phi_eV']
    return all(checks.values()), {'checks': checks, 'errors': err}
