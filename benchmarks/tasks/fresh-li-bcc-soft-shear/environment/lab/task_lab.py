"""Lab logic for fresh-li-bcc-soft-shear: fresh Quantum ESPRESSO 6.7 pw.x calculations for (strained) bcc Li.

validate() checks and charges a request; run() writes the pw.x input for the requested cubic lattice parameter,
Voigt strain, k mesh, smearing width and cutoff, runs it with mpirun on the task's cores (qe_driver.py) and
returns the parsed physical results. Nothing is looked up. grade() compares the terminal submission with
reference values computed at construction time through this same code path (private/hidden.json; see
benchmarks/construction/fresh-li-bcc-soft-shear/).
"""
import json
import math
from pathlib import Path

import qe_driver

HERE = Path(__file__).parent
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
PUBLIC = json.loads((HERE/'public.json').read_text())
LIMITS = PUBLIC['limits']
COSTS = PUBLIC['costs']
CLASS_BY_OPS = {48: 'cubic', 16: 'tetragonal', 12: 'trigonal', 8: 'orthorhombic', 4: 'monoclinic', 2: 'triclinic'}


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def symmetry_class(voigt):
    ops = qe_driver.symmetry_operations(voigt)
    return CLASS_BY_OPS.get(ops, 'triclinic')


def validate(request):
    allowed = {'kind', 'lattice_A', 'strain_voigt', 'kmesh', 'smearing_Ry', 'ecutwfc_Ry'}
    extra = set(request) - allowed
    if extra:
        raise ValueError('unknown fields %s; allowed: %s' % (sorted(extra), sorted(allowed)))
    if request.get('kind', 'scf') != 'scf':
        raise ValueError('kind must be "scf" (the only calculation this lab runs)')
    for key in ('lattice_A', 'kmesh', 'smearing_Ry'):
        if key not in request:
            raise ValueError('missing field %s' % key)
    a = request['lattice_A']
    if not number(a) or not LIMITS['lattice_A'][0] <= a <= LIMITS['lattice_A'][1]:
        raise ValueError('lattice_A must be a number in %s (Angstrom)' % LIMITS['lattice_A'])
    voigt = request.get('strain_voigt', [0.0] * 6)
    if not isinstance(voigt, list) or len(voigt) != 6 or not all(number(e) for e in voigt):
        raise ValueError('strain_voigt must be a list of 6 numbers [e1, e2, e3, e4, e5, e6]')
    if max(abs(e) for e in voigt) > LIMITS['strain_component_max_abs']:
        raise ValueError('each strain component must satisfy |e| <= %g' % LIMITS['strain_component_max_abs'])
    sigma = request['smearing_Ry']
    if not number(sigma) or sigma not in LIMITS['smearing_Ry']:
        raise ValueError('smearing_Ry must be one of %s' % LIMITS['smearing_Ry'])
    ecut = request.get('ecutwfc_Ry', 40)
    if isinstance(ecut, bool) or ecut not in LIMITS['ecutwfc_Ry']:
        raise ValueError('ecutwfc_Ry must be one of %s' % LIMITS['ecutwfc_Ry'])
    cell = symmetry_class(voigt)
    table = COSTS[cell]
    kmesh = request['kmesh']
    if isinstance(kmesh, bool) or not isinstance(kmesh, int) or str(kmesh) not in table:
        raise ValueError('the requested strain leaves a lattice of %s symmetry; for it kmesh must be one of %s'
                         % (cell, sorted(map(int, table))))
    cost = table[str(kmesh)]
    if ecut != 40:
        if cost > LIMITS['max_cost_at_ecut60']:
            raise ValueError('ecutwfc_Ry 60 is available only for requests that cost at most %d units at 40 Ry'
                             % LIMITS['max_cost_at_ecut60'])
        cost *= 2
    return int(cost)


def run(request, workdir):
    result = qe_driver.calculate(request, workdir)
    result['symmetry_class_charged'] = symmetry_class(result['strain_voigt'])
    result['source'] = 'new Quantum ESPRESSO 6.7 pw.x calculation (PBE, SSSP 1.3 precision Li pseudopotential)'
    return result


def grade(payload, trace, validation):
    ref, tol = PRIVATE['reference'], PRIVATE['tolerances']
    numeric = list(tol['absolute']) + list(tol['relative'])
    booleans = list(PRIVATE.get('booleans', {}))
    missing = [f for f in numeric + booleans + ['evidence_ids'] if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    bad = [f for f in numeric if not number(payload[f])]
    if bad:
        return False, {'error': 'fields must be finite numbers: %s' % bad}
    bad = [f for f in booleans if not isinstance(payload[f], bool)]
    if bad:
        return False, {'error': 'fields must be true or false: %s' % bad}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    checks, errors = {'evidence_executed': evidence}, {}
    for field in numeric:
        error = payload[field] - ref[field]
        if field in tol['relative']:
            error = error / abs(ref[field])
            checks[field] = abs(error) <= tol['relative'][field]
        else:
            checks[field] = abs(error) <= tol['absolute'][field]
        errors[field] = error
    for field in booleans:
        checks[field] = payload[field] is PRIVATE['booleans'][field]
    report = {'checks': checks, 'errors': errors,
              'note': 'absolute errors in the field units, except relative errors for %s' % sorted(tol['relative'])}
    return all(checks.values()), report
