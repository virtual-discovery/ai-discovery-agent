"""Lab logic for fresh-inp-polar-phonons: fresh Quantum ESPRESSO 6.7 calculations for zincblende InP.

validate() checks and charges a request; run() writes pw.x/ph.x inputs, runs them with mpirun on the task's
cores (qe_driver.py) and returns the parsed physical results. Nothing is looked up. grade() compares the
terminal submission with reference values computed at construction time through this same code path
(private/hidden.json; see benchmarks/construction/fresh-inp-polar-phonons/).
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
FIELDS = ['a0_A', 'omega_TO_cm-1', 'omega_LO_cm-1', 'born_charge_In', 'eps_inf', 'eps_0',
          'dlnZ_dlnV', 'dln_eps_inf_dlnV', 'gruneisen_TO',
          'strained_omega_LO_q_z_cm-1', 'strained_omega_LO_q_x_cm-1']
CHOICES = {'LOTO_splitting_under_compression': ('grows', 'shrinks')}


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def validate(request):
    allowed = {'kind', 'lattice_A', 'strain_zz', 'kmesh', 'kshift', 'ecutwfc_Ry'}
    kind = request.get('kind')
    if kind not in COSTS:
        raise ValueError('kind must be one of %s' % sorted(COSTS))
    extra = set(request) - allowed
    if extra:
        raise ValueError('unknown fields %s; allowed: %s' % (sorted(extra), sorted(allowed)))
    for key in ('lattice_A', 'kmesh', 'kshift'):
        if key not in request:
            raise ValueError('missing field %s' % key)
    a, strain = request['lattice_A'], request.get('strain_zz', 0.0)
    if not number(a) or not LIMITS['lattice_A'][0] <= a <= LIMITS['lattice_A'][1]:
        raise ValueError('lattice_A must be a number in %s (Angstrom)' % LIMITS['lattice_A'])
    if not number(strain) or not LIMITS['strain_zz'][0] <= strain <= LIMITS['strain_zz'][1]:
        raise ValueError('strain_zz must be a number in %s' % LIMITS['strain_zz'])
    if request['kshift'] not in (0, 1) or isinstance(request['kshift'], bool):
        raise ValueError('kshift must be 0 (Gamma-centred) or 1 (half-step offset)')
    ecut = request.get('ecutwfc_Ry', 50)
    if ecut not in LIMITS['ecutwfc_Ry'] or isinstance(ecut, bool):
        raise ValueError('ecutwfc_Ry must be one of %s' % LIMITS['ecutwfc_Ry'])
    cell = 'strained' if strain != 0 else 'cubic'
    table = COSTS[kind][cell]
    key = str(request['kmesh'])
    if isinstance(request['kmesh'], bool) or key not in table:
        raise ValueError('kmesh for kind %s in a %s cell must be one of %s' % (kind, cell, sorted(map(int, table))))
    cost = table[key]
    if ecut != 50:
        if kind != 'scf' and (cell == 'strained' or request['kmesh'] > LIMITS['max_kmesh_at_ecut65'][kind]):
            raise ValueError('ecutwfc_Ry 65 is available for scf, and for %s in cubic cells up to kmesh %d'
                             % (kind, LIMITS['max_kmesh_at_ecut65'][kind]))
        cost = 2 * cost
    return int(cost)


def run(request, workdir):
    result = qe_driver.calculate(request, workdir)
    result['source'] = 'new Quantum ESPRESSO 6.7 calculation (PBE, SSSP 1.3 precision pseudopotentials)'
    return result


def grade(payload, trace, validation):
    ref, tol = PRIVATE['reference'], PRIVATE['tolerances']
    missing = [f for f in FIELDS + list(CHOICES) + ['evidence_ids'] if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    bad = [f for f in FIELDS if not number(payload[f])]
    if bad:
        return False, {'error': 'fields must be finite numbers: %s' % bad}
    bad = [f for f, options in CHOICES.items() if payload[f] not in options]
    if bad:
        return False, {'error': '%s must be one of %s' % (bad, [CHOICES[f] for f in bad])}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    checks, errors = {'evidence_executed': evidence}, {}
    for field in FIELDS:
        error = payload[field] - ref[field]
        if field in tol['relative']:
            error = error / abs(ref[field])
            ok = abs(error) <= tol['relative'][field]
        else:
            ok = abs(error) <= tol['absolute'][field]
        checks[field] = ok
        errors[field] = error
    for field in CHOICES:
        checks[field] = payload[field] == ref[field]
    anisotropy = payload['strained_omega_LO_q_x_cm-1'] - payload['strained_omega_LO_q_z_cm-1']
    errors['strained_LO_anisotropy'] = anisotropy - ref['strained_LO_anisotropy_cm-1']
    checks['strained_LO_anisotropy'] = abs(errors['strained_LO_anisotropy']) <= tol['absolute']['strained_LO_anisotropy']
    report = {'checks': checks, 'errors': errors,
              'note': 'absolute errors in the field units, except relative errors for %s' % sorted(tol['relative'])}
    return all(checks.values()), report
