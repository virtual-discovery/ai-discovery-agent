"""Lab logic for hydrogen-plasma-facing-materials: retrieval of archived Quantum ESPRESSO runs and grading.

A purchase returns one released run of Materials Cloud record gpxwh-dx645 (input settings, structure, energies,
convergence messages, run header) exactly as parsed from its pwscf.in / pwscf.out. Nothing derived from the
runs (formation, surface-binding or cohesive energies, site identities) is ever served; those grade the
submission against references computed from the same archived runs.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
SOURCE = 'retrieved author-released Quantum ESPRESSO run (no new calculation)'


def validate(request):
    if request.get('kind') != 'calculation':
        raise ValueError('kind must be "calculation"')
    if set(request) != {'kind', 'id'} or not isinstance(request['id'], str) or request['id'] not in DATA:
        raise ValueError('Use {"kind": "calculation", "id": <a calculation id listed in /app/catalog.json>}')
    return 1


def run(request, workdir):
    return {'id': request['id'], **DATA[request['id']], 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    ref, t = PRIVATE['reference'], PRIVATE['thresholds']
    fields = ['h_formation_energy_eV', 'distinct_site_energies_eV', 'preferred_site', 'surface_binding_energy_eV',
              'cohesive_energy_eV_per_atom', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    hf, ds, ps = payload['h_formation_energy_eV'], payload['distinct_site_energies_eV'], payload['preferred_site']
    sb, ce = payload['surface_binding_energy_eV'], payload['cohesive_energy_eV_per_atom']
    hmats = set(ref['h_formation_energy_eV'])
    if not isinstance(hf, dict) or set(hf) != hmats or not all(finite(v) for v in hf.values()):
        return False, {'error': 'h_formation_energy_eV must give a finite number for every hydrogen material'}
    if not isinstance(ds, dict) or set(ds) != hmats or not all(
            isinstance(v, list) and len(v) > 0 and all(finite(x) for x in v) for v in ds.values()):
        return False, {'error': 'distinct_site_energies_eV must give a non-empty list of numbers for every hydrogen material'}
    if not isinstance(ps, dict) or set(ps) != hmats or not all(
            isinstance(v, dict) and set(v) == {'nearest_host_element', 'nearest_host_distance_A'}
            and isinstance(v['nearest_host_element'], str) and finite(v['nearest_host_distance_A']) for v in ps.values()):
        return False, {'error': 'preferred_site must give nearest_host_element (string) and nearest_host_distance_A '
                                '(number) for every hydrogen material'}
    sref, cref = ref['surface_binding_energy_eV'], ref['cohesive_energy_eV_per_atom']
    if not isinstance(sb, dict) or set(sb) != set(sref) or not all(
            isinstance(sb[m], dict) and set(sb[m]) == set(sref[m]) and all(v is None or finite(v) for v in sb[m].values())
            for m in sref):
        return False, {'error': 'surface_binding_energy_eV must give a number or null for every listed material and element'}
    if not isinstance(ce, dict) or set(ce) != set(cref) or not all(v is None or finite(v) for v in ce.values()):
        return False, {'error': 'cohesive_energy_eV_per_atom must give a number or null for every listed material'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    report = {'h_formation_error_eV': {}, 'site_count': {}, 'site_energy_max_error_eV': {}, 'preferred_site': {},
              'surface_binding_error_eV': {}, 'cohesive_error_eV': {}}

    lenient = set(PRIVATE.get('null_also_accepted', []))

    def close(value, target, tol, key):
        if value is None:
            return target is None or key in lenient
        return target is not None and abs(value - target) <= tol

    for m in sorted(hmats):
        err = hf[m] - ref['h_formation_energy_eV'][m]
        report['h_formation_error_eV'][m] = err
        checks[f'h_formation_{m}'] = abs(err) <= t['h_formation_tol_eV']
        target, given = ref['distinct_site_energies_eV'][m], sorted(ds[m])
        report['site_count'][m] = [len(given), len(target)]
        same = len(given) == len(target)
        report['site_energy_max_error_eV'][m] = max(abs(a - b) for a, b in zip(given, target)) if same else None
        checks[f'distinct_sites_{m}'] = same and all(abs(a - b) <= t['site_energy_tol_eV'] for a, b in zip(given, target))
        pr, pg = ref['preferred_site'][m], ps[m]
        derr = pg['nearest_host_distance_A'] - pr['nearest_host_distance_A']
        report['preferred_site'][m] = {'element_ok': pg['nearest_host_element'].strip() == pr['nearest_host_element'],
                                       'distance_error_A': derr}
        checks[f'preferred_site_{m}'] = (pg['nearest_host_element'].strip() == pr['nearest_host_element']
                                         and abs(derr) <= t['distance_tol_A'])
    for m in sorted(sref):
        for el in sorted(sref[m]):
            given, target = sb[m][el], sref[m][el]
            report['surface_binding_error_eV'][f'{m}/{el}'] = (
                'null expected' if target is None and given is not None else
                None if given is None else given - target)
            checks[f'surface_binding_{m}_{el}'] = close(given, target, t['sbe_tol_eV'], f'{m}/{el}')
    for m in sorted(cref):
        given, target = ce[m], cref[m]
        report['cohesive_error_eV'][m] = ('null expected' if target is None and given is not None else
                                          None if given is None else given - target)
        checks[f'cohesive_{m}'] = close(given, target, t['cohesive_tol_eV'], m)
    report['checks'] = checks
    return all(checks.values()), report
