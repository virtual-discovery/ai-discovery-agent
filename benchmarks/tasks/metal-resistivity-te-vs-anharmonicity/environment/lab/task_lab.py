"""Lab logic for metal-resistivity-te-vs-anharmonicity: archived EPW transport results plus numpy phonons.

`transport` returns an author-archived EPW ibte calculation (Materials Cloud record zxf0b-rs406) as parsed from its
input and output files; the runs listed in data.json['withheld'] are never served. `phonons` computes harmonic
phonons of an archived DFPT or TDEP force-constant set (phonons.py; identical to matdyn.x with asr = 'simple').
Graded values are archived iterative-BTE resistivities (or ratios of them) held in private/hidden.json.
"""
import json
import math
from pathlib import Path
import phonons

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
COSTS = DATA['costs']
SOURCE = 'retrieved author-archived EPW result (no new calculation)'


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def validate(request):
    kind = request.get('kind')
    if kind == 'transport':
        if set(request) != {'kind', 'run'}:
            raise ValueError('Use {"kind": "transport", "run": "<run from transport_runs with served = true>"}')
        run = request['run']
        if run in DATA['withheld']:
            raise ValueError(f'{run} is withheld by this task')
        if run in DATA['no_results']:
            raise ValueError(f'{run}: the archived output contains no transport results')
        if run not in DATA['runs']:
            raise ValueError('unknown run; see transport_runs in /app/catalog.json')
        return COSTS['transport']
    if kind == 'phonons':
        if not {'kind', 'ifc', 'grid'} <= set(request) <= {'kind', 'ifc', 'grid', 'eigenvectors'}:
            raise ValueError('Use {"kind": "phonons", "ifc": "<ifc id>", "grid": <4-20>, "eigenvectors": <optional bool>}')
        if request['ifc'] not in DATA['ifcs']:
            raise ValueError('unknown ifc; see force_constant_sets in /app/catalog.json')
        if type(request['grid']) is not int or not 4 <= request['grid'] <= 20:
            raise ValueError('grid must be an integer from 4 to 20')
        if not isinstance(request.get('eigenvectors', False), bool):
            raise ValueError('eigenvectors must be true or false')
        return COSTS['phonons']
    raise ValueError('kind must be "transport" or "phonons"')


def run(request, workdir):
    if request['kind'] == 'transport':
        r = DATA['runs'][request['run']]
        keys = ('id', 'metal', 'condition', 'archive_files', 'lattice_parameter_bohr', 'lattice_parameter_A', 'ibte_in',
                'fine_k_mesh', 'fine_q_mesh', 'n_ibte_iterations', 'ibte_max_error_last', 'results')
        out = {k: r[k] for k in keys}
        out['run'] = out.pop('id')
        out['units'] = {'conductivity': '(Ohm cm)^-1', 'fermi': 'eV', 'dos': 'states/Ry as printed by EPW'}
        out['source'] = SOURCE
        return out
    rec = DATA['ifcs'][request['ifc']]
    meta = {'ifc': rec['id'], 'metal': rec['metal'], 'kind': rec['kind'], 'archive_file': rec['archive_file'],
            'lattice_temperature_K': rec['lattice_temperature_K'], 'sampling_temperature_K': rec['sampling_temperature_K'],
            'mass_amu': rec['mass_amu'], 'q2r_mesh': rec['mesh'],
            'ifc_file_header': {'alat_bohr': rec['header_alat_bohr'], 'unit_cell_volume_bohr3': rec['header_volume_bohr3']}}
    if rec['kind'] == 'TDEP':
        meta['tdep_unit_cell_alat_bohr'] = rec['tdep_cell_alat_bohr']
        meta['tdep_iteration'] = rec['tdep_iteration']
    report = phonons.report(rec, request['grid'], eigenvectors=request.get('eigenvectors', False))
    return {**meta, **report, 'source': 'new numpy phonon calculation from the archived force constants'}


def grade(payload, trace, validation):
    t, th = PRIVATE['targets'], PRIVATE['thresholds']
    fields = ['al_te_change_percent_481K', 'pb_te_resistivity_610K_uohm_cm', 'pb_anharmonic_change_percent_0K_lattice',
              'pb_te_anharm_resistivity_510K_uohm_cm', 'nb_te_anharm_resistivity_2000K_uohm_cm', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    anh = payload['pb_anharmonic_change_percent_0K_lattice']
    if not isinstance(anh, dict) or set(anh) != {'410', '610'} or not all(number(v) for v in anh.values()):
        return False, {'error': 'pb_anharmonic_change_percent_0K_lattice must give finite values for "410" and "610"'}
    scalars = fields[:2] + fields[3:5]
    if not all(number(payload[k]) for k in scalars):
        return False, {'error': ', '.join(scalars) + ' must be finite numbers'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {'al_te_change_481K_pp': payload['al_te_change_percent_481K'] - t['al_te_change_percent_481K'],
              'pb_te_610K_rel': payload['pb_te_resistivity_610K_uohm_cm'] / t['pb_te_resistivity_610K'] - 1,
              'pb_te_anharm_510K_rel': payload['pb_te_anharm_resistivity_510K_uohm_cm'] / t['pb_te_anharm_resistivity_510K'] - 1,
              'nb_te_anharm_2000K_rel': payload['nb_te_anharm_resistivity_2000K_uohm_cm'] / t['nb_te_anharm_resistivity_2000K'] - 1}
    checks['al_te_change_481K'] = abs(errors['al_te_change_481K_pp']) <= th['al_te_change_abs_pp']
    checks['pb_te_610K'] = abs(errors['pb_te_610K_rel']) <= th['pb_te_610_rel']
    for T in ('410', '610'):
        errors[f'pb_anharmonic_{T}K_pp'] = anh[T] - t['pb_anharmonic_change_percent_0K_lattice'][T]
        checks[f'pb_anharmonic_{T}K'] = abs(errors[f'pb_anharmonic_{T}K_pp']) <= th['pb_anharm_abs_pp']
    checks['pb_te_anharm_510K'] = abs(errors['pb_te_anharm_510K_rel']) <= th['pb_te_anharm_510_rel']
    checks['nb_te_anharm_2000K'] = abs(errors['nb_te_anharm_2000K_rel']) <= th['nb_te_anharm_2000_rel']
    return all(checks.values()), {'checks': checks, 'errors': errors,
                                  'error_units': 'pp = percentage points; rel = relative error'}
