"""Lab logic for dense-hydrogen-llt-crossover: retrieval of author-released MD/PIMD results and grading.

Purchases return archived data of one released run: the stitched LAMMPS thermodynamic series of a classical
MLIP NVT run (with its input deck and job history), structural observables computed from frames sampled out of
that run's released trajectory, or the PIMD thermodynamic series of every replica at one density. Nothing is
computed from the submission; the graded values were derived at construction time and live in private/.
"""
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).parent
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
META = json.loads((HERE/'data/meta.json').read_text())
THERMO = np.load(HERE/'data/thermo.npz')
PIMD = np.load(HERE/'data/pimd.npz')
STRUCT = json.loads((HERE/'data/structure.json').read_text())
COST = {'thermo': {128: 1, 256: 1, 512: 2, 2048: 2}, 'structure': {128: 2, 256: 2, 512: 3, 2048: 4}, 'pimd': 1}
SOURCE = 'retrieved author-released data (Zenodo 14790037); no new simulation'


def run_key(request, kinds):
    usage = '{"kind": "%s", "N": 128|256|512|2048, "T_K": <int>, "rs": "<string exactly as in the catalog>"}' % request.get('kind')
    if set(request) != {'kind', 'N', 'T_K', 'rs'}:
        raise ValueError('Use ' + usage)
    n, t, rs = request['N'], request['T_K'], request['rs']
    if type(n) is not int or type(t) is not int or not isinstance(rs, str):
        raise ValueError('Use ' + usage)
    key = f'{n}/{t}/{rs}'
    if key not in kinds:
        raise ValueError(f'No released {request["kind"]} data for N={n}, T={t} K, rs="{rs}" (see the catalog)')
    return key, n


def validate(request):
    kind = request.get('kind')
    if kind == 'thermo':
        return COST['thermo'][run_key(request, META['runs'])[1]]
    if kind == 'structure':
        return COST['structure'][run_key(request, STRUCT)[1]]
    if kind == 'pimd':
        if set(request) != {'kind', 'T_label', 'rs'} or f"{request.get('T_label')}/{request.get('rs')}" not in META['pimd_sets']:
            raise ValueError('Use {"kind": "pimd", "T_label": "900K"|"1000K", "rs": "<string exactly as in the catalog>"}')
        return COST['pimd']
    raise ValueError('kind must be "thermo", "structure" or "pimd"')


def run(request, workdir):
    kind = request['kind']
    if kind == 'thermo':
        key, _ = run_key(request, META['runs'])
        m = META['runs'][key]
        a = THERMO[key].astype(float)
        return {'N': m['N'], 'T_K': m['T'], 'rs': m['rs'], 'volume_A3': m['volume_A3'],
                'timestep_ps': 0.0002, 'first_step': m['first_step'], 'step_interval': m['stride_steps'],
                'columns': ['Temp [K]', 'PotEng [eV]', 'TotEng [eV]', 'Press [bar]'],
                'rows': np.c_[a[:, 0].round(2), (a[:, 1] + m['E_ref_eV']).round(4), (a[:, 2] + m['E_ref_eV']).round(4),
                              a[:, 3].round(0)].tolist(),
                'input_deck': m['input_deck'], 'jobs': m['jobs'],
                'note': 'rows are every 5th thermo output (every 500 steps); all job outputs stitched by step',
                'source': SOURCE}
    if kind == 'structure':
        key, _ = run_key(request, STRUCT)
        return {'N': request['N'], 'T_K': request['T_K'], 'rs': request['rs'], **STRUCT[key], 'source': SOURCE}
    s = META['pimd_sets'][f"{request['T_label']}/{request['rs']}"]
    reps = {str(r): PIMD[k].astype(float).round(6).tolist() for r, k in s['replicas'].items()}
    return {'T_label': request['T_label'], 'rs': request['rs'],
            'columns': ['nmove', 'H [Ha]', 'E_pot [Ha]', 'Temp [K]', 'quantum_kin_virial [Ha]', 'pressure [GPa]'],
            'row_interval_moves': 20, 'replicas': reps, 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    tg, th = PRIVATE['targets'], PRIVATE['thresholds']
    fields = ['classification', 'low_density_phase_900K', 'widom_rs_N2048', 'pimd_crossover_rs_1000K',
              'nqe_crossover_pressure_shift_GPa', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    cl, w, nq = payload['classification'], payload['widom_rs_N2048'], payload['nqe_crossover_pressure_shift_GPa']
    if not isinstance(cl, dict) or set(cl) != set(tg['classification']) or not all(v in ('first_order', 'crossover') for v in cl.values()):
        return False, {'error': 'classification must map "900", "950", "1000", "1100" to "first_order" or "crossover"'}
    if payload['low_density_phase_900K'] not in ('crystal', 'liquid'):
        return False, {'error': 'low_density_phase_900K must be "crystal" or "liquid"'}
    if not isinstance(w, dict) or set(w) != set(tg['widom_rs_N2048']) or not all(finite(v) for v in w.values()):
        return False, {'error': 'widom_rs_N2048 must give a finite rs for "1100", "1200" and "1400"'}
    if not finite(payload['pimd_crossover_rs_1000K']):
        return False, {'error': 'pimd_crossover_rs_1000K must be a finite number'}
    if not isinstance(nq, dict) or set(nq) != {'900', '1000'} or not (nq['900'] is None or finite(nq['900'])) \
            or not (nq['1000'] is None or finite(nq['1000'])):
        return False, {'error': 'nqe_crossover_pressure_shift_GPa must map "900" and "1000" to a number or null'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    checks = {'evidence_executed': evidence}
    report = {'widom_error': {}, 'classification': {}}
    for t, ref in tg['classification'].items():
        checks[f'classification_{t}K'] = cl[t] == ref
    checks['low_density_phase_900K'] = payload['low_density_phase_900K'] == tg['low_density_phase_900K']
    for t, ref in tg['widom_rs_N2048'].items():
        report['widom_error'][t] = w[t] - ref
        checks[f'widom_rs_{t}K'] = abs(w[t] - ref) <= th['widom_rs']
    report['pimd_rs_error'] = payload['pimd_crossover_rs_1000K'] - tg['pimd_crossover_rs_1000K']
    checks['pimd_crossover_rs_1000K'] = abs(report['pimd_rs_error']) <= th['pimd_rs']
    ref1000 = tg['nqe_shift_GPa']['1000']
    report['nqe_shift_1000_error_GPa'] = None if nq['1000'] is None else nq['1000'] - ref1000
    checks['nqe_shift_1000K'] = nq['1000'] is not None and abs(nq['1000'] - ref1000) <= th['nqe_shift_GPa']
    checks['nqe_shift_900K_undetermined'] = nq['900'] is None
    report['checks'] = checks
    return all(checks.values()), report
