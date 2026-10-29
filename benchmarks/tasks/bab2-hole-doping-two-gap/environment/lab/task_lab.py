"""Lab logic for bab2-hole-doping-two-gap: author-archived BaB2 results plus a numpy isotropic Eliashberg solver.

Purchases return archived content of Materials Cloud record 0kqen-w8h90 (bab2-analysis.zip) under the archive's
folder names, or a new isotropic Migdal-Eliashberg solution on an archived alpha2F (eliashberg.py). The DFPT
results of sensitivity/3 and sensitivity/5 are never served. Graded values are derived quantities held in
private/hidden.json.
"""
import json
import math
from pathlib import Path
import numpy as np
import eliashberg

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
COSTS = DATA['costs']
RUNS = [f'sensitivity/{i}' for i in range(7)]
SOURCE = 'retrieved author-archived result (no new calculation)'
SOLVER_KEYS = {'kind', 'run', 'mu_star', 'cutoff_meV', 'temperatures_K', 'a2f_column', 'max_iter', 'conv_thr', 'mixing',
               'initial_gap_meV'}


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def run_of(request):
    run = request.get('run')
    if run not in RUNS:
        raise ValueError('run must be one of ' + ', '.join(RUNS))
    return DATA['runs'][run.split('/')[1]]


def validate(request):
    kind = request.get('kind')
    if kind in ('dfpt_phonons', 'a2f', 'iso_log'):
        if set(request) != {'kind', 'run'}:
            raise ValueError(f'Use {{"kind": "{kind}", "run": "sensitivity/<0-6>"}}')
        entry = run_of(request)
        if kind == 'dfpt_phonons' and entry['dfpt'] is None:
            raise ValueError(f"the DFPT results of {request['run']} are withheld")
        return COSTS[kind]
    if kind == 'aniso_gap':
        if set(request) != {'kind', 'run', 'temperature_K', 'axis'}:
            raise ValueError('Use {"kind": "aniso_gap", "run": "sensitivity/<0-6>", "temperature_K": <archived T>, '
                             '"axis": "imaginary" or "real"}')
        entry = run_of(request)
        if request['axis'] not in ('imaginary', 'real'):
            raise ValueError('axis must be "imaginary" or "real"')
        if not number(request['temperature_K']) or f"{float(request['temperature_K']):g}" not in entry['aniso']:
            raise ValueError(f"archived temperatures of {request['run']} (K): " +
                             ', '.join(sorted(entry['aniso'], key=float)))
        return COSTS[kind]
    if kind == 'eliashberg_iso':
        if not {'kind', 'run', 'mu_star', 'cutoff_meV', 'temperatures_K'} <= set(request) <= SOLVER_KEYS:
            raise ValueError('Use {"kind": "eliashberg_iso", "run": "sensitivity/<0-6>", "mu_star": <0-0.5>, '
                             '"cutoff_meV": <100-2000>, "temperatures_K": [up to 8 values, 5-150]} with optional '
                             'a2f_column, max_iter, conv_thr, mixing, initial_gap_meV')
        run_of(request)
        if not number(request['mu_star']) or not 0 <= request['mu_star'] <= 0.5:
            raise ValueError('mu_star must be between 0 and 0.5')
        if not number(request['cutoff_meV']) or not 100 <= request['cutoff_meV'] <= 2000:
            raise ValueError('cutoff_meV must be between 100 and 2000')
        temps = request['temperatures_K']
        if (not isinstance(temps, list) or not 1 <= len(temps) <= 8
                or not all(number(t) and 5 <= t <= 150 for t in temps)):
            raise ValueError('temperatures_K must be a list of 1 to 8 temperatures between 5 and 150 K')
        col = request.get('a2f_column', 1)
        if type(col) is not int or not 1 <= col <= 10:
            raise ValueError('a2f_column must be an integer from 1 to 10')
        it = request.get('max_iter', 500)
        if type(it) is not int or not 1 <= it <= 5000:
            raise ValueError('max_iter must be an integer from 1 to 5000')
        thr = request.get('conv_thr', 1e-3)
        if not number(thr) or not 1e-10 <= thr <= 0.1:
            raise ValueError('conv_thr must be between 1e-10 and 0.1')
        mix = request.get('mixing', 0.4)
        if not number(mix) or not 0 < mix <= 1:
            raise ValueError('mixing must be in (0, 1]')
        gap = request.get('initial_gap_meV', 5.0)
        if not number(gap) or not 0 < gap <= 50:
            raise ValueError('initial_gap_meV must be in (0, 50]')
        return COSTS[kind]
    if kind == 'supercell_phonons':
        if set(request) != {'kind', 'run'} or request['run'] not in DATA['supercells']:
            raise ValueError('Use {"kind": "supercell_phonons", "run": "331" or "331_cs"}')
        return COSTS[kind]
    raise ValueError('kind must be "dfpt_phonons", "a2f", "iso_log", "aniso_gap", "eliashberg_iso" or "supercell_phonons"')


def run(request, workdir):
    kind = request['kind']
    if kind == 'supercell_phonons':
        return {**DATA['supercells'][request['run']], 'source': SOURCE}
    entry = run_of(request)
    if kind == 'dfpt_phonons':
        return {**entry['header'], **entry['dfpt'], 'source': SOURCE}
    if kind == 'a2f':
        return {**entry['header'], **entry['a2f'], 'source': SOURCE}
    if kind == 'iso_log':
        return {**entry['header'], **entry['iso_log'], 'source': SOURCE}
    if kind == 'aniso_gap':
        T = f"{float(request['temperature_K']):g}"
        a = entry['aniso'][T]
        stem = 'imag' if request['axis'] == 'imaginary' else 'pade'
        return {**entry['header'], 'epw_int.in': entry['aniso_in'], 'epw_int.out (head)': entry['aniso_log_head'],
                'temperature_K': float(T), 'axis': request['axis'], 'file': f'epw/bab2.{stem}_aniso_gap0_{float(T):06.2f}',
                'first_line': a[request['axis']]['header'], 'rows': a[request['axis']]['rows'],
                'epw_int.out (lines for this temperature)': a['log'], 'source': SOURCE}
    col = request.get('a2f_column', 1)
    w = np.array(entry['solver_a2f']['w_meV'])
    a2f = np.array(entry['solver_a2f']['columns'][col - 1])
    results = [eliashberg.solve(w, a2f, float(T), float(request['mu_star']), float(request['cutoff_meV']),
                                max_iter=request.get('max_iter', 500), conv_thr=float(request.get('conv_thr', 1e-3)),
                                mixing=float(request.get('mixing', 0.4)),
                                initial_gap=float(request.get('initial_gap_meV', 5.0)))
               for T in request['temperatures_K']]
    return {'run': request['run'], 'a2f_column': col, 'phonon_smearing_meV': entry['solver_a2f']['smearing_meV'][col - 1],
            'lambda': eliashberg.coupling(w, a2f), 'mu_star': request['mu_star'], 'cutoff_meV': request['cutoff_meV'],
            'max_iter': request.get('max_iter', 500), 'conv_thr': request.get('conv_thr', 1e-3),
            'mixing': request.get('mixing', 0.4), 'initial_gap_meV': request.get('initial_gap_meV', 5.0),
            'results': results, 'source': 'new numpy solution of the isotropic Eliashberg equations on the archived alpha2F'}


def grade(payload, trace, validation):
    t, th = PRIVATE['targets'], PRIVATE['thresholds']
    fields = ['critical_hole_doping', 'iso_tc_K', 'aniso_tc_K', 'gap_ratio_25K', 'soft_mode_cm1', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    iso, soft = payload['iso_tc_K'], payload['soft_mode_cm1']
    if not isinstance(iso, dict) or set(iso) != set(t['iso_tc_K']) or not all(number(v) for v in iso.values()):
        return False, {'error': 'iso_tc_K must give a finite Tc for "0.065" and "0.10"'}
    if not isinstance(soft, dict) or set(soft) != set(t['soft_mode_cm1']) or not all(number(v) for v in soft.values()):
        return False, {'error': 'soft_mode_cm1 must give finite values for "homogeneous_one_ninth" and "cs_substituted"'}
    if not all(number(payload[k]) for k in ('critical_hole_doping', 'aniso_tc_K', 'gap_ratio_25K')):
        return False, {'error': 'critical_hole_doping, aniso_tc_K and gap_ratio_25K must be finite numbers'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {'critical_hole_doping': payload['critical_hole_doping'] - t['critical_hole_doping']}
    checks['critical_hole_doping'] = abs(errors['critical_hole_doping']) <= th['critical_doping_abs']
    for k, ref in t['iso_tc_K'].items():
        errors[f'iso_tc:{k}'] = iso[k] / ref - 1
        checks[f'iso_tc:{k}'] = abs(errors[f'iso_tc:{k}']) <= th['iso_tc_rel']
    errors['aniso_tc'] = payload['aniso_tc_K'] / t['aniso_tc_K'] - 1
    checks['aniso_tc'] = abs(errors['aniso_tc']) <= th['aniso_tc_rel']
    errors['gap_ratio'] = payload['gap_ratio_25K'] / t['gap_ratio_25K'] - 1
    checks['gap_ratio'] = abs(errors['gap_ratio']) <= th['gap_ratio_rel']
    for k, tol in (('homogeneous_one_ninth', th['homogeneous_soft_mode_abs_cm1']),
                   ('cs_substituted', th['cs_soft_mode_abs_cm1'])):
        errors[f'soft_mode:{k}'] = soft[k] - t['soft_mode_cm1'][k]
        checks[f'soft_mode:{k}'] = abs(errors[f'soft_mode:{k}']) <= tol
    return all(checks.values()), {'checks': checks, 'errors': errors,
                                  'error_units': 'critical doping and soft modes absolute; Tc and gap ratio relative'}
