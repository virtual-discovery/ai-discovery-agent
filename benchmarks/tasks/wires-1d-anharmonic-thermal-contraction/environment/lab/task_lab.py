"""Lab logic for wires-1d-anharmonic-thermal-contraction: retrieval of author-archived SSCHA/DFT files and grading.

Purchases return archived files (dynamical matrices, DFT inputs/outputs) or verbatim excerpts of the SSCHA
minimization logs (Materials Cloud record k45dq-ebf35). The TaSe3 free-energy Hessians at 650 K and 700 K are
withheld. Graded values are computed from the archive (benchmarks/construction/<task>/calib.py) and held in
private/hidden.json.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
SOURCE = 'retrieved author-archived file (no new calculation)'
WIRES = ('CuC2', 'TaSe3', 'AuSe2')


def run_key(request):
    T, n = request.get('temperature_K'), request.get('configurations')
    if type(T) is not int or type(n) is not int:
        raise ValueError('temperature_K and configurations must be integers (see sscha_runs in /app/catalog.json)')
    return f'{T}K/{n}'


def validate(request):
    kind = request.get('kind')
    wire = request.get('wire')
    if kind in ('sscha_log', 'sscha_matrices', 'dfpt', 'dft_file') and wire not in WIRES:
        raise ValueError('wire must be one of CuC2, TaSe3, AuSe2')
    if kind == 'sscha_log':
        if set(request) != {'kind', 'wire', 'temperature_K', 'configurations'}:
            raise ValueError('Use {"kind": "sscha_log", "wire": ..., "temperature_K": ..., "configurations": ...}')
        r = DATA['runs'][wire].get(run_key(request))
        if r is None or r['log'] is None:
            raise ValueError('no archived log for this run; see sscha_runs in /app/catalog.json')
        return 1
    if kind == 'sscha_matrices':
        if set(request) != {'kind', 'wire', 'temperature_K', 'configurations', 'matrix'} or \
                request['matrix'] not in ('auxiliary', 'hessian'):
            raise ValueError('Use {"kind": "sscha_matrices", "wire": ..., "temperature_K": ..., "configurations": ..., '
                             '"matrix": "auxiliary" | "hessian"}')
        r = DATA['runs'][wire].get(run_key(request))
        if r is None:
            raise ValueError('no archived run with this temperature and configuration count; see sscha_runs')
        if r[request['matrix']] is None:
            raise ValueError(f'{request["matrix"]} matrices of this run are not available (not archived or withheld by this task)')
        return 1
    if kind == 'dfpt':
        if set(request) != {'kind', 'wire', 'asr'} or request['asr'] not in ('none', '1d'):
            raise ValueError('Use {"kind": "dfpt", "wire": ..., "asr": "none" | "1d"}')
        return 1
    if kind == 'dft_file':
        if set(request) != {'kind', 'wire', 'file'} or request['file'] not in DATA['dft_files'][wire]:
            raise ValueError('Use {"kind": "dft_file", "wire": ..., "file": <name from dft_files>}')
        return 1
    if kind == 'strained_calculation':
        if set(request) != {'kind'}:
            raise ValueError('Use {"kind": "strained_calculation"}')
        return 1
    raise ValueError('kind must be sscha_log, sscha_matrices, dfpt, dft_file or strained_calculation')


def run(request, workdir):
    kind = request['kind']
    if kind == 'sscha_log':
        r = DATA['runs'][request['wire']][run_key(request)]
        return {'wire': request['wire'], 'temperature_K': r['temperature_K'], 'configurations': r['configurations'],
                **r['log'], 'source': SOURCE}
    if kind == 'sscha_matrices':
        r = DATA['runs'][request['wire']][run_key(request)]
        return {'wire': request['wire'], 'temperature_K': r['temperature_K'], 'configurations': r['configurations'],
                'matrix': request['matrix'], 'files': r[request['matrix']], 'source': SOURCE}
    if kind == 'dfpt':
        return {'wire': request['wire'], 'asr': request['asr'], 'files': DATA['dfpt'][request['wire']][request['asr']],
                'source': SOURCE}
    if kind == 'dft_file':
        return {'wire': request['wire'], 'file': request['file'], 'text': DATA['dft_files'][request['wire']][request['file']],
                'source': SOURCE}
    s = DATA['strained']
    return {'wire': 'TaSe3', 'pw_input': s['pw_input'], 'pw_output': s['pw_output'], 'ph_input': s['ph_input'],
            'dyn_mat_files': s['dyn_mat'], 'aiida_structure_label': s['structure_label'], 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    t, th = PRIVATE['targets'], PRIVATE['thresholds']
    ta, cu = payload.get('TaSe3'), payload.get('CuC2')
    if not isinstance(ta, dict) or not isinstance(cu, dict) or 'evidence_ids' not in payload:
        return False, {'error': 'submission needs TaSe3, CuC2 and evidence_ids'}
    fields = {'TaSe3.axial_stiffness_eV': ta.get('axial_stiffness_eV'),
              'TaSe3.alpha_L_300K_1e-6_per_K': ta.get('alpha_L_300K_1e-6_per_K'),
              'TaSe3.cv_ratio_300K': ta.get('cv_ratio_300K'),
              'TaSe3.instability_temperature_K': ta.get('instability_temperature_K'),
              'CuC2.cv_ratio_300K': cu.get('cv_ratio_300K')}
    bad = [k for k, v in fields.items() if not finite(v)]
    if bad:
        return False, {'error': 'these fields must be finite numbers', 'fields': bad}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    rel = lambda v, ref: v / ref - 1
    k_err = rel(fields['TaSe3.axial_stiffness_eV'], t['k_eV'])
    a_err = rel(fields['TaSe3.alpha_L_300K_1e-6_per_K'], t['alpha_1e6_per_K'])
    cv_cu = rel(fields['CuC2.cv_ratio_300K'], t['cv_cuc2'])
    cv_ta = rel(fields['TaSe3.cv_ratio_300K'], t['cv_tase3'])
    dT = fields['TaSe3.instability_temperature_K'] - t['T_inst_K']
    checks = {'evidence_executed': evidence,
              'TaSe3_axial_stiffness': abs(k_err) <= th['k_rel'],
              'TaSe3_thermal_expansion_300K': abs(a_err) <= th['alpha_rel'],
              'CuC2_heat_capacity_300K': abs(cv_cu) <= th['cv_cuc2_rel'],
              'TaSe3_heat_capacity_300K': abs(cv_ta) <= th['cv_tase3_rel'],
              'TaSe3_instability_temperature': abs(dT) <= th['T_inst_abs_K']}
    report = {'checks': checks, 'relative_errors': {'axial_stiffness': k_err, 'alpha_L_300K': a_err,
                                                    'CuC2_cv_ratio': cv_cu, 'TaSe3_cv_ratio': cv_ta},
              'instability_temperature_error_K': dT}
    return all(checks.values()), report
