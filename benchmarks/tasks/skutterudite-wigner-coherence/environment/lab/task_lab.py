"""Lab logic for skutterudite-wigner-coherence: retrieval of archived phono3py data and grading.

Purchases return the authors' archived per-mode arrays (mode_data), velocity-operator products
(velocity_operator), archived output files of the released runs (archived_file) and phonopy harmonic
summaries (harmonic). Nothing is recomputed. The archived outputs (logs, kappa tables) of the IrSb3 and
YbIr4Sb12 runs are never served; reference values derived from their archived inputs grade the submission.
"""
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).parent
META = json.loads((HERE/'data/meta.json').read_text())
PUBLIC = json.loads((HERE/'public.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
COSTS = {k: v['cost'] for k, v in PUBLIC['requests'].items()}
RUNS = META['runs']
_CACHE = {}


def arrays(run):
    if run not in _CACHE:
        with np.load(HERE/'data'/(run.replace('/', '__') + '.npz')) as z:
            _CACHE[run] = {k: z[k] for k in z.files}
    return _CACHE[run]


def rounded(a, sig):
    a = np.asarray(a, float)
    out = np.zeros_like(a)
    nz = a != 0
    mag = np.floor(np.log10(np.abs(a[nz])))
    scale = 10.0 ** (sig - 1 - mag)
    out[nz] = np.round(a[nz] * scale) / scale
    return out.tolist()


def validate(request):
    kind = request.get('kind')
    if kind == 'mode_data':
        if set(request) != {'kind', 'run', 'temperature_K'} or request['run'] not in RUNS:
            raise ValueError('Use {"kind": "mode_data", "run": <run from catalog>, "temperature_K": <archived temperature>}')
        temps = arrays(request['run'])['temperature'].astype(int).tolist()
        if not isinstance(request['temperature_K'], (int, float)) or isinstance(request['temperature_K'], bool) \
                or request['temperature_K'] not in temps:
            raise ValueError('temperature_K must be one of %s' % temps)
    elif kind == 'velocity_operator':
        if set(request) != {'kind', 'run'} or request['run'] not in RUNS:
            raise ValueError('Use {"kind": "velocity_operator", "run": <run from catalog>}')
        if not RUNS[request['run']]['velocity_operator']:
            raise ValueError('No velocity operator is archived for this run (phono3py 1.19 file)')
    elif kind == 'archived_file':
        if set(request) != {'kind', 'run', 'file'} or request['run'] not in RUNS:
            raise ValueError('Use {"kind": "archived_file", "run": <run from catalog>, "file": <one of archived_files>}')
        if RUNS[request['run']]['withheld_outputs']:
            raise ValueError('The archived outputs of this run are withheld')
        if request['file'] not in RUNS[request['run']]['files']:
            raise ValueError('file must be one of %s' % sorted(RUNS[request['run']]['files']))
    elif kind == 'harmonic':
        if set(request) != {'kind', 'compound'} or request['compound'] not in META['harmonic']:
            raise ValueError('Use {"kind": "harmonic", "compound": <one of %s>}' % sorted(META['harmonic']))
    else:
        raise ValueError('kind must be "mode_data", "velocity_operator", "archived_file" or "harmonic"')
    return COSTS[kind]


def run(request, workdir):
    kind = request['kind']
    src = 'retrieved author-archived data (no new calculation)'
    if kind == 'mode_data':
        d = arrays(request['run'])
        it = d['temperature'].astype(int).tolist().index(int(request['temperature_K']))
        return {'run': request['run'], 'temperature_K': int(request['temperature_K']), 'q_mesh': d['mesh'].tolist(),
                'qpoints': rounded(d['qpoint'], 8), 'weights': d['weight'].astype(int).tolist(),
                'frequency_THz': rounded(d['frequency'], 7), 'gamma_THz': rounded(d['gamma'][it], 7),
                'gamma_isotope_THz': rounded(d['gamma_isotope'], 7) if 'gamma_isotope' in d else None,
                'heat_capacity_eV_per_K': rounded(d['heat_capacity'][it], 7), 'v2_A2THz2': rounded(d['v2'], 7),
                'source': src}
    if kind == 'velocity_operator':
        d = arrays(request['run'])
        return {'run': request['run'], 'q_mesh': d['mesh'].tolist(), 'qpoints': rounded(d['qpoint'], 8),
                'weights': d['weight'].astype(int).tolist(), 'frequency_THz': rounded(d['frequency'], 7),
                'G_A2THz2': rounded(d['G'], 6), 'source': src}
    if kind == 'archived_file':
        entry = RUNS[request['run']]['files'][request['file']]
        out = {'run': request['run'], 'file': request['file'], 'kind': entry['kind'], 'source': src}
        if entry['kind'] == 'phono3py log':
            out.update(header=entry['header'], thermal_conductivity_xx_W_mK=entry['tables'])
        else:
            out['text'] = entry['text']
        return out
    h = META['harmonic'][request['compound']]
    return {'compound': request['compound'], **h, 'source': 'retrieved author-archived phonopy/DFPT results (no new calculation)'}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def get(payload, *path):
    for key in path:
        if not isinstance(payload, dict) or key not in payload:
            return None
        payload = payload[key]
    return payload


def grade(payload, trace, validation):
    ref, t = PRIVATE['reference'], PRIVATE['thresholds']
    paths = [('YbIr4Sb12', 'kappa_P_W_mK', '300'), ('YbIr4Sb12', 'kappa_P_W_mK', '700'),
             ('YbIr4Sb12', 'kappa_C_W_mK', '100'), ('YbIr4Sb12', 'kappa_C_W_mK', '300'), ('YbIr4Sb12', 'kappa_C_W_mK', '700'),
             ('YbIr4Sb12', 'T_half_K'), ('YbIr4Sb12', 'rattler_frequency_THz'),
             ('IrSb3', 'kappa_P_50K_W_mK'), ('IrSb3', 'kappa_TOT_300K_W_mK')]
    values = {'/'.join(p): get(payload, *p) for p in paths}
    bad = [k for k, v in values.items() if not finite(v)]
    if bad or 'evidence_ids' not in payload:
        return False, {'error': 'missing or non-finite fields', 'fields': bad + ([] if 'evidence_ids' in payload else ['evidence_ids'])}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    rel = {k: values[k] / get(ref, *k.split('/')) - 1 for k in values if not k.endswith('T_half_K')}
    dT = values['YbIr4Sb12/T_half_K'] - ref['YbIr4Sb12']['T_half_K']
    checks = {'evidence_executed': evidence,
              'YbIr4Sb12_kappa_P_300K': abs(rel['YbIr4Sb12/kappa_P_W_mK/300']) <= t['kappa_P_rel'],
              'YbIr4Sb12_kappa_P_700K': abs(rel['YbIr4Sb12/kappa_P_W_mK/700']) <= t['kappa_P_rel'],
              'YbIr4Sb12_kappa_C_100K': abs(rel['YbIr4Sb12/kappa_C_W_mK/100']) <= t['kappa_C_rel'],
              'YbIr4Sb12_kappa_C_300K': abs(rel['YbIr4Sb12/kappa_C_W_mK/300']) <= t['kappa_C_rel'],
              'YbIr4Sb12_kappa_C_700K': abs(rel['YbIr4Sb12/kappa_C_W_mK/700']) <= t['kappa_C_rel'],
              'YbIr4Sb12_T_half': abs(dT) <= t['T_half_K'],
              'YbIr4Sb12_rattler_frequency': abs(rel['YbIr4Sb12/rattler_frequency_THz']) <= t['rattler_rel'],
              'IrSb3_kappa_P_50K': abs(rel['IrSb3/kappa_P_50K_W_mK']) <= t['irsb3_kappa_P_50K_rel'],
              'IrSb3_kappa_TOT_300K': abs(rel['IrSb3/kappa_TOT_300K_W_mK']) <= t['irsb3_kappa_TOT_300K_rel']}
    report = {'checks': checks, 'relative_errors': {k: round(v, 4) for k, v in rel.items()}, 'T_half_error_K': round(dT, 1)}
    return all(checks.values()), report
