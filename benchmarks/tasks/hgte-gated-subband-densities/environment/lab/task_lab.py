"""Lab logic for hgte-gated-subband-densities: retrieval of author-released data and grading.

Purchases return slices of the raw mK magnetotransport records (Zenodo 18429111) and archived kdotpy
self-consistent runs (Zenodo 18709991, FigureS3/calculations). No new calculation is run. The authors'
processed SdH charts and their analysis notebooks are never served; the graded 70 nm subband densities
come from the archived 70 nm chart and the archived runs (private/hidden.json).
"""
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
CALCS = json.loads((HERE/'data/calcs.json').read_text())
POOL = np.load(HERE/'data/lab_pool.npz')
SAMPLES = ('26nm', '70nm', '107nm')
BANDS = ['IDOS(-3)', 'IDOS(-2)', 'IDOS(-1)', 'IDOS(1)', 'IDOS(2)', 'IDOS(3)', 'IDOS(4)', 'IDOS(5)', 'IDOS(6)',
         'IDOS(total)']
WEDGES = {(105, 165): 'kphi_105_165', (165, 285): 'kphi_165_285'}


def g7(x):
    return float('%.7g' % x)


def _sample(request):
    s = request.get('sample')
    if s not in SAMPLES:
        raise ValueError('sample must be one of %s' % list(SAMPLES))
    return s


def _calc_key(s, run, n):
    if not isinstance(n, (int, float)) or isinstance(n, bool):
        raise ValueError('n_nm2 must be a number from the catalog')
    for key, info in CALCS.items():
        if info['sample'] == s and info['run'] == run and abs(info['n_nm2'] - n) < 1e-7:
            return key
    raise ValueError('no archived %s run for %s at n_nm2 = %s (see catalog)' % (run, s, n))


def _centers(s):
    vg = POOL['tr_%s_vg' % s]
    return [round(c, 2) for c in np.arange(-0.8, 0.651, 0.05) if vg.min() - 1e-9 <= c <= vg.max() + 1e-9]


def validate(request):
    kind = request.get('kind')
    if kind == 'hall_scan':
        if set(request) != {'kind', 'sample'}:
            raise ValueError('Use {"kind": "hall_scan", "sample": "<sample>"}')
        _sample(request)
        return 1
    if kind == 'transport':
        if set(request) != {'kind', 'sample', 'vg_center'}:
            raise ValueError('Use {"kind": "transport", "sample": "<sample>", "vg_center": <catalog value>}')
        s = _sample(request)
        c = request['vg_center']
        if not isinstance(c, (int, float)) or isinstance(c, bool) or not any(abs(c - x) < 1e-6 for x in _centers(s)):
            raise ValueError('vg_center must be one of %s for %s' % (_centers(s), s))
        return 1
    if kind == 'calculation':
        if set(request) != {'kind', 'sample', 'n_nm2'}:
            raise ValueError('Use {"kind": "calculation", "sample": "<sample>", "n_nm2": <catalog value>}')
        _calc_key(_sample(request), 'main', request['n_nm2'])
        return 1
    if kind == 'pocket_calculation':
        if set(request) != {'kind', 'sample', 'n_nm2', 'kphi_deg'}:
            raise ValueError('Use {"kind": "pocket_calculation", "sample": "<sample>", "n_nm2": <value>, '
                             '"kphi_deg": [105, 165] or [165, 285]}')
        w = request['kphi_deg']
        if not isinstance(w, list) or tuple(w) not in WEDGES:
            raise ValueError('kphi_deg must be [105, 165] or [165, 285]')
        _calc_key(_sample(request), WEDGES[tuple(w)], request['n_nm2'])
        return 1
    raise ValueError('kind must be "hall_scan", "transport", "calculation" or "pocket_calculation"')


def _calc_result(key):
    info = CALCS[key]
    E, I, P = POOL['E|' + key], POOL['I|' + key], POOL['P|' + key]
    out = {'sample': info['sample'], 'run': info['run'], 'n_nm2': info['n_nm2'],
           'archive_folder': info['archive_folder'], 'kpout': info['kpout'],
           'potential': {'z_nm': [g7(x) for x in P[:, 0]], 'V_meV': [g7(x) for x in P[:, 1]]},
           'integrated_dos_per_nm2': {'E_meV': [g7(x) for x in E]}}
    for j, name in enumerate(BANDS):
        out['integrated_dos_per_nm2'][name] = [g7(x) for x in I[:, j]]
    out['source'] = 'archived kdotpy result (no new calculation)'
    return out


def run(request, workdir):
    kind, s = request['kind'], request['sample']
    if kind == 'hall_scan':
        h = POOL['hall_%s' % s]
        return {'sample': s, 'columns': ['VG_V', 'slope_Vxy_over_VI_per_T', 'Vxx_over_VI_at_B0'],
                'rows': [[g7(x) for x in r] for r in h], 'field_window_T': '|B| <= 0.2',
                'source': 'derived from the archived raw record (no new measurement)'}
    if kind == 'transport':
        vg, B, D = POOL['tr_%s_vg' % s], POOL['tr_%s_B' % s], POOL['tr_%s_data' % s]
        m = np.abs(vg - request['vg_center']) <= 0.025 + 1e-9
        return {'sample': s, 'VG_V': [g7(x) for x in vg[m]], 'B_T': [g7(x) for x in B],
                'VI_V': [[g7(x) for x in row] for row in D[m, :, 0]],
                'Vxx_V': [[g7(x) for x in row] for row in D[m, :, 1]],
                'Vxy_V': [[g7(x) for x in row] for row in D[m, :, 2]],
                'layout': 'arrays indexed [gate][field]', 'source': 'archived raw lock-in record'}
    if kind == 'calculation':
        return _calc_result(_calc_key(s, 'main', request['n_nm2']))
    return _calc_result(_calc_key(s, WEDGES[tuple(request['kphi_deg'])], request['n_nm2']))


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    fields = PRIVATE['fields']
    targets, tol = PRIVATE['targets'], PRIVATE['tolerances']
    required = ['experiment_1e11_cm2', 'model_1e11_cm2', 'hole_onset_gate_V', 'majority_pocket_multiplicity',
                'evidence_ids']
    missing = [f for f in required if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    values = {}
    for key, (section, name) in fields.items():
        container = payload if section is None else payload.get(section)
        if not isinstance(container, dict) or name not in container:
            return False, {'error': 'missing %s%s' % ('' if section is None else section + '.', name)}
        values[key] = container[name]
    mult = values.pop('majority_pocket_multiplicity')
    if not all(finite(v) for v in values.values()):
        return False, {'error': 'all graded values must be finite numbers'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    for key, v in values.items():
        t = targets[key]
        allowed = max(tol[key]['rel'] * abs(t), tol[key]['abs'])
        checks[key] = abs(v - t) <= allowed
        errors[key] = v - t
    checks['majority_pocket_multiplicity'] = (isinstance(mult, (int, float)) and not isinstance(mult, bool)
                                              and mult == PRIVATE['majority_pocket_multiplicity'])
    report = {'checks': checks, 'errors': errors}
    return all(checks.values()), report
