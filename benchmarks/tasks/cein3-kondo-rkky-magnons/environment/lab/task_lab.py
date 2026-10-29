"""Lab logic for cein3-kondo-rkky-magnons: retrieval of author-archived CeIn3 neutron data and grading.

Purchases return archived chi''(Q,E) data of Simeth et al. (Zenodo 10146787): constant-Q energy spectra of the
12-meV maps (paths RG, GX, XM, MG), constant-q energy spectra of the high-resolution maps along RG and RX, the
high-resolution constant-energy cuts along RG and RX, and the band structure. The high-resolution R-M map and
R-M cuts and the Q-integrated chi'' are never served; they grade the submission.
"""
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
COSTS = DATA['costs']
SOURCE = 'retrieved author-archived neutron data (no new measurement or calculation)'
RANGES = ('nearest_neighbour', 'short_range', 'long_range', 'no_spin_waves')


def index_of(request, table):
    i = request.get('index')
    if type(i) is not int or str(i) not in table:
        raise ValueError('index must be one of the integers listed in /app/catalog.json')
    return str(i)


def resolve(request):
    kind = request.get('kind')
    if kind == 'he_spectrum':
        if set(request) != {'kind', 'path', 'index'} or request.get('path') not in DATA['he']:
            raise ValueError('Use {"kind": "he_spectrum", "path": "RG"|"GX"|"XM"|"MG", "index": <integer from the catalog>}')
        return kind, (request['path'], index_of(request, DATA['he'][request['path']]['rows']))
    if kind == 'hr_spectrum':
        if set(request) != {'kind', 'path', 'index'} or request.get('path') not in DATA['hr']:
            raise ValueError('Use {"kind": "hr_spectrum", "path": "RG"|"RX", "index": <integer from the catalog>}; the R-M map is not served')
        return kind, (request['path'], index_of(request, DATA['hr'][request['path']]['rows']))
    if kind == 'hr_cut':
        if set(request) != {'kind', 'path', 'energy_meV'} or request.get('path') not in DATA['cuts']:
            raise ValueError('Use {"kind": "hr_cut", "path": "RG"|"RX", "energy_meV": 0.6|1.4}; R-M cuts are not served')
        e = request['energy_meV']
        if isinstance(e, bool) or not isinstance(e, (int, float)) or '%.1f' % e not in DATA['cuts'][request['path']] or abs(e - round(e, 1)) > 1e-9:
            raise ValueError('energy_meV must be 0.6 or 1.4')
        return kind, (request['path'], '%.1f' % e)
    if kind == 'band_structure':
        if set(request) != {'kind'}:
            raise ValueError('Use {"kind": "band_structure"}')
        return kind, ()
    raise ValueError('kind must be one of: he_spectrum, hr_spectrum, hr_cut, band_structure')


def validate(request):
    kind, _ = resolve(request)
    return COSTS[kind]


def run(request, workdir):
    kind, key = resolve(request)
    if kind == 'he_spectrum':
        path, i = key
        block = DATA['he'][path]
        return {'path': path, **block['rows'][i], 'energies_meV': block['energies_meV'], 'energy_bin_meV': block['energy_bin_meV'],
                'path_halfwidth_rlu': 0.0125, 'perpendicular_halfwidth_rlu': block['perpendicular_halfwidth_rlu'],
                'units': DATA['units'], 'source': SOURCE}
    if kind == 'hr_spectrum':
        path, i = key
        block = DATA['hr'][path]
        return {'path': path, **block['rows'][i], 'energies_meV': block['energies_meV'], 'energy_bin_meV': block['energy_bin_meV'],
                'q_bin_rlu': 0.02, 'units': DATA['units'], 'source': SOURCE}
    if kind == 'hr_cut':
        path, e = key
        return {'path': path, **DATA['cuts'][path][e], 'units': DATA['units'], 'source': SOURCE}
    return {**DATA['bands'], 'source': SOURCE}


def finite_list(x, n):
    return (isinstance(x, list) and len(x) == n and
            all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) for v in x))


def chi2r(m, y, e, offset):
    """Reduced chi^2 after the best non-negative scale (and constant offset when offset=True)."""
    m, y, e = (np.asarray(v, float) for v in (m, y, e))
    w = 1 / e ** 2
    if offset:
        A = np.c_[m, np.ones_like(m)] * np.sqrt(w)[:, None]
        c = np.linalg.lstsq(A, y * np.sqrt(w), rcond=None)[0]
        if c[0] < 0:
            c = np.array([0.0, np.sum(w * y) / np.sum(w)])
        f = c[0] * m + c[1]
        dof = len(y) - 2
    else:
        a = max(float(np.sum(w * m * y) / max(np.sum(w * m * m), 1e-300)), 0.0)
        f = a * m
        dof = len(y) - 1
    return float(np.sum(((f - y) / e) ** 2) / dof)


def grade(payload, trace, validation):
    P, th = PRIVATE, PRIVATE['thresholds']
    fields = ['rm_map_prediction', 'rm_cut_1p4meV_prediction', 'integrated_chi_prediction', 'bandwidth_meV',
              'velocity_bandwidth_ratio', 'J1_nearest_neighbour_meV', 'exchange_range', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    rm = P['rm_map']
    grid = payload['rm_map_prediction']
    if not (isinstance(grid, list) and len(grid) == len(rm['values']) and all(finite_list(r, len(rm['energies_meV'])) for r in grid)):
        return False, {'error': 'rm_map_prediction must be a %d x %d list of finite numbers (rows = rm_map q_rlu, columns = energies_meV)'
                       % (len(rm['values']), len(rm['energies_meV']))}
    if not finite_list(payload['rm_cut_1p4meV_prediction'], len(P['rm_cut_1p4meV']['y'])):
        return False, {'error': 'rm_cut_1p4meV_prediction must list %d finite numbers' % len(P['rm_cut_1p4meV']['y'])}
    if not finite_list(payload['integrated_chi_prediction'], len(P['integrated']['y'])):
        return False, {'error': 'integrated_chi_prediction must list %d finite numbers' % len(P['integrated']['y'])}
    scalars = ['bandwidth_meV', 'velocity_bandwidth_ratio', 'J1_nearest_neighbour_meV']
    if not finite_list([payload[k] for k in scalars], 3):
        return False, {'error': ', '.join(scalars) + ' must be finite numbers'}
    if payload['exchange_range'] not in RANGES:
        return False, {'error': 'exchange_range must be one of ' + ', '.join(RANGES)}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    m, y = [], []
    for i, row in enumerate(rm['values']):
        for j, v in enumerate(row):
            if v is not None and rm['energies_meV'][j] > rm['emin_meV']:
                m.append(grid[i][j]); y.append(v)
    scores = {'rm_map_chi2r': chi2r(m, y, np.full(len(y), rm['sigma']), True),
              'rm_cut_1p4meV_chi2r': chi2r(payload['rm_cut_1p4meV_prediction'], P['rm_cut_1p4meV']['y'], P['rm_cut_1p4meV']['e'], True),
              'integrated_chi2r': chi2r(payload['integrated_chi_prediction'], P['integrated']['y'], P['integrated']['e'], False)}
    lo_w, hi_w = th['bandwidth_meV']
    lo_e, hi_e = th['velocity_bandwidth_ratio']
    lo_j, hi_j = th['J1_nearest_neighbour_meV']
    checks = {'evidence_executed': evidence,
              'rm_map': scores['rm_map_chi2r'] <= th['rm_map_chi2r'],
              'rm_cut_1p4meV': scores['rm_cut_1p4meV_chi2r'] <= th['rm_cut_1p4meV_chi2r'],
              'integrated_chi': scores['integrated_chi2r'] <= th['integrated_chi2r'],
              'bandwidth': lo_w <= payload['bandwidth_meV'] <= hi_w,
              'velocity_bandwidth_ratio': lo_e <= payload['velocity_bandwidth_ratio'] <= hi_e,
              'J1_nearest_neighbour': lo_j <= payload['J1_nearest_neighbour_meV'] <= hi_j,
              'exchange_range': payload['exchange_range'] == P['exchange_range']}
    return all(checks.values()), {'checks': checks, 'scores': scores}
