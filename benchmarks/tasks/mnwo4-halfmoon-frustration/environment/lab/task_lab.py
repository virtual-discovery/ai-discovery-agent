"""Lab logic for mnwo4-halfmoon-frustration: retrieval of author-archived UppASD outputs and grading.

Purchases return archived spin-dynamics outputs of the 12-exchange + single-ion-anisotropy model of Natsch et al.
(Materials Cloud p8297-g9n39): 1.5 K S(q,w) spectra along the archived path, static S(q) maps at 20/30/60 K in the
half-moon plane, two reduced-model maps at 20 K and 11 K full-zone slices. The 14 K map, the other reduced-model
maps and the exchange constants themselves are never served; they grade the submission.
"""
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).parent
DATA = np.load(HERE/'data.npz')
META = json.loads((HERE/'data_meta.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
COSTS = META['costs']
SOURCE = 'retrieved author-archived UppASD output (no new calculation)'
ORIGINS = ('spin_liquid_local_constraint', 'quasi_1d_zigzag_chain', 'competing_3d_exchange', 'single_ion_anisotropy')
JLABELS = ['J%d' % i for i in range(1, 13)]


def sig4(a):
    """archived values rounded to 4 significant digits, as nested lists"""
    a = np.asarray(a, float)
    return np.vectorize(lambda x: float('%.4g' % x))(a).tolist()


def resolve(request):
    kind = request.get('kind')
    if kind == 'spectrum':
        i = request.get('index')
        if set(request) != {'kind', 'index'} or type(i) is not int or not 0 <= i < 101:
            raise ValueError('Use {"kind": "spectrum", "index": <integer 0-100>}')
        return kind, i
    if kind == 'map':
        T = request.get('temperature_K')
        if set(request) != {'kind', 'temperature_K'} or isinstance(T, bool) or T not in (20, 30, 60):
            raise ValueError('Use {"kind": "map", "temperature_K": 20 | 30 | 60}; other temperatures are not sold')
        return kind, int(T)
    if kind == 'reduced_model_map':
        m = request.get('model')
        if set(request) != {'kind', 'model'} or m not in ('5b', '5c'):
            raise ValueError('Use {"kind": "reduced_model_map", "model": "5b" | "5c"}; the other reduced models are not sold')
        return kind, m
    if kind == 'zone_slice':
        box, j = request.get('supercell'), request.get('k_index')
        if set(request) != {'kind', 'supercell', 'k_index'} or box not in ('36x24x36', '48x24x24') \
                or type(j) is not int or not 0 <= j < 24:
            raise ValueError('Use {"kind": "zone_slice", "supercell": "36x24x36" | "48x24x24", "k_index": <integer 0-23>}')
        return kind, (box, j)
    raise ValueError('kind must be one of: spectrum, map, reduced_model_map, zone_slice')


def validate(request):
    kind, _ = resolve(request)
    return COSTS[kind]


def run(request, workdir):
    kind, key = resolve(request)
    if kind == 'spectrum':
        S = DATA['spectra_S'][key]
        return {'index': key, 'q_rlu': DATA['spectra_q'][key].tolist(), 'temperature_K': 1.5,
                'energies_meV': 'E_n = n * %.10f for n = 0 .. %d' % (META['spectra_dE'], S.shape[0] - 1),
                'Sxx': sig4(S[:, 0]), 'Syy': sig4(S[:, 1]), 'Szz': sig4(S[:, 2]), 'source': SOURCE}
    if kind == 'map':
        S = DATA['map_T%d' % key]
        return {'temperature_K': key, 'grid': META['map_grid'], 'rows': 'one row per k_i (i = 0 .. 186), one value per l_j (j = 0 .. 320)',
                'Sxx': sig4(S[:, :, 0]), 'Syy': sig4(S[:, :, 1]), 'Szz': sig4(S[:, :, 2]), 'source': SOURCE}
    if kind == 'reduced_model_map':
        S = DATA['ko_' + key]
        return {'model': key, 'temperature_K': 20, 'grid': META['ko_grid'],
                'rows': 'one row per k_i (i = 0 .. 160), one value per l_j (j = 0 .. 160)',
                'Sxx': sig4(S[:, :, 0]), 'Syy': sig4(S[:, :, 1]), 'Szz': sig4(S[:, :, 2]), 'source': SOURCE}
    box, j = key
    S = DATA['zone_' + box.replace('x', '_')][:, j]
    g = META['zone_grids'][box]
    return {'supercell': box, 'temperature_K': 11, 'k_index': j, 'k_rlu': (j - 11) / 24, 'grid': g,
            'rows': 'one row per h_i, one value per l_m', 'Sxx': sig4(S[:, :, 0]), 'Syy': sig4(S[:, :, 1]),
            'Szz': sig4(S[:, :, 2]), 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grid_ok(x, shape):
    return (isinstance(x, list) and len(x) == shape[0] and
            all(isinstance(r, list) and len(r) == shape[1] and all(finite(v) for v in r) for r in x))


def pearson(a, b):
    a = np.asarray(a, float).ravel(); b = np.asarray(b, float).ravel()
    a = a - a.mean(); b = b - b.mean()
    den = math.sqrt(float(a @ a) * float(b @ b))
    return float(a @ b) / den if den > 0 else 0.0


def k_distance(k1, k2):
    """largest component difference between propagation vectors, modulo sign and reciprocal-lattice vectors"""
    best = 9.0
    for s in (1, -1):
        d = np.asarray(k1, float) - s * np.asarray(k2, float)
        d = d - np.round(d)
        best = min(best, float(np.abs(d).max()))
    return best


def grade(payload, trace, validation):
    P, th = PRIVATE, PRIVATE['thresholds']
    fields = ['exchange_meV', 'anisotropy_D_meV', 'propagation_vector_rlu', 'spiral_minus_af1_meV_per_Mn',
              'map_14K_prediction', 'reduced_model_predictions', 'half_moon_origin', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    J = payload['exchange_meV']
    if not isinstance(J, dict) or set(J) != set(JLABELS) or not all(finite(v) for v in J.values()):
        return False, {'error': 'exchange_meV must give a finite value for each of J1 .. J12'}
    k = payload['propagation_vector_rlu']
    if not (isinstance(k, list) and len(k) == 3 and all(finite(v) for v in k)):
        return False, {'error': 'propagation_vector_rlu must be a list of three finite numbers'}
    if not finite(payload['anisotropy_D_meV']) or not finite(payload['spiral_minus_af1_meV_per_Mn']):
        return False, {'error': 'anisotropy_D_meV and spiral_minus_af1_meV_per_Mn must be finite numbers'}
    m14 = payload['map_14K_prediction']
    if not grid_ok(m14, P['map14_shape']):
        return False, {'error': 'map_14K_prediction must be a %d x %d list of finite numbers' % tuple(P['map14_shape'])}
    red = payload['reduced_model_predictions']
    if not isinstance(red, dict) or set(red) != set(P['reduced_refs']) or \
            not all(grid_ok(red[m], P['reduced_shape']) for m in P['reduced_refs']):
        return False, {'error': 'reduced_model_predictions must map each of %s to a %d x %d list of finite numbers'
                       % (sorted(P['reduced_refs']), *P['reduced_shape'])}
    if payload['half_moon_origin'] not in ORIGINS:
        return False, {'error': 'half_moon_origin must be one of ' + ', '.join(ORIGINS)}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    scores = {'map_14K_r': pearson(m14, P['map14_ref'])}
    for m, ref in P['reduced_refs'].items():
        scores['reduced_%s_r' % m] = pearson(red[m], ref)
    scores['k_distance'] = k_distance(k, P['propagation_vector_rlu'])
    # exchange_meV, anisotropy_D_meV and spiral_minus_af1_meV_per_Mn are reported and reviewed with the report; they
    # are not graded numerically (the archived data do not determine them uniquely without literature priors)
    checks = {'evidence_executed': evidence,
              'propagation_vector': scores['k_distance'] <= th['k_tolerance']}
    checks['map_14K'] = scores['map_14K_r'] >= th['map14_min_r']
    for m in P['reduced_refs']:
        checks['reduced_model_' + m] = scores['reduced_%s_r' % m] >= th['reduced_min_r'][m]
    checks['half_moon_origin'] = payload['half_moon_origin'] == P['half_moon_origin']
    return all(checks.values()), {'checks': checks, 'scores': scores}
