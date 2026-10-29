"""Lab logic for tmd-moire-chern-twist: retrieval of author-archived results for twisted WSe2 / MoTe2, and grading.

Data (data/, written by benchmarks/construction/tmd-moire-chern-twist/make_task.py from Zenodo record 10888769):
SIESTA moire valence bands along gamma-m-kappa-gamma (fig2), C3z eigenvalues of the top moire doublets (transcribed
from the labels of the paper's Fig. 2), and the tWSe2 surface Hartree potentials (fig3: 2D maps, line cuts and the
archived Fourier magnitudes) and relaxed structures. No new calculation is run. Graded values: private/hidden.json.
"""
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data/archive.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
SOURCE = 'retrieved author-archived result (no new calculation)'
ANGLES = {'tWSe2': ('3.15', '1.70', '1.54', '1.25'), 'tMoTe2': ('3.89', '3.15', '2.14', '1.25')}
HARTREE_ANGLES = ('3.15', '1.70', '1.54', '1.25')
COSTS = {'dft_bands': 1, 'c3': 1, 'hartree_map': 3, 'hartree_linecut': 1, 'hartree_fourier': 1, 'structure': 1}


def check_angle(request):
    mat = request.get('material')
    if mat not in ANGLES:
        raise ValueError('material must be "tWSe2" or "tMoTe2"')
    if request.get('angle') not in ANGLES[mat]:
        raise ValueError(f'angle must be one of {list(ANGLES[mat])} (a string) for {mat}')


def validate(request):
    kind = request.get('kind')
    if kind in ('dft_bands', 'c3'):
        keys = {'kind', 'material', 'angle'} | ({'point'} if kind == 'c3' else set())
        if set(request) != keys:
            raise ValueError('Use {"kind": "dft_bands", "material": ..., "angle": ...}' if kind == 'dft_bands' else
                             'Use {"kind": "c3", "material": ..., "angle": ..., "point": "gamma", "kappa" or "kappa_prime"}')
        check_angle(request)
        if kind == 'c3' and request['point'] not in ('gamma', 'kappa', 'kappa_prime'):
            raise ValueError('point must be "gamma", "kappa" or "kappa_prime" (C3z-invariant momenta of the moire BZ)')
        return COSTS[kind]
    if kind in ('hartree_map', 'hartree_linecut', 'hartree_fourier', 'structure'):
        if set(request) != {'kind', 'angle'} or request['angle'] not in HARTREE_ANGLES:
            raise ValueError(f'Use {{"kind": "{kind}", "angle": one of {list(HARTREE_ANGLES)}}} (tWSe2 only)')
        return COSTS[kind]
    raise ValueError('kind must be one of ' + ', '.join(COSTS))


def run(request, workdir):
    kind = request['kind']
    if kind == 'dft_bands':
        return {**DATA['bands'][request['material']][request['angle']], 'material': request['material'],
                'angle_deg': request['angle'], 'source': SOURCE}
    if kind == 'c3':
        rec = DATA['c3'][request['material']][request['angle']][request['point']]
        return {'material': request['material'], 'angle_deg': request['angle'], 'point': request['point'], **rec,
                'source': SOURCE + '; C3z eigenvalues as labelled by the authors in Fig. 2 of the paper'}
    angle = request['angle']
    if kind == 'hartree_map':
        z = np.load(HERE/'data'/f'hartree_{angle}.npz')
        n = int(z['n'])
        return {'material': 'tWSe2', 'angle_deg': angle, 'n': n, **DATA['hartree_map_note'],
                'v_top_surface_eV': z['top'].round(7).tolist(), 'v_bottom_surface_eV': z['bot'].round(7).tolist(),
                'source': SOURCE}
    if kind == 'hartree_linecut':
        return {'material': 'tWSe2', 'angle_deg': angle, **DATA['linecut'][angle], 'source': SOURCE}
    if kind == 'hartree_fourier':
        return {'material': 'tWSe2', 'angle_deg': angle, **DATA['fourier'][angle], 'source': SOURCE}
    if kind == 'structure':
        return {'material': 'tWSe2', 'angle_deg': angle, 'poscar': (HERE/'data'/f'POSCAR_{angle}').read_text(),
                'note': 'relaxed structure (VASP POSCAR format) on which the Hartree potentials of this angle were computed',
                'source': SOURCE}
    raise ValueError('unknown kind')


def integer_pair(v):
    return isinstance(v, list) and len(v) == 2 and all(type(x) is int for x in v)


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    T, th = PRIVATE['targets'], PRIVATE['thresholds']
    fields = ['dft_chern_K_up', 'continuum_chern_3p89', 'continuum_set_matching_dft', 'continuum_P_chern_2p14',
              'hartree_xm_meV', 'hartree_first_star_xm_meV', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    dc, cc = payload['dft_chern_K_up'], payload['continuum_chern_3p89']
    if not isinstance(dc, dict) or set(dc) != set(T['dft_chern_K_up']) or not all(integer_pair(v) for v in dc.values()):
        return False, {'error': 'dft_chern_K_up must give an integer pair [C1, C2] for every key ' +
                                ', '.join(sorted(T['dft_chern_K_up']))}
    if not isinstance(cc, dict) or set(cc) != {'P', 'Q'} or not all(integer_pair(v) for v in cc.values()):
        return False, {'error': 'continuum_chern_3p89 must give integer pairs for "P" and "Q"'}
    if not integer_pair(payload['continuum_P_chern_2p14']):
        return False, {'error': 'continuum_P_chern_2p14 must be an integer pair [C1, C2]'}
    hx, hf = payload['hartree_xm_meV'], payload['hartree_first_star_xm_meV']
    if not isinstance(hx, dict) or set(hx) != set(T['hartree_xm_meV']) or not all(number(v) for v in hx.values()):
        return False, {'error': 'hartree_xm_meV must give a number for ' + ', '.join(sorted(T['hartree_xm_meV']))}
    if not isinstance(hf, dict) or set(hf) != set(T['hartree_first_star_xm_meV']) or not all(number(v) for v in hf.values()):
        return False, {'error': 'hartree_first_star_xm_meV must give a number for ' +
                                ', '.join(sorted(T['hartree_first_star_xm_meV']))}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    for key, ref in T['dft_chern_K_up'].items():
        checks[f'dft_chern:{key}'] = dc[key] == ref
    for key, ref in T['continuum_chern_3p89'].items():
        checks[f'continuum_chern_3p89:{key}'] = cc[key] == ref
    checks['continuum_set_matching_dft'] = payload['continuum_set_matching_dft'] == T['continuum_set_matching_dft']
    checks['continuum_P_chern_2p14'] = payload['continuum_P_chern_2p14'] == T['continuum_P_chern_2p14']
    errors = {}
    for key, ref in T['hartree_xm_meV'].items():
        errors[f'hartree_xm:{key}'] = hx[key] - ref
        checks[f'hartree_xm:{key}'] = abs(hx[key] - ref) <= th['hartree_xm_abs_meV']
    for key, ref in T['hartree_first_star_xm_meV'].items():
        errors[f'first_star_xm:{key}'] = hf[key] - ref
        tol = max(th['first_star_abs_meV'], th['first_star_rel'] * abs(ref))
        checks[f'first_star_xm:{key}'] = abs(hf[key] - ref) <= tol
    return all(checks.values()), {'checks': checks, 'errors_meV': errors}
