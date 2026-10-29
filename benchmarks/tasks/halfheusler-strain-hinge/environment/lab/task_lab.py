"""Lab logic for halfheusler-strain-hinge: a numpy/scipy evaluator of the four archived FPLO Wannier models.

The models are the four +hamdata files of Zenodo record 16927281 (wannier/), converted unchanged to data/*.npz
(hopping vectors in units of each file's own lattice vectors, energies in eV on the file's energy scale).
tb.py builds bulk Bloch Hamiltonians, slabs and prisms from them.  Graded values are derived quantities held in
private/hidden.json and are never served.
"""
import json
import math
from pathlib import Path
import numpy as np
import tb

HERE = Path(__file__).parent
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
MODELS = ('LiSbZn_orig', 'LiSbZn_strained_setup', 'LiBiZn_orig', 'LiBiZn_strained_setup')
REGIONS = ('c00', 'c0N', 'cN0', 'cNN', 'f1lo', 'f1hi', 'f2lo', 'f2hi', 'bulk')
_CACHE = {}
NOTE = ('new evaluation of the archived Wannier Hamiltonian; energies in eV on the energy scale of the archived file '
        '(its zero is the zero of the FPLO calculation); bands numbered from 1 (lowest)')


def model(name):
    if name not in _CACHE:
        _CACHE[name] = tb.Model(HERE/'data'/f'{name}.npz')
    return _CACHE[name]


def number(x, lo=-math.inf, hi=math.inf):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x) and lo <= x <= hi


def integer(x, lo, hi):
    return type(x) is int and lo <= x <= hi


def check_model(request):
    if request.get('model') not in MODELS:
        raise ValueError('model must be one of ' + ', '.join(MODELS))


def kpoints(request, limit):
    given = [k for k in ('k_frac', 'k_cart_invA') if k in request]
    if len(given) != 1:
        raise ValueError('give exactly one of k_frac or k_cart_invA')
    pts = request[given[0]]
    if (not isinstance(pts, list) or not 1 <= len(pts) <= limit or
            not all(isinstance(p, list) and len(p) == 3 and all(number(x, -20, 20) for x in p) for p in pts)):
        raise ValueError(f'{given[0]} must be a list of 1 to {limit} points [x, y, z]')
    pts = np.array(pts, float)
    return pts @ model(request['model']).B if given[0] == 'k_frac' else pts


def check_basis(basis):
    if (not isinstance(basis, list) or len(basis) != 3 or
            not all(isinstance(r, list) and len(r) == 3 and all(integer(x, -6, 6) for x in r) for r in basis)):
        raise ValueError('basis must be three rows of three integers (-6..6), in units of the model\'s lattice vectors')
    if abs(round(float(np.linalg.det(np.array(basis, float))))) != 1:
        raise ValueError('basis must be unimodular (determinant +1 or -1)')


def check_window(request, width):
    w = request.get('window_eV')
    if not (isinstance(w, list) and len(w) == 2 and all(number(x, -8, 4) for x in w) and 0 < w[1] - w[0] <= width):
        raise ValueError(f'window_eV must be [emin, emax] with 0 < emax - emin <= {width} eV')


def check_phases(request, limit, ndim):
    ks = request.get('k_phases')
    if not (isinstance(ks, list) and 1 <= len(ks) <= limit):
        raise ValueError(f'k_phases must be a list of 1 to {limit} entries')
    for k in ks:
        ok = number(k, -7, 7) if ndim == 1 else (isinstance(k, list) and len(k) == 2 and all(number(x, -7, 7) for x in k))
        if not ok:
            raise ValueError('each k_phases entry must be ' + ('a number (radians)' if ndim == 1 else '[phi1, phi2] (radians)'))


def prism_cost(n):
    return 1 if n <= 12 else 2 if n <= 20 else 3


def validate(request):
    kind = request.get('kind')
    if kind == 'bulk_bands':
        if not set(request) <= {'kind', 'model', 'k_frac', 'k_cart_invA'}:
            raise ValueError('Use {"kind": "bulk_bands", "model": ..., "k_frac" or "k_cart_invA": [[x, y, z], ...] (<= 300)}')
        check_model(request)
        kpoints(request, 300)
        return 1
    if kind == 'band_edges':
        if set(request) != {'kind', 'model', 'band', 'mesh'}:
            raise ValueError('Use {"kind": "band_edges", "model": ..., "band": n (1-7), "mesh": m (4-48)}')
        check_model(request)
        if not integer(request['band'], 1, 7) or not integer(request['mesh'], 4, 48):
            raise ValueError('band must be an integer 1-7 and mesh an integer 4-48')
        return 1 if request['mesh'] <= 24 else 2
    if kind == 'fermi_level':
        if not {'kind', 'model', 'electrons_per_cell', 'mesh'} <= set(request) <= {'kind', 'model', 'electrons_per_cell', 'mesh', 'offset'}:
            raise ValueError('Use {"kind": "fermi_level", "model": ..., "electrons_per_cell": x, "mesh": m (4-64), "offset": 0 or 0.5}')
        check_model(request)
        if not number(request['electrons_per_cell'], 0.01, 7.99) or not integer(request['mesh'], 4, 64):
            raise ValueError('electrons_per_cell must be in (0, 8) and mesh an integer 4-64')
        if request.get('offset', 0) not in (0, 0.5):
            raise ValueError('offset must be 0 (Gamma-centred) or 0.5')
        return 1 if request['mesh'] <= 32 else 2
    if kind == 'slab':
        keys = {'kind', 'model', 'basis', 'layers', 'k_phases', 'window_eV'}
        if not keys <= set(request) <= keys | {'region_depth'}:
            raise ValueError('Use {"kind": "slab", "model": ..., "basis": [[p1], [p2], [n]], "layers": L (2-40), '
                             '"k_phases": [[phi1, phi2], ...] (<= 30), "window_eV": [emin, emax], "region_depth": 1-4}')
        check_model(request)
        check_basis(request['basis'])
        if not integer(request['layers'], 2, 40):
            raise ValueError('layers must be an integer 2-40')
        check_phases(request, 30, 2)
        check_window(request, 1.2)
        if not integer(request.get('region_depth', 2), 1, 4):
            raise ValueError('region_depth must be an integer 1-4')
        return 1
    if kind == 'prism':
        keys = {'kind', 'model', 'basis', 'size', 'k_phases', 'window_eV'}
        if not keys <= set(request) <= keys | {'region_depth', 'potential'}:
            raise ValueError('Use {"kind": "prism", "model": ..., "basis": [[axis], [side1], [side2]], "size": N (4-28), '
                             '"k_phases": [k, ...], "window_eV": [emin, emax], "region_depth": 1-4, "potential": [...]}')
        check_model(request)
        check_basis(request['basis'])
        n = request['size']
        if not integer(n, 4, 28):
            raise ValueError('size must be an integer 4-28')
        check_phases(request, 8 if n <= 12 else 4 if n <= 20 else 2, 1)
        check_window(request, 1.0)
        if not integer(request.get('region_depth', 2), 1, 4):
            raise ValueError('region_depth must be an integer 1-4')
        pot = request.get('potential', [])
        if not isinstance(pot, list) or len(pot) > 4:
            raise ValueError('potential must be a list of at most 4 entries')
        for p in pot:
            if (not isinstance(p, dict) or set(p) != {'region', 'depth', 'V_eV', 'species'} or p['region'] not in REGIONS[:-1]
                    or not integer(p['depth'], 1, 6) or not number(p['V_eV'], -1.5, 1.5) or p['species'] not in ('all', 'Zn', 'X')):
                raise ValueError('each potential entry is {"region": c00|c0N|cN0|cNN|f1lo|f1hi|f2lo|f2hi, "depth": 1-6, '
                                 '"V_eV": -1.5..1.5, "species": "all"|"Zn"|"X"}')
        return prism_cost(n)
    raise ValueError('kind must be one of bulk_bands, band_edges, fermi_level, slab, prism')


def grid(m, mesh, offset=0.0):
    f = (np.arange(mesh) + offset) / mesh
    return np.array(np.meshgrid(f, f, f, indexing='ij')).reshape(3, -1).T


def r6(x):
    return np.round(np.asarray(x, float), 6).tolist()


def point(m, kf):
    kf = kf - np.floor(kf + 0.5)
    return {'k_frac': r6(kf), 'k_cart_invA': r6(kf @ m.B)}


def geometry(m, basis):
    a_new, _, g = m.regroup(basis)
    zn = (m.centers[0] - m.centers[2]) - g @ a_new
    return {'basis_vectors_A': r6(a_new), 'x_site_reduced': [0.0, 0.0, 0.0],
            'zn_site_reduced': r6(zn @ np.linalg.inv(a_new)), 'zn_site_offset_A': r6(zn)}


def region_masks(cells, n, h):
    t1, t2 = cells[:, -2], cells[:, -1]
    lo1, hi1, lo2, hi2 = t1 < h, t1 >= n - h, t2 < h, t2 >= n - h
    masks = {'c00': lo1 & lo2, 'c0N': lo1 & hi2, 'cN0': hi1 & lo2, 'cNN': hi1 & hi2}
    corner = masks['c00'] | masks['c0N'] | masks['cN0'] | masks['cNN']
    masks.update(f1lo=lo1 & ~corner, f1hi=hi1 & ~corner, f2lo=lo2 & ~corner, f2hi=hi2 & ~corner)
    masks['bulk'] = ~(corner | lo1 | hi1 | lo2 | hi2)
    return masks


def run_prism(request):
    m = model(request['model'])
    n, h = request['size'], request.get('region_depth', 2)
    emin, emax = request['window_eV']
    out = {'model': request['model'], 'size': n, 'region_depth': h, 'geometry': geometry(m, request['basis']), 'source': NOTE}
    results = []
    t1, t2 = [x.ravel() for x in np.meshgrid(np.arange(n), np.arange(n), indexing='ij')]
    cells = np.stack([t1, t2], 1)
    onsite = None
    if request.get('potential'):
        onsite = np.zeros((len(cells), tb.NW))
        for p in request['potential']:
            d = p['depth']
            lo1, hi1, lo2, hi2 = t1 < d, t1 >= n - d, t2 < d, t2 >= n - d
            sel = {'c00': lo1 & lo2, 'c0N': lo1 & hi2, 'cN0': hi1 & lo2, 'cNN': hi1 & hi2,
                   'f1lo': lo1, 'f1hi': hi1, 'f2lo': lo2, 'f2hi': hi2}[p['region']]
            orb = {'all': np.ones(tb.NW, bool), 'Zn': tb.ZN, 'X': ~tb.ZN}[p['species']]
            onsite[np.ix_(sel, orb)] += p['V_eV']
    masks = region_masks(cells, n, h)
    for k in request['k_phases']:
        M, cells_m = m.finite_matrix(request['basis'], [n, n], [float(k)], onsite=onsite)
        assert np.array_equal(cells_m, cells)
        e, v, truncated = m.finite_states(M, emin, emax, nmax=150)
        w = (np.abs(v) ** 2).reshape(len(cells), tb.NW, -1)
        states = []
        for i in range(len(e)):
            states.append({'E_eV': round(float(e[i]), 6),
                           'weights': {r: round(float(w[masks[r], :, i].sum()), 4) for r in REGIONS},
                           'zn_s_weight': round(float(w[:, :2, i].sum()), 4)})
        results.append({'k_phase': k, 'states': states, 'truncated': bool(truncated)})
    out.update(results=results, note=('prism periodic along basis row 1 (k_phase = Bloch phase per period, radians) and '
                                      f'finite along rows 2 and 3 with {n} x {n} cells; weights are the probability in '
                                      'each region (corners cXY: cells within region_depth of both ends named; facets '
                                      'f1lo/f1hi/f2lo/f2hi: cells within region_depth of one end of side 1 or 2, corners '
                                      'excluded; bulk: the rest); at most 150 states per k (nearest the window centre; truncated = true flags an incomplete window)'))
    return out


def run(request, workdir):
    kind = request['kind']
    m = model(request['model'])
    if kind == 'bulk_bands':
        k = kpoints(request, 300)
        e, v = m.eig(k)
        zn = (np.abs(v[:, :2, :]) ** 2).sum(1)
        return {'model': request['model'], 'k_cart_invA': r6(k), 'energies_eV': r6(e), 'zn_s_weight': np.round(zn, 4).tolist(),
                'note': 'zn_s_weight: weight of each band on the two Zn-4s orbitals (the rest is pnictogen p)', 'source': NOTE}
    if kind == 'band_edges':
        n, mesh = request['band'], request['mesh']
        kf = grid(m, mesh)
        e = m.bands(kf @ m.B)
        lo, hi, gap = e[:, n - 1], e[:, n], e[:, n] - e[:, n - 1]
        top = np.argsort(-lo, kind='stable')[:6]
        bot = np.argsort(hi, kind='stable')[:6]
        i = int(np.argmin(gap))
        return {'model': request['model'], 'band': n, 'mesh': mesh,
                'band_n_highest': [{'E_eV': round(float(lo[j]), 6), **point(m, kf[j])} for j in top],
                'band_n_plus_1_lowest': [{'E_eV': round(float(hi[j]), 6), **point(m, kf[j])} for j in bot],
                'smallest_direct_gap': {'gap_eV': round(float(gap[i]), 6), **point(m, kf[i])},
                'note': f'values on the Gamma-centred {mesh}^3 grid of k_frac (no refinement)', 'source': NOTE}
    if kind == 'fermi_level':
        mesh, off, ne = request['mesh'], request.get('offset', 0), request['electrons_per_cell']
        e = m.bands(grid(m, mesh, off) @ m.B)
        flat = np.sort(e.ravel())
        nocc = int(round(ne * mesh ** 3))
        ef = 0.5 * (flat[nocc - 1] + flat[nocc])
        return {'model': request['model'], 'mesh': mesh, 'offset': off, 'electrons_per_cell': ne,
                'fermi_level_eV': round(float(ef), 6),
                'band_filling': np.round((e < ef).mean(0), 5).tolist(),
                'note': ('each band holds one electron per cell; E_F is midway between the highest occupied and lowest empty '
                         'state of the mesh; band_filling = fraction of k-points of each band below E_F'), 'source': NOTE}
    if kind == 'slab':
        L, h = request['layers'], request.get('region_depth', 2)
        emin, emax = request['window_eV']
        results = []
        for ph in request['k_phases']:
            M, cells = m.finite_matrix(request['basis'], [L], [float(ph[0]), float(ph[1])])
            e, v, truncated = m.finite_states(M, emin, emax, nmax=100)
            w = (np.abs(v) ** 2).reshape(len(cells), tb.NW, -1)
            t = cells[:, 0]
            results.append({'k_phases': ph, 'truncated': bool(truncated), 'states': [
                {'E_eV': round(float(e[i]), 6), 'weight_low': round(float(w[t < h, :, i].sum()), 4),
                 'weight_high': round(float(w[t >= L - h, :, i].sum()), 4), 'zn_s_weight': round(float(w[:, :2, i].sum()), 4)}
                for i in range(len(e))]})
        return {'model': request['model'], 'layers': L, 'region_depth': h, 'geometry': geometry(m, request['basis']),
                'results': results, 'source': NOTE,
                'note': ('slab periodic along basis rows 1 and 2 (k_phases = Bloch phases per period) and finite along row 3; '
                         'weight_low/high: probability in the cells within region_depth of the low/high end of row 3')}
    if kind == 'prism':
        return run_prism(request)
    raise ValueError('unknown kind')


# ---------------------------------------------------------------- grading
def grade(payload, trace, validation):
    t, th = PRIVATE['targets'], PRIVATE['thresholds']
    fields = ['strain_percent', 'band_inverted', 'fundamental_gap_meV', 'libizn_fermi_level_meV', 'hinge', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    sp, bi, gap, hinge = payload['strain_percent'], payload['band_inverted'], payload['fundamental_gap_meV'], payload['hinge']
    if not isinstance(sp, dict) or set(sp) != set(MODELS) or not all(number(v) for v in sp.values()):
        return False, {'error': 'strain_percent must give a number for each of ' + ', '.join(MODELS)}
    if not isinstance(bi, dict) or set(bi) != {'LiSbZn', 'LiBiZn'} or not all(isinstance(v, bool) for v in bi.values()):
        return False, {'error': 'band_inverted must be {"LiSbZn": true|false, "LiBiZn": true|false}'}
    if not isinstance(gap, dict) or set(gap) != set(t['fundamental_gap_meV']) or not all(number(v) for v in gap.values()):
        return False, {'error': 'fundamental_gap_meV must give a number for ' + ', '.join(t['fundamental_gap_meV'])}
    if not number(payload['libizn_fermi_level_meV']):
        return False, {'error': 'libizn_fermi_level_meV must be a finite number'}
    hk = {'termination', 'energy_k0_meV', 'band_window_meV', 'zn_s_weight_k0', 'removal_potential_eV'}
    if not isinstance(hinge, dict) or set(hinge) != hk:
        return False, {'error': 'hinge must have exactly the keys ' + ', '.join(sorted(hk))}
    ek, win = hinge['energy_k0_meV'], hinge['band_window_meV']
    if not isinstance(ek, dict) or set(ek) != set(t['hinge_energy_k0_meV']) or not all(number(v) for v in ek.values()):
        return False, {'error': 'hinge.energy_k0_meV must give a number for ' + ', '.join(t['hinge_energy_k0_meV'])}
    if not (isinstance(win, list) and len(win) == 2 and all(number(v) for v in win)):
        return False, {'error': 'hinge.band_window_meV must be [lowest, highest]'}
    if not number(hinge['zn_s_weight_k0']):
        return False, {'error': 'hinge.zn_s_weight_k0 must be a number'}
    rp = hinge['removal_potential_eV']
    if not (rp is None or number(rp)):
        return False, {'error': 'hinge.removal_potential_eV must be a number or null'}
    if not isinstance(hinge['termination'], str):
        return False, {'error': 'hinge.termination must be "Zn/Zn", "X/X" or "Zn/X"'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    for k in MODELS:
        errors[f'strain_percent:{k}'] = sp[k] - t['strain_percent'][k]
        checks[f'strain_percent:{k}'] = abs(errors[f'strain_percent:{k}']) <= th['strain_percent_abs']
    for k in ('LiSbZn', 'LiBiZn'):
        checks[f'band_inverted:{k}'] = bi[k] == t['band_inverted'][k]
    for k, ref in t['fundamental_gap_meV'].items():
        errors[f'gap:{k}'] = gap[k] - ref
        checks[f'gap:{k}'] = abs(errors[f'gap:{k}']) <= th['gap_abs_meV'][k]
    errors['libizn_fermi_level_meV'] = payload['libizn_fermi_level_meV'] - t['libizn_fermi_level_meV']
    checks['libizn_fermi_level'] = abs(errors['libizn_fermi_level_meV']) <= th['fermi_level_abs_meV']
    term = hinge['termination'].replace(' ', '')
    term = '/'.join(sorted(term.split('/'))) if term.count('/') == 1 else term
    checks['hinge_termination'] = term == '/'.join(sorted(t['hinge_termination'].split('/')))
    for k, ref in t['hinge_energy_k0_meV'].items():
        errors[f'hinge_energy_k0:{k}'] = ek[k] - ref
        checks[f'hinge_energy_k0:{k}'] = abs(errors[f'hinge_energy_k0:{k}']) <= th['hinge_energy_abs_meV']
    lo, hi = sorted(win)
    errors['hinge_band_bottom_meV'] = lo - t['hinge_band_window_meV'][0]
    errors['hinge_band_top_meV'] = hi - t['hinge_band_window_meV'][1]
    checks['hinge_band_bottom'] = abs(errors['hinge_band_bottom_meV']) <= th['band_bottom_abs_meV']
    checks['hinge_band_top'] = abs(errors['hinge_band_top_meV']) <= th['band_top_abs_meV']
    errors['hinge_zn_s_weight'] = hinge['zn_s_weight_k0'] - t['hinge_zn_s_weight_k0']
    checks['hinge_zn_s_weight'] = abs(errors['hinge_zn_s_weight']) <= th['zn_weight_abs']
    if rp is None:
        checks['hinge_removal_potential'] = False
        errors['hinge_removal_potential_eV'] = 'null'
    else:
        errors['hinge_removal_potential_eV'] = rp - t['hinge_removal_potential_eV']
        checks['hinge_removal_potential'] = abs(errors['hinge_removal_potential_eV']) <= th['removal_potential_abs_eV']
    return all(checks.values()), {'checks': checks, 'errors': errors, 'error_units': 'percent, meV, eV or weight as named'}
