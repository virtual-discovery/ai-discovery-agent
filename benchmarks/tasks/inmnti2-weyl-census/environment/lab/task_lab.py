"""Lab logic for inmnti2-weyl-census: archived InMnTi2 results plus a numpy Wannier tight-binding evaluator.

Purchases return archived content of Materials Cloud record hkprp-xw234 (upload.zip): DFT band energies of
four non-self-consistent runs, the Kubo outputs of hubbard/wannier/kubo/<folder>, and WannierTools node and
chirality outputs. Computing requests evaluate the two archived Wannier90 Hamiltonians (hubbard and standard
InMnTi2_hr.dat, stored as data/*.npz) with tbmodel.py: band energies, gap scans, gap minimisation from seeds
and Berry flux through spheres. Graded values are derived quantities held in private/hidden.json.
"""
import json
import math
from pathlib import Path
import numpy as np
import tbmodel

HERE = Path(__file__).parent
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
ARCHIVE = json.loads((HERE/'data/archive.json').read_text())
MODEL_NAMES = ('hubbard', 'standard')
DFT_RUNS = tuple(ARCHIVE['dft'])
KUBO_RUNS = tuple(ARCHIVE['kubo'])
WT_RUNS = tuple(ARCHIVE['wannier_tools'])
_MODELS = {}
SOURCE = 'retrieved author-archived result (no new calculation)'
TB_NOTE = ('new evaluation of the archived Wannier90 Hamiltonian (InMnTi2_hr.dat of this model); energies in eV on the '
           'scale of the hr.dat; Wannier bands numbered from 1; k_cart in 1/Angstrom including 2 pi')


def model(name):
    if name not in _MODELS:
        _MODELS[name] = tbmodel.Model(HERE/'data'/f'{name}.npz')
    return _MODELS[name]


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def integer(x):
    return type(x) is int


def kpoints(request, keys, limit):
    given = [k for k in keys if k in request]
    if len(given) != 1:
        raise ValueError(f'give exactly one of {keys[0]} or {keys[1]}')
    pts = request[given[0]]
    if (not isinstance(pts, list) or not 1 <= len(pts) <= limit or
            not all(isinstance(p, list) and len(p) == 3 and all(number(x) and abs(x) < 20 for x in p) for p in pts)):
        raise ValueError(f'{given[0]} must be a list of 1 to {limit} points [x, y, z]')
    pts = np.array(pts, float)
    return pts if 'cart' in given[0] else pts @ tbmodel.B


def check_model(request):
    if request.get('model') not in MODEL_NAMES:
        raise ValueError('model must be "hubbard" or "standard"')
    return ARCHIVE['models'][request['model']]['num_wann']


def band(request, key, nw):
    if not integer(request.get(key)) or not 1 <= request[key] <= nw - 1:
        raise ValueError(f'{key} must be an integer from 1 to {nw - 1} (Wannier band numbering of this model)')


def validate(request):
    kind = request.get('kind')
    if kind == 'wannier_bands':
        if not set(request) <= {'kind', 'model', 'k_cart_invA', 'k_frac'}:
            raise ValueError('Use {"kind": "wannier_bands", "model": "hubbard"|"standard", "k_cart_invA" or "k_frac": [[x, y, z], ...] (<= 500)}')
        check_model(request); kpoints(request, ('k_cart_invA', 'k_frac'), 500)
        return 1
    if kind == 'gap_scan':
        if set(request) != {'kind', 'model', 'lower_band', 'mesh'}:
            raise ValueError('Use {"kind": "gap_scan", "model": ..., "lower_band": n, "mesh": N (6-36)}')
        nw = check_model(request); band(request, 'lower_band', nw)
        if not integer(request['mesh']) or not 6 <= request['mesh'] <= 36:
            raise ValueError('mesh must be an integer from 6 to 36')
        return 1 if request['mesh'] <= 24 else 2
    if kind == 'find_nodes':
        if not set(request) <= {'kind', 'model', 'lower_band', 'seeds_cart_invA', 'seeds_frac'} or 'lower_band' not in request:
            raise ValueError('Use {"kind": "find_nodes", "model": ..., "lower_band": n, "seeds_cart_invA" or "seeds_frac": [[x, y, z], ...] (<= 8)}')
        nw = check_model(request); band(request, 'lower_band', nw)
        kpoints(request, ('seeds_cart_invA', 'seeds_frac'), 8)
        return 1
    if kind == 'berry_flux':
        if set(request) != {'kind', 'model', 'occupied_bands', 'spheres'}:
            raise ValueError('Use {"kind": "berry_flux", "model": ..., "occupied_bands": n, "spheres": [{"center_cart_invA": [x, y, z], '
                             '"radius_invA": r, "n_theta": 16, "n_phi": 32}, ...] (<= 8)}')
        nw = check_model(request); band(request, 'occupied_bands', nw)
        sph = request['spheres']
        if not isinstance(sph, list) or not 1 <= len(sph) <= 8:
            raise ValueError('spheres must be a list of 1 to 8 spheres')
        for s in sph:
            if not isinstance(s, dict) or not {'center_cart_invA', 'radius_invA'} <= set(s) <= {'center_cart_invA', 'radius_invA', 'n_theta', 'n_phi'}:
                raise ValueError('each sphere needs center_cart_invA and radius_invA (optional n_theta, n_phi)')
            c = s['center_cart_invA']
            if not isinstance(c, list) or len(c) != 3 or not all(number(x) and abs(x) < 20 for x in c):
                raise ValueError('center_cart_invA must be [kx, ky, kz] in 1/Angstrom')
            if not number(s['radius_invA']) or not 0.001 <= s['radius_invA'] <= 0.5:
                raise ValueError('radius_invA must be between 0.001 and 0.5')
            if not (integer(s.get('n_theta', 16)) and 4 <= s.get('n_theta', 16) <= 40 and
                    integer(s.get('n_phi', 32)) and 8 <= s.get('n_phi', 32) <= 80):
                raise ValueError('n_theta must be an integer 4-40 and n_phi an integer 8-80')
        return 1
    if kind == 'dft_bands':
        if set(request) != {'kind', 'run'} or request['run'] not in DFT_RUNS:
            raise ValueError('Use {"kind": "dft_bands", "run": one of ' + ', '.join(DFT_RUNS) + '}')
        return 3
    if kind == 'kubo':
        if set(request) != {'kind', 'run', 'components'} or request['run'] not in KUBO_RUNS:
            raise ValueError('Use {"kind": "kubo", "run": one of ' + ', '.join(KUBO_RUNS) + ', "components": [...]}')
        comps = request['components']
        avail = ARCHIVE['kubo'][request['run']]
        if not isinstance(comps, list) or not comps or not all(isinstance(c, str) and c in avail for c in comps):
            raise ValueError(f"components of folder {request['run']}: " + ', '.join(sorted(avail)))
        return 2
    if kind == 'wannier_tools':
        if set(request) != {'kind', 'run'} or request['run'] not in WT_RUNS:
            raise ValueError('Use {"kind": "wannier_tools", "run": one of ' + ', '.join(WT_RUNS) + '}')
        return 2
    raise ValueError('kind must be one of wannier_bands, gap_scan, find_nodes, berry_flux, dft_bands, kubo, wannier_tools')


def run(request, workdir):
    kind = request['kind']
    if kind == 'dft_bands':
        return {**ARCHIVE['dft'][request['run']], 'source': SOURCE}
    if kind == 'kubo':
        data = ARCHIVE['kubo'][request['run']]
        return {'folder': f"hubbard/wannier/kubo/{request['run']}",
                'files': {c: data[c]['file'] for c in request['components']},
                'rows': {c: data[c]['rows'] for c in request['components']},
                'note': 'rows of the archived file with omega <= 1.5 eV, columns exactly as in the file', 'source': SOURCE}
    if kind == 'wannier_tools':
        return {'run': request['run'], 'files': ARCHIVE['wannier_tools'][request['run']], 'source': SOURCE}
    m = model(request['model'])
    out = {'model': request['model'], 'num_wann': m.nw, 'source': TB_NOTE}
    if kind == 'wannier_bands':
        k = kpoints(request, ('k_cart_invA', 'k_frac'), 500)
        out.update(k_cart_invA=k.round(8).tolist(), energies_eV=m.bands(k).round(6).tolist())
        return out
    if kind == 'gap_scan':
        n = request['lower_band']
        rows, nmin, gmin = m.gap_scan(n, request['mesh'])
        out.update(lower_band=n, mesh=request['mesh'], grid_local_minima=nmin, smallest_grid_gap_eV=round(gmin, 6),
                   minima=rows, note=f'local minima of E_{n + 1} - E_{n} on the Gamma-centred mesh^3 grid of k_frac, '
                                      'lowest 60 listed (no refinement)')
        return out
    if kind == 'find_nodes':
        n = request['lower_band']
        seeds = kpoints(request, ('seeds_cart_invA', 'seeds_frac'), 8)
        res = []
        for s in seeds:
            x, g, e, it = m.minimize_gap(s, n)
            kf = x @ np.linalg.inv(tbmodel.B)
            res.append({'seed_cart_invA': s.round(6).tolist(), 'k_cart_invA': x.round(6).tolist(),
                        'k_frac': (kf - np.floor(kf)).round(6).tolist(), 'gap_eV': float('%.3e' % g),
                        'mean_energy_eV': round(e, 6), 'iterations': it})
        out.update(lower_band=n, results=res,
                   note=f'BFGS minimisation of (E_{n + 1} - E_{n})^2 from each seed; k_frac wrapped to [0, 1)')
        return out
    if kind == 'berry_flux':
        n = request['occupied_bands']
        res = []
        for s in request['spheres']:
            nt, nph = s.get('n_theta', 16), s.get('n_phi', 32)
            flux, gmin = m.berry_flux(np.array(s['center_cart_invA'], float), float(s['radius_invA']), n, nt, nph)
            res.append({'center_cart_invA': s['center_cart_invA'], 'radius_invA': s['radius_invA'], 'n_theta': nt, 'n_phi': nph,
                        'flux_over_2pi': round(flux, 6), 'min_gap_on_sphere_eV': round(gmin, 6)})
        out.update(occupied_bands=n, results=res,
                   note=f'outward flux of curl A, A = i<u|grad u>, of Wannier bands 1..{n} through each sphere, divided by 2 pi '
                        '(link-variable discretisation on an n_theta x n_phi grid)')
        return out
    raise ValueError('unknown kind')


def grade(payload, trace, validation):
    t, th = PRIVATE['targets'], PRIVATE['thresholds']
    fields = ['node_count', 'orbit_sizes', 'node_energy_meV', 'pair_separation_invA', 'kz_offset_plus_invA',
              'pbe_node_energy_meV', 'ahc_S_per_cm', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    nc, orb, ahc = payload['node_count'], payload['orbit_sizes'], payload['ahc_S_per_cm']
    if not isinstance(nc, dict) or set(nc) != {'dftu', 'pbe'} or not all(number(v) for v in nc.values()):
        return False, {'error': 'node_count must be {"dftu": n, "pbe": n}'}
    if not isinstance(orb, list) or not orb or not all(number(v) for v in orb):
        return False, {'error': 'orbit_sizes must be a non-empty list of numbers'}
    if not isinstance(ahc, dict) or set(ahc) != {'xy', 'yz', 'zx'} or not all(number(v) for v in ahc.values()):
        return False, {'error': 'ahc_S_per_cm must be {"xy": x, "yz": x, "zx": x} with finite numbers'}
    scalars = ['node_energy_meV', 'pair_separation_invA', 'kz_offset_plus_invA', 'pbe_node_energy_meV']
    if not all(number(payload[k]) for k in scalars):
        return False, {'error': ', '.join(scalars) + ' must be finite numbers'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    for k in ('dftu', 'pbe'):
        checks[f'node_count:{k}'] = nc[k] == t['node_count'][k]
    checks['orbit_sizes'] = sorted(orb) == sorted(t['orbit_sizes'])
    errors['node_energy_meV'] = payload['node_energy_meV'] - t['node_energy_meV']
    checks['node_energy'] = abs(errors['node_energy_meV']) <= th['node_energy_abs_meV']
    errors['pair_separation_rel'] = payload['pair_separation_invA'] / t['pair_separation_invA'] - 1
    checks['pair_separation'] = abs(errors['pair_separation_rel']) <= th['pair_separation_rel']
    errors['kz_offset_invA'] = payload['kz_offset_plus_invA'] - t['kz_offset_plus_invA']
    checks['kz_offset'] = abs(errors['kz_offset_invA']) <= th['kz_offset_abs_invA']
    errors['pbe_node_energy_meV'] = payload['pbe_node_energy_meV'] - t['pbe_node_energy_meV']
    checks['pbe_node_energy'] = abs(errors['pbe_node_energy_meV']) <= th['pbe_node_energy_abs_meV']
    errors['ahc_xy_rel'] = ahc['xy'] / t['ahc_S_per_cm']['xy'] - 1
    checks['ahc_xy'] = abs(errors['ahc_xy_rel']) <= th['ahc_xy_rel']
    for c in ('yz', 'zx'):
        errors[f'ahc_{c}_S_cm'] = ahc[c] - t['ahc_S_per_cm'][c]
        checks[f'ahc_{c}'] = abs(errors[f'ahc_{c}_S_cm']) <= th['ahc_forbidden_abs_S_cm']
    return all(checks.values()), {'checks': checks, 'errors': errors,
                                  'error_units': 'meV, relative, 1/Angstrom or S/cm as named'}
