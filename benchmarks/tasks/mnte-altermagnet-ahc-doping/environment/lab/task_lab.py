"""Lab logic for mnte-altermagnet-ahc-doping: a numpy evaluator of the archived MnTe Wannier Hamiltonian.

The archived wannier90_hr.dat (Zenodo 18863304, dft_data/2_wannier90; 32 spinor Wannier functions, Neel vector
along +y) is stored as data/archived_hr.npz and evaluated by mnte_tb.py: band energies and spin expectation
values, band-resolved Berry curvature at k-lists, and the intrinsic anomalous Hall conductivity and electron
count on uniform meshes. The mesh requests are evaluated from Berry-curvature histograms of the same evaluator,
precomputed on each mesh of the menu (data/mesh_N.npz); the Neel-reversed model (all moments reversed) is the
time reverse of the archived one. The archived WannierTools output (sigma_ahc.txt) is withheld; graded values
are held in private/hidden.json.
"""
import json
import math
from pathlib import Path
import numpy as np
import mnte_tb as T

HERE = Path(__file__).parent
DATA = HERE/'data'
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
MODEL_FILE = DATA/'archived_hr.npz'
MESH_COST = {48: 1, 72: 1, 100: 1, 150: 2, 200: 3, 300: 5, 400: 8}
MU_RANGE = (5.05, 8.55)
KT_MAX = 0.025
_MODELS, _CACHES = {}, {}
NOTE = ('new evaluation of the archived Wannier90 Hamiltonian (dft_data/2_wannier90/wannier90_hr.dat); energies in eV '
        'on the absolute scale of the hr.dat (the VASP scale); Wannier bands numbered 1-32 from the bottom')


def model(neel):
    if neel not in _MODELS:
        _MODELS[neel] = T.Model(MODEL_FILE, neel=+1 if neel == '+y' else -1)
    return _MODELS[neel]


def cache(n):
    if n not in _CACHES:
        z = np.load(DATA/f'mesh_{n}.npz')
        _CACHES[n] = {k: z[k] for k in ('nk', 'hist', 'below', 'count', 'nlow')}
    return _CACHES[n]


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def kpoints(request, limit):
    given = [k for k in ('k_frac', 'k_cart_invA') if k in request]
    if len(given) != 1:
        raise ValueError('give exactly one of k_frac or k_cart_invA')
    pts = request[given[0]]
    if (not isinstance(pts, list) or not 1 <= len(pts) <= limit or
            not all(isinstance(p, list) and len(p) == 3 and all(number(x) and abs(x) < 50 for x in p) for p in pts)):
        raise ValueError(f'{given[0]} must be a list of 1 to {limit} points [x, y, z]')
    pts = np.array(pts, float)
    return pts @ T.B if given[0] == 'k_frac' else pts


def neel_of(request):
    neel = request.get('neel', '+y')
    if neel not in ('+y', '-y'):
        raise ValueError('neel must be "+y" (archived state) or "-y" (all moments reversed)')
    return neel


def mesh_args(request):
    n = request.get('mesh')
    if type(n) is not int or n not in MESH_COST:
        raise ValueError('mesh must be one of ' + ', '.join(map(str, MESH_COST)) + ' (an N x N x N/2 Gamma-centred mesh)')
    mu = request.get('mu_eV')
    if (not isinstance(mu, list) or not 1 <= len(mu) <= 500 or
            not all(number(x) and MU_RANGE[0] <= x <= MU_RANGE[1] for x in mu)):
        raise ValueError(f'mu_eV must be a list of 1 to 500 chemical potentials between {MU_RANGE[0]} and {MU_RANGE[1]} eV')
    kT = request.get('kT_eV')
    if not number(kT) or not 0 <= kT <= KT_MAX:
        raise ValueError(f'kT_eV must be a number from 0 to {KT_MAX} (Fermi-Dirac smearing; 0 = step occupation)')
    return n, np.array(mu, float), float(kT)


def validate(request):
    kind = request.get('kind')
    if kind == 'bands':
        if not set(request) <= {'kind', 'k_frac', 'k_cart_invA', 'neel', 'spin'}:
            raise ValueError('Use {"kind": "bands", "k_frac" or "k_cart_invA": [[x, y, z], ...] (<= 1000, or <= 500 with spin), '
                             '"neel": "+y"|"-y", "spin": false}')
        spin = request.get('spin', False)
        if not isinstance(spin, bool):
            raise ValueError('spin must be true or false')
        kpoints(request, 500 if spin else 1000); neel_of(request)
        return 1
    if kind == 'berry_curvature':
        if not set(request) <= {'kind', 'k_frac', 'k_cart_invA', 'neel'}:
            raise ValueError('Use {"kind": "berry_curvature", "k_frac" or "k_cart_invA": [[x, y, z], ...] (<= 500), "neel": "+y"|"-y"}')
        kpoints(request, 500); neel_of(request)
        return 1
    if kind == 'ahc':
        if not set(request) <= {'kind', 'mesh', 'mu_eV', 'kT_eV', 'neel'}:
            raise ValueError('Use {"kind": "ahc", "mesh": N, "mu_eV": [...], "kT_eV": x, "neel": "+y"|"-y"}')
        n, _, _ = mesh_args(request); neel_of(request)
        return MESH_COST[n]
    if kind == 'electron_count':
        if not set(request) <= {'kind', 'mesh', 'mu_eV', 'kT_eV'}:
            raise ValueError('Use {"kind": "electron_count", "mesh": N, "mu_eV": [...], "kT_eV": x}')
        mesh_args(request)
        return 1
    raise ValueError('kind must be one of bands, berry_curvature, ahc, electron_count')


def run(request, workdir):
    kind = request['kind']
    if kind in ('bands', 'berry_curvature'):
        neel = neel_of(request)
        m = model(neel)
        k = kpoints(request, 1000)
        out = {'neel': neel, 'k_cart_invA': k.round(8).tolist(), 'source': NOTE}
        if kind == 'bands':
            if request.get('spin', False):
                E, s = m.spin(k)
                out.update(energies_eV=E.round(6).tolist(),
                           spin_expectation={c: s[i].round(4).tolist() for i, c in enumerate('xyz')},
                           note='<sigma_x>, <sigma_y>, <sigma_z> of each band (Pauli matrices, -1 to 1; spin quantised along '
                                'Cartesian z, SAXIS = 0 0 1)')
            else:
                out.update(energies_eV=m.bands(k).round(6).tolist())
            return out
        E, om = m.curvature(k)
        out.update(energies_eV=E.round(6).tolist(),
                   omega_A2={c: om[i].round(5).tolist() for i, c in enumerate(('x', 'y', 'z'))},
                   note='band-resolved Berry curvature Omega_n = curl_k i<u_n|grad_k u_n> (Angstrom^2) from the Kubo formula '
                        'with the velocity operator in the atomic (Wannier-centre) gauge, as in WannierTools; '
                        'Omega_z enters sigma_xy, Omega_x sigma_yz, Omega_y sigma_zx')
        return out
    n, mu, kT = mesh_args(request)
    c = cache(n)
    ne = T.electrons_from_cache(c, mu, kT)
    out = {'mesh': [n, n, n // 2], 'k_points': int(c['nk']), 'kT_eV': kT, 'mu_eV': mu.round(6).tolist(),
           'electrons_per_cell': ne.round(7).tolist(),
           'note': 'Gamma-centred mesh of fractional coordinates; electrons per cell in the 32 Wannier bands with '
                   'Fermi-Dirac occupations (step function for kT = 0)', 'source': NOTE}
    if kind == 'ahc':
        neel = neel_of(request)
        s = T.sigma_from_cache(c, mu, kT) * (1 if neel == '+y' else -1)
        out.update(neel=neel, sigma_S_per_cm={'xy': s[2].round(4).tolist(), 'yz': s[0].round(4).tolist(),
                                              'zx': s[1].round(4).tolist()},
                   note=out['note'] + '; intrinsic anomalous Hall conductivity sigma_ab = (e^2/hbar) (1/(N_k V)) '
                        'sum_k sum_n f_n Omega^c_n in S/cm, (a, b, c) cyclic, sign of the archived WannierTools output')
    return out


def grade(payload, trace, validation):
    t, th = PRIVATE['targets'], PRIVATE['thresholds']
    fields = ['sigma_xy_S_per_cm', 'sigma_xy_sign_changes', 'sign_change_hole_density_cm3',
              'hole_density_0p30_below_vbm_cm3', 'forbidden_components_max_abs_S_per_cm', 'allowed_hall_components',
              'electron_doped_sigma_xy_S_per_cm', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    sx, fb, al = payload['sigma_xy_S_per_cm'], payload['forbidden_components_max_abs_S_per_cm'], payload['allowed_hall_components']
    if not isinstance(sx, dict) or set(sx) != set(t['sigma_xy_S_per_cm']) or not all(number(v) for v in sx.values()):
        return False, {'error': 'sigma_xy_S_per_cm must give a finite number for each hole density key: '
                                + ', '.join(t['sigma_xy_S_per_cm'])}
    if not isinstance(fb, dict) or set(fb) != {'yz', 'zx'} or not all(number(v) for v in fb.values()):
        return False, {'error': 'forbidden_components_max_abs_S_per_cm must be {"yz": x, "zx": x}'}
    comps = {'xy', 'yz', 'zx'}
    if (not isinstance(al, dict) or set(al) != {'neel_y', 'neel_x'} or
            not all(isinstance(v, list) and all(isinstance(c, str) and c in comps for c in v) for v in al.values())):
        return False, {'error': 'allowed_hall_components must be {"neel_y": [...], "neel_x": [...]} with entries "xy", "yz", "zx"'}
    scalars = ['sign_change_hole_density_cm3', 'hole_density_0p30_below_vbm_cm3']
    if not all(number(payload[k]) and payload[k] > 0 for k in scalars):
        return False, {'error': ', '.join(scalars) + ' must be positive numbers'}
    nsc, ed = payload['sigma_xy_sign_changes'], payload['electron_doped_sigma_xy_S_per_cm']
    if not number(nsc) or ed is not None and not number(ed):
        return False, {'error': 'sigma_xy_sign_changes must be a number and electron_doped_sigma_xy_S_per_cm a number or null'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    for key, ref in t['sigma_xy_S_per_cm'].items():
        errors[f'sigma_xy_rel:{key}'] = sx[key] / ref - 1
        checks[f'sigma_xy:{key}'] = abs(errors[f'sigma_xy_rel:{key}']) <= th['sigma_xy_rel']
    checks['sign_changes'] = nsc == t['sigma_xy_sign_changes']
    errors['sign_change_density_rel'] = payload['sign_change_hole_density_cm3'] / t['sign_change_hole_density_cm3'] - 1
    checks['sign_change_density'] = abs(errors['sign_change_density_rel']) <= th['sign_change_density_rel']
    errors['hole_density_rel'] = payload['hole_density_0p30_below_vbm_cm3'] / t['hole_density_0p30_below_vbm_cm3'] - 1
    checks['hole_density_0p30'] = abs(errors['hole_density_rel']) <= th['hole_density_rel']
    for c in ('yz', 'zx'):
        errors[f'forbidden_{c}_S_cm'] = fb[c]
        checks[f'forbidden_{c}'] = abs(fb[c]) <= th['forbidden_abs_S_cm']
    for k in ('neel_y', 'neel_x'):
        checks[f'allowed:{k}'] = sorted(set(al[k])) == sorted(t['allowed_hall_components'][k])
    checks['electron_doping_not_supported'] = ed is None
    return all(checks.values()), {'checks': checks, 'errors': errors}
