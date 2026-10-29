"""Lab logic for cu2wse4-orbital-edelstein: a numpy evaluator of the archived Cu2WSe4 Wannier Hamiltonians.

The two OpenMX Wannier Hamiltonians of Zenodo record 20131063 (Cu2WSe4_soc.HWR: 108 spinor functions;
Cu2WSe4_wosoc.HWR: 54 functions without spin-orbit coupling) are stored as data/model_<model>.npz and evaluated
by cws_tb.py: band energies, band velocities, Berry curvature and the spin, Berry-phase orbital and intra-atomic
orbital moments at k-lists. The Fermi-surface requests are evaluated from energy histograms of the same
evaluator, precomputed on each mesh of the menu (data/mesh_<model>_<N>.npz). The archived figure data of the
Edelstein and nonlinear Hall responses (Fig2-Fig4) are withheld; graded values are held in private/hidden.json.
"""
import json
import math
from pathlib import Path
import numpy as np
import cws_tb as T

HERE = Path(__file__).parent
DATA = HERE/'data'
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
MODELS = {'soc': 'model_soc.npz', 'nosoc': 'model_nosoc.npz'}
MESH_COST = {'soc': {16: 1, 24: 1, 32: 1, 48: 2, 64: 3, 96: 6, 128: 10},
             'nosoc': {16: 1, 24: 1, 32: 1, 48: 1, 64: 2, 96: 3, 128: 5}}
MU_RANGE = (-1.9, 1.9)
KT_RANGE = (0.005, 0.05)
_MODELS, _CACHES, _FERMI = {}, {}, {}


def model(name):
    if name not in _MODELS:
        _MODELS[name] = T.Model(DATA/MODELS[name])
        _FERMI[name] = float(np.load(DATA/MODELS[name])['fermi_Ha']) * T.HARTREE
    return _MODELS[name]


def fermi(name):
    model(name)
    return _FERMI[name]


def cache(name, n):
    key = (name, n)
    if key not in _CACHES:
        z = np.load(DATA/f'mesh_{name}_{n}.npz')
        c = {k: z[k].astype(float) for k in z.files}
        c['centres'] = c['e0'] + (np.arange(len(c['count'])) + 0.5) * T.EBIN
        _CACHES[key] = c
    return _CACHES[key]


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def model_of(request):
    name = request.get('model')
    if name not in MODELS:
        raise ValueError('model must be "soc" (Cu2WSe4_soc.HWR) or "nosoc" (Cu2WSe4_wosoc.HWR)')
    return name


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


def mesh_args(request, name):
    n = request.get('mesh')
    if type(n) is not int or n not in MESH_COST[name]:
        raise ValueError('mesh must be one of ' + ', '.join(map(str, MESH_COST[name])) + ' (an N x N x N/2 Gamma-centred mesh)')
    mu = request.get('mu_eV')
    if (not isinstance(mu, list) or not 1 <= len(mu) <= 400 or
            not all(number(x) and MU_RANGE[0] <= x <= MU_RANGE[1] for x in mu)):
        raise ValueError(f'mu_eV must be a list of 1 to 400 chemical potentials between {MU_RANGE[0]} and {MU_RANGE[1]} eV')
    kT = request.get('kT_eV')
    if not number(kT) or not KT_RANGE[0] <= kT <= KT_RANGE[1]:
        raise ValueError(f'kT_eV must be a number from {KT_RANGE[0]} to {KT_RANGE[1]} (Fermi-Dirac temperature)')
    return n, np.array(mu, float), float(kT)


def validate(request):
    kind = request.get('kind')
    if kind == 'bands':
        if not set(request) <= {'kind', 'model', 'k_frac', 'k_cart_invA'}:
            raise ValueError('Use {"kind": "bands", "model": "soc"|"nosoc", "k_frac" or "k_cart_invA": [[x, y, z], ...] (<= 1000)}')
        model_of(request); kpoints(request, 1000)
        return 1
    if kind == 'band_properties':
        if not set(request) <= {'kind', 'model', 'k_frac', 'k_cart_invA'}:
            raise ValueError('Use {"kind": "band_properties", "model": "soc"|"nosoc", "k_frac" or "k_cart_invA": [[x, y, z], ...] (<= 200)}')
        model_of(request); kpoints(request, 200)
        return 1
    if kind == 'fermi_surface':
        if not set(request) <= {'kind', 'model', 'mesh', 'mu_eV', 'kT_eV'}:
            raise ValueError('Use {"kind": "fermi_surface", "model": "soc"|"nosoc", "mesh": N, "mu_eV": [...], "kT_eV": x}')
        name = model_of(request)
        n, _, _ = mesh_args(request, name)
        return MESH_COST[name][n]
    raise ValueError('kind must be one of bands, band_properties, fermi_surface')


NOTE = ('new evaluation of the archived OpenMX Wannier Hamiltonian ({file}); energies in eV relative to the Fermi '
        'level printed in its header; bands numbered from 1 at the bottom of the model')
FILES = {'soc': 'Cu2WSe4_soc.HWR', 'nosoc': 'Cu2WSe4_wosoc.HWR'}


def run(request, workdir):
    kind, name = request['kind'], request['model']
    m, ef = model(name), fermi(name)
    out = {'model': name, 'source': NOTE.format(file=FILES[name])}
    if kind == 'bands':
        k = kpoints(request, 1000)
        out.update(k_cart_invA=k.round(8).tolist(), energies_eV=(m.bands(k) - ef).round(6).tolist())
        return out
    if kind == 'band_properties':
        k = kpoints(request, 200)
        p = m.properties(k)
        vec = lambda a, nd: {c: a[i].round(nd).tolist() for i, c in enumerate('xyz')}
        out.update(k_cart_invA=k.round(8).tolist(), energies_eV=(p['E'] - ef).round(6).tolist(),
                   velocity_eV_A=vec(p['v'], 6), berry_curvature_A2=vec(p['omega'], 5),
                   spin_moment_muB=None if p['m_spin'] is None else vec(p['m_spin'], 6),
                   orbital_moment_berry_phase_muB=vec(p['m_orb'], 6),
                   orbital_moment_intra_atomic_muB=vec(p['m_atomic'], 6),
                   note='per band and k-point: velocity = dE/dk (eV Angstrom, divide by hbar for m/s); Berry curvature '
                        'Omega = curl_k i<u|grad_k u> (Angstrom^2); spin moment -mu_B <sigma>; Berry-phase orbital moment '
                        '(|e|/2 hbar) Im <grad_k u| x (H - E)|grad_k u>; intra-atomic orbital moment -mu_B <L> with L '
                        'the angular momentum of the d and p Wannier functions about their own atom. Wannier centres at '
                        'the atomic positions; band pairs closer than 5 meV are treated as degenerate.')
        return out
    n, mu, kT = mesh_args(request, name)
    c = cache(name, n)
    mu_abs = mu + ef
    tensors = {}
    for key, label in (('spin', 'spin'), ('orb', 'orbital_berry_phase'), ('atomic', 'orbital_intra_atomic'),
                       ('bcd', 'berry_curvature_dipole')):
        if key == 'spin' and name == 'nosoc':
            continue
        S = T.surface_integral(c, key, mu_abs, kT) * T.VOL       # per cell
        tensors[label] = [s.round(9).tolist() for s in S]
    w = 1.0 / (np.exp(np.clip((c['centres'][None, :] - mu_abs[:, None]) / kT, -300, 300)) + 1.0)
    electrons = (c['nlow'] + w @ c['count']) / c['nk']
    out.update(mesh=[n, n, n // 2], k_points=int(c['nk']), kT_eV=kT, mu_eV=mu.round(6).tolist(),
               electrons_per_cell=electrons.round(6).tolist(), tensors_per_cell=tensors,
               note='Gamma-centred mesh of fractional coordinates. electrons_per_cell: Fermi-Dirac occupation of the '
                    'model bands (one state per band and k-point). For each mu, tensors_per_cell[q][i][j] = '
                    '(1/N_k) sum_nk (-df/dE) X_ij with X_ij = m_i dE/dk_j for the moments (spin, Berry-phase orbital, '
                    'intra-atomic orbital; mu_B eV Angstrom per eV, i.e. mu_B Angstrom) and X_ij = (dE/dk_i) Omega_j for '
                    'berry_curvature_dipole (Angstrom^3); sums over the bands of the model, per unit cell of the model.')
    return out


def grade(payload, trace, validation):
    t, th = PRIVATE['targets'], PRIVATE['thresholds']
    fields = ['orbital_to_spin_ratio', 'orbital_magnetization_A_per_m', 'soc_to_nosoc_orbital_ratio',
              'edelstein_allowed_components', 'alpha_yy_over_alpha_xx', 'chirality_indicator',
              'nonlinear_hall_allowed_components', 'sigma_zxy_over_sigma_xyz', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    r = payload['orbital_to_spin_ratio']
    if not isinstance(r, dict) or set(r) != set(t['orbital_to_spin_ratio']) or not all(number(v) for v in r.values()):
        return False, {'error': 'orbital_to_spin_ratio must give a finite number for each key: ' + ', '.join(t['orbital_to_spin_ratio'])}
    scalars = ['orbital_magnetization_A_per_m', 'soc_to_nosoc_orbital_ratio', 'alpha_yy_over_alpha_xx',
               'chirality_indicator', 'sigma_zxy_over_sigma_xyz']
    if not all(number(payload[k]) for k in scalars):
        return False, {'error': ', '.join(scalars) + ' must be finite numbers'}
    ea, na = payload['edelstein_allowed_components'], payload['nonlinear_hall_allowed_components']
    if not (isinstance(ea, list) and all(isinstance(c, str) and len(c) == 2 and set(c) <= set('xyz') for c in ea)):
        return False, {'error': 'edelstein_allowed_components must be a list of strings "ij" with i, j in x, y, z'}
    if not (isinstance(na, list) and all(isinstance(c, str) and len(c) == 3 and set(c) <= set('xyz') for c in na)):
        return False, {'error': 'nonlinear_hall_allowed_components must be a list of strings "ijl" with i, j, l in x, y, z'}
    executed = {rec['id'] for rec in trace if rec['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    for key, ref in t['orbital_to_spin_ratio'].items():
        errors[f'orbital_to_spin_ratio_rel:{key}'] = r[key] / ref - 1
        checks[f'orbital_to_spin_ratio:{key}'] = abs(r[key] / ref - 1) <= th['ratio_rel']
    errors['orbital_magnetization_rel'] = payload['orbital_magnetization_A_per_m'] / t['orbital_magnetization_A_per_m'] - 1
    checks['orbital_magnetization'] = abs(errors['orbital_magnetization_rel']) <= th['magnetization_rel']
    errors['soc_to_nosoc_rel'] = payload['soc_to_nosoc_orbital_ratio'] / t['soc_to_nosoc_orbital_ratio'] - 1
    checks['soc_to_nosoc_orbital_ratio'] = abs(errors['soc_to_nosoc_rel']) <= th['soc_ratio_rel']
    checks['edelstein_allowed_components'] = sorted(set(ea)) == sorted(t['edelstein_allowed_components'])
    errors['alpha_yy_over_alpha_xx'] = payload['alpha_yy_over_alpha_xx']
    checks['alpha_yy_over_alpha_xx'] = abs(payload['alpha_yy_over_alpha_xx'] - t['alpha_yy_over_alpha_xx']) <= th['symmetry_abs']
    errors['chirality_indicator'] = payload['chirality_indicator']
    checks['chirality_indicator'] = abs(payload['chirality_indicator'] - t['chirality_indicator']) <= th['symmetry_abs']
    checks['nonlinear_hall_allowed_components'] = sorted(set(na)) == sorted(t['nonlinear_hall_allowed_components'])
    errors['sigma_zxy_over_sigma_xyz'] = payload['sigma_zxy_over_sigma_xyz']
    checks['sigma_zxy_over_sigma_xyz'] = abs(payload['sigma_zxy_over_sigma_xyz'] - t['sigma_zxy_over_sigma_xyz']) <= th['nlh_ratio_abs']
    return all(checks.values()), {'checks': checks, 'errors': errors}
