"""Lab logic for zrsix-butterfly-amr: numpy evaluators of the archived ZrSiX Wannier models and grading.

The three Wannier Hamiltonians of Materials Cloud record tyan3-q3j67 (ZrSiS and ZrSiSe symmetrized
'wannier90_hr.dat_nsymm16', ZrSiTe 'wannier90_hr.dat'; 24 orbitals, no spin-orbit coupling) are stored in
data/model_<compound>.npz and evaluated exactly for band energies, velocities, Fermi-surface slices and k-space orbits
(zrsix_lab.py). Electron counts use band-energy histograms of the uniform meshes (data/counts_<compound>.npz).
Magnetotransport uses precomputed Chambers tables (data/mr_<compound>.npz): orbit-averaged <v_i vbar_j> of every
start point of the 160 x 160 x 80 mesh in an energy window, binned by energy and band, for the 40, 80 and 160 meshes
and the menu of field directions. Graded values are held in private/hidden.json and never served.
"""
import json
import math
from pathlib import Path
import numpy as np
import zrsix_lab as Z

HERE = Path(__file__).parent
DATA = HERE/'data'
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
PUBLIC = json.loads((HERE/'public.json').read_text())
SPEC = json.loads((DATA/'spec.json').read_text())
COMPOUNDS = list(SPEC['compounds'])
TRANSPORT_COST = {int(k): v for k, v in SPEC['transport_cost'].items()}
_M, _C, _T = {}, {}, {}


def model(c):
    if c not in _M:
        _M[c] = Z.Model(DATA/f'model_{c}.npz')
    return _M[c]


def counts(c):
    if c not in _C:
        z = np.load(DATA/f'counts_{c}.npz')
        _C[c] = {k: z[k] for k in z.files}
    return _C[c]


def table(c):
    if c not in _T:
        z = np.load(DATA/f'mr_{c}.npz')
        _T[c] = {k: z[k] for k in z.files}
    return _T[c]


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def vec3(x, name):
    if not (isinstance(x, list) and len(x) == 3 and all(number(v) and abs(v) < 100 for v in x)):
        raise ValueError(f'{name} must be a list of three numbers')
    return np.array(x, float)


def compound_of(req):
    c = req.get('compound')
    if c not in COMPOUNDS:
        raise ValueError('compound must be one of ' + ', '.join(COMPOUNDS))
    return c


def klist(req, limit):
    pts = req.get('k_frac')
    if not (isinstance(pts, list) and 1 <= len(pts) <= limit and
            all(isinstance(p, list) and len(p) == 3 and all(number(x) and abs(x) < 50 for x in p) for p in pts)):
        raise ValueError(f'k_frac must be a list of 1 to {limit} points [k1, k2, k3]')
    return np.array(pts, float)


def mu_list(req, c, lo, hi, limit):
    mu = req.get('mu_eV')
    if not (isinstance(mu, list) and 1 <= len(mu) <= limit and all(number(x) and lo <= x <= hi for x in mu)):
        raise ValueError(f'mu_eV must be a list of 1 to {limit} chemical potentials between {lo} and {hi} eV for {c}')
    return np.array(mu, float)


def temperature(req, lo, hi):
    T = req.get('temperature_K')
    if not (number(T) and lo <= T <= hi):
        raise ValueError(f'temperature_K must be a number from {lo} to {hi}')
    return float(T)


def check_keys(req, allowed, usage):
    if not set(req) <= allowed | {'kind'}:
        raise ValueError('Use ' + usage)


def field_of(req, c):
    f = req.get('field')
    planes = SPEC['field_menu'][c]
    if not (isinstance(f, dict) and set(f) == {'plane', 'angle_deg'} and f['plane'] in planes and number(f['angle_deg'])
            and any(abs(f['angle_deg'] - a) < 1e-6 for a in planes[f['plane']])):
        raise ValueError('field must be {"plane": "bc"|"ac"|"ab", "angle_deg": a} with a from the menu in /app/catalog.json')
    return SPEC['dir_index'][c][f"{f['plane']}:{float(f['angle_deg']):.1f}"]


def validate(req):
    kind = req.get('kind')
    if kind == 'bands':
        check_keys(req, {'compound', 'k_frac', 'velocities'}, '{"kind": "bands", "compound": c, "k_frac": [[k1, k2, k3], ...] (<= 500), "velocities": false}')
        compound_of(req); klist(req, 500)
        if not isinstance(req.get('velocities', False), bool):
            raise ValueError('velocities must be true or false')
        return 1
    if kind == 'fermi_slice':
        check_keys(req, {'compound', 'origin', 'u', 'v', 'n'}, '{"kind": "fermi_slice", "compound": c, "origin": [..], "u": [..], "v": [..], "n": [nu, nv]}')
        compound_of(req); vec3(req.get('origin'), 'origin'); vec3(req.get('u'), 'u'); vec3(req.get('v'), 'v')
        n = req.get('n')
        if not (isinstance(n, list) and len(n) == 2 and all(type(x) is int and 2 <= x <= 121 for x in n) and n[0] * n[1] <= 121 * 121):
            raise ValueError('n must be [nu, nv] with 2 <= nu, nv <= 121')
        return 1
    if kind == 'electron_count':
        check_keys(req, {'compound', 'mesh', 'mu_eV', 'temperature_K'}, '{"kind": "electron_count", "compound": c, "mesh": N, "mu_eV": [...], "temperature_K": T}')
        c = compound_of(req)
        if req.get('mesh') not in SPEC['count_meshes'] or type(req.get('mesh')) is not int:
            raise ValueError('mesh must be one of ' + ', '.join(map(str, SPEC['count_meshes'])))
        mu_list(req, c, *SPEC['count_range'][c], 200); temperature(req, 0, 300)
        return 1
    if kind == 'orbit':
        check_keys(req, {'compound', 'band', 'k_frac', 'field_direction', 'max_length_invA'},
                   '{"kind": "orbit", "compound": c, "band": n, "k_frac": [[...]] (<= 8), "field_direction": [bx, by, bz], "max_length_invA": L}')
        c = compound_of(req)
        if type(req.get('band')) is not int or not 1 <= req['band'] <= 24:
            raise ValueError('band must be an integer from 1 to 24')
        klist(req, 8)
        b = vec3(req.get('field_direction'), 'field_direction')
        if np.linalg.norm(b) < 1e-6:
            raise ValueError('field_direction must be nonzero')
        L = req.get('max_length_invA', 20)
        if not (number(L) and 0.1 <= L <= 40):
            raise ValueError('max_length_invA must be between 0.1 and 40 (1/Angstrom)')
        return 1
    if kind == 'magnetotransport':
        check_keys(req, {'compound', 'field', 'mesh', 'mu_eV', 'temperature_K', 'bands'},
                   '{"kind": "magnetotransport", "compound": c, "field": {"plane": p, "angle_deg": a}, "mesh": N, "mu_eV": [...], "temperature_K": T, "bands": [...]}')
        c = compound_of(req)
        field_of(req, c)
        if type(req.get('mesh')) is not int or req['mesh'] not in TRANSPORT_COST:
            raise ValueError('mesh must be one of ' + ', '.join(map(str, TRANSPORT_COST)))
        mu_list(req, c, *SPEC['transport_mu_range'][c], 4); temperature(req, 10, 40)
        bands = req.get('bands')
        allowed = SPEC['transport_bands'][c]
        if not (isinstance(bands, list) and bands and all(type(b) is int for b in bands) and len(set(bands)) == len(bands)):
            raise ValueError('bands must be a non-empty list of band numbers')
        bad = [b for b in bands if b not in allowed]
        if bad:
            raise ValueError(f'bands {bad} have no states in the transport energy window of {c}; available: {allowed}')
        return TRANSPORT_COST[req['mesh']]
    raise ValueError('kind must be one of bands, fermi_slice, electron_count, orbit, magnetotransport')


def run(req, workdir):
    kind, c = req['kind'], req['compound']
    m = model(c)
    out = {'compound': c, 'source': f"new evaluation of the archived Wannier Hamiltonian {SPEC['files'][c]}"}
    if kind == 'bands':
        k = klist(req, 500)
        if req.get('velocities'):
            e, v = m.bands_vel(k)
            out.update(energies_eV=e.round(6).tolist(), velocity_eV_A=v.round(5).tolist(),
                       note='velocity = dE/dk (eV Angstrom, Cartesian; divide by hbar for m/s), per k, band, component')
        else:
            out.update(energies_eV=m.bands(k).round(6).tolist())
        out['k_cart_invA'] = (k @ m.B).round(8).tolist()
        return out
    if kind == 'fermi_slice':
        o, u, v = (np.array(req[x], float) for x in ('origin', 'u', 'v'))
        nu, nv = req['n']
        a, b = np.meshgrid(np.arange(nu) / (nu - 1), np.arange(nv) / (nv - 1), indexing='ij')
        k = o + a.reshape(-1, 1) * u + b.reshape(-1, 1) * v
        e = m.bands(k)
        lo, hi = SPEC['slice_bands']
        out.update(bands=list(range(lo, hi + 1)), energies_eV=e[:, lo - 1:hi].reshape(nu, nv, -1).round(5).tolist(),
                   note='energies_eV[i][j][b]: k = origin + i/(nu-1) u + j/(nv-1) v (fractional), band = bands[b]')
        return out
    if kind == 'electron_count':
        z = counts(c)
        mi = list(z['meshes']).index(req['mesh'])
        mu = np.array(req['mu_eV'], float); T = float(req['temperature_K'])
        centres = float(z['lo']) + (np.arange(z['hist'].shape[2]) + 0.5) * float(z['ebin'])
        n = req['mesh']; nk = n * n * (n // 2)
        occ = np.zeros((len(mu), z['hist'].shape[1]))
        for i, x in enumerate(mu):
            _, f = Z.fermi_weights(centres, x, Z.KB * T)
            occ[i] = (z['nlow'][mi] + z['hist'][mi] @ f) / nk
        out.update(mesh=[n, n, n // 2], k_points=nk, temperature_K=T, mu_eV=mu.tolist(),
                   electrons_per_cell=occ.sum(1).round(7).tolist(),
                   band_occupation=occ.round(7).tolist(),
                   note='Fermi-Dirac occupation of the model bands on the Gamma-centred mesh, one state per band and k-point '
                        '(the model has no spin index). band_occupation[i][b]: mean occupation of band b+1 at mu_eV[i]; '
                        'electrons_per_cell is their sum over the 24 bands.')
        return out
    if kind == 'orbit':
        k0 = klist(req, 8)
        r = Z.trace_orbits(m, req['band'] - 1, k0, np.array(req['field_direction'], float), float(req.get('max_length_invA', 20)))
        orbits = []
        for i in range(len(k0)):
            p = np.array(r['paths'][i]); step = max(1, len(p) // 300)
            orbits.append({'energy_eV': round(float(r['energy'][i]), 6), 'closed': bool(r['closed'][i]),
                           'G_per_period': r['G'][i].tolist() if r['closed'][i] else None,
                           'period_Tps': None if not r['closed'][i] else round(float(r['period'][i]), 5),
                           'cyclotron_mass_me': None if not r['closed'][i] else round(float(r['mass'][i]), 5),
                           'arc_length_invA': round(float(r['length'][i]), 4),
                           'path_k_frac': p[::step].round(5).tolist()})
        out.update(band=req['band'], orbits=orbits,
                   note='Semiclassical orbit hbar dk/dt = -e v x B on the constant-energy surface of the start point, '
                        'traced backward in time (the history that enters vbar) up to max_length_invA. closed: returned to '
                        'the start modulo a reciprocal-lattice vector G_per_period (fractional; G != 0 is a periodic open '
                        'orbit). period_Tps: B*T_c in T ps; cyclotron_mass_me = e*period/(2 pi m_e).')
        return out
    # magnetotransport
    t = table(c)
    di = field_of(req, c)
    mi = list(t['mesh']).index(req['mesh'])
    bands = req['bands']; group = list(t['bands'])
    mu = np.array(req['mu_eV'], float); T = float(req['temperature_K']); kT = Z.KB * T
    centres = float(t['lo']) + np.arange(t['K'].shape[3]) * float(t['ebin'])
    n = int(req['mesh']); nk = n * n * (n // 2)
    to_ms = 1e-10 / Z.HBAR
    pref = 2 * Z.E ** 2 / (float(t['volume']) * 1e-30 * nk) / Z.E * to_ms ** 2 * 1e-12     # S/(m ps)
    res = []
    for x in mu:
        w, _ = Z.fermi_weights(centres, x, kT)
        per_band = {}
        total = np.zeros((len(t['X']), 3, 3))
        openw = {}
        for b in bands:
            ib = group.index(b)
            s = pref * np.einsum('e,exij->xij', w, t['K'][di, mi, ib].astype(float))
            per_band[str(b)] = s.tolist()
            total += s
            ww = (w * t['W'][mi, ib]).sum()
            openw[str(b)] = float((w * t['O'][di, mi, ib]).sum() / ww) if ww > 0 else 0.0
        rho = np.linalg.inv(total)
        res.append({'mu_eV': float(x), 'rho_tau_Ohm_m_ps': rho.tolist(), 'sigma_over_tau_S_per_m_ps': total.tolist(),
                    'sigma_over_tau_by_band': per_band, 'open_orbit_fraction_by_band': openw})
    f = req['field']
    out.update(field=f, field_direction=t['bvec'][di].round(6).tolist(), mesh=[n, n, n // 2], temperature_K=T, bands=bands,
               Btau_Tps=t['X'].tolist(), results=res,
               note='Chambers solution of the Boltzmann equation with one constant relaxation time tau for all bands. '
                    'Tensors are [Btau][i][j] with i, j = x, y, z; rho_tau = (sum over the requested bands of sigma/tau)^-1, '
                    'so rho = rho_tau / tau. E_i = rho_ij j_j. open_orbit_fraction_by_band: thermal-weighted fraction of '
                    'states whose orbit for this field direction is open.')
    return out


FIELDS = ['neutral_mu_eV', 'carrier_density_cm3', 'butterfly_peak_angle_deg', 'mr_b_over_c', 'mr_exponent_b',
          'zrsite_mr_c_over_b', 'mr_ratio_se_over_te', 'evidence_ids']


def grade(payload, trace, validation):
    t, th = PRIVATE['targets'], PRIVATE['thresholds']
    missing = [f for f in FIELDS if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    for key in ('neutral_mu_eV', 'carrier_density_cm3'):
        d = payload[key]
        if not isinstance(d, dict) or set(d) != set(t[key]) or not all(number(v) for v in d.values()):
            return False, {'error': f'{key} must give a finite number for each of ' + ', '.join(t[key])}
    scalars = ['butterfly_peak_angle_deg', 'mr_b_over_c', 'mr_exponent_b', 'zrsite_mr_c_over_b', 'mr_ratio_se_over_te']
    if not all(number(payload[k]) for k in scalars):
        return False, {'error': ', '.join(scalars) + ' must be finite numbers'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    for c, ref in t['neutral_mu_eV'].items():
        errors[f'neutral_mu_eV:{c}'] = payload['neutral_mu_eV'][c] - ref
        checks[f'neutral_mu:{c}'] = abs(errors[f'neutral_mu_eV:{c}']) <= th['mu_abs_eV']
    for c, ref in t['carrier_density_cm3'].items():
        errors[f'carrier_density_rel:{c}'] = payload['carrier_density_cm3'][c] / ref - 1
        checks[f'carrier_density:{c}'] = abs(errors[f'carrier_density_rel:{c}']) <= th['density_rel']
    errors['butterfly_peak_angle_deg'] = payload['butterfly_peak_angle_deg'] - t['butterfly_peak_angle_deg']
    checks['butterfly_peak_angle'] = abs(errors['butterfly_peak_angle_deg']) <= th['angle_abs_deg']
    errors['mr_exponent_b'] = payload['mr_exponent_b'] - t['mr_exponent_b']
    checks['mr_exponent_b'] = abs(errors['mr_exponent_b']) <= th['exponent_abs']
    for key in ('mr_b_over_c', 'zrsite_mr_c_over_b', 'mr_ratio_se_over_te'):
        errors[key + '_rel'] = payload[key] / t[key] - 1
        checks[key] = abs(errors[key + '_rel']) <= th[key + '_rel']
    return all(checks.values()), {'checks': checks, 'errors': errors}
