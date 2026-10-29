"""Lab logic for hydride-coulomb-mustar: archived Mantel/IsoME data of Bozier & Pickard (Zenodo 23235290)
and a numpy isotropic Eliashberg solver (eliashberg.py) that reproduces the archived IsoME runs.

Retrieval requests return archived files. `eliashberg` runs eliashberg.py on an archived alpha2F and DOS,
with a fixed mu* (constant or energy-dependent DOS) or with the archived screened interaction W(e,e').
At the withheld pressures (Mg2IrH6 15/25 GPa, H3S 180/190 GPa) no Coulomb data or archived runs are served
and the W solver is unavailable; the authors' archived full-W Tc there, an equivalent mu* and an N_F
grade the submission.
"""
import json
import math
from pathlib import Path
import numpy as np
import eliashberg as E

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
COSTS = DATA['costs']
SOURCE = 'retrieved author-released data (Zenodo 23235290, Bozier & Pickard 2026); no new calculation'
KEYS = {'inputs': ('system', 'pressure_GPa'), 'coulomb_table': ('system', 'pressure_GPa'),
        'w_kernel': ('system', 'pressure_GPa'), 'archived_run': ('system', 'pressure_GPa', 'run')}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def entry(request, need):
    system, p = request.get('system'), request.get('pressure_GPa')
    if not isinstance(system, str) or not finite(p):
        raise ValueError('system must be a string and pressure_GPa a number (see /app/catalog.json)')
    key = f'{system}/{int(p)}' if float(p).is_integer() else None
    if key not in DATA['entries']:
        raise ValueError(f'no archived calculation for {system} at {p} GPa (see /app/catalog.json)')
    rec = DATA['entries'][key]
    if need not in rec:
        raise ValueError(f'{need} is not available for {system} at {int(p)} GPa'
                         + (' (withheld: this pressure is graded)' if rec.get('withheld') else ''))
    return key, rec


def mu_list(request):
    m = request.get('mu_star')
    values = m if isinstance(m, list) else [m]
    if not 1 <= len(values) <= 6 or not all(finite(v) and 0 <= v <= 3 for v in values):
        raise ValueError('mu_star must be a number or a list of 1-6 numbers in [0, 3]')
    return [float(v) for v in values]


def validate(request):
    kind = request.get('kind')
    if kind in KEYS:
        if set(request) != {'kind', *KEYS[kind]}:
            raise ValueError('Use {"kind": "%s", %s}' % (kind, ', '.join(f'"{k}": ...' for k in KEYS[kind])))
        need = {'inputs': 'inputs', 'coulomb_table': 'coulomb', 'w_kernel': 'w_kernel', 'archived_run': 'runs'}[kind]
        _, rec = entry(request, need)
        if kind == 'archived_run' and request['run'] not in rec['runs']:
            raise ValueError('run must be one of the archived_runs listed for this system and pressure')
        return COSTS[kind]
    if kind == 'eliashberg':
        solver = request.get('solver')
        allowed = {'kind', 'system', 'pressure_GPa', 'solver', 'omega_c_meV'} | ({'mu_star'} if solver in ('cdos', 'vdos') else set())
        if solver not in ('cdos', 'vdos', 'vdos_W'):
            raise ValueError('solver must be "cdos", "vdos" or "vdos_W"')
        if not set(request) <= allowed or not {'kind', 'system', 'pressure_GPa', 'solver'} <= set(request) \
                or (solver != 'vdos_W' and 'mu_star' not in request):
            raise ValueError('Use {"kind": "eliashberg", "system": ..., "pressure_GPa": ..., "solver": "cdos"|"vdos", '
                             '"mu_star": x or [x1, ...], "omega_c_meV": 7000} or {..., "solver": "vdos_W", "omega_c_meV": 7000}')
        entry(request, 'solver_W' if solver == 'vdos_W' else 'eliashberg')
        if solver != 'vdos_W':
            mu_list(request)
        wc = request.get('omega_c_meV', 7000)
        if not finite(wc) or not 300 <= wc <= 7000:
            raise ValueError('omega_c_meV must be a number in [300, 7000]')
        return COSTS[kind]
    raise ValueError('kind must be one of: ' + ', '.join(COSTS))


def spectrum(rec):
    a = rec['inputs']['alpha2F']
    return E.A2F(np.array(a['omega_THz']) * E.THZ2MEV, a['columns'][str(a['solver_column'])])


def electrons(rec, with_W=False):
    d = rec['inputs']['dos_file']
    e = (np.array(d['E_eV']) - d['EF_eV']) * 1000
    n = np.array(d['dos']) / 2 / 1000          # spin-summed states/eV -> per spin per meV
    if not with_W:
        return E.Electrons(e, n)
    w = rec['solver_W']
    return E.Electrons(e, n, np.array(w['E_eV']) * 1000, np.array(w['W_eV']) * 1000)


def tc_guess(spec, mustar):
    return max(E.allen_dynes(spec.lam, spec.wlog, spec.w2, min(mustar, 0.3)), 20.0)


def solve(rec, solver, mu_stars, wc):
    spec = spectrum(rec)
    base = {'lambda': round(spec.lam, 4), 'omega_log_meV': round(spec.wlog, 3), 'omega_2_meV': round(spec.w2, 3),
            'omega_c_meV': wc, 'alpha2F_column': rec['inputs']['alpha2F']['solver_column']}
    rows = []
    if solver == 'cdos':
        for m in mu_stars:
            tc = E.critical_temperature(lambda T: E.eig_cdos(spec, m, T, wc), tc_guess(spec, m))
            info = E.eig_cdos(spec, m, tc if tc else E.T_FLOOR, wc)[1]
            rows.append({'mu_star': m, 'Tc_K': None if tc is None else round(tc, 2), 'Z0_at_Tc': round(info['Z0'], 4),
                         'matsubara_frequencies': info['n_matsubara']})
        return {**base, 'solver': 'constant DOS, fixed mu*', 'results': rows}
    el = electrons(rec, with_W=(solver == 'vdos_W'))
    base.update({'N_EF_states_per_eV_per_spin': round(el.NF * 1000, 5), 'energy_window_meV': [float(el.e[0]), float(el.e[-1])],
                 'energy_grid_points': len(el.e)})
    if solver == 'vdos':
        for m in mu_stars:
            tc = E.critical_temperature(lambda T: E.eig_vdos(spec, el, m, T, wc), tc_guess(spec, m))
            info = E.eig_vdos(spec, el, m, tc if tc else E.T_FLOOR, wc)[1]
            rows.append({'mu_star': m, 'Tc_K': None if tc is None else round(tc, 2), 'Z0_at_Tc': round(info['Z0'], 4),
                         'chi0_at_Tc_meV': round(info['chi0_meV'], 3), 'chemical_potential_shift_meV': round(info['mu_shift_meV'], 3),
                         'matsubara_frequencies': info['n_matsubara']})
        return {**base, 'solver': 'energy-dependent DOS (full bandwidth), fixed mu*', 'results': rows}
    tc = E.critical_temperature(lambda T: E.eig_vdos_W(spec, el, T, wc), tc_guess(spec, min(el.mu, 0.3)))
    info = E.eig_vdos_W(spec, el, tc if tc else E.T_FLOOR, wc)[1]
    rows.append({'Tc_K': None if tc is None else round(tc, 2), 'mu_from_W': round(el.mu, 5),
                 'W_EF_EF_eV': round(el.mu / el.NF / 1000, 5), 'Z0_at_Tc': round(info['Z0'], 4),
                 'chi0_at_Tc_meV': round(info['chi0_meV'], 3), 'chemical_potential_shift_meV': round(info['mu_shift_meV'], 3),
                 'matsubara_frequencies': info['n_matsubara']})
    return {**base, 'solver': 'energy-dependent DOS (full bandwidth), static screened W(e,e\')', 'results': rows}


def run(request, workdir):
    kind = request['kind']
    if kind == 'eliashberg':
        _, rec = entry(request, 'inputs')
        mus = mu_list(request) if request['solver'] != 'vdos_W' else []
        out = solve(rec, request['solver'], mus, float(request.get('omega_c_meV', 7000)))
        for r in out['results']:
            if r['Tc_K'] is None:
                r['note'] = f'Tc below {E.T_FLOOR} K'
        return {'request': request, **out, 'source': 'new numpy Eliashberg solution on archived inputs'}
    key, rec = entry(request, {'inputs': 'inputs', 'coulomb_table': 'coulomb', 'w_kernel': 'w_kernel',
                               'archived_run': 'runs'}[kind])
    if kind == 'archived_run':
        return {'system': request['system'], 'pressure_GPa': request['pressure_GPa'], 'run': request['run'],
                'folder': rec['folder'] + '/IsoME/' + request['run'], **rec['runs'][request['run']], 'source': SOURCE}
    body = {'inputs': 'inputs', 'coulomb_table': 'coulomb', 'w_kernel': 'w_kernel'}[kind]
    return {'system': request['system'], 'pressure_GPa': request['pressure_GPa'], 'folder': rec['folder'],
            **rec[body], 'source': SOURCE}


def grade(payload, trace, validation):
    ref, tol = PRIVATE['reference'], PRIVATE['tolerance']
    fields = ['tc_full_W_K', 'mu_star_cdos_1eV_Mg2IrH6_20GPa', 'N_EF_H3S_180GPa_per_fu_per_spin', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    tc = payload['tc_full_W_K']
    if not isinstance(tc, dict) or set(tc) != set(ref['tc_full_W_K']) or not all(finite(v) for v in tc.values()):
        return False, {'error': 'tc_full_W_K must give a finite Tc for each of: ' + ', '.join(sorted(ref['tc_full_W_K']))}
    if not finite(payload['mu_star_cdos_1eV_Mg2IrH6_20GPa']) or not finite(payload['N_EF_H3S_180GPa_per_fu_per_spin']):
        return False, {'error': 'mu_star_cdos_1eV_Mg2IrH6_20GPa and N_EF_H3S_180GPa_per_fu_per_spin must be finite numbers'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    rel = {k: tc[k] / v - 1 for k, v in ref['tc_full_W_K'].items()}
    dmu = payload['mu_star_cdos_1eV_Mg2IrH6_20GPa'] - ref['mu_star_cdos_1eV_Mg2IrH6_20GPa']
    nrel = payload['N_EF_H3S_180GPa_per_fu_per_spin'] / ref['N_EF_H3S_180GPa_per_fu_per_spin'] - 1
    checks = {'evidence_executed': evidence}
    for k, r in rel.items():
        checks[f'tc_{k}'] = abs(r) <= tol['tc_rel']
    checks['mu_star_equivalent'] = abs(dmu) <= tol['mu_star_abs']
    checks['N_EF'] = abs(nrel) <= tol['nef_rel']
    report = {'checks': checks, 'tc_relative_errors': rel, 'mu_star_error': dmu, 'N_EF_relative_error': nrel}
    return all(checks.values()), report
