"""Lab logic for mos2-doping-dome-cdw: author-released MoS2 data and an isotropic Eliashberg solver.

Retrieval requests return archived data of Girotto Erhardt et al. (Materials Cloud 51yyy-ksf72). The
`eliashberg` and `eliashberg_custom` requests run eliashberg.py, which reproduces the authors' archived
Eliashberg tables. The 1x1 H spectra, DOS and Tc of the 0.15 and 0.17 calculations are never served;
their archived Tc, the critical doping from the PES, the SOC Tc at that doping and the high-doping CDW Tc
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
SOURCE = 'retrieved author-released data (Materials Cloud 51yyy-ksf72); no new calculation'
WITHHELD = {'0.15', '0.17'}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def keys(request, allowed, usage):
    if set(request) != set(allowed):
        raise ValueError('Use ' + usage)


def spectrum_key(request):
    key = f"{request.get('phase')}/{request.get('label')}/{request.get('flavor')}"
    if request.get('phase') == '1x1H' and request.get('label') in WITHHELD:
        raise ValueError(f"1x1H {request.get('label')} is withheld for grading")
    if key not in DATA['spectra']:
        raise ValueError('Not available: see spectra_available in /app/catalog.json (phase, label, flavor A/NA/FA)')
    return key


def dos_key(spec):
    if not isinstance(spec, dict) or set(spec) != {'phase', 'label'}:
        raise ValueError('dos must be "constant" or {"phase": ..., "label": ...} from dos_available')
    if spec['phase'] == '1x1H' and spec['label'] in WITHHELD:
        raise ValueError(f"the 1x1H {spec['label']} DOS is withheld for grading")
    key = f"{spec['phase']}/{spec['label']}"
    if key not in DATA['dos']:
        raise ValueError('Not available: see dos_available in /app/catalog.json')
    return key


def mu_star(request):
    m = request.get('mu_star')
    if not finite(m) or not 0 <= m <= 0.3:
        raise ValueError('mu_star must be a number in [0, 0.3]')
    return float(m)


def custom_spectrum(request):
    rows = request.get('spectrum_meV')
    if not isinstance(rows, list) or not 3 <= len(rows) <= 600:
        raise ValueError('spectrum_meV must be a list of 3-600 [omega_meV, alpha2F] pairs')
    for r in rows:
        if not isinstance(r, list) or len(r) != 2 or not all(finite(v) for v in r):
            raise ValueError('each spectrum_meV row must be [omega_meV, alpha2F] with finite numbers')
    w = [r[0] for r in rows]
    if not all(0 < a < b <= 100 for a, b in zip(w, w[1:])) or w[0] <= 0 or any(r[1] < 0 for r in rows):
        raise ValueError('omega must increase strictly within (0, 100] meV and alpha2F must be >= 0')
    s = E.Spectrum(w, [r[1] for r in rows])
    if s.lam > 10:
        raise ValueError('lambda of the supplied spectrum exceeds 10')
    return s


def validate(request):
    kind = request.get('kind')
    if kind == 'spectrum':
        keys(request, ('kind', 'phase', 'label', 'flavor'), '{"kind": "spectrum", "phase": ..., "label": ..., "flavor": "A"|"NA"|"FA"}')
        spectrum_key(request)
        return COSTS[kind]
    if kind == 'dos':
        keys(request, ('kind', 'phase', 'label'), '{"kind": "dos", "phase": "1x1H"|"1x1H-SOC", "label": ...}')
        dos_key({'phase': request['phase'], 'label': request['label']})
        return COSTS[kind]
    if kind == 'eliashberg':
        keys(request, ('kind', 'phase', 'label', 'flavor', 'dos', 'mu_star'),
             '{"kind": "eliashberg", "phase": ..., "label": ..., "flavor": ..., "dos": "full"|"constant", "mu_star": 0.13}')
        spectrum_key(request)
        mu_star(request)
        if request['dos'] not in ('full', 'constant'):
            raise ValueError('dos must be "full" or "constant"')
        if request['dos'] == 'full' and f"{request['phase']}/{request['label']}" not in DATA['dos']:
            raise ValueError('no DOS is archived for this phase; use "dos": "constant"')
        return COSTS[kind]
    if kind == 'eliashberg_custom':
        if request.get('dos') == 'constant':
            keys(request, ('kind', 'spectrum_meV', 'dos', 'mu_star'),
                 '{"kind": "eliashberg_custom", "spectrum_meV": [[w, a2F], ...], "dos": "constant", "mu_star": 0.13}')
        else:
            keys(request, ('kind', 'spectrum_meV', 'dos', 'electrons_per_cell', 'mu_star'),
                 '{"kind": "eliashberg_custom", "spectrum_meV": [[w, a2F], ...], "dos": {"phase": ..., "label": ...}, '
                 '"electrons_per_cell": x, "mu_star": 0.13}')
            dos_key(request['dos'])
            n = request['electrons_per_cell']
            if not finite(n) or not 0.02 <= n <= 0.6:
                raise ValueError('electrons_per_cell must be a number in [0.02, 0.6]')
        mu_star(request)
        custom_spectrum(request)
        return COSTS[kind]
    if kind in ('moments_table', 'a1g_raman', 'experiment_dome'):
        keys(request, ('kind',), '{"kind": "%s"}' % kind)
        return COSTS[kind]
    if kind == 'pes':
        keys(request, ('kind', 'label'), '{"kind": "pes", "label": <one of pes_available>}')
        if request['label'] not in DATA['pes']:
            raise ValueError('label must be one of pes_available in /app/catalog.json')
        return COSTS[kind]
    if kind in ('soft_mode', 'mo_displacement'):
        keys(request, ('kind', 'branch'), '{"kind": "%s", "branch": "1x1H"|"2x2CDW"}' % kind)
        if request['branch'] not in DATA[kind]:
            raise ValueError('branch must be "1x1H" or "2x2CDW"')
        return COSTS[kind]
    if kind == 'archived_tc_table':
        keys(request, ('kind', 'phase'), '{"kind": "archived_tc_table", "phase": "1x1H"|"2x2CDW"|"1Tprime"}')
        if request['phase'] not in DATA['archived_tc_tables']:
            raise ValueError('phase must be "1x1H", "2x2CDW" or "1Tprime"')
        return COSTS[kind]
    raise ValueError('kind must be one of: ' + ', '.join(COSTS))


def solve(spec, mu, dos_entry=None, n=None):
    dos = None if dos_entry is None else E.DOS(dos_entry['energy_eV'], dos_entry['total'], n)
    tc, diag = E.critical_temperature(spec, mu, dos)
    out = {'Tc_K': tc, 'lambda': spec.lam, 'omega_log_meV': spec.wlog * 1e3, 'omega_2_meV': spec.w2 * 1e3,
           'mu_star': mu, 'mu_star_at_matsubara_cutoff': diag['mu_star_C'], 'matsubara_frequencies': diag['n_matsubara'],
           'dos_treatment': 'constant' if dos is None else 'full'}
    if tc is None:
        out['Tc_K'] = None
        out['note'] = 'Tc below %.1f K (the largest eigenvalue at 1 K is below 1)' % E.T_FLOOR
    if dos is not None:
        out.update(electrons_per_cell=n, chemical_potential_noninteracting_meV=diag['mu0_eV'] * 1e3,
                   chemical_potential_meV=diag['mu_eV'] * 1e3, dos_at_mu0_per_eV_spin=diag['N_at_mu0'], Z_first_matsubara=diag['Z0'],
                   energy_zero='DOS energy zero = Fermi level of the doped cell at the 800 K smearing')
    return out


def run(request, workdir):
    kind = request['kind']
    if kind == 'spectrum':
        key = spectrum_key(request)
        s = DATA['spectra'][key]
        return {'phase': request['phase'], 'label': request['label'], 'flavor': request['flavor'],
                'source_file': s['source_file'], 'header': s['header'],
                'columns': ['omega_meV', 'alpha2F', 'lambda_cumulative'], 'rows': s['rows'], 'source': SOURCE}
    if kind == 'dos':
        key = dos_key({'phase': request['phase'], 'label': request['label']})
        return {'phase': request['phase'], 'label': request['label'], **DATA['dos'][key],
                'units': 'states per eV, per spin and 1x1 cell; energy in eV relative to the 800 K Fermi level', 'source': SOURCE}
    if kind == 'eliashberg':
        key = spectrum_key(request)
        rows = np.array(DATA['spectra'][key]['rows'])
        spec = E.Spectrum(rows[:, 0], rows[:, 1])
        dkey = f"{request['phase']}/{request['label']}"
        out = solve(spec, mu_star(request), DATA['dos'][dkey] if request['dos'] == 'full' else None,
                    float(request['label']) if request['dos'] == 'full' else None)
        return {'phase': request['phase'], 'label': request['label'], 'flavor': request['flavor'], **out,
                'source': 'eliashberg.py on the archived spectrum (new calculation)'}
    if kind == 'eliashberg_custom':
        spec = custom_spectrum(request)
        if request['dos'] == 'constant':
            out = solve(spec, mu_star(request))
        else:
            out = solve(spec, mu_star(request), DATA['dos'][dos_key(request['dos'])], float(request['electrons_per_cell']))
            out['dos_source'] = request['dos']
        return {**out, 'source': 'eliashberg.py on the supplied spectrum (new calculation)'}
    if kind == 'moments_table':
        return {**DATA['moments'], 'note': 'rows: tot_charge per 1x1 cell, then adiabatic, nonadiabatic, fully nonadiabatic',
                'source': SOURCE}
    if kind == 'pes':
        return {'label': request['label'], **DATA['pes'][request['label']],
                'path': '11 equally spaced points from the 1x1 H structure (first) to the 2x2 CDW structure (last), 2x2 supercell',
                'source': SOURCE}
    if kind in ('soft_mode', 'mo_displacement'):
        units = ('n (10^14 cm^-2), frequency (meV; negative = imaginary)' if kind == 'soft_mode'
                 else 'n (10^14 cm^-2), Mo displacement (A)')
        return {'branch': request['branch'], **DATA[kind][request['branch']], 'columns': units, 'source': SOURCE}
    if kind == 'archived_tc_table':
        return {'phase': request['phase'], **DATA['archived_tc_tables'][request['phase']], 'source': SOURCE}
    if kind == 'a1g_raman':
        return {**DATA['a1g'], 'columns': {'theory': 'tot_charge per 1x1 cell, adiabatic, nonadiabatic A1g (meV)',
                                           'experiment': 'Li per MoS2 in LixMoS2, A1g (cm^-1)'}, 'source': SOURCE}
    return {**DATA['experiment_dome'], 'columns': 'carrier density (10^14 cm^-2), Tc (K)', 'source': SOURCE}


def grade(payload, trace, validation):
    tg, th = PRIVATE['targets'], PRIVATE['thresholds']
    fields = ['tc_1x1_K', 'critical_doping_e_per_cell', 'tc_1x1_soc_at_critical_doping_K', 'tc_cdw_high_K', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    t1, tcdw = payload['tc_1x1_K'], payload['tc_cdw_high_K']
    if not isinstance(t1, dict) or set(t1) != set(tg['tc_1x1_K']) or not all(finite(v) and v > 0 for v in t1.values()):
        return False, {'error': 'tc_1x1_K must give a positive Tc for "0.15" and "0.17"'}
    if not isinstance(tcdw, dict) or set(tcdw) != set(tg['tc_cdw_high_K']) or not all(finite(v) and v > 0 for v in tcdw.values()):
        return False, {'error': 'tc_cdw_high_K must give a positive Tc for "0.35" and "0.40"'}
    for f in ('critical_doping_e_per_cell', 'tc_1x1_soc_at_critical_doping_K'):
        if not finite(payload[f]):
            return False, {'error': f + ' must be a finite number'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    errors = {'tc_1x1_rel': {x: t1[x] / v - 1 for x, v in tg['tc_1x1_K'].items()},
              'critical_doping_abs': payload['critical_doping_e_per_cell'] - tg['critical_doping_e_per_cell'],
              'tc_soc_rel': payload['tc_1x1_soc_at_critical_doping_K'] / tg['tc_1x1_soc_at_critical_doping_K'] - 1,
              'tc_cdw_rel': {x: tcdw[x] / v - 1 for x, v in tg['tc_cdw_high_K'].items()}}
    checks = {'evidence_executed': evidence,
              **{f'tc_1x1_{x}': abs(e) <= th['tc_1x1_rel'] for x, e in errors['tc_1x1_rel'].items()},
              'critical_doping': abs(errors['critical_doping_abs']) <= th['critical_doping_abs'],
              'tc_1x1_soc_at_critical_doping': abs(errors['tc_soc_rel']) <= th['tc_soc_rel'],
              **{f'tc_cdw_high_{x}': abs(e) <= th['tc_cdw_rel'] for x, e in errors['tc_cdw_rel'].items()}}
    return all(checks.values()), {'checks': checks, 'errors': errors}
