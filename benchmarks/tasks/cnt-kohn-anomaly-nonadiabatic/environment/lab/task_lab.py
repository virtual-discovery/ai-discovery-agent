"""Lab logic for cnt-kohn-anomaly-nonadiabatic: archived electron-phonon data of the (3,3) nanotube,
a numpy phonon self-energy evaluator, and grading.

Purchases return author-archived files and tables (archive_file, bands, phonons, coupling, dfpt_scan,
self_energy_table at the released temperatures) or a new self-energy evaluation from the archived
coupling tables (self_energy). The archived self-energy tables at the other temperatures are never
served; reference values derived from them (and from the archived threshold model) grade the submission.
"""
import json
import math
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import selfenergy as SE  # noqa: E402

PUBLIC = json.loads((HERE/'public.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
ARCHIVE = json.loads((HERE/'data/archive.json').read_text())
TABLES = json.loads((HERE/'data/tables.json').read_text())
RELEASED = {float(k): k for k in TABLES}
MODES = (1, 20, 23)
SRC = 'retrieved author-archived data (no new calculation)'
_Z = None
_COUPLING = {}


def coupling_arrays():
    global _Z
    if _Z is None:
        with np.load(HERE/'data/coupling.npz') as z:
            _Z = {k: z[k] for k in z.files}
    return _Z


def coupling(nu, process):
    key = (nu, process)
    if key not in _COUPLING:
        z, tag = coupling_arrays(), SE.PROCESSES[process].replace('+', 'p').replace('-', 'm')
        _COUPLING[key] = SE.Coupling(z['k'], z[f'e_{nu}_{tag}'], z[f'eq_{nu}_{tag}'], z[f'g_{nu}_{tag}'], z[f'wq_{nu}_{tag}'])
    return _COUPLING[key]


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def check_mode(request):
    if request.get('mode') not in MODES or isinstance(request.get('mode'), bool):
        raise ValueError('mode must be 1, 20 or 23')


def validate(request):
    kind = request.get('kind')
    keys = set(request)
    if kind == 'archive_file':
        if keys != {'kind', 'file'} or request['file'] not in ARCHIVE['files']:
            raise ValueError('Use {"kind": "archive_file", "file": <one of %s>}' % sorted(ARCHIVE['files']))
        return 1
    if kind in ('bands', 'phonons', 'dfpt_scan'):
        if keys != {'kind'}:
            raise ValueError('Use {"kind": "%s"}' % kind)
        return 1
    if kind == 'coupling':
        if keys != {'kind', 'mode'}:
            raise ValueError('Use {"kind": "coupling", "mode": 1 | 20 | 23}')
        check_mode(request)
        return 1
    if kind == 'self_energy_table':
        if keys != {'kind', 'mode', 'temperature_K'}:
            raise ValueError('Use {"kind": "self_energy_table", "mode": 1 | 20 | 23, "temperature_K": <released temperature>}')
        check_mode(request)
        t = request['temperature_K']
        if not number(t) or not any(abs(t - r) < 1e-6 for r in RELEASED):
            raise ValueError('temperature_K must be one of the released temperatures %s' % sorted(RELEASED))
        return 1
    if kind == 'self_energy':
        need = {'kind', 'mode', 'temperature_K', 'omega_meV', 'nk', 'eta_meV', 'processes', 'band_pairs'}
        if keys != need:
            raise ValueError('self_energy needs exactly the fields %s' % sorted(need))
        check_mode(request)
        t, w, nk, eta = request['temperature_K'], request['omega_meV'], request['nk'], request['eta_meV']
        if not number(t) or not 30 <= t <= 2000:
            raise ValueError('temperature_K must be a number between 30 and 2000')
        if not isinstance(w, list) or not 1 <= len(w) <= 41 or not all(number(x) and -500 <= x <= 500 for x in w):
            raise ValueError('omega_meV must be a list of 1-41 numbers in [-500, 500]')
        if type(nk) is not int or not 2000 <= nk <= 4000000:
            raise ValueError('nk must be an integer between 2000 and 4000000')
        if not number(eta) or not 0 < eta <= 50:
            raise ValueError('eta_meV must be a number in (0, 50]')
        p = request['processes']
        if not isinstance(p, list) or not p or len(set(map(str, p))) != len(p) or not all(x in SE.PROCESSES for x in p):
            raise ValueError('processes must be a non-empty list of distinct entries from ["k+q", "k-q"]')
        bp = request['band_pairs']
        ok = isinstance(bp, list) and bp and all(isinstance(x, list) and len(x) == 2 and all(type(v) is int and v in (1, 2) for v in x)
                                                  for x in bp)
        if not ok or len({tuple(x) for x in bp}) != len(bp):
            raise ValueError('band_pairs must be a non-empty list of distinct [n, m] with n, m in {1, 2}')
        return 1 if nk <= 100000 else (2 if nk <= 1000000 else 3)
    raise ValueError('kind must be one of %s' % sorted(PUBLIC['requests']))


def run(request, workdir):
    kind = request['kind']
    if kind == 'archive_file':
        return {'file': request['file'], 'text': ARCHIVE['files'][request['file']], 'source': SRC}
    if kind == 'bands':
        return {'file': 'Fig2/e_band/band.eig', 'k_units': '2 pi / alat along z (alat = 13.0 A)', 'energy_units': 'eV',
                'k': ARCHIVE['bands']['k_alat_2pi'], 'bands': ARCHIVE['bands']['energies_eV'],
                'band_plot_script_line': ARCHIVE['fermi_shift_comment'] + '   (fig2_a_b.py; the plot shows e - FERMI_SHIFT)',
                'source': SRC}
    if kind == 'phonons':
        return {'file': 'Fig2/ph_band/matdyn.freq', 'q_units': '2 pi / alat along z (alat = 13.0 A)',
                'frequency_units': 'cm-1', 'q': ARCHIVE['phonons']['q_alat_2pi'],
                'branches': ARCHIVE['phonons']['frequencies_cm1'], 'source': SRC}
    if kind == 'dfpt_scan':
        s = ARCHIVE['dfpt_scan']
        return {'description': 'DFPT phonon frequencies at Q with Fermi-Dirac smearing sigma (authors\' Fig. 2 script, "Eq. (1)")',
                'sigma_mRy': s['sigma_mRy'], 'frequency_cm1': {str(nu): row for nu, row in zip(s['modes'], s['frequency_cm1'])},
                'source': SRC}
    if kind == 'coupling':
        nu, z = request['mode'], coupling_arrays()
        out = {'mode': nu, 'k_crystal': z['k'].round(6).tolist(), 'omega_Q_meV': float(z[f'wq_{nu}_p']),
               'columns': 'e_k[n] (eV) for n = 1, 2; e_kq[m] (eV) for m = 1, 2; g[n][m] = |g_mn(k)| (meV)',
               'authors_note': ARCHIVE['k_plus_minus_q_comment'] + '   (comment in fig2_c_e_g.py, the coupling-plot script)',
               'source': SRC}
        for proc, tag in (('k+q', 'p'), ('k-q', 'm')):
            out[proc] = {'files': f'Fig2/g/g{"+" if tag == "p" else "-"}q.<n><m>{nu}',
                         'e_k': z[f'e_{nu}_{tag}'].round(4).tolist(), 'e_kq': z[f'eq_{nu}_{tag}'].round(4).tolist(),
                         'g': z[f'g_{nu}_{tag}'].tolist()}
        return out
    if kind == 'self_energy_table':
        name = RELEASED[min(RELEASED, key=lambda r: abs(r - request['temperature_K']))]
        table = TABLES[name]
        return {'file': f'Fig2/PiT/{name}', 'mode': request['mode'], 'lines': table['header'] + table[str(request['mode'])],
                'note': 'verbatim header and the data lines of this mode (one line per frequency)', 'source': SRC}
    nu = request['mode']
    procs = [coupling(nu, p) for p in request['processes']]
    re, im = SE.self_energy(procs, float(request['temperature_K']), request['omega_meV'], request['nk'],
                            float(request['eta_meV']), [tuple(x) for x in request['band_pairs']])
    return {'mode': nu, 'temperature_K': request['temperature_K'], 'nk': request['nk'], 'eta_meV': request['eta_meV'],
            'processes': request['processes'], 'band_pairs': request['band_pairs'],
            'omega_Q_meV': float(coupling_arrays()[f'wq_{nu}_p']), 'chemical_potential_eV': SE.MU_EV,
            'omega_meV': request['omega_meV'], 're_pi_meV': np.round(re, 6).tolist(), 'im_pi_meV': np.round(im, 6).tolist(),
            'source': 'new evaluation from the archived EPW coupling tables (numpy evaluator)'}


def get(payload, *path):
    for key in path:
        if not isinstance(payload, dict) or key not in payload:
            return None
        payload = payload[key]
    return payload


def grade(payload, trace, validation):
    ref, t = PRIVATE['reference'], PRIVATE['thresholds']
    paths = [('omega_nonadiabatic_meV', '100'), ('omega_nonadiabatic_meV', '300'), ('linewidth_meV', '100'),
             ('linewidth_meV', '300'), ('omega_adiabatic_100K_meV',), ('threshold_coupling_meV', '1'),
             ('threshold_coupling_meV', '20')]
    values = {'/'.join(p): get(payload, *p) for p in paths}
    bad = [k for k, v in values.items() if not number(v)]
    tstar = payload.get('rbm_instability_temperature_K', 'missing') if isinstance(payload, dict) else 'missing'
    if tstar == 'missing' or not (tstar is None or number(tstar)):
        bad.append('rbm_instability_temperature_K')
    if bad or 'evidence_ids' not in payload:
        return False, {'error': 'missing or non-finite fields', 'fields': bad + ([] if 'evidence_ids' in payload else ['evidence_ids'])}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    err = {}
    for key in ('100', '300'):
        err[f'omega_nonadiabatic_{key}K'] = values[f'omega_nonadiabatic_meV/{key}'] - ref['omega_nonadiabatic_meV'][key]
        err[f'linewidth_{key}K'] = values[f'linewidth_meV/{key}'] / ref['linewidth_meV'][key] - 1
    err['omega_adiabatic_100K'] = values['omega_adiabatic_100K_meV'] - ref['omega_adiabatic_100K_meV']
    for nu in ('1', '20'):
        err[f'threshold_coupling_{nu}'] = values[f'threshold_coupling_meV/{nu}'] / ref['threshold_coupling_meV'][nu] - 1
    err['rbm_instability_temperature'] = (None if tstar is None else
                                          tstar / ref['rbm_instability_temperature_K'] - 1)
    checks = {'evidence_executed': evidence,
              'omega_nonadiabatic_100K': abs(err['omega_nonadiabatic_100K']) <= t['omega_nonadiabatic_abs_meV'],
              'omega_nonadiabatic_300K': abs(err['omega_nonadiabatic_300K']) <= t['omega_nonadiabatic_abs_meV'],
              'linewidth_100K': abs(err['linewidth_100K']) <= t['linewidth_rel'],
              'linewidth_300K': abs(err['linewidth_300K']) <= t['linewidth_rel'],
              'omega_adiabatic_100K': abs(err['omega_adiabatic_100K']) <= t['omega_adiabatic_abs_meV'],
              'rbm_instability_temperature': (err['rbm_instability_temperature'] is not None
                                              and abs(err['rbm_instability_temperature']) <= t['instability_temperature_rel']),
              'threshold_coupling_RBM': abs(err['threshold_coupling_1']) <= t['threshold_coupling_rel'],
              'threshold_coupling_nu20': abs(err['threshold_coupling_20']) <= t['threshold_coupling_rel']}
    report = {'checks': checks, 'errors': {k: (None if v is None else round(v, 4)) for k, v in err.items()}}
    return all(checks.values()), report
