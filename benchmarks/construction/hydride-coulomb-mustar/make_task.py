"""Build the staged, private and ablation files of hydride-coulomb-mustar from DATA/pool.json.gz.

    bash run_in_image.sh make_task.py        (about 15 min on 4 cores; solver results cached in DATA/tc_cache.json)

Writes, in benchmarks/tasks/hydride-coulomb-mustar:
  environment/agent/catalog.json, environment/lab/data.json, environment/lab/eliashberg.py (copied),
  environment/lab/private/hidden.json, private/hidden_values.json, solution/ablation_values.json,
  the budget in both public.json files; and calibration.json here.
All Tc values below come from the lab's own code path (task_lab.solve / eliashberg.py).
"""
import copy
import gzip
import importlib.util
import json
import math
import os
import re
import shutil
import sys
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
TASK = os.path.join(HERE, '../../tasks/hydride-coulomb-mustar')
LAB = os.path.join(TASK, 'environment/lab')
DATA_DIR = '/data'
POOL = json.load(gzip.open(os.path.join(DATA_DIR, 'pool.json.gz'), 'rt'))
BUDGET = 38
COSTS = {'inputs': 1, 'coulomb_table': 1, 'w_kernel': 1, 'archived_run': 1, 'eliashberg': 2}
WITHHELD = {'Mg2IrH6/15', 'Mg2IrH6/25', 'H3S/180', 'H3S/190'}
TARGETS = {'Mg2IrH6_15GPa': 'Mg2IrH6/15', 'Mg2IrH6_25GPa': 'Mg2IrH6/25', 'H3S_180GPa': 'H3S/180', 'H3S_190GPa': 'H3S/190'}
WINDOW_EV = 5.6
TOL = {'tc_rel': 0.035, 'mu_star_abs': 0.025, 'nef_rel': 0.03}


def info_field(info, name):
    m = re.search(r'^%s: (.*)$' % re.escape(name), info, re.M)
    return m.group(1).strip() if m else None


def summary_tc(text):
    m = re.search(r'Tc = ([\d.]+)', text)
    return float(m.group(1))


def w_matrix(w):
    dos = np.array(w['dos'])
    nnw = np.array(w['NNW'])
    with np.errstate(divide='ignore', invalid='ignore'):
        W = nnw / np.outer(dos, dos)
    W[~np.isfinite(W)] = 0.0
    W[np.outer(dos <= 1e-10, np.ones_like(dos, bool)) | np.outer(np.ones_like(dos, bool), dos <= 1e-10)] = 0.0
    return np.array(w['E_eV']), W


def solver_column(rec):
    for run in sorted(rec.get('runs', {})):
        info = rec['runs'][run]['Info.txt']
        if info_field(info, 'include_Weep') == '1':
            return int(info_field(info, 'ind_smear'))
    raise RuntimeError('no W run')


# ---------------------------------------------------------------- staged data
entries, catalog_systems = {}, {}
for key, rec in POOL.items():
    if key.startswith('_'):
        continue
    system, p = key.split('/')
    withheld = key in WITHHELD
    e = {'folder': rec['folder'], 'withheld': withheld}
    inputs = {'cell_scf_in': rec['scf_in']}
    if rec.get('relax_in'):
        m = re.search(r'press\s*=\s*([-\d.eEdD+]+)', rec['relax_in'])
        inputs['relax_target_press_kbar'] = float(m.group(1).replace('d', 'e').replace('D', 'e')) if m else None
    d = rec['dos_file']
    inputs['dos_file'] = {'header': d['header'], 'EF_eV': d['EF_eV'], 'E_eV': d['E_eV'], 'dos': d['dos'],
                          'int_dos': d['int_dos'],
                          'columns': 'E (eV, absolute), dos(E) (states/eV, as written by QE dos.x), integrated dos'}
    if 'a2f' in rec:
        col = solver_column(rec)
        table = np.array(rec['a2f']['table'])
        cols = sorted({1, col, 10, 20} & set(range(1, table.shape[1])))
        smear = rec['a2f']['header'].split()[2:]
        inputs['alpha2F'] = {'file': 'IsoME/rescaled-alpha2F.dat', 'header': rec['a2f']['header'],
                             'omega_THz': table[:, 0].tolist(),
                             'columns': {str(c): table[:, c].tolist() for c in cols},
                             'column_smearings': {str(c): smear[c - 1] for c in cols},
                             'solver_column': col,
                             'note': 'subset of the 50 phonon-smearing columns; the archived IsoME runs and the lab solver use solver_column'}
        e['eliashberg'] = True
    e['inputs'] = inputs
    if not withheld:
        if 'mu_results' in rec:
            e['coulomb'] = {'mu_results.dat': rec['mu_results'], 'mantel_log_header': rec.get('mantel_log_header')}
        if 'W' in rec:
            w = rec['W']
            Egrid = np.array(w['E_eV'])
            if 'NNW' in w:
                Eg, W = w_matrix(w)
                win = np.abs(Eg) <= WINDOW_EV
                kern = {'E_eV': Eg.tolist(), 'mantel_dos': w['dos'], 'W_diagonal_eV': np.diag(W).tolist(),
                        'window_E_eV': Eg[win].tolist(), 'W_window_eV': np.round(W[np.ix_(win, win)], 7).tolist()}
                if 'eliashberg' in e:
                    e['solver_W'] = {'E_eV': Eg[win].tolist(), 'W_eV': W[np.ix_(win, win)].tolist()}
            else:
                dos = np.array(w['dos'])
                with np.errstate(divide='ignore', invalid='ignore'):
                    diag = np.where(dos > 1e-10, np.array(w['NNW_diagonal']) / dos ** 2, 0.0)
                kern = {'E_eV': w['E_eV'], 'mantel_dos': w['dos'], 'W_diagonal_eV': diag.tolist(),
                        'note': '2000-point grid: only the diagonal W(e,e) is staged'}
            kern['description'] = ('Mantel W(e,e\') = N(e)N(e\')W(e,e\') / N(e)N(e\') from W_ee_raw.npy and dos_raw.npy; '
                                   'energies relative to E_F; mantel_dos is dos_raw.npy (states/eV, as written by Mantel)')
            e['w_kernel'] = kern
        if rec.get('runs'):
            e['runs'] = rec['runs']
    entries[key] = e
    av = {'inputs': True, 'alpha2F': 'eliashberg' in e, 'coulomb_table': 'coulomb' in e, 'w_kernel': 'w_kernel' in e,
          'eliashberg_solvers': (['cdos', 'vdos'] + (['vdos_W'] if 'solver_W' in e else [])) if 'eliashberg' in e else [],
          'archived_runs': sorted(e.get('runs', {}))}
    if withheld:
        av['note'] = 'graded pressure: Coulomb data, W and archived runs withheld'
    catalog_systems.setdefault(system, {})[int(p)] = av

data = {'costs': COSTS, 'entries': entries}
os.makedirs(os.path.join(LAB, 'private'), exist_ok=True)
with open(os.path.join(LAB, 'data.json'), 'w') as f:
    json.dump(data, f, separators=(',', ':'))
shutil.copy(os.path.join(HERE, 'eliashberg.py'), os.path.join(LAB, 'eliashberg.py'))
for name in ('public.json',):
    for path in (os.path.join(LAB, name), os.path.join(TASK, 'tests', name)):
        pub = json.load(open(path))
        pub['solver_budget'] = BUDGET
        if 'tests' in path:
            pub['min_report_chars'] = 400
        else:
            pub['requests'] = {k: {'cost': v} for k, v in COSTS.items()}
            pub['catalog'] = '/app/catalog.json'
        json.dump(pub, open(path, 'w'), indent=1)

# private placeholder so task_lab imports; replaced below
json.dump({'reference': {}, 'tolerance': TOL}, open(os.path.join(LAB, 'private/hidden.json'), 'w'))
sys.path.insert(0, LAB)
spec_ = importlib.util.spec_from_file_location('task_lab', os.path.join(LAB, 'task_lab.py'))
TL = importlib.util.module_from_spec(spec_)
spec_.loader.exec_module(TL)
E = TL.E

CACHE_PATH = os.path.join(DATA_DIR, 'tc_cache.json')
CACHE = json.load(open(CACHE_PATH)) if os.path.exists(CACHE_PATH) else {}


def cached(tag, fn):
    if tag not in CACHE:
        CACHE[tag] = fn()
        json.dump(CACHE, open(CACHE_PATH, 'w'), indent=1)
        print(tag, CACHE[tag], flush=True)
    return CACHE[tag]


def lab_tc(key, solver, mu=None, wc=7000.0):
    rec = entries[key]
    return cached(f'{key}|{solver}|{mu}|{wc}', lambda: TL.solve(rec, solver, [mu] if mu is not None else [], wc)['results'][0]['Tc_K'])


def full_rec(key):
    """Entry including the withheld W (construction checks only, never staged)."""
    rec = copy.deepcopy(entries[key])
    w = POOL[key]['W']
    Eg, W = w_matrix(w)
    win = np.abs(Eg) <= WINDOW_EV
    rec['solver_W'] = {'E_eV': Eg[win].tolist(), 'W_eV': W[np.ix_(win, win)].tolist()}
    return rec


def el_W(key, rec=None):
    rec = rec or entries[key]
    return TL.electrons(rec, with_W=True)


def el_plain(key):
    return TL.electrons(entries[key])


def archived_W_tc(key):
    runs = POOL[key]['runs']
    names = [r for r in runs if info_field(runs[r]['Info.txt'], 'include_Weep') == '1'
             and info_field(runs[r]['Info.txt'], 'omega_c') == '7000.0']
    best = sorted(names, key=lambda r: ('±0.05' not in runs[r]['Summary.dat'], r))[0]
    return summary_tc(runs[best]['Summary.dat']), best


# ---------------------------------------------------------------- references
ref_tc = {t: archived_W_tc(k)[0] for t, k in TARGETS.items()}
print('archived references', ref_tc, {t: archived_W_tc(k)[1] for t, k in TARGETS.items()})

checks = {}
# solver validation against every archived full-W run with staged or withheld W (omega_c = 7 eV)
for key in ['Mg2IrH6/0', 'Mg2IrH6/5', 'Mg2IrH6/10', 'Mg2IrH6/15', 'Mg2IrH6/20', 'Mg2IrH6/25', 'Mg2IrH6/30',
            'H3S/140', 'H3S/160', 'H3S/190', 'H3S/200', 'H3S/210', 'H3S/220', 'H3S/230', 'H3S/240']:
    rec = full_rec(key)
    t = cached(f'{key}|vdos_W_full|7000', lambda: TL.solve(rec, 'vdos_W', [], 7000.0)['results'][0]['Tc_K'])
    arch = archived_W_tc(key)[0]
    checks[key] = {'lab_vdos_W': t, 'archived': arch, 'rel': t / arch - 1}
    assert abs(t / arch - 1) < 0.02, (key, t, arch)
# fixed-mu* archived runs
for key, run in [('Mg2IrH6/20', 'vdos_muStar0.1'), ('Mg2IrH6/20', 'vdos_mustar_0.16'), ('Mg2IrH6/20', 'vdos_muStar0.1_1'),
                 ('Mg2IrH6/20', 'cdos_mustar_0.16'), ('H3S/200', 'vdos_muStar0.1'), ('H3S/220', 'vdos_muStar0.2'),
                 ('H3S/220', 'vdos_muStar0.2_1')]:
    info = POOL[key]['runs'][run]['Info.txt']
    m = float(info_field(info, 'muc_ME'))
    solver = 'cdos' if info_field(info, 'cDOS_flag') == '1' else 'vdos'
    t = lab_tc(key, solver, m)
    arch = summary_tc(POOL[key]['runs'][run]['Summary.dat'])
    checks[f'{key}/{run}'] = {'mu_star_in_Info': m, 'lab': t, 'archived': arch, 'rel': t / arch - 1}
    assert abs(t / arch - 1) < 0.02, (key, run, t, arch)
# H3S 180: W transferred from neighbouring pressures (the archive lacks W at 180 GPa)
for q in (170, 190):
    rec = copy.deepcopy(entries['H3S/180'])
    w = POOL[f'H3S/{q}']['W']
    Eg, W = w_matrix(w)
    win = np.abs(Eg) <= WINDOW_EV
    rec['solver_W'] = {'E_eV': Eg[win].tolist(), 'W_eV': W[np.ix_(win, win)].tolist()}
    t = cached(f'H3S/180|vdos_W_of_{q}|7000', lambda: TL.solve(rec, 'vdos_W', [], 7000.0)['results'][0]['Tc_K'])
    checks[f'H3S/180 with W of {q} GPa'] = {'lab': t, 'archived': ref_tc['H3S_180GPa'], 'rel': t / ref_tc['H3S_180GPa'] - 1}
    assert abs(t / ref_tc['H3S_180GPa'] - 1) < 0.015
print(json.dumps(checks, indent=1))

# equivalent mu*: constant DOS, omega_c = 1 eV, reproducing the archived full-W Tc at Mg2IrH6 20 GPa
TcW20 = archived_W_tc('Mg2IrH6/20')[0]


def mu_eq(target, solver='cdos', wc=1000.0, lo=0.3, hi=1.2):
    return brentq(lambda m: lab_tc('Mg2IrH6/20', solver, round(m, 6), wc) - target, lo, hi, xtol=2e-4)


mu_star_ref = cached('mu_eq_cdos_1eV_archived', lambda: mu_eq(TcW20))
# N_F of H3S at 180 GPa per formula unit per spin (the cell holds two formula units)
scf = POOL['H3S/180']['scf_in']
nat = int(re.search(r'nat\s*=\s*(\d+)', scf).group(1))
assert nat == 8                      # H6S2: two H3S formula units
d = POOL['H3S/180']['dos_file']
nef_cell_spin = float(np.interp(d['EF_eV'], d['E_eV'], d['dos'])) / 2
nef_ref = nef_cell_spin / 2
reference = {'tc_full_W_K': ref_tc, 'mu_star_cdos_1eV_Mg2IrH6_20GPa': mu_star_ref,
             'N_EF_H3S_180GPa_per_fu_per_spin': nef_ref}

# ---------------------------------------------------------------- oracle route and alternatives
NF = {k: el_plain(k).NF * 1000 for k in entries if entries[k].get('inputs', {}).get('dos_file')}


def W00(key):
    el = el_W(key) if 'solver_W' in entries[key] else None
    if el is None:      # H3S 170: DOS file + W (no alpha2F)
        rec = entries[key]
        dd = rec['inputs']['dos_file']
        w = POOL[key]['W']
        Eg, W = w_matrix(w)
        el = E.Electrons((np.array(dd['E_eV']) - dd['EF_eV']) * 1000, np.array(dd['dos']) / 2000, Eg * 1000, W * 1000)
    return el.mu / el.NF / 1000, el.mu


W0 = {k: W00(k) for k in ['Mg2IrH6/0', 'Mg2IrH6/10', 'Mg2IrH6/20', 'Mg2IrH6/30', 'H3S/160', 'H3S/170', 'H3S/200']}
print('W00, mu', W0)
NEIGH = {'Mg2IrH6/15': ('Mg2IrH6/10', 'Mg2IrH6/20'), 'Mg2IrH6/25': ('Mg2IrH6/20', 'Mg2IrH6/30'),
         'H3S/180': ('H3S/170', 'H3S/200'), 'H3S/190': ('H3S/170', 'H3S/200')}
P = lambda k: int(k.split('/')[1])


def mu_pred(key, nearest=False):
    a, b = NEIGH[key]
    if nearest:
        w = W0[min((a, b), key=lambda q: abs(P(q) - P(key)))][0]
    else:
        w = W0[a][0] + (W0[b][0] - W0[a][0]) * (P(key) - P(a)) / (P(b) - P(a))
    return NF[key] / 1000 * w * 1000


def kappa(key):
    TcW = lab_tc(key, 'vdos_W')
    m = W0[key][1]
    eff = cached(f'mu_eff_vdos_{key}', lambda: brentq(lambda x: lab_tc(key, 'vdos', round(x, 6)) - TcW, 0.5 * m, 1.3 * m, xtol=2e-4))
    return eff / m


KAPPA = {'Mg2IrH6': kappa('Mg2IrH6/20'), 'H3S': kappa('H3S/200')}
print('kappa', KAPPA)


def route(nearest=False, use_kappa=True):
    out = {}
    for t, k in TARGETS.items():
        m = mu_pred(k, nearest) * (KAPPA[k.split('/')[0]] if use_kappa else 1.0)
        out[t] = lab_tc(k, 'vdos', round(m, 5))
    return out


oracle_tc = route()
rob = {'tc_bare_mu_route': route(use_kappa=False), 'tc_nearest_W00': route(nearest=True),
       'mu_star_vs_lab_TcW': cached('mu_eq_cdos_1eV_labTcW', lambda: mu_eq(lab_tc('Mg2IrH6/20', 'vdos_W'))),
       'nef_nearest_grid': float(np.array(d['dos'])[int(np.argmin(np.abs(np.array(d['E_eV']) - d['EF_eV'])))]) / 4}

# ---------------------------------------------------------------- ablations (realistic wrong paths)
ab = {}
arch = {k: archived_W_tc(k)[0] for k in ['H3S/160', 'H3S/200', 'Mg2IrH6/10', 'Mg2IrH6/20', 'Mg2IrH6/30']}
lin = lambda k, a, b: arch[a] + (arch[b] - arch[a]) * (P(k) - P(a)) / (P(b) - P(a))
ab['tc_interpolated_in_pressure'] = {'Mg2IrH6_15GPa': lin('Mg2IrH6/15', 'Mg2IrH6/10', 'Mg2IrH6/20'),
                                     'Mg2IrH6_25GPa': lin('Mg2IrH6/25', 'Mg2IrH6/20', 'Mg2IrH6/30'),
                                     'H3S_180GPa': lin('H3S/180', 'H3S/160', 'H3S/200'),
                                     'H3S_190GPa': lin('H3S/190', 'H3S/160', 'H3S/200')}
ab['tc_conventional_mustar_0.16'] = {t: lab_tc(k, 'vdos', 0.16) for t, k in TARGETS.items()}
ab['tc_constant_mu_from_0GPa'] = {**oracle_tc, **{t: lab_tc(k, 'vdos', round(W0['Mg2IrH6/0'][1] * KAPPA['Mg2IrH6'], 5))
                                                  for t, k in TARGETS.items() if k.startswith('Mg')}}


def mu_table_col(key, sigma, col):
    for line in POOL[key]['mu_results'].splitlines():
        if not line.startswith('#') and line.split() and abs(float(line.split()[0]) - sigma) < 1e-6:
            return float(line.split()[col])
    raise KeyError


def mg_from_table(sigma, col):
    out = dict(oracle_tc)
    for t, k in TARGETS.items():
        if k.startswith('Mg'):
            a, b = NEIGH[k]
            m = (mu_table_col(a, sigma, col) + mu_table_col(b, sigma, col)) / 2
            out[t] = lab_tc(k, 'vdos', round(m, 5))
    return out


ab['tc_smeared_NF_mu_sigma0.3'] = mg_from_table(0.3, 4)
ab['tc_mu_rescaled_sigma0.01'] = mg_from_table(0.01, 5)
ab['tc_constant_dos_solver'] = {t: lab_tc(k, 'cdos', round(mu_pred(k), 5)) for t, k in TARGETS.items()}


def ad_tc(key):
    spec = TL.spectrum(entries[key])
    a, b = NEIGH[key]
    muc = []
    for q in (a, b):
        for r, v in POOL[q].get('runs', {}).items():
            if info_field(v['Info.txt'], 'include_Weep') == '1':
                muc.append(float(info_field(v['Info.txt'], 'muc_AD')))
                break
    if not muc:   # H3S 170 has no runs: use the other neighbour's value
        muc = [float(info_field(v['Info.txt'], 'muc_AD')) for r, v in POOL[b]['runs'].items()
               if info_field(v['Info.txt'], 'include_Weep') == '1'][:1]
    return E.allen_dynes(spec.lam, spec.wlog, spec.w2, float(np.mean(muc))), float(np.mean(muc))


ab['tc_allen_dynes_archived_mucAD'] = {t: ad_tc(k)[0] for t, k in TARGETS.items()}
ab['mu_star_bare_mu'] = W0['Mg2IrH6/20'][1]
EF20 = float(info_field(POOL['Mg2IrH6/20']['runs']['vdos_W']['Info.txt'], 'typEl'))
ab['mu_star_morel_anderson'] = W0['Mg2IrH6/20'][1] / (1 + W0['Mg2IrH6/20'][1] * math.log(EF20 / 1000.0))
ab['mu_star_vdos_7eV'] = cached('mu_eq_vdos_7eV', lambda: mu_eq(TcW20, 'vdos', 7000.0, 0.3, 0.9))
ab['nef_per_cell'] = nef_cell_spin
ab['nef_spin_summed_per_fu'] = nef_cell_spin

# ---------------------------------------------------------------- write private, ablations, calibration
hidden = {'reference': reference, 'tolerance': TOL}
json.dump(hidden, open(os.path.join(LAB, 'private/hidden.json'), 'w'), indent=1)
hv = list(ref_tc.values()) + [mu_star_ref, nef_ref]
json.dump(hv, open(os.path.join(LAB, 'private/hidden_values.json'), 'w'))
json.dump({'oracle_expected': {'tc': oracle_tc}, 'ablations': ab, 'robustness': rob, 'kappa': KAPPA},
          open(os.path.join(TASK, 'solution/ablation_values.json'), 'w'), indent=1)
rel = lambda d_: {t: round(d_[t] / ref_tc[t] - 1, 4) for t in d_}
calib = {'reference': reference, 'tolerance': TOL, 'solver_checks': checks, 'oracle_route_tc': oracle_tc,
         'oracle_route_rel': rel(oracle_tc), 'robustness': rob,
         'robustness_rel': {k: rel(v) for k, v in rob.items() if isinstance(v, dict)},
         'ablations': ab, 'ablations_rel': {k: rel(v) for k, v in ab.items() if isinstance(v, dict)},
         'W00_mu': W0, 'kappa': KAPPA, 'N_EF_states_per_eV_per_spin_per_cell': {k: v for k, v in NF.items()}}
json.dump(calib, open(os.path.join(HERE, 'calibration.json'), 'w'), indent=1)

# ---------------------------------------------------------------- catalog
catalog = {
    'study': 'Archived first-principles calculations of Bozier and Pickard (2026): static RPA-screened Coulomb '
             'interaction W(e,e\') from Mantel (Quantum ESPRESSO + Yambo), the Coulomb parameter mu, and isotropic '
             'Migdal-Eliashberg Tc from the IsoME solver, for Mg2IrH6 and H3S versus pressure and for other hydrides.',
    'systems': {s: {str(p): v for p, v in sorted(ps.items())} for s, ps in catalog_systems.items()},
    'requests': {
        'inputs': {'cost': COSTS['inputs'], 'example': {'kind': 'inputs', 'system': 'Mg2IrH6', 'pressure_GPa': 20},
                   'returns': 'alpha2F(omega) (omega in THz; selected phonon-smearing columns), the dense-k QE DOS file '
                              'used by IsoME, the scf input (cell and atoms) and the relaxation target pressure'},
        'coulomb_table': {'cost': COSTS['coulomb_table'], 'example': {'kind': 'coulomb_table', 'system': 'H3S', 'pressure_GPa': 200},
                          'returns': 'Mantel mu_results.dat (mu versus electronic smearing sigma) and the Mantel log header'},
        'w_kernel': {'cost': COSTS['w_kernel'], 'example': {'kind': 'w_kernel', 'system': 'H3S', 'pressure_GPa': 200},
                     'returns': 'Mantel W(e,e\') in eV: diagonal on the full grid and the matrix for |e|,|e\'| <= 5.6 eV, with the Mantel DOS'},
        'archived_run': {'cost': COSTS['archived_run'],
                         'example': {'kind': 'archived_run', 'system': 'Mg2IrH6', 'pressure_GPa': 20, 'run': 'vdos_W'},
                         'returns': 'Summary.dat and Info.txt of one archived IsoME run folder (folder names as archived)'},
        'eliashberg': {'cost': COSTS['eliashberg'],
                       'examples': [{'kind': 'eliashberg', 'system': 'Mg2IrH6', 'pressure_GPa': 20, 'solver': 'cdos', 'mu_star': [0.1, 0.16], 'omega_c_meV': 7000},
                                    {'kind': 'eliashberg', 'system': 'Mg2IrH6', 'pressure_GPa': 20, 'solver': 'vdos', 'mu_star': 0.16},
                                    {'kind': 'eliashberg', 'system': 'H3S', 'pressure_GPa': 200, 'solver': 'vdos_W'}],
                       'returns': 'Tc (K) of the linearized isotropic Migdal-Eliashberg equations for each mu*, with lambda, '
                                  'omega_log, N(E_F), Z and the chemical-potential shift at Tc'},
    },
    'solver': {
        'description': 'numpy re-implementation of the IsoME imaginary-axis equations (linearized at Tc): '
                       'cdos = constant DOS with a fixed Morel-Anderson mu* applied up to the Matsubara cutoff omega_c; '
                       'vdos = full-bandwidth equations with the energy-dependent DOS N(e), the energy shift chi and a '
                       'self-consistent chemical potential, fixed mu*; vdos_W = the same with the archived static W(e,e\') '
                       'in place of mu* (phi = phi_ph(i w_n) + phi_c(e)).',
        'conventions': 'alpha2F from the solver_column listed in each inputs result (the column used by the archived '
                       'IsoME runs); DOS from the inputs DOS file divided by 2 for spin; energy grid 1/5/50 meV inside '
                       '+-100/+-500/+-5000 meV of E_F; default omega_c = 7000 meV (allowed 300-7000); Tc to 0.01 K, '
                       'search range 10-500 K.',
        'validation': 'Where archived IsoME runs exist, the lab reproduces their Tc (fixed mu* and full W) to within 2 %.',
    },
    'purchase_costs': COSTS,
    'budget': BUDGET,
}
json.dump(catalog, open(os.path.join(TASK, 'environment/agent/catalog.json'), 'w'), indent=1)
print('reference', reference)
print('oracle route', oracle_tc, rel(oracle_tc))
print('ablations', json.dumps(calib['ablations_rel'], indent=1), ab['mu_star_bare_mu'], ab['mu_star_morel_anderson'], ab['mu_star_vdos_7eV'])
print('robustness', json.dumps(calib['robustness_rel'], indent=1), rob['mu_star_vs_lab_TcW'])
print('data.json MB', os.path.getsize(os.path.join(LAB, 'data.json')) / 2 ** 20)
