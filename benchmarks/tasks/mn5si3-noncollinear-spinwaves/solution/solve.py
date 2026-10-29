"""Reference solution for mn5si3-noncollinear-spinwaves; runs as the Harbor oracle and in check_task.

Route: buy the KKR exchange tables and the Mn1 couplings, relax the classical ground state of the 16-site
orthorhombic cell (8 Mn2a + 8 Mn1) from several starts, fit one uniaxial anisotropy to the 10 K gap, compute
noncollinear LSWT neutron weights, calibrate the momentum resolution on purchased cuts and the overall energy
scale (Hamiltonian normalisation) on the purchased h and k dispersion tables, then predict the withheld
(0,3,Ql) branch positions (archive's two-Gaussian extraction), the (2,1,0) zone-centre excitation and the canting.
No paper parameters are used.
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np
import mnsi_lswt as L

APP = Path(os.environ.get('APP_DIR', '/app'))
S_SPIN = 2.44 / 2.0            # S = m / (g muB) with the DFT moment, g = 2
SIGMA_E = 0.8                  # meV, energy resolution used for constant-energy scans
SCALES = (0.5, 0.7, 0.85, 1.0, 1.2, 1.4, 2.0)
SIGMA_Q = (0.04, 0.06, 0.08)


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(request, evidence):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(1)
        r = lab('job/%d' % r['id'])
    evidence.append(r['id'])
    return r['result']


# ---------------------------------------------------------------- model construction
def model_sites(catalog, with_mn1=True):
    labels = [s['label'] for s in catalog['sites'] if s['label'].startswith('Mn2a')]
    if with_mn1:
        labels += [s['label'] for s in catalog['sites'] if s['label'].startswith('Mn1')]
    frac = {s['label']: np.array(s['frac'], float) for s in catalog['sites']}
    return labels, np.array([frac[l] for l in labels])


def bonds_to_pairs(bonds, labels, pos, scale=1.0):
    idx = {l: k for k, l in enumerate(labels)}
    rows = []
    for b in bonds:
        if b['i'] not in idx or b['j'] not in idx:
            continue
        i, j = idx[b['i']], idx[b['j']]
        d = np.array(b['vector_A']) / L.LAT
        R = d + pos[i] - pos[j]
        assert np.abs(R - np.round(R)).max() < 1e-3, b
        rows.append([i, j, *np.round(R), scale * b['J_meV']])
    return np.array(rows)


def build(catalog, tables, with_mn1=True, aniso=(), seed=0):
    labels, pos = model_sites(catalog, with_mn1)
    bonds = list(tables['mn2']) + (list(tables['mn1']) if with_mn1 else [])
    pairs = bonds_to_pairs(bonds, labels, pos)
    e, E, torque = L.ground_state(pairs, aniso, n=len(labels), starts=8, seed=seed)
    return {'labels': labels, 'pos': pos, 'pairs': pairs, 'e': e, 'E': E, 'aniso': list(aniso), 'S': S_SPIN}


def as_tuple(m, weights=True):
    return (m['pos'], m['pairs'], m['e'], m['S'], m['aniso'])


def canting_deg(m):
    n = len([l for l in m['labels'] if l.startswith('Mn2a')])
    if n == len(m['labels']):
        return 0.0
    e = m['e']
    mn1 = e[n:]
    ax = mn1[0] / np.linalg.norm(mn1[0])
    return float(np.degrees(np.arcsin(np.abs(e[:n] @ ax).mean())))


def gap_at_120(m):
    w, i = L.sqw(np.array([1.0, 2.0, 0.0]), *as_tuple(m))
    sel = (w < 3.0) & (i > 1e-6)
    return float(np.sum(w[sel] * i[sel]) / np.sum(i[sel])) if sel.any() else 0.0


def fit_anisotropy(catalog, tables, gap_target, with_mn1=True):
    best = None
    for K in (0.005, 0.01, 0.02, 0.04, 0.08):
        m = build(catalog, tables, with_mn1, [(K, [0, 1, 0])])
        g = gap_at_120(m)
        if best is None or abs(g - gap_target) < best[0]:
            best = (abs(g - gap_target), K, g)
    return best[1], best[2]


# ---------------------------------------------------------------- comparison with purchased data
def chi2_cut(model_profile, cut):
    y, e = np.array(cut['counts'], float), np.array(cut['error'], float)
    A = np.c_[model_profile, np.ones_like(y), np.linspace(0, 1, len(y))]
    c = np.linalg.lstsq(A / e[:, None], y / e, rcond=None)[0]
    r = (A @ c - y) / e
    return float(r @ r / (len(y) - 3))


def windowed(cut):
    lo, hi = L.WIN[cut['direction']]
    x = np.array(cut['x_rlu'])
    m = (x >= lo) & (x <= hi)
    return {'direction': cut['direction'], 'energy_meV': cut['energy_meV'], 'x': x[m],
            'counts': np.array(cut['counts'])[m], 'error': np.array(cut['error'])[m]}


def table_rms(spec, tables, sigQ, scale, emax=None):
    res = []
    for key, pts in tables.items():
        c = L.DIRS[key][1]
        pts = [p for p in pts if p['q_err_rlu'] > 0.001 and (emax is None or p['E_meV'] < emax)]
        Es = sorted({p['E_meV'] for p in pts})
        guesses = [np.mean([abs(p['q_rlu'] - c) for p in pts if p['E_meV'] == E]) for E in Es]
        pk = L.predicted_peaks(spec[key], key, Es, guesses, SIGMA_E, sigQ, scale)
        for p in pts:
            lo, hi = pk[p['E_meV']]
            res.append(p['q_rlu'] - (lo if p['q_rlu'] < c else hi))
    return float(np.sqrt(np.mean(np.square(res))))


def gaussian_peak(scan, emin=5.0, emax=20.0):
    from scipy.optimize import curve_fit
    E, y, e = (np.array(scan[k], float) for k in ('E_meV', 'counts', 'error'))
    m = (E >= emin) & (E <= emax)
    E, y, e = E[m], y[m], e[m]
    f = lambda x, a, x0, w, b: a * np.exp(-0.5 * ((x - x0) / w) ** 2) + b
    p, cov = curve_fit(f, E, y, p0=[y.max() - y.min(), E[np.argmax(y)], 2.0, y.min()], sigma=e, maxfev=20000)
    flat = np.sum(((y - np.average(y, weights=1 / e ** 2)) / e) ** 2)
    peak = np.sum(((f(E, *p) - y) / e) ** 2)
    return float(p[1]), float(p[0] / np.sqrt(cov[0, 0])), float(flat - peak)


def strongest_mode(m, Q, emin, emax, scale):
    w, i = L.sqw(np.array(Q, float), *as_tuple(m))
    w = scale * w
    sel = (w >= emin) & (w <= emax)
    if not sel.any():
        return None, 0.0
    tot = np.sum(i[w < 40])
    k = np.argmax(np.where(sel, i, -1))
    return float(w[k]), float(i[k] / tot)


# ---------------------------------------------------------------- main analysis (shared with calibrate.py)
def zone_centre_mode(m, scale):
    E, frac = strongest_mode(m, (2, 1, 0), 5, 20, scale)
    if E is None:      # no mode in the window: report the lowest mode above 5 meV
        w, i = L.sqw(np.array([2.0, 1.0, 0.0]), *as_tuple(m))
        w = scale * w
        E = float(w[w > 5][0]) if (w > 5).any() else 0.0
    return E


def analyse(catalog, data, variant=None):
    """data: purchased results. variant: switches used by the construction scripts for robustness checks and ablations."""
    v = dict(variant or {})
    tables = {'mn2': data['exchange'][v.get('exchange', 'noncollinear_r4.5')]['bonds'],
              'mn1': data['exchange']['mn1_couplings']['bonds']}
    gap_target = data['gap_table']['10K']['E_meV'][0]
    with_mn1 = v.get('mn1', True)    # KKR gives ~1.2 muB on Mn1 in the noncollinear state; see report
    K, gap = fit_anisotropy(catalog, tables, gap_target, with_mn1) if v.get('aniso', 'fit') == 'fit' else (0.0, None)
    aniso = [(K, [0, 1, 0])] if K > 0 else []
    m = build(catalog, tables, with_mn1, aniso)
    weights = v.get('weights', True)
    spec = {k: L.line_spectrum(k, as_tuple(m), weights=weights) for k in 'hkl'}
    cuts = [windowed(c) for c in data['cuts']]
    sq = min(SIGMA_Q, key=lambda s: np.mean([chi2_cut(L.cut_profile(spec[c['direction']], c['direction'], c['energy_meV'],
                                                                    c['x'], SIGMA_E, s), c) for c in cuts]))
    sq = v.get('sigQ', sq)
    tabs = {k: data['dispersion'][k]['points'] for k in data['dispersion']}
    scan = {s: table_rms(spec, tabs, sq, s) for s in SCALES}
    scale = v.get('scale', min(scan, key=scan.get))
    energies = catalog['to_predict']['l_dispersion']['energies_meV']
    x = np.linspace(0.5, 1.5, 401)
    guesses = []
    for E in energies:
        prof = L.cut_profile(spec['l'], 'l', E, x, SIGMA_E, sq, scale)
        guesses.append(max(0.01, abs(x[x < 1][np.argmax(prof[x < 1])] - 1)))
    pk = L.predicted_peaks(spec['l'], 'l', energies, guesses, SIGMA_E, sq, scale)
    pred = {('%g' % E): [round(float(pk[E][0]), 4), round(float(pk[E][1]), 4)] for E in energies}
    return {'pred': pred, 'mn1': with_mn1, 'canting': canting_deg(m), 'optical': zone_centre_mode(m, scale),
            'scale': scale, 'scale_scan': scan, 'sigQ': sq, 'K_y': K, 'gap_model': gap, 'E': m['E'] / len(m['labels'])}


PLAN = {'exchange': ['noncollinear_r4.5', 'mn1_couplings', 'collinear_AFMz_r4.5', 'noncollinear_r17'],
        'dispersion': ['h', 'k'],
        'cuts': [('h', 4), ('h', 5), ('h', 8), ('k', 5), ('k', 7), ('k', 9)],
        'escans': ['gap_10K_dh0.00', 'gap_10K_dh0.02', '030_10K']}


def purchase(evidence):
    data = {'exchange': {}, 'dispersion': {}, 'cuts': []}
    for name in PLAN['exchange']:
        data['exchange'][name] = buy({'kind': 'exchange_set', 'name': name}, evidence)
    for d in PLAN['dispersion']:
        data['dispersion'][d] = buy({'kind': 'dispersion_table', 'direction': d}, evidence)
    for d, E in PLAN['cuts']:
        data['cuts'].append(buy({'kind': 'constant_energy_cut', 'direction': d, 'energy_meV': E}, evidence))
    data['gap_scans'] = [buy({'kind': 'energy_scan', 'name': n}, evidence) for n in PLAN['escans']]
    data['gap_table'] = buy({'kind': 'gap_table'}, evidence)
    return data


def main():
    catalog = json.loads((APP/'catalog.json').read_text())
    evidence = []
    data = purchase(evidence)
    r = analyse(catalog, data)
    alt = {n: analyse(catalog, data, {'exchange': n}) for n in ('collinear_AFMz_r4.5', 'noncollinear_r17')}
    m2 = analyse(catalog, data, {'mn1': False})
    spread = {n: round(max(abs(a - b) for E in r['pred'] for a, b in zip(r['pred'][E], a_['pred'][E])), 3)
              for n, a_ in alt.items()}
    submission = {'l_branch_positions_rlu': r['pred'], 'zone_centre_excitation_meV': round(r['optical'], 2),
                  'mn1_carries_moment': r['mn1'], 'mn2_canting_deg': round(r['canting'], 2), 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# Mn5Si3 AFM1 spin waves (reference solution)

Model: classical spins on the 16 magnetic sites of the orthorhombic cell (8 Mn2a, 8 Mn1; Mn2b nonmagnetic), unit vectors,
E = -sum over all listed table entries J e_i.e_j - K_y sum (e_i.y)^2, Landau-Lifshitz dynamics with S = 2.44/2.
Exchange: KKR noncollinear-reference table (r <= 4.5 A) for Mn2a-Mn2a plus the authors' Mn1 couplings; the bond vectors fix
the lattice translations. Mn1 is kept magnetic because KKR gives it ~1.2 muB in the noncollinear state and the field-temperature
phase diagram shows several field-induced transitions inside AFM1, which a collinear Mn2-only structure would not have.
The ground state was relaxed from several random starts in the k = 0 orthorhombic cell (the magnon cones are centred on integer
(H,K,L) of this cell): E = {r['E']:.4f} meV/site, Mn2a moments canted by {r['canting']:.1f} deg out of the plane perpendicular to
the Mn1 moments (noncollinear AFM1). A Mn2-only model relaxes to a collinear state (canting 0) and has no magnon between 5 and
20 meV at (2,1,0) (lowest above 5 meV: {m2['optical']:.1f} meV).

Anisotropy: one easy axis K_y = {r['K_y']} meV, fitted to the 10 K gap at (1,2,0) (model {r['gap_model']}); it fixes the moment
directions and hence the polarisation factors. Neutron weights: noncollinear LSWT with local frames, symplectic normalisation,
perpendicular projection and the Mn form factor; constant-energy scans simulated with sigma_E = {SIGMA_E} meV and Gaussian
momentum resolution sigma_Q = {r['sigQ']} 1/A chosen on the purchased cuts; branch positions extracted with the archive's
two-Gaussian fit.

Energy scale: the table normalisation (each bond listed from both ends; spin length) was checked on the h and k dispersion tables:
rms by scale factor {json.dumps({k: round(v, 3) for k, v in r['scale_scan'].items()})}; scale {r['scale']} retained. A doubled
Hamiltonian (bonds counted twice) or the free-ion spin 5/2 are excluded. The h and k directions individually prefer somewhat
different scales (about 0.8-0.9 and 1.0-1.1), a sign that the DFT model is only semi-quantitative.

Predictions: (0,3,Ql) branch positions from the LSWT intensities (at high energy the intensity sits on the steep upper branch, not
on the acoustic branch); the zone-centre excitation at (2,1,0) between 5 and 20 meV is the strongest Mn1-Mn2 optical mode,
{r['optical']:.1f} meV. Alternative exchange tables change the l predictions by at most {spread} r.l.u.

Not supported by these data: individual Mn1 couplings (not from DFT), the full anisotropy tensor, the long-range exchange tail,
the soft AFM2 spectra (80 K gaps of 0.18 and 0.41 meV are far below any of these models) and the k-direction plateau near
2.75 r.l.u. at 9-17 meV, which no table reproduces. A phenomenological cone fitted to h and k predicts the l dispersion about as well
as LSWT; the microscopic model is needed for the optical mode and the moment arrangement.
""")
    lab('submit', submission)
    for _ in range(120):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
