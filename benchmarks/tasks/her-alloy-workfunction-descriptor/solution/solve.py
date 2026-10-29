"""Reference solution for her-alloy-workfunction-descriptor; runs as the Harbor oracle and in check_task."""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

APP = Path(os.environ.get('APP_DIR', '/app'))
CAT = json.loads((APP/'catalog.json').read_text())
DESC = json.loads((APP/'descriptors.json').read_text())
KT = 8.617333e-5 * 298.15
PHI_REF = 5.49888                       # pure-Pt RPBE work function of the released DFT set
POTS = (100, 200, 300, 400)
ROWS = [3, 7, 12, 16]
EL = ('Ag', 'Au', 'Pd', 'Pt')


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    while result.get('status') == 'running':
        time.sleep(1)
        result = lab('job/%d' % result['id'])
    return result


class Spots:
    """Descriptor arrays for a list of MDIs of one library; j (n, 4) may be filled in later."""
    def __init__(self, lib, mdis):
        d = DESC['libraries'][lib]
        self.lib, self.mdi, self.n = lib, list(mdis), len(mdis)
        e0, w = DESC['bin_start_eV'], DESC['bin_width_eV']
        self.hist = []
        for m in mdis:
            h = d[str(m)]['gnn']['hads_histogram']
            c = e0 + w * (h['first_bin'] + np.arange(len(h['counts'])) + 0.5)
            self.hist.append((c, np.array(h['counts'], float)))
        self.emean = np.array([(c * n).sum() / n.sum() for c, n in self.hist])
        self.wf = np.array([d[str(m)]['gnn']['wf_mean_eV'] for m in mdis])
        self.j = np.full((self.n, 4), np.nan)

    def f(self, a, e):
        return np.array([(n * np.exp(-a * np.abs(c - e) / KT)).sum() / n.sum() for c, n in self.hist])


def predict(s, p, use_wf):
    U = np.array(POTS) / 1000.0
    wf = np.exp(p['alpha_wf'] * (s.wf - PHI_REF) / KT) if use_wf else 1.0
    jk = np.exp(p['logk']) * np.exp(p['alpha_c'] * U / KT)[None, :] * (s.f(p['alpha_ads'], p['e_opt']) * wf)[:, None]
    return jk / (1 + jk / np.exp(p['logjlim']))


def resid_for(s, p, use_wf):
    obs = -s.j
    ok = np.isfinite(obs) & (obs > 1e-9)          # positive (anodic) readings carry no HER information
    return np.log10(predict(s, p, use_wf)[ok]) - np.log10(obs[ok])


def fit(libs, names, fixed, use_wf):
    """Shared physics, library-specific prefactor k and transport limit j_lim (Koutecky-Levich)."""
    lo = {'alpha_c': 0.0, 'alpha_ads': 0.0, 'e_opt': -1.0, 'alpha_wf': -1.0}
    hi = {'alpha_c': 1.0, 'alpha_ads': 1.0, 'e_opt': 0.5, 'alpha_wf': 1.0}
    starts = {'alpha_c': [0.3, 0.5], 'alpha_ads': [0.5, 1.0], 'e_opt': [-0.3, -0.45], 'alpha_wf': [0.0, 0.1]}
    nl = len(libs)

    def unpack(x):
        ps = []
        for i in range(nl):
            p = dict(fixed, **{k: x[j] for j, k in enumerate(names)})
            p['logk'], p['logjlim'] = x[len(names) + i], x[len(names) + nl + i]
            ps.append(p)
        return ps

    def resid(x):
        return np.concatenate([resid_for(s, p, use_wf) for s, p in zip(libs, unpack(x))])

    lb = [lo[k] for k in names] + [-60] * nl + [-10] * nl
    ub = [hi[k] for k in names] + [60] * nl + [10] * nl
    best = None
    for st in range(2):
        x0 = [starts[k][st] for k in names]
        x0 += [np.log(np.nanmean(-s.j[:, 0][s.j[:, 0] < 0])) for s in libs] + [np.log(np.nanmax(-s.j)) for s in libs]
        x0 = np.clip(x0, np.array(lb) + 1e-9, np.array(ub) - 1e-9)
        r = least_squares(resid, x0, bounds=(lb, ub), loss='huber', f_scale=0.3, max_nfev=4000)
        if best is None or r.cost < best.cost:
            best = r
    return unpack(best.x)


def calibrate(p, s, use_wf):
    """Lib-3 prefactor and transport limit from the calibration spots; physics frozen."""
    def resid(x):
        return resid_for(s, dict(p, logk=x[0], logjlim=x[1]), use_wf)
    r = least_squares(resid, [p['logk'], p['logjlim']], loss='huber', f_scale=0.3)
    return dict(p, logk=r.x[0], logjlim=r.x[1])


def r2(s, p, use_wf):
    obs = -s.j
    ok = np.isfinite(obs) & (obs > 1e-9)
    a, b = np.log10(obs[ok]), np.log10(predict(s, p, use_wf)[ok])
    return 1 - ((a - b) ** 2).sum() / ((a - a.mean()) ** 2).sum()


def main():
    evidence = []
    # 1. SECCM rows of Lib-1 and Lib-2 at all four potentials (saturation and Tafel slope both needed).
    train = {}
    for lib in ('Lib-1', 'Lib-2'):
        mdis = [sp['mdi'] for sp in CAT['libraries'][lib]['spots'] if sp['row'] in ROWS]
        s = Spots(lib, mdis)
        pos = {m: i for i, m in enumerate(mdis)}
        for row in ROWS:
            for k, mv in enumerate(POTS):
                r = lab('evaluate', {'kind': 'seccm_row', 'library': lib, 'row': row, 'potential_mV': -mv})
                evidence.append(r['id'])
                for sp in r['result']['spots']:
                    s.j[pos[sp['mdi']], k] = sp['j_A_cm2']
        train[lib] = s
    # 2. Lib-3 calibration spots spread over composition (farthest-point sampling).
    l3 = {sp['mdi']: sp for sp in CAT['libraries']['Lib-3']['spots']}
    elig = CAT['libraries']['Lib-3']['calibration_spots_available']
    comp = np.array([[l3[m]['composition_at_frac'][e] for e in EL] for m in elig])
    idx = [int(np.argmax(comp[:, 2]))]
    while len(idx) < CAT['max_lib3_spot_purchases']:
        d = np.min(np.linalg.norm(comp[:, None, :] - comp[None, idx, :], axis=2), axis=1)
        idx.append(int(np.argmax(d)))
    cal_mdi = sorted(elig[i] for i in idx)
    cal = Spots('Lib-3', cal_mdi)
    for i, m in enumerate(cal_mdi):
        r = lab('evaluate', {'kind': 'lib3_spot', 'mdi': m})
        evidence.append(r['id'])
        cal.j[i] = [r['result']['j_A_cm2_by_potential_mV'][str(-mv)] for mv in POTS]
    # 3. Coverage: H-H interaction slope from the equimolar paths, with each site's own dilute energy removed.
    omegas, alphas = [], []
    for pth in CAT['coverage_paths']:
        if not pth['composition_id'].startswith('ag025'):
            continue
        r = lab('evaluate', {'kind': 'coverage_path', 'id': pth['id']})
        evidence.append(r['id'])
        res = r['result']
        th = np.array([x['theta_ML'] for x in res['filling_path']])
        diff = np.array([x['differential_adsorption_energy_eV'] for x in res['filling_path']])
        own = np.array([res['single_site_adsorption_energies_eV'][x['site_added']] for x in res['filling_path']])
        omegas.append(np.polyfit(th, diff - own, 1)[0])
        alphas.append(np.polyfit(th, diff, 1)[0])
    shift = float(np.mean(omegas) * (1 - 1 / 9))
    # 4. Pure-metal anchors from DFT (for the report).
    anchors = {}
    for cid in ('pd100', 'pt100'):
        r = lab('evaluate', {'kind': 'dft_composition', 'id': cid})
        evidence.append(r['id'])
        anchors[cid] = float(np.mean([e for v in r['result']['variants'].values() for e in v['sites'].values()]))
    # 5. Activity models fitted to Lib-1 + Lib-2.
    libs = [train['Lib-1'], train['Lib-2']]
    base = {'alpha_c': 0.4, 'alpha_ads': 1.0, 'e_opt': -0.24, 'alpha_wf': 0.0}
    pA = fit(libs, ['alpha_c', 'alpha_ads', 'e_opt'], base, False)
    pC = fit(libs, ['alpha_c', 'alpha_ads', 'e_opt', 'alpha_wf'], base, True)
    pB = fit(libs, ['alpha_c', 'alpha_wf'], dict(base, alpha_ads=0.0), True)
    p024 = fit(libs, ['alpha_c', 'alpha_ads'], base, False)
    r2s = {name: [round(r2(s, p, w), 3) for s, p in zip(libs, ps)]
           for name, ps, w in (('A', pA, False), ('C', pC, True), ('B', pB, True), ('A_e024', p024, False))}
    p0 = dict(pA[0], logk=np.mean([p['logk'] for p in pA]), logjlim=np.mean([p['logjlim'] for p in pA]))
    p3 = calibrate(p0, cal, False)
    graded = CAT['libraries']['Lib-3']['spots_to_predict']
    G = Spots('Lib-3', graded)
    pred = np.log10(predict(G, p3, False))
    submission = {'lib3_log10_abs_current': {str(m): [float(v) for v in row] for m, row in zip(graded, pred)},
                  'volcano_centre_eV': float(pA[0]['e_opt']), 'coverage_shift_eV': shift,
                  'alpha_wf': float(pC[0]['alpha_wf']), 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    ks = [round(float(np.exp(p['logk'])), 3) for p in pA] + [round(float(np.exp(p3['logk'])), 3)]
    jl = [round(float(np.exp(p['logjlim'])), 3) for p in pA] + [round(float(np.exp(p3['logjlim'])), 3)]
    (APP/'report.md').write_text(f"""# Ag-Au-Pd-Pt HER: reference analysis

Data: SECCM rows {ROWS} of Lib-1 and Lib-2 at all four potentials; Lib-3 calibration spots {cal_mdi}
(spread over composition); all ten equimolar coverage paths; DFT pure Pd and Pt.

Coverage. On the alloy the differential adsorption energy along the filling path rises steeply
(slope {np.mean(alphas):.2f} eV/ML) mostly because the strongest hollows fill first: that is site
heterogeneity, which the dilute site-energy distribution already contains. Subtracting each added
site's own dilute energy leaves the H-H interaction, slope {np.mean(omegas):.3f} eV/ML (pure metals: the two
slopes coincide). From 1/9 to 1 ML this shifts every site by {shift:.3f} eV, so the dilute optimum
Delta E = -0.24 eV (Delta G = 0 with the +0.24 eV ZPE/entropy term) moves to about {-0.24 - shift:.2f} eV.

Activity model. j = k exp(alpha_c e|U|/kT) <exp(-alpha_ads |dE - dE*|/kT)>_sites, capped by a
transport limit j_lim (Koutecky-Levich): the 300-400 mV currents saturate near 2 A/cm^2. The average runs
over the whole site distribution of each spot: the spot-mean energies all lie on the weak-binding side of
the optimum, so a mean-energy volcano cannot even locate the peak. k and j_lim differ strongly between
libraries (k = {ks[0]}, {ks[1]}; j_lim = {jl[0]}, {jl[1]} A/cm^2), so the physics (alpha_c = {pA[0]['alpha_c']:.3f},
alpha_ads = {pA[0]['alpha_ads']:.2f}, dE* = {pA[0]['e_opt']:.3f} eV) was fitted jointly with library-specific k
and j_lim, and Lib-3's k = {ks[2]} and j_lim = {jl[2]} A/cm^2 were calibrated on the purchased spots only.
The fitted centre agrees with the DFT coverage estimate; a fixed literature centre of -0.24 eV fits far
worse (R2_log {r2s['A_e024']} vs {r2s['A']}). DFT anchors: Pd {anchors['pd100']:.3f} eV and Pt {anchors['pt100']:.3f} eV, so at
realistic coverage Pd-rich sites sit near the optimum and the most active Lib-3 region is Pd-rich with some Pt.

Work function. A work-function-only model fits well (R2_log {r2s['B']}, alpha_wf = {pB[0]['alpha_wf']:.3f}) because
Phi tracks the Pd/Pt content, i.e. the binding term. With the coverage-corrected binding term present,
alpha_wf drops to {pC[0]['alpha_wf']:.3f} (R2_log {r2s['C']} vs {r2s['A']}), and it is not stable between libraries.
Pooling libraries without library-specific k inflates it further (library offsets masquerade as a descriptor
effect). So the work function carries little information beyond coverage-corrected adsorption for these
four metals; its physical role (interfacial field / charge transfer) cannot be separated from composition here.

Limits: the coverage slope varies strongly between decorations (single-decoration values 0.01-0.65 eV/ML);
alpha_c differs between libraries; absolute Lib-3 currents depend on the calibration spots.
""")
    lab('submit', submission)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
