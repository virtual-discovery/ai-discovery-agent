"""Activity-model utilities (numpy/scipy only) used by the construction scripts.

j_kin,p(U) = k_lib * exp(alpha_c*e|U|/kT) * f_p * exp(alpha_wf*(Phi_p - Phi_ref)/kT)
f_p        = <exp(-alpha_ads*|dE - e_opt|/kT)> over the spot's dilute site-energy distribution
j_p(U)     = j_kin / (1 + j_kin/j_lim_lib)                       (Koutecky-Levich cap)
Residuals are log10|j| with a Huber loss (delta 0.3), as in the authors' production fits.
"""
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
KT = 8.617333e-5 * 298.15
PHI_REF = 5.49888
POTS = (100, 200, 300, 400)
POOL = json.loads(gzip.open(HERE/'pool.json.gz', 'rt').read())
E0, DE = POOL['bin_edges_start_eV'], POOL['bin_width_eV']


def centres_counts(h):
    c = E0 + DE * (h['first_bin'] + np.arange(len(h['counts'])) + 0.5)
    return c, np.array(h['counts'], float)


class Lib:
    """Arrays for one library's spots (descriptor surrogate `sur`)."""
    def __init__(self, name, spots=None, sur='dimenet'):
        allspots = POOL['libraries'][name]['spots']
        self.spots = [s for s in allspots if spots is None or s['mdi'] in spots]
        self.name, self.n = name, len(self.spots)
        self.mdi = np.array([s['mdi'] for s in self.spots])
        self.comp = np.array([[s['composition'][e] for e in ('Ag', 'Au', 'Pd', 'Pt')] for s in self.spots])
        self.j = np.array([[s['j_A_cm2'][str(mv)] for mv in POTS] for s in self.spots])
        self.wf = np.array([s[sur]['wf_mean_eV'] for s in self.spots])
        self.hist = [centres_counts(s[sur]['hads']) for s in self.spots]
        self.emean = np.array([s[sur]['hads']['mean'] for s in self.spots])

    def f(self, alpha_ads, e_opt, mode='dist'):
        if mode == 'mean':
            return np.exp(-alpha_ads * np.abs(self.emean - e_opt) / KT)
        return np.array([(n * np.exp(-alpha_ads * np.abs(c - e_opt) / KT)).sum() / n.sum() for c, n in self.hist])


def predict(lib, p, cap=True, mode='dist', use_wf=True):
    """p: dict with logk, alpha_c, alpha_ads, e_opt, alpha_wf, logjlim. Returns |j| (n, 4)."""
    U = np.array(POTS) / 1000.0
    fp = lib.f(p['alpha_ads'], p['e_opt'], mode)
    wf = np.exp(p.get('alpha_wf', 0.0) * (lib.wf - PHI_REF) / KT) if use_wf else 1.0
    jk = np.exp(p['logk']) * np.exp(p['alpha_c'] * U / KT)[None, :] * (fp * wf)[:, None]
    if cap:
        jk = jk / (1 + jk / np.exp(p['logjlim']))
    return jk


def fit(libs, free=('alpha_c', 'alpha_ads', 'e_opt'), fixed=None, cap=True, mode='dist', use_wf=False,
        pots=POTS, huber=0.3, x0=None, shared_k=False):
    """Joint fit: per-library logk and logjlim, shared physics parameters."""
    fixed = dict({'alpha_c': 0.4, 'alpha_ads': 1.0, 'e_opt': -0.24, 'alpha_wf': 0.0}, **(fixed or {}))
    names = list(free) + (['alpha_wf'] if use_wf and 'alpha_wf' not in free else [])
    pidx = [POTS.index(mv) for mv in pots]
    nk = 1 if shared_k else len(libs)

    def unpack(x):
        phys = dict(fixed)
        for i, nme in enumerate(names):
            phys[nme] = x[i]
        out = []
        for li in range(len(libs)):
            q = dict(phys)
            q['logk'] = x[len(names) + (0 if shared_k else li)]
            q['logjlim'] = x[len(names) + nk + (0 if shared_k else li)] if cap else 50.0
            out.append(q)
        return out

    def resid(x):
        r = []
        for lib, q in zip(libs, unpack(x)):
            pr = predict(lib, q, cap, mode, use_wf)[:, pidx]
            obs = -lib.j[:, pidx]
            ok = obs > 1e-9
            r.append(np.log10(np.maximum(pr[ok], 1e-30)) - np.log10(obs[ok]))
        return np.concatenate(r)

    starts = {'alpha_c': [0.3, 0.5], 'alpha_ads': [0.5, 1.0], 'e_opt': [-0.3, -0.45], 'alpha_wf': [0.0, 0.1]}
    lo = {'alpha_c': 0.0, 'alpha_ads': 0.0, 'e_opt': -1.0, 'alpha_wf': -1.0}
    hi = {'alpha_c': 1.0, 'alpha_ads': 1.0, 'e_opt': 0.5, 'alpha_wf': 1.0}
    best = None
    for s in range(2):
        xp = [starts[n][s] for n in names]
        xk = [np.log(np.mean(-l.j[:, 0][l.j[:, 0] < 0])) for l in libs][:nk]
        xj = [np.log(np.max(-l.j)) for l in libs][:nk] if cap else []
        x = np.array(x0 if x0 is not None else xp + xk + xj, float)
        lb = [lo[n] for n in names] + [-60] * nk + ([-10] * nk if cap else [])
        ub = [hi[n] for n in names] + [60] * nk + ([10] * nk if cap else [])
        x = np.clip(x, np.array(lb) + 1e-9, np.array(ub) - 1e-9)
        res = least_squares(resid, x, bounds=(lb, ub), loss='huber', f_scale=huber, max_nfev=4000)
        if best is None or res.cost < best.cost:
            best = res
    return unpack(best.x), best


def r2log(lib, q, cap=True, mode='dist', use_wf=False, pots=POTS):
    pidx = [POTS.index(mv) for mv in pots]
    pr = predict(lib, q, cap, mode, use_wf)[:, pidx]
    obs = -lib.j[:, pidx]
    ok = obs > 1e-9
    a, b = np.log10(obs[ok]), np.log10(pr[ok])
    return 1 - ((a - b) ** 2).sum() / ((a - a.mean()) ** 2).sum(), np.sqrt(((a - b) ** 2).mean())
