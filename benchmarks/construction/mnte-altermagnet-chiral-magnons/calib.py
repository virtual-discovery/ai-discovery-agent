"""Reference analysis, defensible variants and wrong paths for mnte-altermagnet-chiral-magnons.

Runs in ai-theorist-benchmark-runtime (numpy 1.21 / scipy 1.8):
    MNTE_CACHE=/cache python3 -B calib.py   -> calibration.json
"""
import json
import numpy as np
from scipy.optimize import brentq
import forward as F
import lswt
import reduce as R

PRIOR = {1: 3.99, 2: -0.120, 3: 0.472}
D = 0.0482
JSUM = 0.0460     # J10 + J11 of the Hamiltonian the authors used (0.0681 - 0.0221); reproduces Fig1c.mat to 0.012 meV
L_NOM = 1.25
L_TRUE = 1.25 * 6.711 / 6.87
QH_SCALE = 4.17 / 4.148
FR = lswt.frames(30.0)
X = R.X + '/Figure2/'
CUTS = {'g': dict(file='Fig2g.dat', along='u', fixed=0.45, rng=(0.80, 0.95), guess=[0.5, 0.7, 0.845, 0.45, 0.88, 0.03]),
        'h': dict(file='Fig2h.dat', along='v', fixed=1.00, rng=(0.52, 0.72), guess=[0.7, 0.5, 0.60, 0.6, 0.635, 0.03])}


def jset(dJ, jsum=JSUM):
    J = dict(PRIOR)
    J[10] = jsum / 2 + dJ / 2
    J[11] = jsum / 2 - dJ / 2
    return J


def data_peaks():
    out = {}
    for k, c in CUTS.items():
        d = np.loadtxt(X + c['file'])
        p, s, chi = R.fit_two_gauss(d[:, 2], d[:, 0] * 1e4, d[:, 1] * 1e4, c['guess'], slope=True)
        out[k] = dict(pos=[p[2], p[4]], err=[s[2], s[4]], width=p[5], sep=p[4] - p[2])
    return out


def cut_crossings(J, cut, n=31):
    c = CUTS[cut]
    xs = np.linspace(c['rng'][0], c['rng'][1], n)
    out = []
    for m in range(2):
        def f(x):
            u, v = (x, c['fixed']) if c['along'] == 'u' else (c['fixed'], x)
            return lswt.modes(F.uv_to_hkl(u, v), J, D, FR)[0][m] - F.hyspec_energy(u, v)
        vals = [f(x) for x in xs]
        out += [brentq(f, xs[i], xs[i + 1]) for i in range(n - 1) if vals[i] * vals[i + 1] < 0]
    out = sorted(out)
    return out[0], out[-1]


def in20_crossings(J, l, E=32.0):
    hs = np.linspace(1.2, 1.5, 31)
    out = []
    for m in range(2):
        f = lambda h: lswt.modes([h, 0, l], J, D, FR)[0][m] - E
        v = [f(h) for h in hs]
        r = [brentq(f, hs[i], hs[i + 1]) for i in range(30) if v[i] * v[i + 1] < 0]
        out.append(r[0])
    return sorted(out)


def in20_truth():
    res = {}
    allrows = []
    for st, folder in (('A', 'Fig3d'), ('B', 'Fig3c')):
        rows = R.in20_points(folder)
        allrows += rows
        x1, y1, e1 = R.binned(rows, (1, -1))
        x2, y2, e2 = R.binned(rows, (-1, 1))
        m = (x1 > 1.2) & (x1 < 1.5)
        ys, es = (y1 + y2) / 2, np.hypot(e1, e2) / 2
        p, s, c = R.fit_two_gauss(x1[m], ys[m], es[m], [3, 55, 1.30, 60, 1.395, 0.05])
        yd = (y2 - y1) / 2
        sg = p[5] / 2.3548
        G = np.c_[np.exp(-0.5 * ((x1[m] - p[2]) / sg) ** 2), np.exp(-0.5 * ((x1[m] - p[4]) / sg) ** 2)]
        A = G / es[m, None]
        amp = np.linalg.lstsq(A, yd[m] / es[m], rcond=None)[0]
        err = np.sqrt(np.diag(np.linalg.inv(A.T @ A)))
        res[st] = dict(centres=[p[2], p[4]], width=p[5], chiral_amp=amp.tolist(), chiral_err=err.tolist(),
                       signs=['+' if a > 0 else '-' for a in amp], significance=(np.abs(amp) / err).tolist())
    x1, y1, e1 = R.binned(allrows, (1, -1))
    x2, y2, e2 = R.binned(allrows, (-1, 1))
    m = (x1 > 1.2) & (x1 < 1.5)
    p, s, c = R.fit_two_gauss(x1[m], ((y1 + y2) / 2)[m], (np.hypot(e1, e2) / 2)[m], [3, 55, 1.30, 60, 1.395, 0.05])
    res['combined'] = dict(centres=[p[2], p[4]], centre_err=[s[2], s[4]], mean=(p[2] + p[4]) / 2, separation=p[4] - p[2])
    return res


def dj_from_cuts(dat, grid=np.arange(0.03, 0.161, 0.01)):
    seps = {k: [] for k in CUTS}
    for dJ in grid:
        for k in CUTS:
            a, b = cut_crossings(jset(dJ), k)
            seps[k].append(b - a)
    est = {k: float(np.interp(dat[k]['sep'], seps[k], grid)) for k in CUTS}
    return est, {k: [float(x) for x in v] for k, v in seps.items()}


if __name__ == '__main__':
    dat = data_peaks()
    truth = in20_truth()
    est, table = dj_from_cuts(dat)
    dj_ref = (est['g'] + est['h']) / 2
    out = {'data_peaks': dat, 'in20_truth': truth, 'dj_estimates': est, 'dj_ref': dj_ref, 'sep_table': table}
    preds = {}
    for name, dJ, l, scale in [('oracle', dj_ref, L_TRUE, 1.0), ('h_scaled', dj_ref, L_TRUE, QH_SCALE),
                               ('dj_g_only', est['g'], L_TRUE, 1.0), ('dj_h_only', est['h'], L_TRUE, 1.0),
                               ('dj_published', 0.0902, L_TRUE, 1.0), ('nominal_l', dj_ref, L_NOM, 1.0),
                               ('nominal_l_h_scaled', dj_ref, L_NOM, QH_SCALE), ('no_split', 0.0, L_TRUE, 1.0),
                               ('half_dj', dj_ref / 2, L_TRUE, 1.0)]:
        c = in20_crossings(jset(dJ), l)
        c = [x * scale for x in c]
        preds[name] = dict(dJ=dJ, centres=c, mean=(c[0] + c[1]) / 2, separation=c[1] - c[0])
    out['in20_predictions'] = preds
    # Wrong path seen in both dev pilots: read the chiral cut i as the antisymmetric difference of two fully chiral
    # branches, M = S(u - du/2)/2 - S(u + du/2)/2 ~ -(du/2) dS/du (no dilution by the domain imbalance).
    gd, idat = np.loadtxt(X + 'Fig2g.dat'), np.loadtxt(X + 'Fig2i.dat')
    xg, S, M, eM = gd[:, 2], gd[:, 0], idat[:, 0], idat[:, 1]
    dS = np.gradient(S, xg)
    a = np.sum(M * dS / eM ** 2) / np.sum(dS ** 2 / eM ** 2)
    du = 2 * abs(a)
    dj_dip = float(np.interp(du, [0.0] + table['g'], [0.0] + list(np.arange(0.03, 0.161, 0.01))))
    out['dipole'] = dict(slope=a, du=du, dJ=dj_dip, peaks_nominal_l=in20_crossings(jset(dj_dip), L_NOM),
                         peaks_true_l=in20_crossings(jset(dj_dip), L_TRUE))
    json.dump(out, open('calibration.json', 'w'), indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))
