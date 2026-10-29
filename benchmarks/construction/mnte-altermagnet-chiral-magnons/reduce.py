"""Reduction of the archived MnTe data (IN20 raw scans, HYSPEC cuts) and the reference peak fits."""
import glob
import os
import numpy as np
from scipy.optimize import least_squares
import ill

CACHE = os.environ.get('MNTE_CACHE', os.path.expanduser('~/.cache/ai-theorist/mnte-altermagnet-chiral-magnons'))
X = os.path.join(CACHE, 'x')


def in20_points(folder):
    """All points of the constant-E scans in Figure3/<folder>: (QH_nominal, I9 sign, I10 sign, counts, monitor)."""
    rows = []
    for f in sorted(glob.glob(os.path.join(X, 'Figure3', folder, '*'))):
        d = ill.parse(f)
        c, a = d['cols'], d['data']
        for r in a:
            rows.append((r[c.index('QH')], np.sign(r[c.index('I9')]), np.sign(r[c.index('I10')]), r[c.index('CNTS')], r[c.index('M1')],
                         os.path.basename(f)))
    return rows


def binned(rows, chan, step=0.005, lo=1.0, hi=1.55):
    """Monitor-weighted rebinning on the nominal QH grid; chan = (s9, s10)."""
    edges = np.arange(lo - step / 2, hi + step, step)
    cnt = np.zeros(len(edges) - 1)
    mon = np.zeros_like(cnt)
    hs = np.zeros_like(cnt)
    for (qh, s9, s10, c, m, f) in rows:
        if (s9, s10) != chan:
            continue
        k = np.searchsorted(edges, qh) - 1
        cnt[k] += c
        mon[k] += m
        hs[k] += qh * m
    ok = mon > 0
    x = hs[ok] / mon[ok]
    y = cnt[ok] / mon[ok] * 1e6
    e = np.sqrt(np.maximum(cnt[ok], 1)) / mon[ok] * 1e6
    return x, y, e


def two_gauss(p, x):
    b, a1, x1, a2, x2, w = p
    s = w / 2.3548
    return b + a1 * np.exp(-0.5 * ((x - x1) / s) ** 2) + a2 * np.exp(-0.5 * ((x - x2) / s) ** 2)


def fit_two_gauss(x, y, e, guess, slope=False):
    f = (lambda p, x: two_gauss(p[:6], x) + p[6] * (x - x.mean())) if slope else two_gauss
    p0 = list(guess) + ([0.0] if slope else [])
    r = least_squares(lambda p: (f(p, x) - y) / e, p0)
    J = r.jac
    cov = np.linalg.pinv(J.T @ J) * max(1.0, 2 * r.cost / max(1, len(x) - len(p0)))
    return r.x, np.sqrt(np.diag(cov)), 2 * r.cost / (len(x) - len(p0))
