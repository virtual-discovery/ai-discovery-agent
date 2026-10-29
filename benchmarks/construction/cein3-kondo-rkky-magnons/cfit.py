"""Data sets, variable-projection fits and the grading metric for the CeIn3 construction (numpy/scipy only).

Every data set is compared with a model after a per-set non-negative scale and a constant offset, both solved
by weighted linear least squares; this is also how task_lab.grade scores the withheld predictions.
"""
import os
import numpy as np
from scipy.optimize import least_squares
import cmodel as M

RAW = os.environ.get('CEIN3_RAW', os.path.expanduser('~/.cache/ai-theorist/cein3-kondo-rkky-magnons/raw')) + '/'
L = lambda n: np.loadtxt(RAW + n)
HEF = {'RG': 'HighEnergy_Int_RG.dat', 'GX': 'HighEnergy_Int_GX.dat', 'XM': 'HighEnergy_IntXM.dat', 'MG': 'HighEnergy_Int_MG.dat'}
HR_EMIN, HE_EMIN = 0.25, 0.3          # energies below these are dominated by the (subtracted) elastic line


def sym_sigma(I, Q, emask):
    """Noise of a high-resolution map from the difference of the +q and -q rows (equivalent by symmetry)."""
    pairs = [(i, int(np.argmin(abs(Q + Q[i])))) for i in range(len(Q)) if Q[i] > 0]
    d = np.array([I[i] - I[j] for i, j in pairs if abs(Q[j] + Q[i]) < 1e-6])[:, emask]
    return float(np.nanstd(d) / np.sqrt(2))


class DS:
    pass


def he_ds(path, rows=None, perp_scale=1.0):
    d = DS(); d.kind, d.name = 'he', 'he' + path
    I, Q, E = L(HEF[path]), L('HighEnergy_Q_%s.dat' % path), L('HighEnergy_E_%s.dat' % path)
    idx = np.arange(len(Q)) if rows is None else np.array(rows)
    y2 = I[idx]; mask = np.isfinite(y2) & (E[None, :] > HE_EMIN)
    keep = mask.any(1)
    d.idx, d.mask = idx[keep], mask[keep]
    d.obj = M.HERows(Q[d.idx], E, perp_scale=perp_scale)
    d.y = y2[keep][d.mask]; d.e = np.sqrt(0.10 ** 2 + (0.15 * np.abs(d.y)) ** 2)
    return d


def hr_ds(path, rows=None, half_perp=0.08, emin=HR_EMIN):
    d = DS(); d.kind, d.name = 'hr', 'hr' + path
    I, Q, E = L('HighResolution_Int_%s.dat' % path), L('HighResolution_Q_%s.dat' % path), L('HighResolution_E_%s.dat' % path)
    sig = sym_sigma(I, Q, E > emin)
    d.idx = np.arange(len(Q)) if rows is None else np.array(rows)
    d.mask = np.isfinite(I[d.idx]) & (E[None, :] > emin)
    d.obj = M.HRRows(path, Q[d.idx], E, half_perp=half_perp)
    d.y = I[d.idx][d.mask]; d.e = np.full(len(d.y), sig); d.sigma = sig
    return d


def cut_ds(path, tag, half_perp=0.08):
    d = DS(); d.kind, d.name = 'cut', 'cut' + path + tag
    q, d.y, d.e = L('HighResolution_%s_ConstEcut_%s.dat' % (path, tag)).T
    d.q, d.mask = q, None
    d.obj = M.HRCut(path, q, float(tag.replace('p', '.').replace('meV', '')), half_perp=half_perp)
    return d


def integrated_ds(perp_scale=1.0):
    d = DS(); d.kind, d.name, d.mask = 'int', 'integrated', None
    d.s, d.y, d.e = L('Experiment_IntegralOverchi.dat').T
    d.obj = M.Integrated(d.s, perp_scale=perp_scale)
    return d


def predict(d, p):
    m = d.obj.model(p)
    return m[d.mask] if d.mask is not None else m


def project(m, y, e, offset=True):
    w = 1 / e ** 2
    if offset:
        A = np.c_[m, np.ones_like(m)] * np.sqrt(w)[:, None]
        c = np.linalg.lstsq(A, y * np.sqrt(w), rcond=None)[0]
        if c[0] < 0:
            c = np.array([0.0, np.sum(w * y) / np.sum(w)])
        return c[0] * m + c[1], c
    a = max(np.sum(w * m * y) / max(np.sum(w * m * m), 1e-300), 0.0)
    return a * m, np.array([a, 0.0])


def chi2r(m, y, e, offset=True):
    f, _ = project(m, y, e, offset)
    return float(np.sum(((f - y) / e) ** 2) / (len(y) - (2 if offset else 1)))


NAMES = ['A', 'kappa', 'gamma_hr', 'sigma_he', 'gamma_he', 'gap']
BND = {'A': (0, 30), 'kappa': (0.015, 0.12), 'gamma_hr': (0, 1.5), 'sigma_he': (0.05, 0.6), 'gamma_he': (0, 1.5), 'gap': (0, 1.0)}


def to_p(x, free, fixed, N):
    p = dict(fixed); p.update(dict(zip(free, x)))
    p['SJ'] = [float(p.pop(f'SJ{i+1}')) for i in range(N)]
    return {k: (float(v) if k != 'SJ' else v) for k, v in p.items()}


def fit(dsets, N, x0, fixed, nfev=200):
    """Least squares over the non-linear parameters; each data set enters with equal total weight."""
    names = [f'SJ{i+1}' for i in range(N)] + NAMES
    free = [n for n in names if n not in fixed]
    lo = np.array([BND.get(n, (-3, 3))[0] for n in free]); hi = np.array([BND.get(n, (-3, 3))[1] for n in free])
    xs = np.clip(np.array([x0.get(n, 0.0) for n in free], float), lo + 1e-6, hi - 1e-6)

    def res(x):
        p = to_p(x, free, fixed, N); out = []
        for d in dsets:
            f, _ = project(predict(d, p), d.y, d.e)
            out.append((f - d.y) / d.e / np.sqrt(len(d.y)))
        return np.concatenate(out)
    r = least_squares(res, xs, bounds=(lo, hi), x_scale='jac', max_nfev=nfev, diff_step=2e-2)
    return to_p(r.x, free, fixed, N), float(2 * r.cost)


def score_all(dsets, p, offset=True):
    return {d.name: round(chi2r(predict(d, p), d.y, d.e, offset), 4) for d in dsets}
