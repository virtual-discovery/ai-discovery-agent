"""Shared thermodynamics for sngeo2-coherent-miscibility (numpy only).

Energies: eV per rutile primitive cell (2 cations) unless stated; x = Ge fraction.
"""
import gzip, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
KB = 8.617333e-5
_POOL = None


def pool():
    global _POOL
    if _POOL is None:
        with gzip.open(HERE/'pool.json.gz', 'rt') as f:
            _POOL = json.load(f)
    return _POOL


# ---------------- SQS (disordered alloy) ----------------
def sqs_rows(fn):
    out = []
    for r in pool()['sqs'][fn]:
        if 0 < r['x_Ge'] < 1:
            A = np.array(r['lattice_A'])
            out.append({'x': r['x_Ge'], 'index': r['index'], 'e_cell': r['energy_as_released'],
                        'vol_per_cation': abs(np.linalg.det(A)) / 24})
    return out


def end_volumes(fn):
    v = {}
    for r in pool()['sqs'][fn]:
        if r['x_Ge'] in (0.0, 1.0):
            v[r['x_Ge']] = abs(np.linalg.det(np.array(r['lattice_A']))) / 2
    return v


def vegard_outliers(rows, fn, tol=0.3):
    """flag SQS whose volume per cation deviates from the end-member Vegard line by > tol A^3"""
    v = end_volumes(fn)
    return [abs(r['vol_per_cation'] - (r['x'] * v[1.0] + (1 - r['x']) * v[0.0])) > tol for r in rows]


def sqs_points(fn, mode='clean', indices=None):
    rows = sqs_rows(fn)
    bad = vegard_outliers(rows, fn)
    keep = [(r, b) for r, b in zip(rows, bad) if indices is None or (r['x'], r['index']) in indices]
    if mode == 'all':
        sel = [r for r, b in keep]
    else:
        sel = [r for r, b in keep if not b]
    x = np.array([r['x'] for r in sel]); y = np.array([r['e_cell'] for r in sel])
    if mode in ('clean', 'all'):
        return x, y
    u = np.unique(x)
    if mode == 'cleanmean':
        return u, np.array([y[x == q].mean() for q in u])
    if mode == 'cleanmin':
        return u, np.array([y[x == q].min() for q in u])
    raise ValueError(mode)


def subregular(x, y):
    M = np.c_[x * (1 - x), x * (1 - x) * (1 - 2 * x)]
    return np.linalg.lstsq(M, y, rcond=None)[0]


# ---------------- coherency strain ----------------
def scan_key(plane, A):
    A = np.array(A); c = A[0, 2]; v2, v3 = A[1], A[2]
    if plane == '001':
        return (round(float(np.linalg.norm(v2)), 4),), float(c)
    if plane == '100':
        return (round(float(np.linalg.norm(v2)), 4), round(float(c), 4)), float(np.linalg.norm(v3))
    return (round(float(np.linalg.norm(v2 - v3)), 4), round(float(c), 4)), float(np.linalg.norm(v2 + v3))


_SCANS = {}


def scans(fn, plane, member):
    if (fn, plane, member) in _SCANS:
        return _SCANS[(fn, plane, member)]
    g = {}
    for r in pool()['strain'][f'{fn}|{plane}|{member}']:
        k, oop = scan_key(plane, r['lattice_A'])
        g.setdefault(k, []).append((oop, r['energy_eV']))
    _SCANS[(fn, plane, member)] = {k: np.array(sorted(v)) for k, v in g.items()}
    return _SCANS[(fn, plane, member)]


def scan_min(arr, method='poly5'):
    i = int(np.argmin(arr[:, 1]))
    if method == 'raw' or i in (0, len(arr) - 1):
        return float(arr[i, 1])
    npt = int(method[4:])
    lo = max(0, min(i - npt // 2, len(arr) - npt)); s = arr[lo:lo + npt]
    p = np.polyfit(s[:, 0], s[:, 1], 2)
    return float(np.polyval(p, -p[1] / (2 * p[0])))


def strain_curve(fn, plane, method='poly5', xsel=None, mode='ok'):
    """Delta e_cs(x) per cell on the Vegard line of the scan grid (end-member minima as anchors)."""
    S = {m: scans(fn, plane, m) for m in ('GeO2', 'SnO2')}
    loc = {m: min(((v[:, 1].min(), k) for k, v in S[m].items()))[1] for m in S}
    E0 = {m: min(v[:, 1].min() for v in S[m].values()) for m in S}
    kG, kS = np.array(loc['GeO2']), np.array(loc['SnO2'])
    keys = sorted(S['GeO2'])
    xs, ec = [], []
    for k in keys:
        kk = np.array(k)
        x = float(np.dot(kk - kS, kG - kS) / np.dot(kG - kS, kG - kS))
        if np.linalg.norm(kk - (x * kG + (1 - x) * kS)) > 2e-3:
            continue
        if xsel is not None and not any(abs(x - q) < 1e-3 for q in xsel):
            continue
        e = 0.0
        for m, w in (('GeO2', x), ('SnO2', 1 - x)):
            kk2 = k
            if mode == 'aonly' and len(k) == 2:      # wrong path: only the first in-plane parameter follows the alloy
                kk2 = min(S[m], key=lambda q: abs(q[0] - k[0]) + abs(q[1] - loc[m][1]))
            arr = S[m][kk2]
            if mode == 'norelax':                    # wrong path: out-of-plane fixed at the member's bulk value
                own = S[m][loc[m]]; oop0 = own[np.argmin(own[:, 1]), 0]
                en = float(arr[np.argmin(abs(arr[:, 0] - oop0)), 1])
            else:
                en = scan_min(arr, method)
            e += w * (en - E0[m])
        xs.append(x); ec.append(e)
    return np.array(xs), np.array(ec)


def strain_fit(fn, plane, **kw):
    xs, ec = strain_curve(fn, plane, **kw)
    m = (xs > 1e-6) & (xs < 1 - 1e-6)
    return subregular(xs[m], ec[m])


def strain_at(fn, plane, x0=0.5, **kw):
    xs, ec = strain_curve(fn, plane, **kw)
    i = int(np.argmin(abs(xs - x0))); assert abs(xs[i] - x0) < 1e-3
    return float(ec[i])


# ---------------- free-energy model (per cell, 2 cations) ----------------
def spinodal_Tc(a, b):
    x = np.linspace(1e-5, 1 - 1e-5, 400001)
    hpp = -2 * a + b * (-6 + 12 * x)
    T = -hpp * x * (1 - x) / (2 * KB)
    i = int(np.argmax(T))
    return float(max(T[i], 0.0)), float(x[i])


def spinodal_at(a, b, T):
    x = np.linspace(1e-6, 1 - 1e-6, 400001)
    hpp = -2 * a + b * (-6 + 12 * x)
    inside = x[-hpp * x * (1 - x) / (2 * KB) > T]
    return (float(inside.min()), float(inside.max())) if len(inside) else None


def g_cell(x, a, b, T):
    x = np.clip(x, 1e-12, 1 - 1e-12)
    return a * x * (1 - x) + b * x * (1 - x) * (1 - 2 * x) + 2 * KB * T * (x * np.log(x) + (1 - x) * np.log(1 - x))


def dg_cell(x, a, b, T):
    return a * (1 - 2 * x) + b * (1 - 6 * x + 6 * x * x) + 2 * KB * T * np.log(x / (1 - x))


def binodal_at(a, b, T):
    """common tangent: lower convex hull of g(x) on a grid, refined by Newton (scipy fsolve)"""
    from scipy.spatial import ConvexHull
    from scipy.optimize import fsolve
    x = np.unique(np.concatenate([np.logspace(-15, -2, 400), np.linspace(1e-2, 1 - 1e-2, 20001), 1 - np.logspace(-2, -15, 400)]))
    g = g_cell(x, a, b, T)
    pts = np.c_[x, g]
    hull = ConvexHull(pts)
    lower = set()
    for simplex, eq in zip(hull.simplices, hull.equations):
        if eq[1] < 0:            # outward normal points down: lower hull facet
            lower.add(tuple(sorted(simplex)))
    best = None
    for i, j in lower:
        if x[j] - x[i] > 2e-3 and (best is None or x[j] - x[i] > best[1] - best[0]):
            best = (x[i], x[j])
    if best is None:
        return None
    def eqs(z):
        x1, x2 = 1 / (1 + np.exp(-z[0])), 1 / (1 + np.exp(-z[1]))
        d1, d2 = dg_cell(x1, a, b, T), dg_cell(x2, a, b, T)
        return [d1 - d2, (g_cell(x2, a, b, T) - g_cell(x1, a, b, T)) - d1 * (x2 - x1)]
    z = fsolve(eqs, [np.log(best[0] / (1 - best[0])), np.log(best[1] / (1 - best[1]))], xtol=1e-12)
    x1, x2 = sorted(1 / (1 + np.exp(-z)))
    if x2 - x1 < 1e-4:
        return None
    return float(x1), float(x2)


def binodal_Tc(a, b, Tlo=50, Thi=6000):
    if binodal_at(a, b, Tlo) is None:
        return 0.0
    for _ in range(25):
        T = (Tlo + Thi) / 2
        if binodal_at(a, b, T) is None:
            Thi = T
        else:
            Tlo = T
    return (Tlo + Thi) / 2


# ---------------- full analysis (one configuration of choices) ----------------
DEFAULT = {'sqs_mode': {'PBE': 'clean', 'r2SCAN': 'clean'}, 'sqs_indices': {'PBE': None, 'r2SCAN': None},
           'sqs_factor': 1.0,          # 2.0 = released per-cell values read as per-cation values
           'strain_factor': 1.0,       # 2.0 = per-6-atom-cell strain energies mixed with per-cation enthalpies
           'strain_method': 'poly5', 'strain_xsel': None, 'strain_mode': 'ok',
           'coherent_plane': 'min',    # or a plane label (wrong path: forced plane)
           'strain_functional': {'PBE': 'PBE', 'r2SCAN': 'r2SCAN'},
           'no_coherency': False,      # wrong path: incoherent model reported as coherent
           'swap_700K': False, 'incoherent_Tc': None, 'dh05_from': 'fit'}


def analysis(**changes):
    cfg = {**DEFAULT, **changes}
    out = {}
    fits = {}
    for fn in ('PBE', 'r2SCAN'):
        x, y = sqs_points(fn, cfg['sqs_mode'][fn], cfg['sqs_indices'][fn])
        a, b = subregular(x, cfg['sqs_factor'] * y)
        fits[fn] = (a, b, x, cfg['sqs_factor'] * y)
    a, b, x, y = fits['PBE']
    if cfg['dh05_from'] == 'fit':
        dh05 = a / 4 / 2 * 1000
    else:                                   # mean of the SQS at x = 0.5
        dh05 = float(y[np.abs(x - 0.5) < 1e-9].mean()) / 2 * 1000
    out['mixing_enthalpy_x050_meV_per_cation'] = dh05
    sf = cfg['strain_factor']
    kw = dict(method=cfg['strain_method'], mode=cfg['strain_mode'])
    strain05 = {pl: sf * strain_at('PBE', pl, 0.5, **kw) / 2 * 1000 for pl in ('001', '100', '101')}
    out['coherency_strain_energy_x050_meV_per_cation'] = strain05
    habit = min(strain05, key=strain05.get) if cfg['coherent_plane'] == 'min' else cfg['coherent_plane']
    out['habit_plane'] = habit
    out['incoherent_Tc_K'] = cfg['incoherent_Tc'] if cfg['incoherent_Tc'] is not None else binodal_Tc(a, b)
    tcs, coh = {}, {}
    for fn in ('PBE', 'r2SCAN'):
        a, b = fits[fn][:2]
        if cfg['no_coherency']:
            ac = bc = 0.0
        else:
            # habit plane of each functional: the plane with the smallest strain energy for that functional's data
            sfn = cfg['strain_functional'][fn]
            if cfg['coherent_plane'] == 'min':
                pl = min(('001', '100', '101'), key=lambda p: strain_at(sfn, p, 0.5, **kw))
            else:
                pl = cfg['coherent_plane']
            ac, bc = sf * strain_fit(sfn, pl, xsel=cfg['strain_xsel'], **kw)
        coh[fn] = (a - ac, b - bc)
        tcs[fn] = spinodal_Tc(*coh[fn])
    out['coherent_spinodal_Tc_K'] = {fn: tcs[fn][0] for fn in tcs}
    bi = binodal_at(*coh['PBE'], 700.0); sp = spinodal_at(*coh['PBE'], 700.0)
    bi = list(bi) if bi else None; sp = list(sp) if sp else None
    if cfg['swap_700K']:
        bi, sp = sp, bi
    out['coherent_boundaries_700K_PBE'] = {'binodal': bi, 'spinodal': sp}
    return out
