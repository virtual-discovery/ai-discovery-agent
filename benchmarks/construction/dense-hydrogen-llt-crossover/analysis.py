"""Reference protocol for dense-hydrogen-llt-crossover, applied to the served (pooled) data.

c_v/k_B = var(E_tot)/(N k_B^2 T^2) over t >= t0 of a run (served series: every 0.1 ps).
Widom line: rs of the c_v maximum along an isotherm (quadratic fit through the 5 points around the largest c_v).
Crossover pressure: pressure at the flattest point of P(rs) (largest dP/drs, i.e. maximal compressibility), from a
quadratic fit of the finite-difference slopes around their maximum.
"""
import collections
import json
import os
from pathlib import Path
import numpy as np

KB = 8.617333262e-5
CACHE = Path(os.environ.get('HCACHE', os.path.expanduser('~/.cache/ai-theorist/dense-hydrogen-llt-crossover')))
TH = np.load(CACHE/'pool_thermo.npz')
PI = np.load(CACHE/'pool_pimd.npz')
META = json.loads((CACHE/'pool_meta.json').read_text())


def series(n, t, rs):
    m = META['runs'][f'{n}/{t}/{rs}']
    a = TH[f'{n}/{t}/{rs}'].astype(float)
    time = (m['first_step'] + m['stride_steps'] * np.arange(len(a))) * 2e-4
    return time, a[:, 0], a[:, 1] + m['E_ref_eV'], a[:, 2] + m['E_ref_eV'], a[:, 3] / 1e4   # ps, K, eV, eV, GPa


def rs_list(n, t):
    return sorted([m['rs'] for m in META['runs'].values() if m['N'] == n and m['T'] == t], key=float)


def cv(n, t, rs, t0=10.0, t1=None):
    time, _, _, e, _ = series(n, t, rs)
    sel = (time >= t0) & ((time <= t1) if t1 is not None else True)
    return float(np.var(e[sel]) / (n * KB * KB * t * t)) if sel.sum() >= 50 else float('nan')


def pressure(n, t, rs, t0=10.0):
    time, _, _, _, p = series(n, t, rs)
    return float(p[time >= t0].mean())


def peak(xs, ys, npts=5):
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    ok = np.isfinite(ys)
    xs, ys = xs[ok], ys[ok]
    i = int(np.argmax(ys))
    lo = max(0, min(i - npts // 2, len(xs) - npts))
    hi = min(len(xs), lo + npts)
    c = np.polyfit(xs[lo:hi], ys[lo:hi], 2)
    if c[0] < 0 and xs[lo] <= -c[1] / (2 * c[0]) <= xs[hi - 1]:
        return float(-c[1] / (2 * c[0]))
    return float(xs[i])


def widom_rs(t, n=2048, t0=10.0, npts=5):
    rr = rs_list(n, t)
    return peak([float(r) for r in rr], [cv(n, t, r, t0) for r in rr], npts)


def flattest(x, p, npts=5):
    x, p = np.asarray(x, float), np.asarray(p, float)
    s = np.diff(p) / np.diff(x)
    xm = 0.5 * (x[1:] + x[:-1])
    xs = peak(xm, s, npts)
    return xs, float(np.interp(xs, x, p))


def classical_crossover(t, n=2048, t0=10.0):
    rr = rs_list(n, t)
    return flattest([float(r) for r in rr], [pressure(n, t, r, t0) for r in rr])


def pimd_eos(label='1000K', skip=0.2):
    d = collections.defaultdict(list)
    for m in META['pimd'].values():
        if m['T_label'] == label:
            a = PI[m['data_key']]
            d[float(m['rs'])].append(float(a[int(len(a) * skip):, 5].mean()))
    rs = sorted(d)
    return np.array(rs), np.array([np.mean(d[r]) for r in rs])


def pimd_crossover(label='1000K', skip=0.2):
    return flattest(*pimd_eos(label, skip))


def pimd_tanh_crossover(label='1000K'):
    from scipy.optimize import curve_fit
    x, p = pimd_eos(label)

    def f(r, a, b, c, d, r0, w):
        return a + b * (r - 1.47) + c * (r - 1.47) ** 2 + d * np.tanh((r - r0) / w)
    par = curve_fit(f, x, p, p0=[145, -500, 0, -3, 1.49, 0.005], maxfev=20000)[0]
    return float(par[4]), float(f(par[4], *par))
