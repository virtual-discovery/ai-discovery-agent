"""Readers and estimators for the BaB2 task construction (numpy only; shared by make_task.py).

The isotropic Eliashberg solver is the lab's own module (environment/lab/eliashberg.py), loaded without
writing bytecode, so that the reference values come from the same code path the agent pays for.
"""
import importlib.util
import json
import re
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/bab2-hole-doping-two-gap'
POOL = json.loads((HERE/'pool.json').read_text())
RUNS, SUPER = POOL['runs'], POOL['supercells']
README_ORDER = [0.0, 0.02, 0.04, 0.06, 0.065, 0.08, 0.1]      # doping list in the archive README, in its order
WITHHELD_DFPT = ['3', '5']                                    # folders whose DFPT results the lab never serves
KB = 0.08617333262

sys.dont_write_bytecode = True
_spec = importlib.util.spec_from_file_location('eliashberg', TASK/'environment/lab/eliashberg.py')
EL = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(EL)


def doping(folder):
    return RUNS[folder]['tot_charge']


def folder_of(x):
    return next(f for f, r in RUNS.items() if abs(r['tot_charge'] - x) < 1e-9)


def soft_mode(entry):
    """Lowest harmonic frequency over the DFPT q-points, excluding the three acoustic modes at Gamma."""
    best = (np.inf, None)
    for d in entry['dfpt']:
        f = sorted(d['frequencies_cm1'])
        if np.allclose(d['q_cartesian_2pi_over_alat'], 0):
            f = f[3:]
        if f[0] < best[0]:
            best = (f[0], d['q_cartesian_2pi_over_alat'])
    return best


def signed_square(w):
    return np.sign(w) * np.asarray(w) ** 2


def crossing(xs, ys):
    """First zero crossing of piecewise-linear y(x) (xs sorted)."""
    for (x0, y0), (x1, y1) in zip(zip(xs, ys), zip(xs[1:], ys[1:])):
        if y0 < 0 <= y1:
            return x0 + (x1 - x0) * (-y0) / (y1 - y0)
    return float('nan')


def poly_root(xs, ys, deg):
    p = np.polyfit(xs, ys, deg)
    roots = [r.real for r in np.roots(p) if abs(r.imag) < 1e-12 and 0 < r.real < 0.2]
    return min(roots), p


def a2f_arrays(folder, column=1):
    rows = np.array(RUNS[folder]['a2f']['rows'])
    return rows[:, 0], rows[:, column]


def allen_dynes(w, a2f, mu, f1f2=True):
    dw = w[1] - w[0]
    lam = 2 * np.sum(a2f / w) * dw
    wlog = np.exp(2 / lam * np.sum(a2f / w * np.log(w)) * dw)
    w2 = np.sqrt(2 / lam * np.sum(a2f * w) * dw)
    tc = wlog / 1.2 * np.exp(-1.04 * (1 + lam) / (lam - mu * (1 + 0.62 * lam))) / KB
    if f1f2:
        f1 = (1 + (lam / (2.46 * (1 + 3.8 * mu))) ** 1.5) ** (1 / 3)
        f2 = 1 + (w2 / wlog - 1) * lam ** 2 / (lam ** 2 + (1.82 * (1 + 6.3 * mu) * w2 / wlog) ** 2)
        tc *= f1 * f2
    return float(tc), float(lam), float(wlog), float(w2)


def iso_archived(folder):
    """(T, gap at w_0, converged) from the archived iso.out: the last iteration line before each convergence line."""
    out, last, T = [], None, None
    for line in RUNS[folder]['iso_log'].splitlines():
        m = re.match(r'\s*temp\(\s*\d+\)\s*=\s*([\d.]+) K', line)
        if m:
            T = float(m.group(1))
        p = line.split()
        if len(p) == 4 and p[0].isdigit():
            last = p
        if 'Convergence was' in line and 'Pade' not in line and 'approximants' not in line:
            ok = 'not reached' not in line
            gap = float(last[3]) if last[3] != 'NaN' else float('nan')
            out.append((T, gap, ok))
    return out


def aniso_series(folder, axis='imaginary'):
    """Per archived temperature: max gap (log), free energy (log), weighted mean of the small- and large-gap clusters."""
    res = []
    for key, a in sorted(RUNS[folder]['aniso'].items(), key=lambda kv: float(kv[0])):
        T = float(key)
        mm = [l for l in a['log'] if 'Min. / Max.' in l]
        line = mm[0] if axis == 'imaginary' else mm[1]
        dmax = float(line.split('=')[1].split()[1])
        F = float(re.search(r'Free energy =\s*([-\d.E+]+)', next(l for l in a['log'] if 'Free energy' in l)).group(1))
        small, large = clusters(a[axis]['rows'], T)
        res.append({'T': T, 'dmax': dmax, 'F': F, 'small': small, 'large': large})
    return res


def clusters(rows, T, estimator='mean'):
    """Split the archived gap histogram (col 1 = T + weight, col 2 = gap in meV) at its widest empty interval."""
    a = np.array(rows)
    wgt, gap = a[:, 0] - T, a[:, 1]
    keep = wgt > 0
    gap, wgt = gap[keep], wgt[keep]
    o = np.argsort(gap)
    gap, wgt = gap[o], wgt[o]
    k = int(np.argmax(np.diff(gap)))
    parts = [(gap[:k + 1], wgt[:k + 1]), (gap[k + 1:], wgt[k + 1:])]
    if estimator == 'mean':
        return [float(np.sum(g * w) / np.sum(w)) for g, w in parts]
    if estimator == 'median':
        out = []
        for g, w in parts:
            c = np.cumsum(w) / np.sum(w)
            out.append(float(g[np.searchsorted(c, 0.5)]))
        return out
    if estimator == 'unweighted':
        return [float(np.mean(g)) for g, w in parts]
    if estimator == 'mode':
        return [float(g[np.argmax(w)]) for g, w in parts]
    raise ValueError(estimator)


def zero_of_square(T, y):
    """Two-point linear extrapolation of y^2 to zero from the last two temperatures."""
    y2 = np.asarray(y[-2:]) ** 2
    return float(T[-1] + (T[-1] - T[-2]) * y2[-1] / (y2[-2] - y2[-1]))


def zero_linear(T, y):
    return float(T[-1] + (T[-1] - T[-2]) * y[-1] / (y[-2] - y[-1]))
