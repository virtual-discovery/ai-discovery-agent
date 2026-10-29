"""Shared construction utilities: frame selection and equation-of-state / enthalpy tools."""
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, curve_fit

POOL = json.load(gzip.open(Path.home()/'.cache/ai-theorist/ge-pressure-polymorph-pathways/pool.json.gz'))
PHASES = ['CD', 'bSn', 'ST12', 'BC8', 'R8']
METASTABLE = ['ST12', 'BC8', 'R8']
EV_A3_GPA = 160.21766208
# Near-equilibrium (hydrostatic, symmetric) frames: deviatoric normal stress < 0.5 GPa, shear stress < 0.05 GPa,
# max force < 0.6 eV/A. Rattled training configurations break the symmetry (shear stress, forces 1-3 eV/A);
# isotropically scaled metastable structures keep zero shear but carry internal forces up to ~0.5 eV/A.
DEV, SHEAR, FMAX = 0.5, 0.05, 0.6


def all_frames(phase):
    out = []
    for key, v in POOL.items():
        parts = key.split('/')
        if len(parts) != 3 or parts[0] != 'training_set' or parts[1] != phase:
            continue
        s = v['settings']
        n = s['NIONS']
        for j, f in enumerate(v['frames']):
            st = np.array(f['stress_kB']) / 10
            p = st[:3].mean()
            out.append({'key': key, 'j': j, 'tag': f'{key}#{j}', 'nf': len(v['frames']), 'n': n, 'V': f['volume'] / n,
                        'E': f['energy'] / n, 'P': float(p), 'dev': float(np.abs(st[:3] - p).max()),
                        'sh': float(np.abs(st[3:]).max()), 'F': f['max_force'], 'pst': s['PSTRESS'] / 10,
                        'isif': s['ISIF'], 'acc': v['reached_accuracy']})
    return out


def clean_frames(phase, dev=DEV, shear=SHEAR, fmax=FMAX):
    return [o for o in all_frames(phase) if o['dev'] < dev and o['sh'] < shear and o['F'] < fmax]


def bm3(V, E0, V0, B0, Bp):
    x = (V0 / V) ** (2 / 3)
    return E0 + 9 * V0 * B0 / 16 * ((x - 1) ** 3 * Bp + (x - 1) ** 2 * (6 - 4 * x))


def p_bm3(V, E0, V0, B0, Bp):
    x = (V0 / V) ** (2 / 3)
    return 1.5 * B0 * (x ** 3.5 - x ** 2.5) * (1 + 0.75 * (Bp - 4) * (x - 1))


def vinet(V, E0, V0, B0, Bp):
    x = (V / V0) ** (1 / 3)
    eta = 1.5 * (Bp - 1)
    return E0 + 2 * B0 * V0 / (Bp - 1) ** 2 * (2 - (5 + 3 * Bp * (x - 1) - 3 * x) * np.exp(-eta * (x - 1)))


def p_vinet(V, E0, V0, B0, Bp):
    x = (V / V0) ** (1 / 3)
    eta = 1.5 * (Bp - 1)
    return 3 * B0 * (1 - x) / x ** 2 * np.exp(eta * (1 - x))


FORMS = {'bm3': (bm3, p_bm3), 'vinet': (vinet, p_vinet)}


def fit(V, E, form='bm3'):
    V, E = np.asarray(V, float), np.asarray(E, float)
    i = int(np.argmin(E))
    p, _ = curve_fit(FORMS[form][0], V, E, p0=[E[i], V[i], 0.4, 4.5], maxfev=20000)
    return p


def enthalpy(p, P, form='bm3'):
    f, pf = FORMS[form]
    Pe = P / EV_A3_GPA
    V = brentq(lambda v: pf(v, *p) - Pe, 0.6 * p[1], 1.4 * p[1])
    return f(V, *p) + Pe * V, V


def crossing(pa, pb, form='bm3', lo=-2.0, hi=20.0):
    def g(P):
        try:
            return enthalpy(pa, P, form)[0] - enthalpy(pb, P, form)[0]
        except ValueError:
            return np.nan
    Ps = np.linspace(lo, hi, 221)
    vals = [g(P) for P in Ps]
    for i in range(len(Ps) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] * vals[i + 1] < 0:
            return float(brentq(g, Ps[i], Ps[i + 1]))
    return float('nan')
