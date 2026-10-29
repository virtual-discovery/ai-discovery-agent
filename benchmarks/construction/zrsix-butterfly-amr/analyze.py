"""Exploration helpers reading the raw caches (same formulas as the lab)."""
import numpy as np, sys, glob, os
sys.path.insert(0, '.')
import tb
from common import WINDOW, XGRID, SUBMESH, dkey
HBAR = 6.582119569e-16; E = 1.602176634e-19; KB = 8.617333262e-5
_vol = {}

def vol(c):
    if c not in _vol:
        _vol[c] = abs(np.linalg.det(tb.Model(f'/data/model/{c}.npz', c).A))
    return _vol[c]

def sigma(c, plane, ang, mesh, mu, T=30, bands=None):
    z = np.load(f'/data/cache/{c}/{dkey(plane, ang)}.npz')
    mi = list(z['mesh']).index(mesh)
    cent = float(z['lo']) + np.arange(z['K'].shape[2]) * float(z['ebin'])
    kT = KB * T
    x = np.clip((cent - mu) / kT, -60, 60); f = 1 / (np.exp(x) + 1); w = f * (1 - f) / kT
    nk = mesh * mesh * (mesh // 2)
    pref = 2 * E ** 2 / (vol(c) * 1e-30 * nk) / E * (1e-10 / HBAR) ** 2 * 1e-12
    nb = z['K'].shape[1]
    bands = range(nb) if bands is None else bands
    s = sum(pref * np.einsum('e,exij->xij', w, z['K'][mi, b].astype(float)) for b in bands)
    O = z['O'][mi]; Wt = z['W'][mi]
    return s, float((w * O.sum(0)).sum() / (w * Wt.sum(0)).sum())

def rho(*a, **k):
    s, o = sigma(*a, **k)
    return np.linalg.inv(s), o

def mr(r, i=0, j=0):
    return r[:, i, j] / r[0, i, j] - 1
