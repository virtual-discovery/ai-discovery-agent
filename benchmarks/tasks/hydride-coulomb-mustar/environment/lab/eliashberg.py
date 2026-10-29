"""Linearized isotropic Migdal-Eliashberg Tc on the Matsubara axis (numpy), IsoME conventions.

Modes (all energies in meV, DOS per spin per meV, W in meV):
  cdos + mu*   constant DOS, Morel-Anderson mu*_ME referenced to the Matsubara cutoff omega_c;
  vdos + mu*   full-bandwidth with energy-dependent DOS N(eps), energy shift chi, chemical-potential update;
  vdos + W     full-bandwidth with the static screened interaction W(eps, eps') (phi = phi_ph(iw) + phi_c(eps)).
The equations follow the IsoME formulation (Kogler et al., CPC 2025; Pickett PRB 26, 1186; Lee et al.
npj Comput. Mater. 9, 156): Matsubara frequencies up to omega_c, energy grid of 1/5/50 meV steps inside
+-100/+-500/+-encut meV, normal-state Z, chi and mu solved self-consistently, then the largest eigenvalue
of the linearized gap equation; Tc is where it equals 1.
"""
import math
import numpy as np
from scipy.optimize import brentq

KB = 0.08617333262          # meV/K
THZ2MEV = 4.13566553853599


def trapz_w(x):
    w = np.empty_like(x)
    w[0] = (x[1] - x[0]) / 2
    w[-1] = (x[-1] - x[-2]) / 2
    w[1:-1] = (x[2:] - x[:-2]) / 2
    return w


class A2F:
    """alpha^2F on a uniform grid from 0.01 meV to omega_max (as IsoME resamples it)."""

    def __init__(self, omega_meV, a2f):
        w = np.asarray(omega_meV, float)
        a = np.asarray(a2f, float)
        grid = np.linspace(1e-2, w[-1], len(w))
        f = np.interp(grid, w, a, left=0.0, right=0.0)
        f[f < 0] = 0.0
        self.w, self.a = grid, f
        self.tw = trapz_w(grid)
        self.lam = float(np.sum(self.tw * 2 * f / grid))
        self.wlog = math.exp(2 / self.lam * np.sum(self.tw * f * np.log(grid) / grid))
        self.w2 = math.sqrt(2 / self.lam * np.sum(self.tw * f * grid))
        support = grid[f > 0.01]
        self.wmax = float(support.max())

    def lam_n(self, n, T):
        nu = 2 * np.pi * np.asarray(n, float)[..., None] * KB * T
        return np.sum(self.tw * 2 * self.a * self.w / (self.w ** 2 + nu ** 2), axis=-1)


def isome_grid(emin, emax, encut=5000.0, steps=(1.0, 5.0, 50.0), bounds=(100.0, 500.0)):
    lo, hi = max(emin, -encut), min(emax, encut)
    pts = [np.arange(0.0, min(bounds[0], hi) + 1e-9, steps[0])]
    edges = list(bounds) + [hi]
    prev = bounds[0]
    for b, s in zip(edges[1:], steps[1:]):
        if prev >= hi:
            break
        pts.append(np.arange(prev, min(b, hi) + 1e-9, s))
        prev = b
    pos = np.unique(np.concatenate(pts))
    pts = [np.arange(0.0, max(-bounds[0], lo) - 1e-9, -steps[0])]
    prev = -bounds[0]
    for b, s in zip([-x for x in edges[1:-1]] + [lo], steps[1:]):
        if prev <= lo:
            break
        pts.append(np.arange(prev, max(b, lo) - 1e-9, -s))
        prev = b
    neg = np.unique(np.concatenate(pts))
    return np.unique(np.concatenate([neg, pos]))


class Electrons:
    """DOS N(eps) per spin per meV on the solver grid, optionally with W(eps, eps') in meV."""

    def __init__(self, e_meV, dos_per_spin_meV, W_e_meV=None, W_meV=None, encut=5000.0):
        e = np.asarray(e_meV, float)
        n = np.asarray(dos_per_spin_meV, float).copy()
        n[n < 0] = 0
        nz = np.nonzero(n > 0)[0]
        e, n = e[nz[0]:nz[-1] + 1], n[nz[0]:nz[-1] + 1]
        emin, emax = e[0], e[-1]
        if W_e_meV is not None:
            emin, emax = max(emin, W_e_meV[0]), min(emax, W_e_meV[-1])
        self.e = isome_grid(emin, emax, encut)
        self.N = np.interp(self.e, e, n)
        self.wq = trapz_w(self.e)
        self.idx_ef = int(np.argmin(np.abs(self.e)))
        self.NF = float(self.N[self.idx_ef])
        self.W = None
        if W_meV is not None:
            We = np.asarray(W_e_meV, float)
            Wm = np.asarray(W_meV, float)
            rows = np.array([np.interp(self.e, We, Wm[i]) for i in range(len(We))])
            self.W = np.array([np.interp(self.e, We, rows[:, j]) for j in range(len(self.e))]).T
            self.mu = float(self.NF * self.W[self.idx_ef, self.idx_ef])


def fermi(x, T):
    return 0.5 * (1 - np.tanh(x / (2 * KB * T)))


def matsubara(T, omega_c):
    M = int(math.ceil((omega_c / (math.pi * KB * T) - 1) / 2))
    return (2 * np.arange(M + 1) + 1) * math.pi * KB * T


def lambda_pm(a2f, T, nw):
    lam = a2f.lam_n(np.arange(2 * nw + 1), T)
    n = np.arange(nw)
    minus = lam[np.abs(n[:, None] - n[None, :])]
    plus = lam[n[:, None] + n[None, :] + 1]
    return minus - plus, minus + plus


def normal_state(a2f, el, T, omega_c, tol=1e-6, maxit=500, mix=0.5):
    """Self-consistent normal-state Z(iw), chi(iw) and chemical potential (vDOS, Delta = 0)."""
    wsi = matsubara(T, omega_c)
    nw = len(wsi)
    Lm, Lp = lambda_pm(a2f, T, nw)
    e, N, wq, NF = el.e, el.N, el.wq, el.NF
    kT = KB * T
    Ne0 = np.sum(wq * 2 * fermi(e, T) * N)
    Z, chi, mu = np.ones(nw), np.zeros(nw), 0.0
    for it in range(maxit):
        xi = e[:, None] - mu + chi[None, :]
        th = 1 / ((Z * wsi)[None, :] ** 2 + xi ** 2)
        zi = np.sum((wq * N)[:, None] * th * (wsi * Z)[None, :], axis=0)
        si = np.sum((wq * N)[:, None] * th * xi, axis=0)
        Zn = 1 + kT * (Lm @ zi) / wsi / NF
        chin = -kT * (Lp @ si) / NF

        def dNe(m):
            d = e - m
            t = (wsi[None, :] * Zn[None, :]) ** 2 + (d[:, None] + chin[None, :]) ** 2
            s = np.sum((d[:, None] + chin[None, :]) / t - d[:, None] / (wsi[None, :] ** 2 + d[:, None] ** 2), axis=1)
            return Ne0 - np.sum(wq * N * (2 * fermi(d, T) - 4 * kT * s))
        a, b = mu - 50, mu + 50
        fa, fb = dNe(a), dNe(b)
        k = 0
        while fa * fb > 0 and k < 60:
            a, b = a - 100, b + 100
            fa, fb = dNe(a), dNe(b)
            k += 1
        mun = brentq(dNe, a, b, xtol=1e-7)
        change = np.max(np.abs(Zn - Z)) + np.max(np.abs(chin - chi)) + abs(mun - mu)
        Z = (1 - mix) * Z + mix * Zn
        chi = (1 - mix) * chi + mix * chin
        mu = (1 - mix) * mu + mix * mun
        if change < tol:
            break
    return wsi, Z, chi, mu


def eig_cdos(a2f, mustar, T, omega_c):
    wsi = matsubara(T, omega_c)
    nw = len(wsi)
    Lm, Lp = lambda_pm(a2f, T, nw)
    kT = KB * T
    Z = 1 + math.pi * kT * (Lm @ np.ones(nw)) / wsi
    # Z_n D_n = pi kT sum_m (L+ - 2 mu*)_nm D_m / w_m (linearized); symmetrized with sqrt(Z w)
    d = np.sqrt(Z * wsi)
    S = math.pi * kT * (Lp - 2 * mustar)
    ev = np.linalg.eigvalsh(S / d[:, None] / d[None, :])[-1]
    return float(ev), {'Z0': float(Z[0]), 'n_matsubara': nw}


def eig_vdos(a2f, el, mustar, T, omega_c, ns=None):
    wsi, Z, chi, mu = ns if ns is not None else normal_state(a2f, el, T, omega_c)
    nw = len(wsi)
    Lm, Lp = lambda_pm(a2f, T, nw)
    kT = KB * T
    e, N, wq, NF = el.e, el.N, el.wq, el.NF
    xi = e[:, None] - mu + chi[None, :]
    th = 1 / ((Z * wsi)[None, :] ** 2 + xi ** 2)
    a = np.sum((wq * N)[:, None] * th, axis=0) * Z          # integral N Z / Theta
    # Z_n Delta_n = (kT/NF) sum_m (Lp - 2 mu*) a_m Delta_m
    S = kT / NF * (Lp - 2 * mustar)
    sa = np.sqrt(a / Z)
    M = sa[:, None] * S * sa[None, :]
    ev = np.linalg.eigvalsh(M)[-1]
    return float(ev), {'Z0': float(Z[0]), 'chi0_meV': float(chi[0]), 'mu_shift_meV': float(mu), 'n_matsubara': nw}


def eig_vdos_W(a2f, el, T, omega_c, ns=None):
    wsi, Z, chi, mu = ns if ns is not None else normal_state(a2f, el, T, omega_c)
    nw = len(wsi)
    Lm, Lp = lambda_pm(a2f, T, nw)
    kT = KB * T
    e, N, wq, NF = el.e, el.N, el.wq, el.NF
    xi = e[:, None] - mu + chi[None, :]
    th = 1 / ((Z * wsi)[None, :] ** 2 + xi ** 2)                 # ndos x nw
    wN = wq * N
    B = wN[:, None] * th                                           # ndos x nw
    K11 = kT / NF * Lp * B.sum(axis=0)[None, :]
    K12 = kT / NF * Lp @ B.T
    x0 = e - mu
    with np.errstate(divide='ignore', invalid='ignore'):
        tail = np.where(np.abs(x0) < 1e-9, 1 / (4 * kT), np.tanh(x0 / (2 * kT)) / (2 * x0))
    c = 2 * kT * np.sum(th - 1 / (wsi[None, :] ** 2 + x0[:, None] ** 2), axis=1) + tail
    WW = el.W * wN[None, :]
    K21 = -2 * kT * WW @ th
    K22 = -WW * c[None, :]
    K = np.block([[K11, K12], [K21, K22]])
    ev = np.linalg.eigvals(K)
    lead = float(np.max(ev.real[np.abs(ev.imag) < 1e-8 * max(1, np.max(np.abs(ev.real)))]))
    return lead, {'Z0': float(Z[0]), 'chi0_meV': float(chi[0]), 'mu_shift_meV': float(mu), 'n_matsubara': nw}


def allen_dynes(lam, wlog, w2, mustar):
    """Allen-Dynes Tc (K) with strong-coupling and shape corrections; wlog, w2 in meV."""
    den = lam - mustar * (1 + 0.62 * lam)
    if den <= 0:
        return 0.0
    f1 = (1 + (lam / (2.46 * (1 + 3.8 * mustar))) ** 1.5) ** (1 / 3)
    L2 = 1.82 * (1 + 6.3 * mustar) * w2 / wlog
    f2 = 1 + (w2 / wlog - 1) * lam ** 2 / (lam ** 2 + L2 ** 2)
    return f1 * f2 * wlog / (1.2 * KB) * math.exp(-1.04 * (1 + lam) / den)


T_FLOOR, T_CAP = 10.0, 500.0


def critical_temperature(f_eig, guess, xtol=0.01):
    """Tc (K) where the leading eigenvalue of the linearized gap equation crosses 1.

    Returns None if Tc lies below T_FLOOR; raises ValueError above T_CAP."""
    g = lambda T: f_eig(T)[0] - 1
    a = min(max(guess, T_FLOOR), T_CAP)
    ga = g(a)
    b = a
    for _ in range(80):
        b = min(max(b * (1.12 if ga > 0 else 1 / 1.12), T_FLOOR), T_CAP)
        gb = g(b)
        if (ga > 0) != (gb > 0):
            break
        if b == T_FLOOR:
            return None
        if b == T_CAP:
            raise ValueError('critical temperature above the solver range (500 K)')
        a, ga = b, gb
    return brentq(g, min(a, b), max(a, b), xtol=xtol)
