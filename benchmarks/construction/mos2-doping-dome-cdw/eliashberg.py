"""Isotropic Migdal-Eliashberg critical temperature on the Matsubara axis (numpy implementation).

Written for the mos2-doping-dome-cdw task from the equations of Girotto Erhardt, Berges, Ponce and Novko,
npj 2D Mater. Appl. 9, 44 (2025), Eqs. (1)-(9): full-bandwidth equations with an energy-dependent DOS
N(eps), the energy shift chi and particle-number conservation; and their constant-DOS limit. mu* is
referenced to the second-moment frequency omega_2 and rescaled to the Matsubara cutoff (Eqs. 7-9); in the
constant-DOS limit mu*_C = mu* / (1 + mu* ln(omega_2 / omega_C)). The Matsubara cutoff is 20 omega_2.
The numerical conventions (trapezoidal weights, cutoff, spectrum rounding) follow the authors' archived
tc.py workflow, and construction checks reproduce their archived Tc tables (fig2e/tc_*.dat) to <= 0.2 %.
This file is copied unchanged into the lab (environment/lab/eliashberg.py) by make_task.py.
"""
import math
import numpy as np
from scipy.optimize import brentq

KB = 8.61733e-5      # eV/K
CUTOFF = 20.0        # Matsubara cutoff in units of omega_2
T_FLOOR = 1.0        # K; lower critical temperatures are reported as below the floor
T_CAP = 150.0


def trap_weights(x):
    dx = np.empty_like(x)
    dx[0] = x[1] - x[0]
    dx[1:-1] = x[2:] - x[:-2]
    dx[-1] = x[-1] - x[-2]
    return dx / 2


class Spectrum:
    """alpha^2F(omega) on a grid in meV; rows with omega = 0 are dropped (as in the authors' workflow)."""

    def __init__(self, omega_meV, a2f):
        w = np.round(np.asarray(omega_meV, float) * 1e-3, 6)
        a = np.round(np.asarray(a2f, float), 6)
        keep = w != 0
        self.w, self.a = w[keep], a[keep]
        if len(self.w) < 3 or np.any(np.diff(self.w) <= 0):
            raise ValueError('spectrum needs at least 3 strictly increasing positive frequencies')
        self.wt = trap_weights(self.w) * self.a
        self.lam = float(self.lam_n(0, 1.0))
        if not self.lam > 0:
            raise ValueError('spectrum has no coupling (lambda <= 0)')
        self.wlog = math.exp(2 / self.lam * np.sum(self.wt * np.log(self.w) / self.w))
        self.w2 = math.sqrt(2 / self.lam * np.sum(self.wt * self.w))

    def lam_n(self, n, T):
        nu = 2 * np.pi * np.asarray(n, float)[..., None] * KB * T
        return 2 * np.sum(self.wt * self.w / (self.w ** 2 + nu ** 2), axis=-1)


class Einstein:
    def __init__(self, lam, wE):
        self.lam, self.wlog, self.w2 = lam, wE, wE

    def lam_n(self, n, T):
        nE = self.w2 / (2 * np.pi * KB * T)
        return self.lam / (1 + (np.asarray(n, float) / nE) ** 2)


class DOS:
    """Electronic DOS per spin and cell on an energy grid (eV); n = electrons in these bands per cell."""

    def __init__(self, energy, dos, n):
        self.e = np.asarray(energy, float)
        self.N = np.asarray(dos, float)
        self.n = float(n)
        self.w = trap_weights(self.e) * self.N
        self.states = float(self.w.sum())
        if not 0 < self.n < 2 * self.states:
            raise ValueError('electron count outside the range these bands can hold')

    def mu0(self, T):
        kT = KB * T
        mu = (self.e[0] * (2 * self.states - self.n) + self.e[-1] * self.n) / (2 * self.states)
        for _ in range(100000):
            x = self.e - mu
            with np.errstate(divide='ignore', invalid='ignore'):
                ms = np.where(np.abs(x) < 1e-14, 1 / (2 * kT), np.tanh(x / (2 * kT)) / x)
            new = (self.n - self.states + np.sum(self.w * ms * self.e)) / np.sum(self.w * ms)
            if abs(new - mu) < 1e-15:
                return float(new)
            mu = new
        return float(mu)


def allen_dynes(lam, wlog, mustar, w2=None, correct=True):
    """Allen-Dynes (correct=True, with shape correction if w2 is given) or McMillan Tc in K; wlog, w2 in eV."""
    tc = wlog / (1.2 * KB) * math.exp(-1.04 * (1 + lam) / max(1e-3, lam - mustar * (1 + 0.62 * lam)))
    if correct:
        tc *= (1 + (lam / (2.46 * (1 + 3.8 * mustar))) ** 1.5) ** (1 / 3)
        if w2 is not None:
            L2 = 1.82 * (1 + 6.3 * mustar) * w2 / wlog
            tc *= 1 + (w2 / wlog - 1) * lam ** 2 / (lam ** 2 + L2 ** 2)
    return tc


def _lambda_blocks(spec, T, no):
    lam = spec.lam_n(np.arange(1 - no, 2 * no), T)
    m = np.arange(no)
    off = no - 1
    minus = lam[m[:, None] - m[None, :] + off]
    plus = lam[m[:, None] + m[None, :] + 1 + off]
    return minus + plus, minus - plus, lam[off:]


def eigen_constant(spec, mustar, T):
    nE = spec.w2 / (2 * np.pi * KB * T)
    no = int(math.ceil(CUTOFF * nE - 0.5))
    L, _, lam_pos = _lambda_blocks(spec, T, no)
    muc = mustar / (1 + mustar * math.log(nE / (no + 0.5)))
    diag = np.empty(no)
    diag[0] = lam_pos[0]
    diag[1:] = lam_pos[0] + 2 * np.cumsum(lam_pos[1:no])
    d = 1 / np.sqrt(2 * np.arange(no) + 1.0)
    S = L - 2 * muc - np.diag(diag)
    return float(np.linalg.eigvalsh(d[:, None] * S * d[None, :])[-1]), {'mu_star_C': muc, 'n_matsubara': no}


def _R(dos, mu, omega, NF):
    x = dos.e - mu
    with np.errstate(divide='ignore', invalid='ignore'):
        ms = np.where(np.abs(x) < 1e-14, 1 / omega, np.arctan2(x, omega) / x)
    return float(np.sum(dos.w / NF * ms) / np.pi)


def eigen_full(spec, dos, mustar, T, tol=1e-11, maxiter=20000):
    kT = KB * T
    domega = 2 * np.pi * kT
    no = int(math.ceil(CUTOFF * spec.w2 / domega - 0.5))
    mu = dos.mu0(T)
    mu0 = mu
    NF = float(dos.N[int(np.argmin(np.abs(dos.e - mu)))])
    L, Lm, _ = _lambda_blocks(spec, T, no)
    G, Gm = L / NF, Lm / NF
    omega = domega * (np.arange(no) + 0.5)
    if mustar > 0:
        mu_full = 1 / (1 / mustar - _R(dos, mu, spec.w2, NF))
        muC = 1 / (1 / mu_full + _R(dos, mu, domega * (no + 0.5), NF))
    else:
        muC = 0.0
    U = -2 * muC / NF
    e, w = dos.e, dos.w
    Z, chi = np.ones(no), np.zeros(no)
    wcut = domega * (no + 0.5)
    for _ in range(maxiter):
        Theta = (omega * Z)[:, None] ** 2 + (e[None, :] - mu + chi[:, None]) ** 2
        A = np.sum(w / Theta, axis=1)
        B = np.sum(w * e / Theta, axis=1)
        Znew = 1 + (Gm @ (A * Z * omega)) * kT / omega
        chinew = -(G @ (A * (chi - mu) + B)) * kT
        Theta = (omega * Znew)[:, None] ** 2 + (e[None, :] - mu + chinew[:, None]) ** 2
        A = np.sum(w / Theta, axis=1)
        B = np.sum(w * e / Theta, axis=1)
        tail = np.sum(w * 2 / np.pi * np.arctan2(e - mu, wcut))
        munew = ((dos.n - dos.states + tail) / (4 * kT) + np.sum(A * chinew + B)) / np.sum(A)
        change = np.max(np.abs(Znew - Z)) + np.max(np.abs(chinew - chi)) + abs(munew - mu)
        Z, chi, mu = Znew, chinew, munew
        if change < tol:
            break
    Theta = (omega * Z)[:, None] ** 2 + (e[None, :] - mu + chi[:, None]) ** 2
    sa = np.sqrt(np.sum(w / Theta, axis=1))
    ev = float(np.linalg.eigvalsh(sa[:, None] * (kT * (G + U)) * sa[None, :])[-1])
    return ev, {'mu_star_C': muC, 'n_matsubara': no, 'mu0_eV': mu0, 'mu_eV': float(mu), 'N_at_mu0': NF,
                'Z0': float(Z[0])}


def critical_temperature(spec, mustar, dos=None):
    """Tc (K) where the largest eigenvalue of the linearized gap equation reaches 1, with diagnostics."""
    eig = (lambda T: eigen_full(spec, dos, mustar, T)) if dos is not None else (lambda T: eigen_constant(spec, mustar, T))
    f = lambda T: eig(T)[0] - 1
    t0 = min(max(allen_dynes(spec.lam, spec.wlog, mustar, spec.w2), T_FLOOR), T_CAP)
    if f(T_FLOOR) < 0:
        return None, eig(T_FLOOR)[1]
    a, fa = t0, f(t0)
    b = a
    for _ in range(200):
        b = min(max(b * (1.1 if fa > 0 else 1 / 1.1), T_FLOOR), T_CAP)
        fb = f(b)
        if (fa > 0) != (fb > 0):
            break
        if b in (T_FLOOR, T_CAP):
            raise ValueError('critical temperature outside the solver range')
        a, fa = b, fb
    tc = brentq(f, min(a, b), max(a, b), xtol=1e-6)
    return tc, eig(tc)[1]
