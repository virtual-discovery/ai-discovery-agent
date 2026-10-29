"""Fast MnWO4 model utilities (numpy/scipy only).

Hamiltonian H = -sum_{i!=j} J_ij e_i.e_j + D sum_i (e_i.n)^2, unit vectors, ordered pairs, J,D in meV.
Bonds are given as shells: list of (s, s', [n1,n2,n3]) with the neighbour in cell n + offset.
"""
import numpy as np
from scipy.optimize import brentq, minimize

KB = 0.08617333262
CELL = np.array([[4.8236999512, 0, 0], [0, 5.7551999092, 0], [-0.0946104705, 0, 5.0010053420]])
TAU = np.array([[0.5, 0.68490, 0.25], [0.5, 0.31510, 0.75]])
NAX = np.array([0.81915, 0, 0.57358]); NAX = NAX / np.linalg.norm(NAX)


class Model:
    def __init__(self, shells):
        """shells: list (len 12) of lists of [s, s2, n1, n2, n3] (every ordered pair listed)."""
        self.shell = []
        self.s1, self.s2, self.d, self.R = [], [], [], []
        for k, bonds in enumerate(shells):
            for s, s2, n1, n2, n3 in bonds:
                self.shell.append(k); self.s1.append(s); self.s2.append(s2); self.R.append((n1, n2, n3))
                self.d.append(np.array([n1, n2, n3]) + TAU[s2] - TAU[s])
        self.shell = np.array(self.shell); self.s1 = np.array(self.s1); self.s2 = np.array(self.s2)
        self.d = np.array(self.d); self.R = np.array(self.R)
        self.nsh = len(shells)
        self._af1_setup()

    # ---------- J(q) ----------
    def Jq(self, J, q):
        q = np.atleast_2d(q)
        ph = np.exp(2j * np.pi * q @ self.d.T) * np.asarray(J)[self.shell][None, :]
        M = np.zeros((len(q), 2, 2), complex)
        for a in range(2):
            for b in range(2):
                m = (self.s1 == a) & (self.s2 == b)
                M[:, a, b] = ph[:, m].sum(axis=1)
        return M

    def lam_max(self, J, q):
        return np.linalg.eigvalsh(self.Jq(J, q))[:, -1]

    def lt(self, J, n=40):
        g = (np.arange(n) + 0.5) / n - 0.5
        G = np.array(np.meshgrid(g, g, g, indexing='ij')).reshape(3, -1).T
        lm = self.lam_max(J, G)
        best = None
        for i in np.argsort(-lm)[:6]:
            r = minimize(lambda x: -self.lam_max(J, x[None])[0], G[i], method='Nelder-Mead',
                         options={'xatol': 1e-8, 'fatol': 1e-13, 'maxiter': 4000})
            if best is None or r.fun < best.fun:
                best = r
        return canonical_k(best.x), -best.fun

    # ---------- AF1 ----------
    @staticmethod
    def af1_sign(s, n1, n2, n3):
        s0 = (1, 1, -1, -1)[n1 % 4] if s == 0 else (1, -1, -1, 1)[n1 % 4]
        return s0 * (-1) ** ((n2 + n3) % 2)

    def _af1_setup(self):
        # magnetic cell: 4 chemical cells (n1 = 0..3, n2 = n3 = 0) x 2 sublattices; translations (2,1,0),(0,1,1),(0,0,2)
        T = np.array([[2, 1, 0], [0, 1, 1], [0, 0, 2]]).T
        Tinv = np.linalg.inv(T)
        sites = [(s, n1) for n1 in range(4) for s in range(2)]
        def find(s, n):
            for n1 in range(4):
                f = Tinv @ (np.array(n) - np.array([n1, 0, 0]))
                if np.allclose(f, np.round(f), atol=1e-9):
                    return sites.index((s, n1))
            raise KeyError
        A, B, Dd, Sh, Par = [], [], [], [], []
        for a, (s, n1) in enumerate(sites):
            for k in np.where(self.s1 == s)[0]:
                n = np.array([n1, 0, 0]) + self.R[k]
                b = find(self.s2[k], n)
                par = self.af1_sign(s, n1, 0, 0) * self.af1_sign(self.s2[k], *n) > 0
                A.append(a); B.append(b); Dd.append(self.d[k]); Sh.append(self.shell[k]); Par.append(par)
        self.af = dict(a=np.array(A), b=np.array(B), d=np.array(Dd), sh=np.array(Sh), par=np.array(Par), n=len(sites))

    def af1_energy(self, J, D):
        """classical energy per spin of AF1 (collinear along n)"""
        af = self.af
        sgn = np.where(af['par'], 1.0, -1.0)
        return -np.sum(np.asarray(J)[af['sh']] * sgn) / af['n'] + D

    def lswt(self, J, D, q, S=2.5):
        """AF1 magnon energies (meV), shape (nq, 8); classical precession with S = mu/(g muB)."""
        af = self.af; n = af['n']
        J = np.asarray(J, float)
        Jb = J[af['sh']] / S ** 2
        sgn = np.where(af['par'], 1.0, -1.0)
        Adiag = np.zeros(n)
        np.add.at(Adiag, af['a'], 2 * S * Jb * sgn)
        Adiag += -2 * (D / S ** 2) * S
        q = np.atleast_2d(q)
        out = np.zeros((len(q), n))
        eta = np.diag(np.r_[np.ones(n), -np.ones(n)])
        for iq, qq in enumerate(q):
            Am = np.diag(Adiag).astype(complex); Amm = Am.copy(); Bm = np.zeros((n, n), complex)
            ph = np.exp(2j * np.pi * af['d'] @ qq)
            c = -S * Jb  # per listing; each bond listed twice (once per end)
            p = af['par']
            np.add.at(Am, (af['a'][p], af['b'][p]), c[p] * ph[p])
            np.add.at(Am, (af['b'][p], af['a'][p]), c[p] * np.conj(ph[p]))
            np.add.at(Amm, (af['a'][p], af['b'][p]), c[p] * np.conj(ph[p]))
            np.add.at(Amm, (af['b'][p], af['a'][p]), c[p] * ph[p])
            np.add.at(Bm, (af['a'][~p], af['b'][~p]), 2 * c[~p] * ph[~p])
            M = np.block([[Am, Bm], [Bm.conj().T, Amm.T]])
            w = np.linalg.eigvals(eta @ M)
            out[iq] = np.sort(np.real(w))[n:]
        return out

    # ---------- SCGA ----------
    def scga(self, J, D, T, q, nbz=24, lam_fixed=None):
        """per-component structure factors (Sxx, Syy, Szz) per spin, Cartesian, easy axis NAX."""
        beta = 1.0 / (KB * T)
        g = (np.arange(nbz) + 0.5) / nbz - 0.5
        G = np.array(np.meshgrid(g, g, g, indexing='ij')).reshape(3, -1).T
        ev = np.linalg.eigvalsh(self.Jq(J, G))
        if lam_fixed is None:
            def f(lam):
                return 0.5 * np.mean(1 / (lam + beta * D - beta * ev)) + np.mean(1 / (lam - beta * ev)) - 1.0
            lo = beta * ev.max() + max(-beta * D, 0.0)
            lam = brentq(f, lo + 1e-12, lo + 100)
        else:
            lam = lam_fixed
        Jm = self.Jq(J, q)
        I2 = np.eye(2)
        Ke = (lam + beta * D) * I2 - beta * Jm
        Kp = lam * I2 - beta * Jm
        Se = 0.25 * np.real(np.linalg.inv(Ke).sum(axis=(1, 2)))
        Sp = 0.25 * np.real(np.linalg.inv(Kp).sum(axis=(1, 2)))
        n2 = NAX ** 2
        return np.c_[n2[0] * Se + (1 - n2[0]) * Sp, n2[1] * Se + (1 - n2[1]) * Sp, n2[2] * Se + (1 - n2[2]) * Sp]


def canonical_k(k):
    """representative of {+k, -k} + G with k_b = 1/2 sign convention: h <= 0 when possible"""
    k = np.asarray(k, float)
    cands = []
    for s in (1, -1):
        v = s * k
        v = v - np.round(v)
        cands.append(v)
    # prefer k_b = +0.5 and h<0
    def key(v):
        return (abs(abs(v[1]) - 0.5) > 1e-6, v[0] > 0, -v[1])
    v = sorted(cands, key=key)[0]
    if abs(abs(v[1]) - 0.5) < 1e-6:
        v[1] = 0.5
    return v


def k_distance(k1, k2):
    """distance between propagation vectors modulo sign and reciprocal lattice (r.l.u., max over components)"""
    best = 9.0
    for s in (1, -1):
        dv = np.asarray(k1) - s * np.asarray(k2)
        dv = dv - np.round(dv)
        best = min(best, np.abs(dv).max())
    return best


# ---------------- archived grids, prediction blocks and loading ----------------
import json as _json
from pathlib import Path as _Path

MAP_DK, MAP_DL = 4.66 / 372, 0.0111          # Fig. 3 grid: k_i = -2.33 + i dk (i = 0..372), l_j = -1.776 + j dl (j = 0..320)
KO_DK, KO_DL = 0.025, 1.0 / 90                # Fig. 5 grid: k_i = -2 + i dk (i = 0..160), l_j = j dl (j = 0..160)
TABLE_I = np.array([-0.592, -0.0032, -0.272, -0.336, -0.0176, -0.544, -0.176, -0.016, -0.32, -0.192, -0.0672, -0.0256])
TABLE_I_D = -0.192
REDUCED = {'5d': [0, 2, 3, 5, 6, 9, 10], '5e': [0, 3, 5, 8]}   # shells kept (0-based) in the withheld reduced models
SOLD_REDUCED = {'5b': list(range(12)), '5c': [0, 2, 3, 5, 6, 8, 9, 10]}  # 5b: no anisotropy


def map14_blocks(full):
    """27 x 29 block averages (7 x 11 grid points) of a full Fig.-3 map (373 x 321)"""
    out = np.zeros((27, 29))
    for r in range(27):
        i0 = 186 + 7 * r
        for c in range(29):
            j0 = 160 + 11 * (c - 14)
            out[r, c] = full[i0 - 3:i0 + 4, j0 - 5:j0 + 6].mean()
    return out


def map14_centres():
    r, c = np.meshgrid(np.arange(27), np.arange(29), indexing='ij')
    k = 7 * r * MAP_DK; l = 11 * (c - 14) * MAP_DL
    return np.c_[(-l / 2).ravel(), k.ravel(), l.ravel()]


def ko_blocks(full):
    """32 x 32 block averages (5 x 5 grid points) of a Fig.-5 map (161 x 161)"""
    return full[:160, :160].reshape(32, 5, 32, 5).mean(axis=(1, 3))


def ko_centres():
    r, c = np.meshgrid(np.arange(32), np.arange(32), indexing='ij')
    k = -2 + KO_DK * (5 * r + 2); l = KO_DL * (5 * c + 2)
    return np.c_[(-l / 2).ravel(), k.ravel(), l.ravel()]


def half_to_full(half):
    """Fig.-3 map from its k >= 0 half (187 x 321) using S(q) = S(-q) (the grid is symmetric about q = 0)"""
    return np.concatenate([half[1:][::-1, ::-1], half], axis=0)


def norm3(S):
    return np.sqrt((np.asarray(S) ** 2).sum(-1))


def truth(model_json):
    """shells and values of the authors' archived Hamiltonian (meV) from extract.py's model.json"""
    m = _json.loads(_Path(model_json).read_text())
    MRY = 13.605693122994
    rows = m['jfile']
    dist = lambda r: np.linalg.norm((np.array(r[2:5]) + TAU[int(r[1]) - 1] - TAU[int(r[0]) - 1]) @ CELL)
    ds = np.array([dist(r) for r in rows]); u = np.unique(np.round(ds, 3))
    shells = [[] for _ in u]; J = np.zeros(len(u))
    for r, dd in zip(rows, ds):
        k = np.searchsorted(u, np.round(dd, 3))
        shells[k].append([int(r[0]) - 1, int(r[1]) - 1, int(r[2]), int(r[3]), int(r[4])]); J[k] = r[5] * MRY
    return shells, J, m['kfile'][0][2] * MRY, u


def _lswt_fast(self, J, D, q, S=2.5):
    """vectorised AF1 LSWT over q (same result as Model.lswt)"""
    af = self.af; n = af['n']
    if not hasattr(self, '_P'):
        nb = len(af['a'])
        P = np.zeros((nb, n * n)); Pt = np.zeros((nb, n * n))
        P[np.arange(nb), af['a'] * n + af['b']] = 1
        Pt[np.arange(nb), af['b'] * n + af['a']] = 1
        self._P, self._Pt = P, Pt
    J = np.asarray(J, float)
    Jb = J[af['sh']] / S ** 2
    sgn = np.where(af['par'], 1.0, -1.0)
    Adiag = np.zeros(n)
    np.add.at(Adiag, af['a'], 2 * S * Jb * sgn)
    Adiag += -2 * (D / S ** 2) * S
    q = np.atleast_2d(q)
    ph = np.exp(2j * np.pi * q @ af['d'].T)                 # (nq, nb)
    c = -S * Jb
    p = af['par']
    cp = np.where(p, c, 0.0); ca = np.where(p, 0.0, 2 * c)
    A = (ph * cp) @ self._P + (np.conj(ph) * cp) @ self._Pt
    Amm = (np.conj(ph) * cp) @ self._P + (ph * cp) @ self._Pt
    B = (ph * ca) @ self._P
    nq = len(q)
    A = A.reshape(nq, n, n) + np.diag(Adiag)[None]
    Amm = Amm.reshape(nq, n, n) + np.diag(Adiag)[None]
    B = B.reshape(nq, n, n)
    M = np.zeros((nq, 2 * n, 2 * n), complex)
    M[:, :n, :n] = A; M[:, :n, n:] = B
    M[:, n:, :n] = -np.conj(np.transpose(B, (0, 2, 1))); M[:, n:, n:] = -np.transpose(Amm, (0, 2, 1))
    w = np.linalg.eigvals(M)
    return np.sort(np.real(w), axis=1)[:, n:]


Model.lswt_fast = _lswt_fast
