"""Plane-wave continuum model of R-stacked twisted TMD homobilayers (K valley, spin up), the form of Eq. (1) of
Zhang et al. (2024) / Wu et al. PRL 122, 086402 (2019):

    H = [[-hbar^2 (k - kappa_+)^2 / 2m* + Delta_b(r), Delta_T(r)], [Delta_T(r)^*, -hbar^2 (k - kappa_-)^2 / 2m* + Delta_t(r)]]

Potentials are given as Fourier coefficients {(m1, m2): c} on a moire reciprocal basis (b1, b2) (Cartesian, 1/A),
V(r) = sum c e^{i (m1 b1 + m2 b2).r}. kappa_pm = (4 pi / 3 aM)(-sqrt3/2, -+1/2), aM = a0 / (2 sin(theta/2)).
Energies meV, lengths A. Bands are returned in descending energy (valence: band 1 = topmost).
"""
import numpy as np

HB2M = 3809.982  # hbar^2 / (2 m_e) in meV A^2


def moire(theta_deg, a0):
    aM = a0 / (2 * np.sin(np.radians(theta_deg) / 2))
    g = 4 * np.pi / (np.sqrt(3) * aM)
    kk = 4 * np.pi / (3 * aM)
    return aM, g, kk * np.array([-np.sqrt(3) / 2, -0.5]), kk * np.array([-np.sqrt(3) / 2, 0.5])


def wu_basis(theta_deg, a0):
    """b1 = g1 (0 deg), b2 = g3 (120 deg): the basis in which the archived maps take the Eq. (1) form."""
    g = moire(theta_deg, a0)[1]
    return np.array([[g, 0.0], [g * np.cos(2 * np.pi / 3), g * np.sin(2 * np.pi / 3)]])


class Continuum:
    def __init__(self, theta_deg, a0, mstar, coef_b, coef_t, coef_T, nshell=5, basis=None):
        self.aM, _, self.kp, self.km = moire(theta_deg, a0)
        self.b = wu_basis(theta_deg, a0) if basis is None else np.asarray(basis, float)
        self.t = HB2M / mstar
        G = [(i, j) for i in range(-nshell, nshell + 1) for j in range(-nshell, nshell + 1)
             if abs(i - j) <= nshell]
        self.G = np.array(G)
        self.Gc = self.G @ self.b
        n = self.n = len(G)
        idx = {x: a for a, x in enumerate(G)}
        mats = []
        for coef in (coef_b, coef_t, coef_T):
            M = np.zeros((n, n), complex)
            for a, ga in enumerate(G):
                for dg, c in coef.items():
                    j = idx.get((ga[0] - dg[0], ga[1] - dg[1]))   # <k+ga|V|k+gb> = V_{ga-gb}
                    if j is not None:
                        M[a, j] += c
            mats.append(M)
        self.Vb, self.Vt, self.T = mats
        self.shift = {}
        for d, gv in ((0, (1, 0)), (1, (0, 1))):
            p = np.array([idx.get((ga[0] + gv[0], ga[1] + gv[1]), -1) for ga in G])
            self.shift[d] = p

    def hamiltonian(self, k):
        kb = np.asarray(k, float) + self.Gc - self.kp
        kt = np.asarray(k, float) + self.Gc - self.km
        n = self.n
        H = np.empty((2 * n, 2 * n), complex)
        H[:n, :n] = self.Vb - np.diag(self.t * (kb ** 2).sum(1))
        H[n:, n:] = self.Vt - np.diag(self.t * (kt ** 2).sum(1))
        H[:n, n:] = self.T
        H[n:, :n] = self.T.conj().T
        return H

    def bands(self, k, nb=6):
        e, v = np.linalg.eigh(self.hamiltonian(k))
        return e[::-1][:nb], v[:, ::-1][:, :nb]

    def _shifted(self, u, d):
        """psi_{k+b_d}(G) = psi_k(G + b_d), applied to both layer blocks."""
        p = self.shift[d]
        out = np.zeros_like(u)
        ok = p >= 0
        for L in range(2):
            blk = u[L * self.n:(L + 1) * self.n]
            out[L * self.n:(L + 1) * self.n][ok] = blk[p[ok]]
        return out

    def chern(self, nk=24, nb=3, groups=None):
        """Fukui-Hatsugai-Suzuki Chern numbers of the topmost bands; `groups` (list of band-index tuples) gives
        non-Abelian (multi-band) Chern numbers, default one group per band."""
        groups = groups or [(i,) for i in range(nb)]
        nb = max(max(g) for g in groups) + 1
        U = np.empty((nk + 1, nk + 1, 2 * self.n, nb), complex)
        E = np.empty((nk, nk, nb))
        for i in range(nk):
            for j in range(nk):
                e, U[i, j] = self.bands(i / nk * self.b[0] + j / nk * self.b[1], nb)
                E[i, j] = e
        for j in range(nk):
            U[nk, j] = self._shifted(U[0, j], 0)
        for i in range(nk + 1):
            U[i, nk] = self._shifted(U[i, 0], 1)
        out = []
        for g in groups:
            g = list(g)
            u = U[:, :, :, g]
            def link(a, b):
                return np.linalg.det(np.einsum('ijak,ijal->ijkl', a.conj(), b))
            U1 = link(u[:-1, :-1], u[1:, :-1]); U2 = link(u[1:, :-1], u[1:, 1:])
            U3 = link(u[1:, 1:], u[:-1, 1:]); U4 = link(u[:-1, 1:], u[:-1, :-1])
            F = np.angle(U1 * U2 * U3 * U4)
            out.append(-F.sum() / (2 * np.pi))    # arg<u|u'> = -A.dk, so the plaquette sum is -2 pi C
        return out, E


def coefficients_from_maps(top, bot, tun, ngrid=100, ngrid_T=10):
    """Fourier coefficients of the archived fig1 maps (grids include the periodic endpoint)."""
    def fc(z, n, thresh=1e-6):
        F = np.fft.fft2(np.asarray(z).reshape(n, n)[:-1, :-1]) / (n - 1) ** 2
        out = {}
        for i in range(n - 1):
            for j in range(n - 1):
                if abs(F[i, j]) > thresh:
                    m1 = i if i <= (n - 1) // 2 else i - (n - 1)
                    m2 = j if j <= (n - 1) // 2 else j - (n - 1)
                    out[(m1, m2)] = complex(F[i, j])
        return out
    return fc(bot, ngrid, 1e-3), fc(top, ngrid, 1e-3), fc(tun, ngrid_T, 1e-3)


def wu_coefficients(V, psi_deg, w):
    """Delta_b = 2V sum_{j odd} cos(g_j.r + psi), Delta_t with -psi; Delta_T = w(1 + e^{-i g2.r} + e^{-i g3.r})."""
    p = np.radians(psi_deg)
    odd = [(1, 0), (0, 1), (-1, -1)]          # g1, g3, g5 in the (b1, b2) basis
    cb, ct = {}, {}
    for s in odd:
        cb[s] = V * np.exp(1j * p); cb[(-s[0], -s[1])] = V * np.exp(-1j * p)
        ct[s] = V * np.exp(-1j * p); ct[(-s[0], -s[1])] = V * np.exp(1j * p)
    return cb, ct, {(0, 0): w, (-1, -1): w, (0, -1): w}   # g2 = b1 + b2, g3 = b2
