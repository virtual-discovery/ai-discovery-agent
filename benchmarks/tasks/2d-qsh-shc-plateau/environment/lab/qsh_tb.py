"""Numpy Wannier tight-binding evaluator for 2d-qsh-shc-plateau (copied unchanged into the task lab).

Each model is an archived Wannier90 aiida_hr.dat (fully relativistic, spinor Wannier functions; Wannier90 3.1,
use_ws_distance = .true.) stored as npz: R (lattice vectors after the Wigner-Seitz image selection of Wannier90),
H(R) (complex64, degeneracy weights already applied) and the final Wannier centres (Angstrom).

Conventions
  H(k) = sum_R exp(2 pi i k . R) H(R),  k in reduced coordinates of the reciprocal basis of the archived cell
         (k_cart = k @ B, B = 2 pi inv(A)^T, A = CELL_PARAMETERS of aiida.in in Angstrom).
  Velocity v_a (eV Angstrom, hbar = 1), Cartesian a = x, y, in one of two position conventions:
     'lattice'          v_a = dH/dk_a with every Wannier function placed at the lattice point of its home cell;
     'wannier_centres'  v_a = dH/dk_a + i [H, tau_a], tau = diagonal position operator at the Wannier centres.
  Spin: Pauli matrices acting on the spinor index of the Wannier functions, ordered (1 up, 1 down, 2 up, ...).
  Spin Berry curvature (Qiao et al., PRB 98, 214402 (2018); Wannier90 berry_task = shc, omega = 0):
     Omega^g_n = sum_{m != n} -2 Im[<n|j^g|m><m|v_y|n>] / ((e_n - e_m)^2 + eta^2),  j^g = {sigma^g, v_x}/2.
  Spin Hall conductivity at T = 0 (step occupation):
     sheet  sigma^g_xy(mu) = pi/(S N_k) sum_k sum_{e_n < mu} Omega^g_n     in units of (hbar/e)(e^2/h) = e/2pi,
     Wannier90 output units: sheet * (e^2/h) / c  in (hbar/e) S/cm, with S the cell area and c the cell height
     (vacuum included), as in the archived aiida-shc-fermiscan.dat files.
"""
import numpy as np

E2_OVER_H = 3.874045864e-5          # e^2/h in S
CHUNK = 400


def pauli():
    return np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], complex)


def select_refine(c, refine):
    """Indices of the coarse k-points to refine: max |contribution| over gamma and mu above the threshold,
    at most max_points of them, largest first."""
    f, thr, cap = refine
    score = np.abs(c).max(axis=(1, 2))
    cand = np.where(score > thr)[0]
    return cand[np.argsort(-score[cand], kind='stable')][:cap]


class Model:
    def __init__(self, path, cell_A):
        z = np.load(path)
        self.R = z['R'].astype(int)
        H = z['H'].astype(complex)
        self.nr, self.nw = H.shape[0], H.shape[1]
        self.Hflat = H.reshape(self.nr, -1)
        self.centres = z['centres'].astype(float) if 'centres' in z else None
        self.A = np.array(cell_A, float)
        self.B = 2 * np.pi * np.linalg.inv(self.A).T
        self.Rc = self.R @ self.A
        self.area = abs(np.linalg.det(self.A[:2, :2]))
        self.c = abs(self.A[2, 2])
        s = pauli()
        self.S = np.array([np.kron(np.eye(self.nw // 2), s[g]) for g in range(3)])
        # row decomposition R = (R1, R2): H(k1, k2) = sum_R1 e^{2 pi i k1 R1} sum_R2 e^{2 pi i k2 R2} H(R1, R2)
        self.r1 = np.unique(self.R[:, 0])
        self.i1 = np.searchsorted(self.r1, self.R[:, 0])
        self.sel1 = np.zeros((len(self.r1), self.nr))
        self.sel1[self.i1, np.arange(self.nr)] = 1.0

    # ------------------------------------------------------------------ basic evaluations
    def hk(self, k):
        k = np.atleast_2d(np.asarray(k, float))
        return (np.exp(2j * np.pi * k @ self.R.T) @ self.Hflat).reshape(-1, self.nw, self.nw)

    def bands(self, k):
        k = np.atleast_2d(np.asarray(k, float))
        return np.concatenate([np.linalg.eigvalsh(self.hk(k[s:s + CHUNK])) for s in range(0, len(k), CHUNK)])

    def spin_expectation(self, k):
        k = np.atleast_2d(np.asarray(k, float))
        out_e, out_s = [], []
        for s in range(0, len(k), CHUNK):
            e, U = np.linalg.eigh(self.hk(k[s:s + CHUNK]))
            out_e.append(e)
            out_s.append(np.stack([np.real(np.sum(U.conj() * (self.S[g] @ U), axis=1)) for g in range(3)], -1))
        return np.concatenate(out_e), np.concatenate(out_s)

    def row(self, k2, k1s, derivs=True):
        """H and Cartesian dH/dk (x, y) on the k-points (k1s, k2) of one mesh row."""
        p2 = np.exp(2j * np.pi * k2 * self.R[:, 1])
        H1 = (self.sel1 * p2[None, :]) @ self.Hflat                      # (nr1, nw^2)
        p1 = np.exp(2j * np.pi * np.outer(k1s, self.r1))                # (nk, nr1)
        H = (p1 @ H1).reshape(-1, self.nw, self.nw)
        if not derivs:
            return H, None
        G2 = (self.sel1 * (p2 * self.R[:, 1])[None, :]) @ self.Hflat
        D1 = (p1 * self.r1[None, :]) @ H1                                # sum R1 e H1
        D2 = p1 @ G2                                                     # sum e R2 H
        dH = [(1j * (self.A[0, a] * D1 + self.A[1, a] * D2)).reshape(-1, self.nw, self.nw) for a in range(2)]
        return H, dH

    def mesh_rows(self, mesh, shift=(0.0, 0.0)):
        n1, n2 = mesh
        k1s = (np.arange(n1) + shift[0]) / n1
        for j in range(n2):
            yield (j + shift[1]) / n2, k1s

    # ------------------------------------------------------------------ band edges
    def band_edges(self, n, mesh):
        """Maximum of band n and minimum of band n+1 (and the smallest direct gap) on a Gamma-centred mesh."""
        best = {'vbm': (-np.inf, None), 'cbm': (np.inf, None), 'direct': (np.inf, None)}
        for k2, k1s in self.mesh_rows(mesh):
            H, _ = self.row(k2, k1s, derivs=False)
            e = np.linalg.eigvalsh(H)
            i = int(np.argmax(e[:, n - 1])); j = int(np.argmin(e[:, n])); d = e[:, n] - e[:, n - 1]; l = int(np.argmin(d))
            if e[i, n - 1] > best['vbm'][0]:
                best['vbm'] = (float(e[i, n - 1]), [float(k1s[i]), float(k2)])
            if e[j, n] < best['cbm'][0]:
                best['cbm'] = (float(e[j, n]), [float(k1s[j]), float(k2)])
            if d[l] < best['direct'][0]:
                best['direct'] = (float(d[l]), [float(k1s[l]), float(k2)])
        return best

    # ------------------------------------------------------------------ spin Hall conductivity
    def hk_dhk(self, k):
        """H and Cartesian dH/dk (x, y) at arbitrary reduced k-points (full Fourier sum)."""
        ph = np.exp(2j * np.pi * np.atleast_2d(k) @ self.R.T)
        H = (ph @ self.Hflat).reshape(-1, self.nw, self.nw)
        dH = [((1j * ph * self.Rc[:, a]) @ self.Hflat).reshape(-1, self.nw, self.nw) for a in range(2)]
        return H, dH

    def contributions(self, H, dH, gammas, eta, mus, tau):
        """Per-k occupied spin Berry curvature sum_{e_n < mu} Omega^g_n (Angstrom^2); shape (nk, n_gamma, n_mu)."""
        e, U = np.linalg.eigh(H)
        Uh = np.conj(np.transpose(U, (0, 2, 1)))
        v = []
        for a in range(2):
            va = Uh @ dH[a] @ U
            if tau is not None:
                va = va + 1j * (e[:, :, None] - e[:, None, :]) * (Uh @ (tau[:, a][None, :, None] * U))
            v.append(va)
        de = e[:, :, None] - e[:, None, :]
        inv = 1.0 / (de ** 2 + eta ** 2)
        idx = np.arange(self.nw)
        inv[:, idx, idx] = 0.0
        vyT = np.transpose(v[1], (0, 2, 1))
        occ = (e[:, :, None] < mus[None, None, :]).astype(float)
        out = np.zeros((len(e), len(gammas), len(mus)))
        for i, g in enumerate(gammas):
            Sb = Uh @ (self.S[g] @ U)
            js = 0.5 * (Sb @ v[0] + v[0] @ Sb)
            om = (-2 * np.imag(js * vyT) * inv).sum(2)
            out[:, i, :] = np.einsum('kn,knm->km', om, occ)
        return out

    def shc(self, mesh, gammas, eta, mus, convention, refine=None):
        """sigma^gamma_xy(mu) for each gamma (0, 1, 2 = x, y, z) on a Gamma-centred N1 x N2 mesh.

        refine = (factor, threshold_A2, max_points): coarse points where max |sum_{e_n<mu} Omega^g_n| over the
        requested gamma and mu exceeds threshold_A2 (at most max_points of them, largest first) are replaced by
        the average over a factor x factor sub-mesh of their mesh cell (Wannier90 berry_curv_adpt_kmesh).
        Returns (sheet in e/2pi, Wannier90 units (hbar/e) S/cm, number of refined coarse points)."""
        if convention not in ('lattice', 'wannier_centres'):
            raise ValueError('unknown position convention')
        tau = self.centres if convention == 'wannier_centres' else None
        mus = np.asarray(mus, float)
        n1, n2 = mesh
        coarse, kpts = [], []
        for k2, k1s in self.mesh_rows(mesh):
            for s in range(0, len(k1s), CHUNK):
                H, dH = self.row(k2, k1s[s:s + CHUNK])
                coarse.append(self.contributions(H, dH, gammas, eta, mus, tau))
                kpts.append(np.c_[k1s[s:s + CHUNK], np.full(len(H), k2)])
        c = np.concatenate(coarse)
        kpts = np.concatenate(kpts)
        nref = 0
        if refine is not None:
            f = refine[0]
            cand = select_refine(c, refine)
            off = (np.arange(f) + 0.5) / f - 0.5
            sub = np.array([(a / n1, b / n2) for a in off for b in off])
            for s in range(0, len(cand), max(1, CHUNK // (f * f))):
                block = cand[s:s + max(1, CHUNK // (f * f))]
                kk = (kpts[block][:, None, :] + sub[None]).reshape(-1, 2)
                H, dH = self.hk_dhk(np.c_[kk, np.zeros(len(kk))])
                cs = self.contributions(H, dH, gammas, eta, mus, tau).reshape(len(block), f * f, len(gammas), len(mus))
                c[block] = cs.mean(1)
            nref = len(cand)
        sheet = np.pi * c.sum(0) / (self.area * n1 * n2)
        return sheet, sheet * E2_OVER_H / (self.c * 1e-8), nref

    # ------------------------------------------------------------------ topology
    def wcc(self, nocc, k_list, n_loop, loop_axis):
        """Hybrid Wannier charge centres (fractions of the lattice vector along loop_axis) of bands 1..nocc, from the
        discrete Wilson loop along reciprocal axis loop_axis (0 or 1) at each k of the other in-plane axis."""
        out = []
        t = np.arange(n_loop) / n_loop
        for k in k_list:
            kf = np.zeros((n_loop, 3))
            kf[:, 1 - loop_axis] = k
            kf[:, loop_axis] = t
            _, U = np.linalg.eigh(self.hk(kf))
            u = U[:, :, :nocc]
            M = np.conj(np.transpose(u, (0, 2, 1))) @ np.roll(u, -1, axis=0)
            a, _, bh = np.linalg.svd(M)
            W = np.eye(nocc, dtype=complex)
            for j in range(n_loop):
                W = W @ (a[j] @ bh[j])
            out.append(np.sort(np.angle(np.linalg.eigvals(W)) / (2 * np.pi) % 1.0))
        return np.array(out)

    def ribbon(self, k_list, width, finite_axis, window, edge_cells=3):
        """Ribbon `width` cells wide along lattice vector finite_axis (0 or 1), periodic along the other lattice vector
        (reduced momentum k). Per k: eigenvalues inside window and the weight of each state on the first and on the
        last edge_cells cells."""
        from scipy.linalg import eigh
        p = 1 - finite_axis
        Rp, Rf = self.R[:, p], self.R[:, finite_axis]
        H = self.Hflat.reshape(self.nr, self.nw, self.nw)
        nw, dim = self.nw, width * self.nw
        out = []
        for k in k_list:
            ph = np.exp(2j * np.pi * k * Rp)
            M = np.zeros((dim, dim), complex)
            for d in np.unique(Rf):
                if abs(d) >= width:
                    continue
                sel = Rf == d
                blk = np.tensordot(ph[sel], H[sel], axes=(0, 0))
                for i in range(max(0, -d), min(width, width - d)):
                    j = i + d
                    M[i * nw:(i + 1) * nw, j * nw:(j + 1) * nw] += blk
            M = 0.5 * (M + M.conj().T)
            e, v = eigh(M, subset_by_value=window)
            w = (np.abs(v) ** 2).reshape(width, nw, -1).sum(1)
            out.append((e, w[:edge_cells].sum(0), w[-edge_cells:].sum(0)))
        return out
