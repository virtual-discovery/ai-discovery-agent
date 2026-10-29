"""Numpy/scipy tight-binding evaluator for the archived FPLO Wannier models (copied unchanged into the task lab).

A model is one archived +hamdata file converted to npz: lattice vectors A (rows, Angstrom), Wannier centres (Angstrom),
integer R (in units of the file's own lattice vectors) and H[n, i, j] = <i, 0|H|j, R_n> in eV (FPLO energy zero).
Orbitals 0-1: Zn-4s (spin up, down); 2-7: pnictogen p (m = -1, 0, +1; spin up, down), as named in the file.

Conventions
  bulk:   H(k)_ij = sum_R H_ij(R) exp(i k.(R + w_j - w_i)), k Cartesian in 1/Angstrom (2 pi included), in the file's frame.
  slabs and prisms: a unimodular integer matrix S (rows = new basis vectors in units of the file's lattice vectors)
          defines the new cell.  Sites are kept by position: every Zn and pnictogen site whose reduced coordinates
          along the finite directions, measured from a pnictogen site, lie in [0, N).  k along a periodic direction is
          the Bloch phase per period (radians, -pi..pi).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

ZN = np.array([1, 1, 0, 0, 0, 0, 0, 0], bool)
NW = 8


class Model:
    def __init__(self, path):
        z = np.load(path)
        self.A = z['lattice_A']
        self.centers = z['centers_A']
        self.R = z['R'].astype(int)
        self.H = z['H'].astype(complex)
        self.B = 2 * np.pi * np.linalg.inv(self.A).T
        self._T = self.R @ self.A
        self._dw = self.centers[None, :, :] - self.centers[:, None, :]
        self._Hflat = self.H.reshape(len(self.R), -1)

    # ---------- bulk ----------
    def hk(self, kcart):
        kcart = np.atleast_2d(np.asarray(kcart, float))
        ph = np.exp(1j * kcart @ self._T.T)
        Hk = (ph @ self._Hflat).reshape(-1, NW, NW)
        return Hk * np.exp(1j * np.einsum('ka,ija->kij', kcart, self._dw))

    def eig(self, kcart, chunk=2000):
        kcart = np.atleast_2d(np.asarray(kcart, float))
        es, vs = [], []
        for s in range(0, len(kcart), chunk):
            e, v = np.linalg.eigh(self.hk(kcart[s:s + chunk]))
            es.append(e)
            vs.append(v)
        return np.concatenate(es), np.concatenate(vs)

    def bands(self, kcart, chunk=4000):
        kcart = np.atleast_2d(np.asarray(kcart, float))
        return np.concatenate([np.linalg.eigvalsh(self.hk(kcart[s:s + chunk])) for s in range(0, len(kcart), chunk)])

    # ---------- supercell bases ----------
    def regroup(self, S):
        """Hoppings in the basis S (rows, integer, det +-1) with Zn regrouped so that its reduced coordinates
        relative to the pnictogen site lie in [0, 1).  Returns (A_new, R_new, H) with H unchanged in order."""
        S = np.array(S, int)
        if abs(round(np.linalg.det(S))) != 1:
            raise ValueError('basis matrix must be unimodular (determinant +1 or -1)')
        Sinv = np.linalg.inv(S)
        Rn = self.R @ Sinv
        Rn_i = np.round(Rn).astype(int)
        assert np.abs(Rn - Rn_i).max() < 1e-6
        A_new = S @ self.A
        f = (self.centers[0] - self.centers[2]) @ np.linalg.inv(A_new)
        g = np.floor(f + 1e-6).astype(int)
        return A_new, Rn_i, g

    def _blocks(self, S, finite_axes, kphase):
        """Dictionary {offset along finite axes: 8x8 block summed over periodic axes with Bloch phases}."""
        A_new, Rn, g = self.regroup(S)
        periodic = [a for a in range(3) if a not in finite_axes]
        out = {}
        # orbital-pair specific shift of R: (Zn, X) -> R - g ; (X, Zn) -> R + g
        shift = np.zeros((NW, NW, 3), int)
        shift[np.ix_(ZN, ~ZN)] = -g
        shift[np.ix_(~ZN, ZN)] = g
        for sh in {tuple(x) for x in shift.reshape(-1, 3)}:
            mask = np.all(shift == np.array(sh), axis=2)
            Rs = Rn + np.array(sh)
            phase = np.exp(1j * sum(kphase[i] * Rs[:, a] for i, a in enumerate(periodic)))
            Hm = self.H * mask[None]
            for key in {tuple(r) for r in Rs[:, finite_axes]}:
                sel = np.all(Rs[:, finite_axes] == np.array(key), axis=1)
                blk = np.tensordot(phase[sel], Hm[sel], axes=(0, 0))
                out[key] = out.get(key, 0) + blk
        return out

    def finite_matrix(self, S, sizes, kphase, onsite=None):
        """Sparse Hamiltonian of a slab (one finite axis, S row 2) or prism (two finite axes, S rows 1 and 2).
        sizes: number of cells along each finite axis; kphase: Bloch phases along the periodic axes.
        onsite: optional array (ncells, 8) of extra onsite energies (eV)."""
        finite_axes = list(range(3 - len(sizes), 3))
        blocks = self._blocks(S, finite_axes, kphase)
        grids = np.meshgrid(*[np.arange(n) for n in sizes], indexing='ij')
        cells = np.stack([gr.ravel() for gr in grids], 1)
        strides = np.cumprod([1] + list(sizes[::-1]))[:-1][::-1]
        ncell = len(cells)
        oi, oj = np.meshgrid(np.arange(NW), np.arange(NW), indexing='ij')
        rows, cols, vals = [], [], []
        for key, blk in blocks.items():
            if not np.any(blk):
                continue
            tgt = cells + np.array(key)
            ok = np.all((tgt >= 0) & (tgt < np.array(sizes)), axis=1)
            a = cells[ok] @ strides * NW
            b = tgt[ok] @ strides * NW
            rows.append((a[:, None] + oi.ravel()[None]).ravel())
            cols.append((b[:, None] + oj.ravel()[None]).ravel())
            vals.append(np.tile(blk.ravel(), ok.sum()))
        dim = ncell * NW
        M = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(dim, dim))
        if onsite is not None:
            M = M + sp.diags(np.asarray(onsite, float).ravel())
        return M, cells

    def finite_states(self, M, emin, emax, nmax=150, dense_limit=1600):
        """Eigenpairs with energies in [emin, emax]: from a dense solver when the matrix is small enough, otherwise from
        shift-invert Lanczos around the window centre.  The k eigenvalues nearest the centre contain every eigenvalue
        within their largest distance from it, so k is doubled (up to nmax + 20) until that distance reaches the window's
        half-width.  If more than nmax states fall in the window, the nmax nearest the centre are kept.
        Returns (energies, vectors, truncated)."""
        dim = M.shape[0]
        centre, half = 0.5 * (emin + emax), 0.5 * (emax - emin)
        if dim <= dense_limit:
            e, v = np.linalg.eigh(M.toarray())
            complete = True
        else:
            k = min(40, dim - 2)
            while True:
                e, v = sla.eigsh(M, k=k, sigma=centre, which='LM', tol=1e-10)
                order = np.argsort(e)
                e, v = e[order], v[:, order]
                complete = np.abs(e - centre).max() >= half
                if complete or k >= min(nmax + 20, dim - 2):
                    break
                k = min(2 * k, nmax + 20, dim - 2)
        sel = (e >= emin) & (e <= emax)
        e, v = e[sel], v[:, sel]
        truncated = (not complete) or len(e) > nmax
        if len(e) > nmax:
            keep = np.sort(np.argsort(np.abs(e - centre), kind='stable')[:nmax])
            e, v = e[keep], v[:, keep]
        return e, v, truncated
