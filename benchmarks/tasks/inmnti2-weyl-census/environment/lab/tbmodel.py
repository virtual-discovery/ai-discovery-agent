"""Numpy Wannier tight-binding evaluator for inmnti2-weyl-census (copied unchanged into the task lab).

Each model is an archived Wannier90 hr.dat Hamiltonian stored as npz (R, deg, H as complex64).
Conventions:
  H(k) = sum_R exp(i k_cart . R_cart) H(R) / deg(R)      (Wannier90 tight-binding convention)
  k_cart in 1/Angstrom (2 pi included), Cartesian axes of the QE CELL_PARAMETERS; k_cart = k_frac @ B.
  Bands are numbered from 1 (lowest Wannier band).
Berry flux: outward flux of Omega = curl A, A = i<u|grad u>, of Wannier bands 1..nocc through a sphere,
from the gauge-invariant link-variable (Fukui-Hatsugai-Suzuki) discretisation, divided by 2 pi.
"""
import numpy as np

CELL = np.array([[0.0, 3.1335920710, 3.1335920710],
                 [3.1335920710, 0.0, 3.1335920710],
                 [3.1335920710, 3.1335920710, 0.0]])
B = 2 * np.pi * np.linalg.inv(CELL).T
CHUNK = 1000


class Model:
    def __init__(self, path):
        z = np.load(path)
        self.Rcart = z['R'].astype(float) @ CELL
        H = z['H'].astype(complex) / z['deg'].astype(float)[:, None, None]
        self.nr, self.nw = H.shape[0], H.shape[1]
        self.Hflat = H.reshape(self.nr, -1)

    def hk(self, kcart):
        kcart = np.atleast_2d(np.asarray(kcart, float))
        ph = np.exp(1j * kcart @ self.Rcart.T)
        return (ph @ self.Hflat).reshape(-1, self.nw, self.nw)

    def dhk(self, kcart):
        """dH/dk_a for a = x, y, z at Cartesian k; shape (3, nk, nw, nw)."""
        kcart = np.atleast_2d(np.asarray(kcart, float))
        ph = np.exp(1j * kcart @ self.Rcart.T)
        return np.array([((1j * ph * self.Rcart[:, a]) @ self.Hflat).reshape(-1, self.nw, self.nw) for a in range(3)])

    def bands(self, kcart):
        kcart = np.atleast_2d(np.asarray(kcart, float))
        return np.concatenate([np.linalg.eigvalsh(self.hk(kcart[s:s + CHUNK])) for s in range(0, len(kcart), CHUNK)])

    def eig(self, kcart):
        kcart = np.atleast_2d(np.asarray(kcart, float))
        parts = [np.linalg.eigh(self.hk(kcart[s:s + CHUNK])) for s in range(0, len(kcart), CHUNK)]
        return np.concatenate([p[0] for p in parts]), np.concatenate([p[1] for p in parts])

    def gap_and_grad(self, kcart, n):
        """Gap E_{n+1} - E_n, its Cartesian gradient (Hellmann-Feynman) and the mean energy of the two bands."""
        e, v = np.linalg.eigh(self.hk(kcart))
        d = self.dhk(kcart)
        lo, hi = v[0][:, n - 1], v[0][:, n]
        grad = np.array([np.real(hi.conj() @ d[a, 0] @ hi - lo.conj() @ d[a, 0] @ lo) for a in range(3)])
        return float(e[0][n] - e[0][n - 1]), grad, float(0.5 * (e[0][n] + e[0][n - 1]))

    def minimize_gap(self, k0, n, max_iter=200):
        """BFGS on g(k)^2 (a Weyl cone becomes a smooth paraboloid) from the Cartesian seed k0."""
        from scipy.optimize import minimize

        def f(k):
            g, gr, _ = self.gap_and_grad(k, n)
            return g * g, 2 * g * gr
        res = minimize(f, np.asarray(k0, float), jac=True, method='BFGS', options={'gtol': 1e-12, 'maxiter': max_iter})
        g, _, e = self.gap_and_grad(res.x, n)
        return res.x, g, e, int(res.nit)

    def gap_scan(self, n, mesh, keep=60):
        """Local minima of E_{n+1} - E_n on a Gamma-centred mesh^3 grid of fractional coordinates."""
        f = np.arange(mesh) / mesh
        kf = np.array(np.meshgrid(f, f, f, indexing='ij')).reshape(3, -1).T
        e = self.bands(kf @ B)
        g = (e[:, n] - e[:, n - 1]).reshape(mesh, mesh, mesh)
        mid = (0.5 * (e[:, n] + e[:, n - 1])).reshape(mesh, mesh, mesh)
        loc = np.ones_like(g, bool)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    if dx or dy or dz:
                        loc &= g <= np.roll(np.roll(np.roll(g, dx, 0), dy, 1), dz, 2)
        idx = np.argwhere(loc)
        order = np.argsort(g[loc], kind='stable')[:keep]
        rows = []
        for i, j, k in idx[order]:
            frac = np.array([i, j, k]) / mesh
            frac = frac - np.floor(frac + 0.5)
            rows.append({'k_frac': frac.round(6).tolist(), 'k_cart_invA': (frac @ B).round(6).tolist(),
                         'gap_eV': round(float(g[i, j, k]), 6), 'mean_energy_eV': round(float(mid[i, j, k]), 6)})
        return rows, int(loc.sum()), float(g.min())

    def berry_flux(self, k0, radius, nocc, ntheta=16, nphi=32):
        th = np.linspace(0, np.pi, ntheta + 1)
        ph = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
        T, P = np.meshgrid(th, ph, indexing='ij')
        pts = np.asarray(k0, float) + radius * np.stack([np.sin(T) * np.cos(P), np.sin(T) * np.sin(P), np.cos(T)], -1)
        e, v = self.eig(pts.reshape(-1, 3))
        u = v[:, :, :nocc].reshape(ntheta + 1, nphi, self.nw, nocc)

        def link(a, b):
            d = np.linalg.det(np.einsum('...mi,...mj->...ij', a.conj(), b))
            return d / np.maximum(np.abs(d), 1e-300)
        ut = link(u[:-1], u[1:])
        up = link(u, np.roll(u, -1, axis=1))
        flux = np.angle(ut * up[1:] / np.roll(ut, -1, axis=1) / up[:-1])
        return float(flux.sum() / (2 * np.pi)), float((e[:, nocc] - e[:, nocc - 1]).min())
