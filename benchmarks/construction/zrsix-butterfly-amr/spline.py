"""Smooth band-group interpolation for the Chambers solver.

The bands that cross the Fermi level form an isolated group (bands 12-13 for ZrSiS and ZrSiSe, 12-14 for ZrSiTe).
Inside the group the sorted bands cross along the nodal lines, so each band has kinks there; the elementary
symmetric polynomials of the group's eigenvalues (the characteristic-polynomial coefficients) are analytic wherever
the group is isolated. We interpolate those with periodic cubic B-splines on a uniform grid of exact eigenvalues and
recover energies (polynomial roots) and velocities (implicit differentiation) at arbitrary k.
"""
import numpy as np

VCAP = 7.5    # eV Angstrom


def coeffs_from_eigs(e):
    n = e.shape[-1]
    if n == 2:
        return np.stack([e[..., 0] + e[..., 1], e[..., 0] * e[..., 1]], -1)
    if n == 3:
        a, b, c = e[..., 0], e[..., 1], e[..., 2]
        return np.stack([a + b + c, a * b + b * c + a * c, a * b * c], -1)
    raise ValueError(n)


def prefilter(f):
    """Periodic cubic B-spline coefficients of samples f (n1, n2, n3, nf)."""
    F = np.fft.fftn(f, axes=(0, 1, 2))
    for d in range(3):
        n = f.shape[d]
        w = (4 + 2 * np.cos(2 * np.pi * np.arange(n) / n)) / 6
        shape = [1, 1, 1, 1]; shape[d] = n
        F = F / w.reshape(shape)
    return np.fft.ifftn(F, axes=(0, 1, 2)).real


def _w(t):
    t2 = t * t; t3 = t2 * t; u = 1 - t
    w = np.stack([u * u * u, 3 * t3 - 6 * t2 + 4, -3 * t3 + 3 * t2 + 3 * t + 1, t3], -1) / 6
    dw = np.stack([-3 * u * u, 9 * t2 - 12 * t, -9 * t2 + 6 * t + 3, 3 * t2], -1) / 6
    return w, dw


class GroupSpline:
    def __init__(self, coef, B, eref):
        self.N = np.array(coef.shape[:3])
        self.nf = coef.shape[3]
        self.flat = np.ascontiguousarray(coef.reshape(-1, self.nf))
        self.Binv = np.linalg.inv(B)
        self.eref = eref
        self.off = np.arange(-1, 3)

    def coeff_and_grad(self, kfrac):
        M = len(kfrac)
        x = kfrac * self.N
        i0 = np.floor(x).astype(np.int64)
        t = x - i0
        w0, d0 = _w(t[:, 0]); w1, d1 = _w(t[:, 1]); w2, d2 = _w(t[:, 2])
        i = (i0[:, 0, None] + self.off) % self.N[0]
        j = (i0[:, 1, None] + self.off) % self.N[1]
        l = (i0[:, 2, None] + self.off) % self.N[2]
        lin = ((i[:, :, None] * self.N[1] + j[:, None, :]) * self.N[2])[:, :, :, None] + l[:, None, None, :]
        vals = self.flat[lin.reshape(M, 64)]                                   # (M, 64, nf)
        w01 = (w0[:, :, None] * w1[:, None, :])
        W = np.empty((M, 4, 64))
        W[:, 0] = (w01[:, :, :, None] * w2[:, None, None, :]).reshape(M, 64)
        W[:, 1] = ((d0[:, :, None] * w1[:, None, :])[:, :, :, None] * w2[:, None, None, :]).reshape(M, 64) * self.N[0]
        W[:, 2] = ((w0[:, :, None] * d1[:, None, :])[:, :, :, None] * w2[:, None, None, :]).reshape(M, 64) * self.N[1]
        W[:, 3] = (w01[:, :, :, None] * d2[:, None, None, :]).reshape(M, 64) * self.N[2]
        out = np.matmul(W, vals)                                                # (M, 4, nf)
        f = out[:, 0]
        g = np.einsum('mdf,cd->mfc', out[:, 1:], self.Binv)                     # d/dk_cart, (M, nf, 3)
        return f, g

    def eig_vel(self, kfrac):
        """Sorted group energies (M, n) in eV (absolute) and velocities dE/dk (M, n, 3) in eV*Angstrom."""
        c, g = self.coeff_and_grad(kfrac)
        n = self.nf
        if n == 2:
            s = np.sqrt(np.maximum(c[:, 0] ** 2 - 4 * c[:, 1], 1e-14))
            lam = np.stack([(c[:, 0] - s) / 2, (c[:, 0] + s) / 2], -1)
        else:
            a, b, cc = -c[:, 0], c[:, 1], -c[:, 2]
            p = np.minimum(b - a * a / 3, -1e-14)
            q = 2 * a ** 3 / 27 - a * b / 3 + cc
            r = 2 * np.sqrt(-p / 3)
            phi = np.arccos(np.clip(3 * q / (p * r), -1, 1)) / 3
            lam = np.sort(r[:, None] * np.cos(phi[:, None] - 2 * np.pi * np.arange(3) / 3) - a[:, None] / 3, -1)
        num = np.zeros(lam.shape + (3,))
        for m in range(1, n + 1):
            num += ((-1) ** (m + 1)) * (lam ** (n - m))[..., None] * g[:, None, m - 1, :]
        dp = n * lam ** (n - 1)
        for m in range(1, n):
            dp = dp + ((-1) ** m) * (n - m) * c[:, m - 1, None] * lam ** (n - m - 1)
        dp = np.where(dp >= 0, np.maximum(dp, 1e-12), np.minimum(dp, -1e-12))
        v = num / dp[..., None]
        # interpolation error makes the implicit derivative singular at the in-group band crossings (nodal lines);
        # the exact model's Fermi velocities never exceed 7 eV A, so cap the speed there
        sp = np.linalg.norm(v, axis=-1, keepdims=True)
        v = v * np.minimum(1.0, VCAP / np.maximum(sp, 1e-12))
        return lam + self.eref, v
