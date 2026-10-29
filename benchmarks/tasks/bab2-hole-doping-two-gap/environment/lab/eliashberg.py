"""Isotropic Migdal-Eliashberg equations on the imaginary axis (numpy), in the form EPW 5.8 solves them.

Matsubara frequencies w_n = (2n+1) pi kB T for n = 0..N-1 with N = int((cutoff/(pi kB T) - 1)/2) + 1;
negative frequencies enter through symmetry. With lambda(m) = int 2 w a2F(w) / (w^2 + (2 m pi kB T)^2) dw,

    Z_n = 1 + (pi kB T / w_n) sum_m [lambda(n-m) - lambda(n+m+1)] w_m / sqrt(w_m^2 + D_m^2)
    Z_n D_n = pi kB T sum_m [lambda(n-m) + lambda(n+m+1) - 2 mu*] D_m / sqrt(w_m^2 + D_m^2)

The equations are iterated with linear mixing D <- mixing * D_new + (1 - mixing) * D_old. The residual is
sum |D_new - D_old| / sum |D_new| (EPW's convergence measure). Units: meV and K.
"""
import numpy as np

KB = 0.08617333262  # meV / K


def n_matsubara(T, cutoff):
    return int(0.5 * (cutoff / (np.pi * KB * T) - 1)) + 1


def coupling(w, a2f):
    dw = w[1] - w[0]
    return float(2 * np.sum(a2f / w) * dw)


def _kernels(w, a2f, T, cutoff):
    N = n_matsubara(T, cutoff)
    wn = np.pi * KB * T * (2 * np.arange(N) + 1)
    nu = 2 * np.pi * KB * T * np.arange(2 * N + 1)
    dw = w[1] - w[0]
    lam = (2 * w[None, :] * a2f[None, :] / (w[None, :] ** 2 + nu[:, None] ** 2)).sum(1) * dw
    n = np.arange(N)
    return wn, lam[np.abs(n[:, None] - n[None, :])], lam[n[:, None] + n[None, :] + 1]


def solve(w, a2f, T, mu, cutoff=500.0, max_iter=500, conv_thr=1e-3, mixing=0.4, initial_gap=5.0):
    """Iterate the isotropic equations at temperature T; returns the gap and renormalization at w_0."""
    w, a2f = np.asarray(w, float), np.asarray(a2f, float)
    wn, Lm, Lp = _kernels(w, a2f, T, cutoff)
    D = np.full(len(wn), float(initial_gap))
    Z = np.ones(len(wn))
    residual, it = float('inf'), 0
    for it in range(1, max_iter + 1):
        r = 1 / np.sqrt(wn ** 2 + D ** 2)
        Z = 1 + np.pi * KB * T * ((Lm - Lp) @ (wn * r)) / wn
        Dn = np.pi * KB * T * ((Lm + Lp - 2 * mu) @ (D * r)) / Z
        norm = np.abs(Dn).sum()
        residual = float(np.abs(Dn - D).sum() / norm) if norm > 0 else 0.0
        D = mixing * Dn + (1 - mixing) * D
        if abs(D[0]) < 1e-6:
            return {'temperature_K': T, 'gap_meV': 0.0, 'Z': float(Z[0]), 'n_matsubara': len(wn), 'iterations': it,
                    'converged': True, 'residual': residual, 'state': 'gap collapsed to zero (normal state)'}
        if residual < conv_thr:
            break
    return {'temperature_K': T, 'gap_meV': float(D[0]), 'Z': float(Z[0]), 'n_matsubara': len(wn), 'iterations': it,
            'converged': bool(residual < conv_thr), 'residual': residual,
            'state': 'converged' if residual < conv_thr else 'not converged within max_iter'}


def leading_eigenvalue(w, a2f, T, mu, cutoff=500.0):
    """Largest eigenvalue of the linearized gap equation (construction reference only; not served)."""
    wn, Lm, Lp = _kernels(np.asarray(w, float), np.asarray(a2f, float), T, cutoff)
    Z = 1 + np.pi * KB * T * (Lm - Lp).sum(1) / wn
    K = np.pi * KB * T * (Lm + Lp - 2 * mu) / wn[None, :] / Z[:, None]
    return float(np.max(np.linalg.eigvals(K).real))


def tc_linearized(w, a2f, mu, cutoff=500.0, lo=2.0, hi=150.0):
    """Tc where the leading eigenvalue of the linearized equations crosses 1 (bisection)."""
    f = lambda T: leading_eigenvalue(w, a2f, T, mu, cutoff) - 1
    if f(hi) > 0 or f(lo) < 0:
        raise ValueError('Tc outside the bracket')
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi)
