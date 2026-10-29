"""Noncollinear linear spin-wave theory and neutron-scan helpers used by the reference solution
(and by the construction scripts in benchmarks/construction/mn5si3-noncollinear-spinwaves).

Classical energy (meV) of unit vectors e_i:
    E = - sum_{entries (i, j, R, J)} J e_i . e_{j+R} - sum_i sum_a K_a (n_a . e_i)^2
Dynamics (hbar = 1): S_i de_i/dt = - e_i x dE/de_i, so magnon energies scale as 1/S.
Bond vectors are handled in fractional coordinates of the orthorhombic magnetic cell.
"""
import numpy as np

LAT = np.array([6.9093, 6.9093 * np.sqrt(3), 4.8134])     # a, b = sqrt(3) a, c  (Angstrom)
MN2_FF = (0.4220, 17.684, 0.5948, 6.005, 0.0043, -0.609, -0.0219)   # Mn2+ <j0>


def formfactor(Q):
    s2 = (np.linalg.norm(2 * np.pi * np.asarray(Q, float) / LAT, axis=-1) / (4 * np.pi)) ** 2
    A, a, B, b, C, c, D = MN2_FF
    return A * np.exp(-a * s2) + B * np.exp(-b * s2) + C * np.exp(-c * s2) + D


def frames(e0):
    e0 = np.asarray(e0, float)
    e0 = e0 / np.linalg.norm(e0, axis=1, keepdims=True)
    F = []
    for e in e0:
        ref = np.array([1.0, 0, 0]) if abs(e[0]) < 0.9 else np.array([0, 1.0, 0])
        x = ref - (ref @ e) * e
        x /= np.linalg.norm(x)
        F.append((x, np.cross(e, x)))
    return e0, np.array(F)


def hessian_k(k, pos, pairs, e0, F, aniso=()):
    N = len(e0)
    ii = pairs[:, 0].astype(int); jj = pairs[:, 1].astype(int); J = pairs[:, 5]
    d = pairs[:, 2:5] + pos[jj] - pos[ii]
    ph = np.exp(2j * np.pi * (d @ np.asarray(k, float)))
    blk = -J[:, None, None] * np.einsum('bax,bcx->bac', F[ii], F[jj]) * ph[:, None, None]
    H4 = np.zeros((N, N, 2, 2), complex)
    np.add.at(H4, (ii, jj), blk)
    np.add.at(H4, (jj, ii), np.conj(np.transpose(blk, (0, 2, 1))))
    c = J * np.einsum('bx,bx->b', e0[ii], e0[jj])
    diag = np.zeros(N)
    np.add.at(diag, ii, c)
    np.add.at(diag, jj, c)
    H = H4.transpose(0, 2, 1, 3).reshape(2 * N, 2 * N)
    H[np.arange(2 * N), np.arange(2 * N)] += np.repeat(diag, 2)
    for i in range(N):
        b = slice(2 * i, 2 * i + 2)
        for K, n in aniso:
            n = np.asarray(n, float)
            fn = F[i] @ n
            H[b, b] += -2 * K * (np.outer(fn, fn) - (n @ e0[i]) ** 2 * np.eye(2))
    return H


def sqw(Q, pos, pairs, e0, S, aniso=(), perp=True, weights=True, tol=1e-6):
    """Magnon energies w > 0 and unpolarised neutron weights |dS_perp|^2 at Q (r.l.u.).
    Eigenvectors are normalised with the symplectic form; near-degenerate modes form one cluster."""
    e0n, F = frames(e0)
    N = len(e0n)
    S = np.asarray(S, float) * np.ones(N)
    Q = np.asarray(Q, float)
    H = hessian_k(Q, pos, pairs, e0n, F, aniso)
    H = 0.5 * (H + H.conj().T)
    Om = np.kron(np.diag(1.0 / S), np.array([[0, -1], [1, 0]]))
    lam, V = np.linalg.eig(Om @ H)
    w = (1j * lam).real
    B = np.kron(np.diag(S), -1j * np.array([[0, 1], [-1, 0]]))
    Qc = 2 * np.pi * Q / LAT
    qhat = Qc / max(np.linalg.norm(Qc), 1e-12)
    P = np.eye(3) - np.outer(qhat, qhat) if perp else np.eye(3)
    D = P @ np.einsum('i,iac,iam->cm', S, F, V.reshape(N, 2, -1))
    keep = np.where(w > 1e-7)[0]
    keep = keep[np.argsort(w[keep])]
    out_w, out_i, k = [], [], 0
    while k < len(keep):
        c = [keep[k]]
        while k + len(c) < len(keep) and abs(w[keep[k + len(c)]] - w[c[0]]) < tol * max(1, w[c[0]]):
            c.append(keep[k + len(c)])
        G = V[:, c].conj().T @ B @ V[:, c]
        tot = abs(np.trace(np.linalg.pinv(G) @ (D[:, c].conj().T @ D[:, c])).real)
        for _ in c:
            out_w.append(float(np.mean(w[c])))
            out_i.append(tot / len(c) if weights else 1.0)
        k += len(c)
    return np.array(out_w), np.array(out_i)


def energy_grad(e, pairs, aniso=()):
    g = np.zeros_like(e)
    ii = pairs[:, 0].astype(int); jj = pairs[:, 1].astype(int); J = pairs[:, 5]
    E = -np.sum(J * np.einsum('ij,ij->i', e[ii], e[jj]))
    np.add.at(g, ii, -J[:, None] * e[jj])
    np.add.at(g, jj, -J[:, None] * e[ii])
    for K, n in aniso:
        n = np.asarray(n, float)
        p = e @ n
        E -= K * np.sum(p ** 2)
        g -= 2 * K * p[:, None] * n[None, :]
    return E, g


def relax(e_init, pairs, aniso=()):
    """Local classical energy minimum inside the given cell (L-BFGS over spherical angles)."""
    from scipy.optimize import minimize
    e = np.asarray(e_init, float)
    e = e / np.linalg.norm(e, axis=1, keepdims=True)
    N = len(e)

    def unpack(x):
        th, ph = x[:N], x[N:]
        return th, ph, np.c_[np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]

    def f(x):
        th, ph, u = unpack(x)
        E, g = energy_grad(u, pairs, aniso)
        dth = np.c_[np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)]
        dph = np.c_[-np.sin(th) * np.sin(ph), np.sin(th) * np.cos(ph), np.zeros(N)]
        return E, np.r_[np.sum(g * dth, 1), np.sum(g * dph, 1)]

    u = e
    for _ in range(3):      # re-chart between rounds to avoid the pole singularity
        x = np.r_[np.arccos(np.clip(u[:, 2], -1, 1)), np.arctan2(u[:, 1], u[:, 0])]
        x = minimize(f, x, jac=True, method='L-BFGS-B', options={'maxiter': 20000, 'ftol': 1e-15, 'gtol': 1e-12}).x
        u = unpack(x)[2]
    E, g = energy_grad(u, pairs, aniso)
    return u, E, float(np.abs(np.cross(u, g)).max())


def ground_state(pairs, aniso=(), n=16, starts=8, seed=0, extra=()):
    """Lowest of several relaxations (random starts plus any given starting states)."""
    rng = np.random.default_rng(seed)
    best = None
    for st in list(extra) + [rng.normal(size=(n, 3)) for _ in range(starts)]:
        u, E, t = relax(st, pairs, aniso)
        if best is None or E < best[1] - 1e-9:
            best = (u, E, t)
    return best


# ---------------------------------------------------------------- neutron scans
DIRS = {'h': (lambda x: (x, 1.0, 0.0), 2.0), 'k': (lambda x: (0.0, x, 1.0), 3.0), 'l': (lambda x: (0.0, 3.0, x), 1.0)}
WIN = {'h': (1.45, 2.55), 'k': (2.0, 4.0), 'l': (0.5, 1.5)}


def line_grid(key, step=0.005):
    lo, hi = WIN[key]
    return np.round(np.arange(lo - 0.12, hi + 0.12 + 1e-9, step), 4)


def line_spectrum(key, model, weights=True, perp=True):
    pos, pairs, e0, S, aniso = model
    f = DIRS[key][0]
    return [sqw(np.array(f(x), float), pos, pairs, e0, S, aniso, perp=perp, weights=weights) for x in line_grid(key)]


def cut_profile(spec, key, E, x, sigE, sigQ, scale=1.0, ff=True):
    """Constant-energy scan along `key` at energy E: Gaussian energy resolution sigE (meV), Gaussian momentum
    resolution sigQ (1/A) along the scan, Mn form factor; magnon energies multiplied by `scale`."""
    axis = {'h': 0, 'k': 1, 'l': 2}[key]
    sx = sigQ * LAT[axis] / (2 * np.pi)
    g = line_grid(key)
    I = np.array([np.sum(i * np.exp(-0.5 * ((scale * w - E) / sigE) ** 2)) for w, i in spec])
    if ff:
        f = DIRS[key][0]
        I = I * formfactor(np.array([f(v) for v in g])) ** 2
    K = np.exp(-0.5 * ((np.asarray(x)[:, None] - g[None, :]) / sx) ** 2)
    return (K @ I) / K.sum(1)


def two_gauss_positions(x, y, c, guess):
    """Positions of the two branches from a fit A1 G(c-d1) + A2 G(c+d2) + const with a common width."""
    from scipy.optimize import curve_fit

    def f(x, a1, a2, d1, d2, w, b):
        return a1 * np.exp(-0.5 * ((x - c + d1) / w) ** 2) + a2 * np.exp(-0.5 * ((x - c - d2) / w) ** 2) + b
    amp = max(y.max() - y.min(), 1e-12)
    try:
        p, _ = curve_fit(f, x, y, p0=[amp, amp, guess, guess, max(guess / 2, 0.02), y.min()],
                         bounds=([0, 0, 0, 0, 0.005, -np.inf], [np.inf, np.inf, 0.5, 0.5, 0.3, np.inf]), maxfev=20000)
        return c - p[2], c + p[3]
    except Exception:
        return c - guess, c + guess


def predicted_peaks(spec, key, energies, guesses, sigE, sigQ, scale=1.0):
    """{E: (low, high)} branch positions of simulated constant-energy scans, extracted as in the archive."""
    c = DIRS[key][1]
    x = np.linspace(WIN[key][0], WIN[key][1], 221)
    out = {}
    for E, g in zip(energies, guesses):
        out[E] = two_gauss_positions(x, cut_profile(spec, key, E, x, sigE, sigQ, scale), c, g)
    return out
