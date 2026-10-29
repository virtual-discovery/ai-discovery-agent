"""Fitting and prediction pipeline shared by calibrate.py (construction) and mirrored in solution/solve.py.

Spectra: variable projection. Each archived 1.5 K component spectrum (smoothed, E in 0.3-3.0 meV, normalised) is fitted
by Gaussians of width sig at the distinct AF1 LSWT energies plus a constant, with non-negative amplitudes solved
linearly; the residual is minimised over J1..J12 and D (sigma continuation 0.06 -> 0.03 meV).
Null direction: the archived paths leave J3 - J11 undetermined (exact); it is fixed by a 1-D scan against the
block-averaged diffuse maps (SCGA at the nominal temperature, best scale and offset).
"""
import numpy as np
from scipy.optimize import least_squares, nnls
from scipy.ndimage import gaussian_filter1d
from model import Model, KB, map14_blocks, map14_centres, ko_blocks, ko_centres, norm3, REDUCED


def spectra_arrays(spectra, emin=0.3, emax=3.0, smooth=1.5):
    Q, Y = [], []
    for q, E, S in spectra:
        win = (E > emin) & (E < emax)
        comps = []
        for c in range(3):
            I = gaussian_filter1d(np.asarray(S)[:, c], smooth)[win]
            comps.append(I / I.max())
        Q.append(q); Y.append(comps); Ew = E[win]
    return np.array(Q), np.array(Y), Ew


def spec_resid(M, p, Q, Y, E, sig, S=2.5):
    w = M.lswt(p[:12], p[12], Q, S=S)
    out = []
    for k in range(len(Q)):
        w_ = np.unique(np.round(w[k], 4))
        A = np.c_[np.exp(-0.5 * ((E[:, None] - w_[None, :]) / sig) ** 2), np.ones(len(E))]
        for c in range(3):
            a, _ = nnls(A, Y[k][c]); out.append(A @ a - Y[k][c])
    return np.concatenate(out)


def fit_spectra(M, spectra, p0, sigs=(0.06, 0.03), nfev=400, fixD=None, S=2.5):
    Q, Y, E = spectra_arrays(spectra)
    p = np.array(p0, float)
    for sig in sigs:
        if fixD is None:
            f = lambda x: spec_resid(M, x, Q, Y, E, sig, S)
            r = least_squares(f, p, diff_step=1e-4, max_nfev=nfev, x_scale=0.05)
            p = r.x
        else:
            f = lambda x: spec_resid(M, np.r_[x, fixD], Q, Y, E, sig, S)
            r = least_squares(f, p[:12], diff_step=1e-4, max_nfev=nfev, x_scale=0.05)
            p = np.r_[r.x, fixD]
    return p, float(r.cost)


def map_chi2(M, p, maps, nbz=24, lam_fixed=None):
    """sum over temperatures of chi^2 (best scale + offset) of block-averaged maps; maps: {T: full map}"""
    qc = map14_centres()
    tot = 0.0
    for T, full in maps.items():
        y = map14_blocks(full).ravel()
        m = norm3(M.scga(p[:12], p[12], T, qc, nbz=nbz, lam_fixed=lam_fixed))
        A = np.c_[m, np.ones_like(m)]
        c = np.linalg.lstsq(A, y, rcond=None)[0]
        tot += float(np.sum((A @ c - y) ** 2) / np.var(y))
    return tot


def null_scan(M, p, maps, grid=None):
    """move along (J3 + d, J11 - d) and return the parameter vector minimising the map chi^2"""
    s = p[2] + p[10]
    grid = np.linspace(min(s, 0.0), max(s, 0.0), 41) if grid is None else grid
    best = None
    for J3 in grid:
        q = p.copy(); q[2] = J3; q[10] = s - J3
        try:
            c = map_chi2(M, q, maps)
        except Exception:
            continue
        if best is None or c < best[0]:
            best = (c, J3)
    # refine with a parabola around the minimum
    J3b = best[1]; h = grid[1] - grid[0]
    cs = []
    for J3 in (J3b - h, J3b, J3b + h):
        q = p.copy(); q[2] = J3; q[10] = s - J3
        cs.append(map_chi2(M, q, maps))
    den = cs[0] - 2 * cs[1] + cs[2]
    if den > 0:
        J3b = J3b + 0.5 * h * (cs[0] - cs[2]) / den
    q = p.copy(); q[2] = J3b; q[10] = s - J3b
    return q


def predictions(M, p, T14=14.0, T20=20.0, lam_fixed=None):
    J, D = np.asarray(p[:12]), float(p[12])
    k, lam = M.lt(J)
    dE = (-lam + D / 2) - M.af1_energy(J, D)
    m14 = norm3(M.scga(J, D, T14, map14_centres(), nbz=32, lam_fixed=lam_fixed)).reshape(27, 29)
    red = {}
    for name, keep in REDUCED.items():
        Jr = np.zeros(12); Jr[keep] = J[keep]
        red[name] = norm3(M.scga(Jr, D, T20, ko_centres(), nbz=32, lam_fixed=lam_fixed)).reshape(32, 32)
    return {'J': J.tolist(), 'D': D, 'k': np.asarray(k).tolist(), 'dE': float(dE), 'map14': m14, 'reduced': red}


# ---------------- non-recall global route (no literature values) ----------------
from scipy.optimize import differential_evolution


class Joint:
    """Joint objective: 10 spectra (Gaussian width 0.08 meV) + block-averaged 20/30 K maps (SCGA at tau * T, best scale
    and offset per map) + a penalty if the maximum of lambda_max(J(q)) on the k = 1/2 plane of the 36x24x36 box grid is
    not at the archived 11 K ordering peak."""
    def __init__(self, M, spectra, maps, kbox, wmap=2.0, s0=75.0, m0=420.0):
        self.M = M
        self.Q, self.Y, self.E = spectra_arrays(spectra)
        self.qc = map14_centres()
        self.ys = {T: map14_blocks(f).ravel() for T, f in maps.items()}
        self.wmap, self.s0, self.m0 = wmap, s0, m0
        h, l = np.meshgrid(np.arange(36) / 36 - 17 / 36, np.arange(36) / 36 - 17 / 36, indexing='ij')
        self.grid = np.c_[h.ravel(), np.full(h.size, 0.5), l.ravel()]
        self.kbox = np.array(kbox)

    def spec(self, p, sig, sub=None):
        ii = list(range(len(self.Q))) if sub is None else sub
        w = self.M.lswt_fast(p[:12], p[12], self.Q[ii]); out = []
        for kk, k in enumerate(ii):
            w_ = np.unique(np.round(w[kk], 4))
            A = np.c_[np.exp(-0.5 * ((self.E[:, None] - w_[None, :]) / sig) ** 2), np.ones(len(self.E))]
            for c in range(3):
                a, _ = nnls(A, self.Y[k][c]); out.append(A @ a - self.Y[k][c])
        return np.concatenate(out)

    def maps(self, J, D, tau, nbz=10):
        tot = 0.0
        for T, y in self.ys.items():
            m = norm3(self.M.scga(J, D, T * tau, self.qc, nbz=nbz)); A = np.c_[m, np.ones_like(m)]
            c = np.linalg.lstsq(A, y, rcond=None)[0]; tot += np.sum((A @ c - y) ** 2) / np.var(y)
        return tot

    def kpen(self, J):
        lam = self.M.lam_max(J, self.grid)
        return max(0.0, lam.max() - self.M.lam_max(J, self.kbox[None])[0])

    def f(self, p, sub):
        try:
            mm = self.maps(p[:12], p[12], p[13])
        except Exception:
            return 1e3
        return np.sum(self.spec(p[:13], 0.08, sub) ** 2) / self.s0 + self.wmap * mm / self.m0 + 50.0 * self.kpen(p[:12])

    def search(self, seed, maxiter=150):
        sub = list(range(0, len(self.Q), 2))
        bounds = [(-0.9, 0.3)] * 12 + [(-0.4, 0.0), (0.3, 2.0)]
        r = differential_evolution(lambda p: self.f(p, sub), bounds, popsize=8, maxiter=maxiter, tol=1e-8, seed=seed,
                                   polish=False, init='sobol', mutation=(0.5, 1.0), recombination=0.7)
        p = r.x
        for sig in (0.06, 0.03):
            def g(x):
                try:
                    mm = np.sqrt(self.wmap * self.maps(x[:12], x[12], x[13], nbz=16) / self.m0)
                except Exception:
                    mm = 30.0
                return np.r_[self.spec(x[:13], sig) / np.sqrt(121.0), mm, np.sqrt(50.0 * self.kpen(x[:12]))]
            p = least_squares(g, p, diff_step=1e-4, max_nfev=300, x_scale=np.r_[np.full(13, 0.05), 0.1]).x
        return p, float(r.fun)
