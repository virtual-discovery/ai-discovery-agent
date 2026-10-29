"""Wannier tight-binding tools and multi-orbital slave-spin mean field (SSMF) for the 135 kagome models.

This file is the evaluator of the task: it is copied unchanged into the lab and is also the construction
engine. numpy >= 1.21 and scipy >= 1.8 only.

Model: 33 Wannier orbitals per cell (archive order): 3 transition-metal sites x 5 d (indices 0-14, site-local
frames), 5 pnictogen x 3 p (15-29), 3 transition-metal s (30-32). Spin-degenerate (no spin-orbit coupling).

SSMF (Z2 slave spins in the de' Medici formulation, paramagnetic, T = 0 with a small Fermi smearing):
- correlated orbitals: the 15 d orbitals; density-density Kanamori interaction
  H_int = U sum_m S_m^ S_m_ + (U - 2J) sum_{m != m'} S_m^ S_m'_ + (U - 3J) sum_{m < m', s} S_ms S_m's,
  written with S^z = n - 1/2;
- every one-body term that involves a d orbital, except the local diagonal d energies, is renormalized by
  sqrt(Z) per d end (local inter-orbital terms included);
- gauge c = 1/sqrt(n(1-n)) - 1 and the corresponding lambda0 term, so that Z = 1 at U = 0;
- fully-localized-limit double counting in the S^z convention: the d levels are shifted by
  -[(2M-1)U - 5(M-1)J] (nbar - 1/2) + (U - 3J)/2, with nbar the mean d occupation per spin orbital of a site;
- the three kagome sites are symmetry-equivalent: Z, lambda and h are site-averaged per orbital.
"""
import numpy as np
from scipy import optimize, sparse
from scipy.sparse import linalg as splinalg

ND, NORB, M = 15, 33, 5


def load_model(path):
    d = np.load(path)
    return {'R': d['R'], 'degen': d['degen'], 'H': d['H'], 'cell': d['cell'], 'efermi': float(d['efermi'])}


def kmesh(n1, n2, n3):
    g = np.stack(np.meshgrid(np.arange(n1) / n1, np.arange(n2) / n2, np.arange(n3) / n3, indexing='ij'), -1)
    return g.reshape(-1, 3)


def kmesh_tr(n1, n2, n3):
    """Gamma-centred mesh reduced by time reversal (k ~ -k): points and integer weights (sum = n1 n2 n3)."""
    idx = np.stack(np.meshgrid(np.arange(n1), np.arange(n2), np.arange(n3), indexing='ij'), -1).reshape(-1, 3)
    n = np.array([n1, n2, n3])
    key = lambda i: (i[..., 0] * n2 + i[..., 1]) * n3 + i[..., 2]
    a, b = key(idx), key((-idx) % n)
    keep = a <= b
    return idx[keep] / n, np.where(a[keep] == b[keep], 1.0, 2.0)


def hk(model, ks):
    ph = np.exp(2j * np.pi * ks @ model['R'].T) / model['degen']
    return np.einsum('kr,rmn->kmn', ph, model['H'], optimize=True)


def fermi(x):
    return 0.5 * (1.0 - np.tanh(0.5 * x))


def find_mu(e, nelec, kT, weights=None):
    """Chemical potential giving nelec electrons (both spins) per cell for eigenvalues e (nk, nb)."""
    w = np.ones(e.shape[0]) if weights is None else weights
    nk = w.sum()
    target = nelec / 2.0
    lo, hi = e.min() - 1.0, e.max() + 1.0
    for _ in range(100):
        mu = 0.5 * (lo + hi)
        if w @ fermi((e - mu) / kT).sum(1) / nk > target:
            hi = mu
        else:
            lo = mu
        if hi - lo < 1e-11:
            break
    return 0.5 * (lo + hi)


def electron_count(e, mu, kT, weights=None):
    w = np.ones(e.shape[0]) if weights is None else weights
    return 2.0 * (w @ fermi((e - mu) / kT).sum(1)) / w.sum()


def sigma_exact(model):
    """Second centred moment of each orbital's projected DOS, exact: sigma_a^2 = sum_{(b,R) != (a,0)} |H_ab(R)|^2."""
    H, R = model['H'], model['R']
    w = 1.0 / model['degen']
    s2 = np.einsum('r,rab->a', w, np.abs(H) ** 2)
    i0 = np.where((R == 0).all(1))[0][0]
    s2 -= np.abs(np.diag(H[i0])) ** 2
    return np.sqrt(s2)


# ---------------------------------------------------------------- slave-spin site problem
class SlaveSpinSite:
    """2M spins (orbital m, spin s -> bit 2m+s), basis of 2^(2M) S^z states."""

    def __init__(self, U, J):
        nb = 2 * M
        self.dim = 1 << nb
        states = np.arange(self.dim)
        self.bits = ((states[:, None] >> np.arange(nb)[None, :]) & 1).astype(float)
        sz = self.bits - 0.5
        up, dn = sz[:, 0::2], sz[:, 1::2]
        e = U * (up * dn).sum(1)
        e += (U - 2 * J) * (up.sum(1) * dn.sum(1) - (up * dn).sum(1))
        for s in (up, dn):
            e += (U - 3 * J) * 0.5 * (s.sum(1) ** 2 - (s ** 2).sum(1))
        self.eint = e
        rows, cols, orb = [], [], []
        for b in range(nb):
            rows.append(states)
            cols.append(states ^ (1 << b))
            orb.append(np.full(self.dim, b // 2))
        self.rows, self.cols, self.orb = np.concatenate(rows), np.concatenate(cols), np.concatenate(orb)
        self.v0 = None

    def sx(self, psi):
        """Spin-averaged <S^x_m> for the ground state psi."""
        out = np.zeros(M)
        for b in range(2 * M):
            out[b // 2] += 0.5 * psi @ psi[np.arange(self.dim) ^ (1 << b)]
        return out / 2.0

    def ground(self, a, lam):
        diag = self.eint + self.bits @ np.repeat(lam, 2)
        vals = 0.5 * a[self.orb]
        H = sparse.csr_matrix((vals, (self.rows, self.cols)), shape=(self.dim, self.dim)) + sparse.diags(diag)
        try:
            w, v = splinalg.eigsh(H, k=1, which='SA', v0=self.v0, tol=1e-13, maxiter=20000)
            psi = v[:, 0]
        except splinalg.ArpackNoConvergence:
            psi = np.linalg.eigh(H.toarray())[1][:, 0]
        if psi.sum() < 0:
            psi = -psi
        self.v0 = psi
        n = ((psi ** 2) @ self.bits).reshape(M, 2).mean(1)
        return n, self.sx(psi)


# ---------------------------------------------------------------- SSMF driver
class Anderson:
    """Anderson (Pulay) mixing for the fixed point x = g(x)."""

    def __init__(self, beta=0.3, depth=6):
        self.beta, self.depth, self.X, self.F = beta, depth, [], []

    def step(self, x, g):
        f = g - x
        self.X.append(x.copy())
        self.F.append(f.copy())
        self.X, self.F = self.X[-self.depth - 1:], self.F[-self.depth - 1:]
        if len(self.F) == 1 or self.depth == 0:
            return x + self.beta * f
        dF = np.array([self.F[i + 1] - self.F[i] for i in range(len(self.F) - 1)])
        dX = np.array([self.X[i + 1] - self.X[i] for i in range(len(self.X) - 1)])
        gamma = np.linalg.lstsq(dF.T, f, rcond=None)[0]
        return x + self.beta * f - (dX.T + self.beta * dF.T) @ gamma


class SSMF:
    def __init__(self, model, mesh=(12, 12, 8), kT=0.01):
        self.model = model
        self.ks, self.wk = kmesh_tr(*mesh)
        self.Hk = hk(model, self.ks)
        i0 = np.where((model['R'] == 0).all(1))[0][0]
        self.eloc = np.real(np.diag(model['H'][i0]))[:ND].copy()     # local diagonal d energies
        self.Hoff = self.Hk.copy()
        idx = np.arange(ND)
        self.Hoff[:, idx, idx] -= self.eloc[None, :]
        self.kT = kT
        self.isd = np.arange(NORB) < ND

    def noninteracting(self, nelec):
        e, v = np.linalg.eigh(self.Hk)
        mu = find_mu(e, nelec, self.kT, self.wk)
        f = fermi((e - mu) / self.kT) * self.wk[:, None]
        occ = np.einsum('kan,kn->a', np.abs(v) ** 2, f) / self.wk.sum()
        return mu, occ

    def fermion_step(self, Z, shift, nelec):
        """Solve the renormalized fermion problem. Z, shift: per d orbital (length 5, site-symmetric)."""
        q = np.ones(NORB)
        q[:ND] = np.sqrt(np.tile(Z, 3))
        Hf = self.Hoff * np.outer(q, q)[None]
        idx = np.arange(ND)
        Hf[:, idx, idx] += (self.eloc + np.tile(shift, 3))[None, :]
        e, v = np.linalg.eigh(Hf)
        mu = find_mu(e, nelec, self.kT, self.wk)
        f = fermi((e - mu) / self.kT) * self.wk[:, None]
        nk = self.wk.sum()
        nf = np.einsum('kan,kn->a', np.abs(v) ** 2, f) / nk
        # h_a = Re sum_k sum_b Hoff_ab q_b <f+_a f_b> = Re sum_k sum_n (Hoff q V)_an f_n conj(V_an)
        AV = np.matmul(self.Hoff[:, :ND, :] * q[None, None, :], v)
        h = np.real(np.einsum('kan,kn,kan->a', AV, f, v[:, :ND, :].conj(), optimize=True)) / nk
        return mu, nf, h, e

    def run(self, nelec, U, J, maxiter=300, tol=2e-6, mix=0.3, depth=6, Z0=None, shift0=None, lam0_guess=None,
            max_shift_step=0.3, newton_after=1e-3, warmup=20, conv_tol=1e-4, verbose=False):
        """Self-consistent SSMF at fixed electron count. Unknowns x = (Z_m, d-level shift_m); the map
        x -> g(x) solves the fermions at x, then the slave spins with lambda_m enforcing n_S = n_f.
        Damped Anderson mixing first, then a Powell-hybrid Newton solve of g(x) = x."""
        site = SlaveSpinSite(U, J)
        mu0, occ0 = self.noninteracting(nelec)
        n0 = occ0[:ND].reshape(3, M).mean(0)
        Ueff = (2 * M - 1) * U - 5 * (M - 1) * J
        state = {'lam': np.zeros(M) if lam0_guess is None else np.array(lam0_guess, float), 'calls': 0}

        def g(x):
            Z, tot = np.clip(x[:M], 1e-10, 1.0), x[M:]
            mu, nf, h, e = self.fermion_step(Z, tot, nelec)
            n = np.clip(nf[:ND].reshape(3, M).mean(0), 1e-6, 1 - 1e-6)
            h = h.reshape(3, M).mean(0)
            c = 1.0 / np.sqrt(n * (1 - n)) - 1.0
            a = 2.0 * h * (1.0 + c)
            sol = optimize.root(lambda l: site.ground(a, l)[0] - n, state['lam'], method='hybr',
                                options={'xtol': 1e-12})
            lam = sol.x
            ns, sx = site.ground(a, lam)
            if np.abs(ns - n).max() < 1e-7:
                state['lam'] = lam
            Znew = ((1.0 + c) * np.abs(sx)) ** 2
            lam0 = np.abs(h) * np.sqrt(Znew) * (1 - 2 * n) / (n * (1 - n))
            dc = -Ueff * (n.mean() - 0.5) + 0.5 * (U - 3 * J)
            state.update(calls=state['calls'] + 1, n=n, mu=mu, h=h, e=e, lam0=lam0, Znew=Znew,
                         constraint=float(np.abs(ns - n).max()))
            return np.concatenate([Znew, -lam + lam0 + dc])

        Z = np.ones(M) if Z0 is None else np.array(Z0, float)
        tot = np.full(M, -Ueff * (n0.mean() - 0.5) + 0.5 * (U - 3 * J)) if shift0 is None else np.array(shift0, float)
        x = np.concatenate([Z, tot])
        anderson = Anderson(beta=mix, depth=depth)
        limit = np.concatenate([np.full(M, 0.1), np.full(M, max_shift_step)])
        res = np.inf
        for it in range(warmup):
            gx = g(x)
            f = gx - x
            res = np.abs(f / np.concatenate([np.ones(M), np.full(M, 10.0)])).max()
            if verbose:
                print('mix', it, np.round(gx[:M], 4), np.round(state['n'], 4), 'res %.2e' % res, flush=True)
            if res < max(tol, newton_after):
                break
            step = anderson.step(x, gx) - x
            if state['n'].min() < 0.02 or state['n'].max() > 0.98:
                anderson.X, anderson.F = [], []
                step = mix * f
            scale = np.max(np.abs(step) / limit)
            x = x + (step / scale if scale > 1 else step)
        if res >= tol:
            sol = optimize.root(lambda y: g(y) - y, x, method='hybr', options={'xtol': 1e-10, 'maxfev': maxiter, 'eps': 1e-7})
            x = sol.x
            gx = g(x)
            res = np.abs((gx - x) / np.concatenate([np.ones(M), np.full(M, 10.0)])).max()
            if verbose:
                print('newton', sol.message, np.round(gx[:M], 4), np.round(state['n'], 4), 'res %.2e' % res, flush=True)
        converged = bool(res < conv_tol and state['constraint'] < 1e-6)
        return {'Z': state['Znew'], 'n': state['n'], 'mu': state['mu'], 'lam': state['lam'], 'lam0': state['lam0'],
                'shift': x[M:], 'h': state['h'], 'iterations': state['calls'], 'residual': float(res),
                'converged': converged, 'nelec_check': electron_count(state['e'], state['mu'], self.kT, self.wk)}

    def scan(self, nelec, j_over_u, u_grid, u_start=1.0, step=0.5, min_step=0.02, **kw):
        """Metallic branch continued from weak coupling: U is increased from u_start through every value of
        u_grid at fixed J/U, each stage starting from the converged previous one; a stage that fails is retried
        with half the step. Returns {U: result} for the grid values reached; the scan stops at the first U where
        the metallic solution cannot be continued (reported under 'stopped_at')."""
        kw.setdefault('maxiter', 150)
        targets = sorted(u_grid)
        out, path = {}, []
        u = u_start
        prev = self.run(nelec, u, j_over_u * u, **kw)
        path.append((u, prev['converged']))
        if not prev['converged']:
            return {'results': out, 'path': path, 'stopped_at': u}
        if abs(u - targets[0]) < 1e-9:
            out[targets[0]] = prev
        h = step
        while targets and u < targets[-1] - 1e-9:
            nxt = min([t for t in targets if t > u + 1e-9] + [u + h])
            trial = self.run(nelec, nxt, j_over_u * nxt, Z0=prev['Z'], shift0=prev['shift'], lam0_guess=prev['lam'],
                             warmup=2, **kw)
            path.append((nxt, trial['converged']))
            if trial['converged']:
                prev, u = trial, nxt
                if any(abs(u - t) < 1e-9 for t in targets):
                    out[min(targets, key=lambda t: abs(t - u))] = prev
                h = min(step, h * 1.5)
            else:
                h = (nxt - u) / 2
                if h < min_step:
                    return {'results': out, 'path': path, 'stopped_at': nxt}
        return {'results': out, 'path': path, 'stopped_at': None}

    def solve(self, nelec, U, J, start=0.5, step=0.15, min_step=0.01, **kw):
        """Metallic solution continued from weak coupling. U and J are ramped together from start*U, each
        stage starting from the converged previous one; a stage that fails to converge is retried with half
        the step. Returns the final stage plus the ramp record."""
        kw.setdefault('maxiter', 150)
        record = []
        f = start
        prev = self.run(nelec, f * U, f * J, **kw)
        record.append({'U': f * U, 'converged': prev['converged'], 'iterations': prev['iterations']})
        if not prev['converged']:
            prev['ramp'] = record
            return prev
        while f < 1.0 - 1e-12:
            fn = min(1.0, f + step)
            trial = self.run(nelec, fn * U, fn * J, Z0=prev['Z'], shift0=prev['shift'], lam0_guess=prev['lam'],
                             warmup=2, **kw)
            record.append({'U': fn * U, 'converged': trial['converged'], 'iterations': trial['iterations']})
            if trial['converged']:
                prev, f = trial, fn
            else:
                step /= 2
                if step < min_step:
                    trial['ramp'] = record
                    return trial
        prev['ramp'] = record
        return prev


# ---------------------------------------------------------------- lab entry point
LAB_MESH = (9, 9, 6)
LAB_KT = 0.01


def lab_ssmf(model, nelec, U, j_over_u, solver=None):
    """The lab's SSMF evaluation: metallic solution continued from weak coupling (U = 1 eV, steps of 0.5 eV at
    fixed J/U, halved on failure down to 0.02 eV) on the 9x9x6 mesh. Returns a dict; 'metallic' is False when the
    metallic branch cannot be continued to the requested U."""
    s = solver or SSMF(model, mesh=LAB_MESH, kT=LAB_KT)
    if U <= 1.0:
        r = s.run(nelec, U, j_over_u * U, maxiter=150)
        ok, res = r['converged'], r
    else:
        out = s.scan(nelec, j_over_u, [U], u_start=1.0, step=0.5, min_step=0.02)
        ok = U in out['results']
        res = out['results'].get(U)
    if not ok:
        return {'metallic': False}
    return {'metallic': True, 'Z': [float(z) for z in res['Z']], 'n': [float(x) for x in res['n']],
            'mu': float(res['mu']), 'residual': float(res['residual'])}
