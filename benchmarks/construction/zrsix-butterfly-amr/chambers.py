"""Chambers (Boltzmann, constant relaxation time) magnetoconductivity of the ZrSiX Wannier models.

sigma_ij(B)/tau = (2 e^2 / (V N_k)) sum_{n,k} (-df/dE) <v_i vbar_j>_orbit(k),
vbar_j(k) = int_0^inf (du/tau) exp(-u/tau) v_j(k(-u)),   hbar dk/dt = -e v x B  (electron charge -e; spin factor 2).
The past trajectory k(-u) obeys dk/du = +(e/hbar) v x B. With s = B u (T ps) one trajectory per start point gives
vbar for every B tau = X: vbar(X) = int_0^inf ds/X exp(-s/X) v(s).

Each start point of the k mesh (inside the thermal window) seeds its orbit; its contribution is replaced by the
average of v_i vbar_j over that orbit (the Shockley tube integral), which has the same expectation as the single
start point and is free of the start-point sampling noise. Closed orbits (returning to the start modulo a reciprocal
lattice vector G; G != 0 means a periodic open orbit) are averaged exactly with a Fourier series in s; orbits that do
not close within s_max are averaged over the first half of the trajectory with vbar from a backward recursion.
Orbits are integrated with RK4 in k-space arc length (step dl, 1/Angstrom) with projection back onto the energy
surface.
"""
import numpy as np

HBAR = 6.582119569e-16       # eV s
E = 1.602176634e-19          # C
KB = 8.617333262e-5          # eV/K
S_PER_L = HBAR ** 2 * 1e32   # (T ps) per (1/Angstrom) of arc for |v_perp| = 1 eV*Angstrom


def _vel(spline, kf, bd):
    e, v = spline.eig_vel(kf)
    r = np.arange(len(kf))
    return e[r, bd], v[r, bd]


def integrate(spline, k0, band, bvec, s_max, dl=0.004, max_steps=8000):
    """RK4 orbits. Returns per-orbit histories s (M, T+1), v (M, T+1, 3), the number of steps, closure flags,
    G vectors and energies."""
    M = len(k0)
    B = np.linalg.inv(spline.Binv)
    Binv = spline.Binv
    b = np.asarray(bvec, float) / np.linalg.norm(bvec)
    k = np.array(k0, float)
    e0, v0 = _vel(spline, k, band)
    S = [np.zeros(M)]; V = [v0.copy()]
    nsteps = np.zeros(M, int)
    closed = np.zeros(M, bool); G = np.zeros((M, 3), int)
    active = np.ones(M, bool)
    s = np.zeros(M); v = v0.copy()

    def f(vv):
        w = np.cross(vv, b)
        n = np.maximum(np.linalg.norm(w, axis=1), 1e-3)
        return (w / n[:, None]) @ Binv, n
    step = 0
    dl0 = dl
    while active.any() and step < max_steps:
        dl = dl0 if step < 1500 else 2 * dl0          # long (heavy or open) orbits: coarser arc-length step
        a = np.where(active)[0]
        ka, va, ba = k[a], v[a], band[a]
        d1, n1 = f(va)
        _, v2 = _vel(spline, ka + 0.5 * dl * d1, ba); d2, n2 = f(v2)
        _, v3 = _vel(spline, ka + 0.5 * dl * d2, ba); d3, n3 = f(v3)
        _, v4 = _vel(spline, ka + dl * d3, ba); d4, n4 = f(v4)
        kn = ka + dl * (d1 + 2 * d2 + 2 * d3 + d4) / 6
        en, vn = _vel(spline, kn, ba)
        kn = kn - (((en - e0[a]) / np.maximum((vn ** 2).sum(1), 1e-6))[:, None] * vn) @ Binv
        en, vn = _vel(spline, kn, ba)
        ds = S_PER_L * dl * (1 / n1 + 2 / n2 + 2 / n3 + 1 / n4) / 6
        k[a], v[a], s[a] = kn, vn, s[a] + ds
        step += 1
        snew = np.full(M, np.nan); snew[a] = s[a]
        vnew = np.full((M, 3), np.nan); vnew[a] = vn
        S.append(snew); V.append(vnew)
        nsteps[a] = step
        delta = kn - k0[a]
        Gn = np.round(delta)
        dist = np.linalg.norm((delta - Gn) @ B, axis=1)
        done = (dist < 0.75 * dl) & (step > 8)
        if done.any():
            i = a[done]
            closed[i] = True; G[i] = Gn[done].astype(int); active[i] = False
        over = active[a] & (s[a] >= s_max)
        active[a[over]] = False
    return dict(s=np.array(S).T, v=np.transpose(np.array(V), (1, 0, 2)), nsteps=nsteps, closed=closed, G=G,
                energy=e0, v0=v0)


def orbit_average(tr, X, nfft_min=64):
    """<v_i vbar_j> over each orbit, (M, nX, 3, 3) in (eV Angstrom)^2, for each B tau in X (T ps)."""
    X = np.asarray(X, float)
    M = len(tr['nsteps'])
    out = np.zeros((M, len(X), 3, 3))
    period = np.full(M, np.nan)
    for i in np.where(tr['closed'])[0]:
        n = tr['nsteps'][i]
        s = tr['s'][i, :n + 1]; v = tr['v'][i, :n + 1]
        Sp = s[-1]
        period[i] = Sp
        # the last point closes the loop (~ the start); resample uniformly in s over one period
        N = max(nfft_min, 1 << int(np.ceil(np.log2(2 * n))))
        su = np.arange(N) * Sp / N
        vu = np.stack([np.interp(su, s, v[:, c]) for c in range(3)], -1)
        a = np.fft.fft(vu, axis=0) / N                              # (N, 3)
        q = 2 * np.pi * np.fft.fftfreq(N, d=Sp / N)                   # 1/(T ps)
        corr = np.einsum('ni,nj->nij', a.conj(), a)                   # (N, 3, 3)
        filt = 1.0 / (1.0 - 1j * q[None, :] * X[:, None])             # (nX, N)
        out[i] = np.einsum('xn,nij->xij', filt, corr).real
    op = np.where(~tr['closed'])[0]
    if len(op):
        s = tr['s'][op]; v = tr['v'][op]; T = s.shape[1] - 1
        half = np.zeros(len(op), int)
        for t, i in enumerate(op):
            n = tr['nsteps'][i]
            half[t] = max(1, n // 2)
        for x_i, Xv in enumerate(X):
            if Xv == 0:
                acc = np.zeros((len(op), 3, 3))
                for t, i in enumerate(op):
                    n = half[t]
                    w = np.diff(s[t, :n + 1])
                    vv = v[t, :n]
                    acc[t] = np.einsum('m,mi,mj->ij', w, vv, vv) / w.sum()
                out[op, x_i] = acc
                continue
            # backward recursion vbar(s_m) = int_m^{m+1} + exp(-ds/X) vbar(s_{m+1}), started at each orbit's end
            vbar = np.zeros((len(op), T + 1, 3))
            n_all = np.array([tr['nsteps'][i] for i in op])
            for m in range(T - 1, -1, -1):
                live = n_all > m
                if not live.any():
                    continue
                sa, sb = s[live, m], s[live, m + 1]
                va, vb = v[live, m], v[live, m + 1]
                dd = (sb - sa) / Xv
                dE = -np.expm1(-dd)
                lin = -(1 - dE) + dE / np.maximum(dd, 1e-300)
                nxt = np.where((n_all[live] > m + 1)[:, None], vbar[live, m + 1], vb)
                vbar[live, m] = va * dE[:, None] + (vb - va) * lin[:, None] + (1 - dE)[:, None] * nxt
            acc = np.zeros((len(op), 3, 3))
            for t in range(len(op)):
                n = half[t]
                w = np.diff(s[t, :n + 1])
                acc[t] = np.einsum('m,mi,mj->ij', w, v[t, :n], vbar[t, :n]) / w.sum()
            out[op, x_i] = acc
    return out, period


def sigma_over_tau(avg, weights, volume_A3, nk):
    """sigma/tau (nX, 3, 3) in S/(m s) from orbit averages (M, nX, 3, 3) in (eV A)^2 and -df/dE weights (1/eV)."""
    to_ms = 1e-10 / HBAR
    pref = 2 * E ** 2 / (volume_A3 * 1e-30 * nk) / E * to_ms ** 2
    return pref * np.einsum('m,mxij->xij', weights, avg)
