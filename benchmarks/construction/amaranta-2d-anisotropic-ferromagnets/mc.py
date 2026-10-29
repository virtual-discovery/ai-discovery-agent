"""Classical Monte Carlo of the teacher model (unit-vector spins), numpy only.

    python3 -B mc.py COMPOUND L T1,T2,... SWEEPS [seed]      (temperatures in K)

Metropolis with adaptive cone moves plus over-relaxation (reflection about the local exchange
field, accepted with the Metropolis probability of the single-ion energy change); sites are
updated colour class by colour class (no two coupled sites share a colour), so each class is
updated in one vectorised step. Prints, per
temperature, <|m|>, <m_z^2>, <m_z^4>, the Binder cumulant U = 1 - <m_z^4>/(3 <m_z^2>^2),
energy and specific heat per site.
"""
import json
import sys
import numpy as np
import model as M

K_PER_MEV = 11.604518


class Lattice:
    def __init__(self, comp, L1, L2, table=None, K=None):
        self.c = comp
        tab = comp.bond_table() if table is None else table
        self.K = comp.K if K is None else K
        nb = len(comp.basis)
        self.N = N = nb * L1 * L2
        idx = lambda k, i, j: (k * L1 + i % L1) * L2 + j % L2
        nbr = [[] for _ in range(N)]
        Jm = [[] for _ in range(N)]
        for k0, k, m1, m2, J in tab:
            for i in range(L1):
                for j in range(L2):
                    a, b = idx(k0, i, j), idx(k, i + m1, j + m2)
                    if a == b:
                        raise ValueError('self-coupling: supercell too small')
                    nbr[a].append(b)
                    Jm[a].append(J)
        z = max(len(x) for x in nbr)
        self.nbr = np.zeros((N, z), int)
        self.J = np.zeros((N, z, 3, 3))
        for a in range(N):
            for t, (b, J) in enumerate(zip(nbr[a], Jm[a])):
                self.nbr[a, t] = b
                self.J[a, t] = J
            for t in range(len(nbr[a]), z):
                self.nbr[a, t] = a           # padding with zero coupling
        # greedy colouring
        color = -np.ones(N, int)
        for a in range(N):
            used = {color[b] for b in nbr[a] if color[b] >= 0}
            cidx = 0
            while cidx in used:
                cidx += 1
            color[a] = cidx
        self.classes = [np.where(color == cc)[0] for cc in range(color.max() + 1)]

    def field(self, S, sites):
        return np.einsum('nzab,nzb->na', self.J[sites], S[self.nbr[sites]])

    def energy(self, S):
        h = self.field(S, np.arange(self.N))
        return 0.5 * np.sum(S * h) + self.K * np.sum(S[:, 2] ** 2)


def random_unit(rng, n):
    v = rng.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def run(lat, T_K, sweeps, rng, S=None, therm=None, axis=2, n_or=2, nbin=20):
    T = T_K / K_PER_MEV
    N = lat.N
    if S is None:
        S = np.tile([0, 0, 1.0], (N, 1)) if axis == 2 else random_unit(rng, N)
    width = 0.5
    therm = sweeps // 5 if therm is None else therm
    series = []
    acc_tot = 0
    for sw in range(therm + sweeps):
        acc = 0
        for _ in range(n_or):
            for sites in lat.classes:
                old = S[sites]
                h = lat.field(S, sites)
                hn = h / np.maximum(np.linalg.norm(h, axis=1, keepdims=True), 1e-12)
                prop = 2 * np.sum(old * hn, axis=1, keepdims=True) * hn - old
                dE = lat.K * (prop[:, 2] ** 2 - old[:, 2] ** 2)
                ok = (dE <= 0) | (rng.random(len(sites)) < np.exp(-np.clip(dE, 0, None) / T))
                S[sites[ok]] = prop[ok]
        for sites in lat.classes:
            old = S[sites]
            prop = old + width * random_unit(rng, len(sites))
            prop /= np.linalg.norm(prop, axis=1, keepdims=True)
            h = lat.field(S, sites)
            dE = np.sum((prop - old) * h, axis=1) + lat.K * (prop[:, 2] ** 2 - old[:, 2] ** 2)
            ok = (dE <= 0) | (rng.random(len(sites)) < np.exp(-np.clip(dE, 0, None) / T))
            S[sites[ok]] = prop[ok]
            acc += ok.sum()
        rate = acc / N
        if sw < therm:
            width *= 1.05 if rate > 0.5 else 0.95
            width = min(max(width, 0.02), 2.0)
            continue
        m = S.mean(axis=0)
        series.append((m[axis], np.linalg.norm(m), lat.energy(S) / N))
        acc_tot += rate
    x = np.array(series)
    mz, mn, e = x[:, 0], x[:, 1], x[:, 2]
    def stats(f):
        """mean and binning error of a function of the series arrays"""
        val = f(mz, mn, e)
        bins = np.array_split(np.arange(len(mz)), nbin)
        vb = np.array([f(mz[b], mn[b], e[b]) for b in bins])
        return float(val), float(vb.std(ddof=1) / np.sqrt(nbin))
    U = stats(lambda a, b, c: 1 - np.mean(a ** 4) / (3 * np.mean(a ** 2) ** 2))
    chi = stats(lambda a, b, c: N * (np.mean(a ** 2) - np.mean(np.abs(a)) ** 2) / T)
    C = stats(lambda a, b, c: N * np.var(c) / T ** 2)
    return S, {'T_K': T_K, 'N': N, 'mabs_z': float(np.mean(np.abs(mz))), 'mnorm': float(mn.mean()), 'U': U[0], 'U_err': U[1],
               'chi': chi[0], 'chi_err': chi[1], 'e': float(e.mean()), 'C': C[0], 'C_err': C[1], 'acc': acc_tot / len(mz)}


def main():
    name, L = sys.argv[1], int(sys.argv[2])
    temps = [float(x) for x in sys.argv[3].split(',')]
    sweeps = int(sys.argv[4])
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    comp = M.Compound(name, M.load_pool()[name])
    lat = Lattice(comp, L, L)
    rng = np.random.default_rng(seed)
    S = None
    out = []
    for T in temps:
        S, r = run(lat, T, sweeps, rng, S)
        out.append(r)
        print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()}), flush=True)


if __name__ == '__main__':
    main()
