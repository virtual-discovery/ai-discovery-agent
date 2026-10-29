"""Construction library for mpea-vacancy-diffusion-kmc: bcc geometry, barrier-model training, kMC driver.

Barrier model (the builders' reference, also the class the lab kMC service accepts):
  KRA(X, env) = K0[X] + sum_o mean_{sites in orbit o} W[X][o][species] + sum_ij Q[X][i][j] c_i c_j,
      orbits = classes of sites around the hop pair by their (sorted) distances to the two hop sites,
      c = composition of orbits 1-3 (the 14 nearest neighbours of the pair);
  end states: pair Hamiltonian E = sum_shells sum_bonds V[s][a][b] (species V Cr Nb Mo Ta W Va);
  hop barrier = KRA + (E_final - E_initial)/2.
Trained per hopping species by ridge regression on the archived KRA (E_TS - (E1+E2)/2); the pair Hamiltonian by
least squares on archived end-state total energies (with element point terms).
"""
import json
import os
import subprocess
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/mpea-vacancy-diffusion-kmc'
EL = ['V', 'Cr', 'Nb', 'Mo', 'Ta', 'W']
EI = {e: i for i, e in enumerate(EL)}
VA = 6
KB = 8.617333262e-5
EXE = os.environ.get('KMC_EXE', str(HERE/'kmc'))

# ---------------------------------------------------------------- geometry (positions in units of a)
NN_DIRS = np.array([(sx * .5, sy * .5, sz * .5) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)])
PTS = np.array([(i + b * .5, j + b * .5, k + b * .5) for i in range(-5, 6) for j in range(-5, 6) for k in range(-5, 6)
                for b in (0, 1)])
SHELL_R = np.array(sorted(set(np.round(np.linalg.norm(PTS, axis=1), 3))))[1:]


def hop_orbits(maxd=1.732):
    """Orbit classes for a hop from the origin (vacancy) to (1/2,1/2,1/2): sites keyed by sorted distances."""
    B = np.array([.5, .5, .5])
    dA, dB = np.round(np.linalg.norm(PTS, axis=1), 3), np.round(np.linalg.norm(PTS - B, axis=1), 3)
    cls = {}
    for p, a, b in zip(PTS, dA, dB):
        if a < 1e-6 or b < 1e-6:
            continue
        key = tuple(sorted((a, b)))
        if key[0] <= maxd + 1e-6:
            cls.setdefault(key, []).append(p)
    keys = sorted(cls, key=lambda k: (k[0], k[1]))
    return [(k, np.array(sorted(cls[k], key=tuple))) for k in keys]


def shell_offsets(n):
    d = np.round(np.linalg.norm(PTS, axis=1), 3)
    return [np.array(sorted(PTS[d == SHELL_R[s]], key=tuple)) for s in range(n)]


ORBITS = hop_orbits()             # 12 orbits within sqrt(3) a of the nearer hop site
SHELLS = shell_offsets(4)


def half_units(arr):
    return np.round(np.asarray(arr) * 2).astype(int)


class Cell:
    def __init__(self, L):
        self.L = L
        self.pos = np.array([(i + b * .5, j + b * .5, k + b * .5) for i in range(L) for j in range(L) for k in range(L)
                             for b in (0, 1)])
        self.N = len(self.pos)

    def idx(self, p):
        L = self.L
        q = np.round(2 * np.asarray(p)).astype(int) % (2 * L)
        b = q[..., 0] % 2
        return (((q[..., 0] - b) // 2 * L + (q[..., 1] - b) // 2) * L + (q[..., 2] - b) // 2) * 2 + b

    def nn_table(self):
        return np.stack([self.idx(self.pos + d) for d in NN_DIRS], 1)

    def orbit_table(self, orbits):
        out = []
        for _, offs in orbits:
            arr = np.zeros((self.N, 8, len(offs)), int)
            for d, delta in enumerate(NN_DIRS):
                o = offs * np.sign(delta)
                arr[:, d, :] = np.stack([self.idx(self.pos + x) for x in o], 1)
            out.append(arr)
        return out

    def shell_table(self, shells):
        return [np.stack([self.idx(self.pos + o) for o in offs], 1) for offs in shells]


# ---------------------------------------------------------------- data
POOL = json.loads((HERE/'pool.json').read_text())
NH = len(POOL)
OCC = np.array([[EI.get(x, VA) for x in h['occupation_initial']] for h in POOL])
V1 = np.array([h['vacancy_site_initial'] for h in POOL])
V2 = np.array([h['vacancy_site_final'] for h in POOL])
SP = np.array([EI[h['hopping_species']] for h in POOL])
E1 = np.array([h['energy_eV']['end_state_1'] for h in POOL])
E2 = np.array([h['energy_eV']['end_state_2'] for h in POOL])
ET = np.array([h['energy_eV']['transition_state'] for h in POOL])
KRA = ET - (E1 + E2) / 2
PATHO = np.minimum(ET - E1, ET - E2) < 0          # saddle below an end state
CELL4 = Cell(4)
DIRS = np.array([list(CELL4.nn_table()[a]).index(b) for a, b in zip(V1, V2)])
OCC2 = OCC.copy()
OCC2[np.arange(NH), V1] = SP
OCC2[np.arange(NH), V2] = VA
PURE_IDX = [int(np.where((OCC == s).sum(1) == 127)[0][0]) for s in range(6)]
PURE_KRA = KRA[PURE_IDX]
_cache = {}


def orbit_frac(n_orb):
    if ('o', n_orb) not in _cache:
        T = CELL4.orbit_table(ORBITS[:n_orb])
        _cache[('o', n_orb)] = np.stack([np.stack([(OCC[np.arange(NH)[:, None], o[V1, DIRS]] == e).mean(1)
                                                   for e in range(6)], 1) for o in T], 1)
    return _cache[('o', n_orb)]


IJ = [(i, j) for i in range(6) for j in range(i, 6)]


def nnfrac():
    C = orbit_frac(3)
    return (C[:, 0] * 6 + C[:, 1] * 6 + C[:, 2] * 2) / 14


def train_kra(mask, n_orb=6, lam=0.03, quad=True, target=None):
    C = orbit_frac(n_orb)
    y = KRA if target is None else target
    K0, W, Q = np.zeros(6), np.zeros((6, n_orb, 6)), np.zeros((6, 6, 6))
    for s in range(6):
        idx = np.where(mask & (SP == s))[0]
        X = [np.ones(len(idx)), C[idx].reshape(len(idx), -1)]
        if quad:
            nn = nnfrac()[idx]
            X.append(np.stack([nn[:, i] * nn[:, j] for i, j in IJ], 1))
        X = np.column_stack(X)
        reg = lam * np.eye(X.shape[1])
        reg[0, 0] = 0
        w = np.linalg.solve(X.T @ X + reg, X.T @ y[idx])
        K0[s], W[s] = w[0], w[1:1 + 6 * n_orb].reshape(n_orb, 6)
        if quad:
            for q, (i, j) in enumerate(IJ):
                Q[s, i, j] = w[1 + 6 * n_orb + q]
    return K0, W, Q


def predict_kra(K0, W, Q, idx):
    C = orbit_frac(W.shape[1])[idx]
    nn = nnfrac()[idx]
    return K0[SP[idx]] + np.einsum('noe,noe->n', W[SP[idx]], C) + np.einsum('nij,ni,nj->n', Q[SP[idx]], nn, nn)


PAIRS = [(a, b) for a in range(7) for b in range(a, 7)]


def _bond_counts(O, nsh):
    out = []
    for t in CELL4.shell_table(SHELLS[:nsh]):
        B = O[:, t]
        cnt = np.zeros((len(O), 7, 7))
        for a in range(7):
            ma = (O == a)[:, :, None]
            for b in range(7):
                cnt[:, a, b] = ((B == b) & ma).sum((1, 2))
        out.append(np.stack([(cnt[:, a, b] + cnt[:, b, a]) / 2 if a != b else cnt[:, a, a] / 2 for a, b in PAIRS], 1))
    return np.stack(out, 1)


def ce_feats(nsh):
    if ('c', nsh) not in _cache:
        _cache[('c', nsh)] = (np.stack([(OCC == e).sum(1) for e in range(6)], 1), _bond_counts(OCC, nsh),
                              np.stack([(OCC2 == e).sum(1) for e in range(6)], 1), _bond_counts(OCC2, nsh))
    return _cache[('c', nsh)]


def train_ce(mask, nsh=3, lam=0.1, mode='total'):
    pt1, c1, pt2, c2 = ce_feats(nsh)
    if mode == 'total':
        X = np.r_[np.c_[np.ones(NH), pt1, c1.reshape(NH, -1)], np.c_[np.ones(NH), pt2, c2.reshape(NH, -1)]]
        y, m, npt = np.r_[E1, E2], np.r_[mask, mask], 7
        X, y = X[m], y[m]
    else:
        X, y, npt = (c2 - c1).reshape(NH, -1)[mask], (E2 - E1)[mask], 0
    keep = X.std(0) > 1e-12
    keep[:npt] = True
    reg = lam * np.eye(keep.sum())
    reg[:npt, :npt] = 0
    w = np.zeros(X.shape[1])
    w[keep] = np.linalg.solve(X[:, keep].T @ X[:, keep] + reg, X[:, keep].T @ y)
    V = np.zeros((nsh, 7, 7))
    for s in range(nsh):
        for q, (a, b) in enumerate(PAIRS):
            V[s, a, b] = V[s, b, a] = w[npt + s * len(PAIRS) + q]
    return V


def predict_de(V, idx):
    pt1, c1, pt2, c2 = ce_feats(V.shape[0])
    w = np.array([[V[s, a, b] for a, b in PAIRS] for s in range(V.shape[0])])
    return np.einsum('nsq,sq->n', (c2 - c1)[idx], w)


# ---------------------------------------------------------------- kMC
def counts_for(x, L):
    n = 2 * L ** 3 - 1
    x = np.asarray(x, float) / np.sum(x)
    raw = x * n
    c = np.floor(raw).astype(int)
    order = sorted(range(6), key=lambda i: (-(raw[i] - c[i]), i))
    for i in order[:n - c.sum()]:
        c[i] += 1
    return c


def kmc_input(L, counts, T, nsteps, ntraj, seed, model, nmc=0, nck=0):
    K0, W, Q, V = model['K0'], model['W'], model['Q'], model['V']
    n_orb = W.shape[1]
    out = [f"{L} {nsteps} {ntraj} {seed} {float(T)} 1.0 {nmc} {nck} {0 if V is None else 1}",
           ' '.join(str(int(c)) for c in counts), ' '.join(repr(float(x)) for x in K0), str(n_orb)]
    for _, off in ORBITS[:n_orb]:
        out.append(f"{len(off)} " + ' '.join(str(int(x)) for x in half_units(off).ravel()))
    out.append(' '.join(repr(float(x)) for x in W.ravel()))
    q = Q is not None and np.any(Q)
    out.append(str(3 if q else 0))
    out.append(' '.join(repr(float(x)) for x in (Q if q else np.zeros((6, 6, 6))).ravel()))
    nsh = 0 if V is None else V.shape[0]
    out.append(str(nsh))
    for off in SHELLS[:nsh]:
        out.append(f"{len(off)} " + ' '.join(str(int(x)) for x in half_units(off).ravel()))
    out.append(' '.join(repr(float(x)) for x in V.ravel()) if nsh else '')
    return '\n'.join(out) + '\n'


def run_kmc(text):
    p = subprocess.run([EXE], input=text, capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError(p.stderr)
    T = [[float(x) for x in line.split()[1:]] for line in p.stdout.splitlines() if line.startswith('T ')]
    return np.array(T)


def summarize(rows, counts, N, x, T):
    """rows: trajectory lines (tr t vacR2 r2[6] hops[6] enc acc mcacc); R^2 in (a/2)^2, t in ps."""
    t, vac, r2, hops = rows[:, 1], rows[:, 2], rows[:, 3:9], rows[:, 9:15]
    Dv = vac.sum() / (6 * t.sum()) / 4.0                          # a^2/ps
    cnt = np.asarray(counts, float)
    with np.errstate(invalid='ignore', divide='ignore'):
        Dt = r2.sum(0) / (6 * t.sum()) / 4.0 / cnt * N            # a^2/ps per unit vacancy site fraction
        Dj = hops.sum(0) * 3.0 / (6 * t.sum()) / 4.0 / cnt * N     # uncorrelated (jump-count) estimate
        f = r2.sum(0) / (hops.sum(0) * 3.0)
    xx = np.asarray(x, float) / np.sum(x)
    rom = np.exp(np.sum(xx * (-PURE_KRA / (KB * T))))              # a^2 nu exp(-E/kT) with nu = 1/ps
    return {'Dv': Dv, 'Dt': Dt, 'Dj': Dj, 'f': f, 'fv': vac.sum() / (hops.sum() * 3.0), 'log10_eta': np.log10(Dv / rom),
            'enc': rows[:, 15].mean(), 'acc': rows[:, 16].mean(), 't': t.mean()}


def simulate(x, T, model, L=10, nsteps=25000, ntraj=80, nmc=100, seed=1, workers=4):
    counts = counts_for(x, L)
    chunks = [ntraj // workers + (1 if i < ntraj % workers else 0) for i in range(workers)]
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(workers) as pool:
        parts = list(pool.map(lambda i: run_kmc(kmc_input(L, counts, T, nsteps, chunks[i], 7919 * seed + 104729 * i,
                                                          model, nmc=nmc)), [i for i in range(workers) if chunks[i]]))
    return summarize(np.concatenate(parts), counts, 2 * L ** 3, x, T)
