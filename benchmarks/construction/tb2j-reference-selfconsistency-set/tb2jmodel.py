"""Parsing and isotropic Heisenberg/LSWT utilities for tb2j-reference-selfconsistency-set (numpy/scipy; runtime image).

Convention (TB2J, as archived): E = -sum_{i != j} J_ij e_i.e_j over ordered pairs (each bond listed twice), unit vectors
along the local moments, J > 0 ferromagnetic, meV. Wave vectors are in reciprocal-lattice units (r.l.u.) of the cell
named in each function.

`Model.h` is the authors' bosonic matrix h(k) (paper Eq. 9, Toth-Lake, unit spin length, meV) around an arbitrary
(collinear or noncollinear) moment configuration of the calculation's own cell. For the FM reference its lowest
eigenvalue as a function of k defines the paper's k_min.
"""
import copy
import re
from collections import defaultdict
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize

RAW = Path.home()/'.cache/ai-theorist/tb2j-reference-selfconsistency-set/raw'
CALCS = {  # task label -> archive member
    'NiO/FM': 'NiO/FM/exchange.out', 'NiO/AFM': 'NiO/AFM/exchange.out',
    'FeCl2/FM': 'FeCl2/FM/exchange.out', 'FeCl2/AFM': 'FeCl2/AFM/exchange.out',
    'FePS3/DFT/FM': 'FePS3/DFT/FM/exchange.out', 'FePS3/DFT/AFM': 'FePS3/DFT/AFM/exchange.out',
    'FePS3/DFT+U/FM': 'FePS3/DFT+U/FM/exchange.out', 'FePS3/DFT+U/AFM': 'FePS3/DFT+U/AFM/exchange.out',
    'FeP/FM': 'FeP/FM/exchange.out', 'FeP/AFM': 'FeP/AFM/exchange.out'}
MATERIALS = {'NiO': ('NiO/FM', 'NiO/AFM'), 'FeCl2': ('FeCl2/FM', 'FeCl2/AFM'),
             'FePS3_DFT': ('FePS3/DFT/FM', 'FePS3/DFT/AFM'), 'FePS3_DFT+U': ('FePS3/DFT+U/FM', 'FePS3/DFT+U/AFM'),
             'FeP': ('FeP/FM', 'FeP/AFM')}
MAGNETIC = ('Ni', 'Fe')
PAIR = re.compile(r'\s+(\w+)\s+(\w+)\s+\(\s*(-?\d+),\s*(-?\d+),\s*(-?\d+)\)\s+(-?[\d.]+)\s+'
                  r'\(\s*(-?[\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+)\)\s+([\d.]+)')


def parse_exchange(path):
    """Header lines, cell, atoms and isotropic pairs of a TB2J exchange.out (DMI/J_ani lines are ignored)."""
    lines = Path(path).read_text().split('\n')
    header, cell, atoms, pairs, columns = [], None, [], [], None
    i = 0
    while i < len(lines):
        line = lines[i]
        if cell is None and not line.startswith('Cell (Angstrom)'):
            if line.strip() and not line.startswith('====') and line.strip() != 'Information:':
                header.append(line.rstrip())
        if line.startswith('Cell (Angstrom)'):
            cell = [[float(x) for x in lines[i + k].split()] for k in (1, 2, 3)]
            i += 4
            continue
        if line.startswith('Atom') and 'w_charge' in line:
            tokens = line.split()
            columns = tokens[tokens.index('z') + 1:]
            i += 1
            while not lines[i].startswith('Total') and not lines[i].startswith('==='):
                p = lines[i].split()
                if p:
                    atoms.append({'name': p[0], 'position_A': [float(x) for x in p[1:4]],
                                  **{c: float(v) for c, v in zip(columns, p[4:])}})
                i += 1
            continue
        m = PAIR.match(line)
        if m:
            pairs.append([m.group(1), m.group(2), int(m.group(3)), int(m.group(4)), int(m.group(5)), float(m.group(6)),
                          float(m.group(7)), float(m.group(8)), float(m.group(9)), float(m.group(10))])
        i += 1
    return {'header': header, 'cell_A': cell, 'atom_columns': columns, 'atoms': atoms, 'pairs': pairs}


def load(label):
    return parse_exchange(RAW/CALCS[label])


def subset(d, keep):
    d2 = dict(d)
    d2['pairs'] = [p for p in d['pairs'] if keep(p)]
    return d2


def minimage(d, axis=2):
    """Re-assign each pair to its minimum image along one supercell axis (pairs listed at R_axis = 0 only)."""
    d2 = copy.deepcopy(d)
    cell = np.array(d['cell_A'])
    inv = np.linalg.inv(cell)
    for p in d2['pairs']:
        v = np.array(p[6:9])
        n = -int(np.round((v @ inv)[axis]))
        if n:
            v = v + n*cell[axis]
            p[2 + axis] += n
            p[6:9] = list(v)
            p[9] = float(np.linalg.norm(v))
    return d2


class Model:
    """Isotropic model of one calculation in its own cell and with its own pair list."""

    def __init__(self, d):
        mag = [a for a in d['atoms'] if a['name'].rstrip('0123456789') in MAGNETIC]
        self.names = [a['name'] for a in mag]
        self.n = len(mag)
        self.cell = np.array(d['cell_A'])
        self.pos = np.array([a['position_A'] for a in mag])
        m = np.array([[a['M(x)'], a['M(y)'], a['M(z)']] if 'M(z)' in a else [0, 0, a['w_magmom']] for a in mag])
        self.m = m
        self.e = m/np.linalg.norm(m, axis=1, keepdims=True)
        idx = {n: k for k, n in enumerate(self.names)}
        P = d['pairs']
        self.I = np.array([idx[p[0]] for p in P], int)
        self.J = np.array([idx[p[1]] for p in P], int)
        self.R = np.array([p[2:5] for p in P], float).reshape(-1, 3)
        self.V = np.array([p[5] for p in P], float)
        self._A = sp.csr_matrix((np.ones(len(P)), (np.arange(len(P)), self.I*self.n + self.J)),
                                shape=(len(P), self.n*self.n))

    def K(self, k):
        """K_ij(k) = sum_R J_ij(R) exp(2 pi i k.R), k in this cell's r.l.u."""
        k = np.atleast_2d(np.asarray(k, float))
        out = []
        for c in range(0, len(k), 256):
            ph = np.exp(2j*np.pi*(k[c:c + 256] @ self.R.T))*self.V
            out.append(np.asarray((self._A.T @ ph.T).T).reshape(-1, self.n, self.n))
        return np.concatenate(out)

    def K0(self):
        return self.K([0, 0, 0])[0].real

    def h(self, k, e=None):
        """Bosonic matrix (unit spin length, meV) around directions e (default: the archived moments)."""
        e = self.e if e is None else np.asarray(e, float)
        u, v = frames(e)
        K0 = self.K0()
        C = np.diag(-(K0*(v @ v.T)).sum(1))
        uu_c, uu = u @ u.conj().T, u @ u.T
        k = np.atleast_2d(np.asarray(k, float))
        out = []
        for Kq, Kmq in zip(self.K(k), self.K(-k)):
            A, B, Am = -0.5*Kq*uu_c, -0.5*Kq*uu, -0.5*Kmq*uu_c
            out.append(2*np.block([[A - C, B], [B.conj().T, Am.conj() - C]]))
        return np.array(out)

    def energy(self, e=None):
        e = self.e if e is None else np.asarray(e, float)
        return float(-np.einsum('ij,ij->', self.K0(), e @ e.T)/self.n)


def frames(e):
    u, v = [], []
    for z in e:
        z = z/np.linalg.norm(z)
        t = np.array([1.0, 0, 0]) if abs(z[0]) < 0.9 else np.array([0, 1.0, 0])
        x = t - (t @ z)*z
        x /= np.linalg.norm(x)
        u.append(x + 1j*np.cross(z, x))
        v.append(z)
    return np.array(u), np.array(v)


def grid(n, shift=0.0):
    g = (np.arange(n) + shift)/n
    return np.array(np.meshgrid(g, g, g, indexing='ij')).reshape(3, -1).T


def omegas(H):
    n = H.shape[1]//2
    g = np.diag(np.r_[np.ones(n), -np.ones(n)])
    ev = np.linalg.eigvals(g @ H)
    return np.sort(ev.real, axis=1)[:, n:], float(np.abs(ev.imag).max())


def energy(K0, e):
    return float(-np.einsum('ij,ij->', K0, e @ e.T)/len(e))


def relax(K0, e0, iters=20000, tol=1e-13):
    """Classical zero-temperature relaxation of unit vectors in their exchange fields (self-images excluded)."""
    e = np.array(e0, float)
    Kx = K0 - np.diag(np.diag(K0))
    old = None
    for _ in range(iters):
        for i in range(len(e)):
            f = Kx[i] @ e
            e[i] = f/np.linalg.norm(f)
        E = energy(K0, e)
        if old is not None and abs(E - old) < tol:
            break
        old = E
    return e, E


class FoldedModel(Model):
    """A supercell pair list mapped onto the primitive (FM-reference) cell; equivalent bonds averaged."""

    def __init__(self, acc, n, average='mean'):
        items = list(acc.items())
        self.n = n
        self.I = np.array([k[0] for k, v in items], int)
        self.J = np.array([k[1] for k, v in items], int)
        self.R = np.array([k[2] for k, v in items], float)
        self.V = np.array([np.mean(v) if average == 'mean' else np.sum(v)/average for k, v in items])
        self.e = np.tile([0, 0, 1.0], (n, 1))
        self._A = sp.csr_matrix((np.ones(len(items)), (np.arange(len(items)), self.I*self.n + self.J)),
                                shape=(len(items), self.n*self.n))


def site_map(fm, af):
    """Primitive-site index and FM-lattice translation of every supercell site; supercell matrix."""
    Msc = np.round(af.cell @ np.linalg.inv(fm.cell)).astype(int)
    inv_eff = np.linalg.inv(np.linalg.inv(Msc) @ af.cell)
    pf = fm.pos @ np.linalg.inv(fm.cell)
    pmap, tvec = [], []
    for fi in af.pos @ inv_eff:
        dd = fi - pf
        a = int(np.argmin(np.linalg.norm(dd - np.round(dd), axis=1)))
        assert np.linalg.norm(dd[a] - np.round(dd[a])) < 0.05
        pmap.append(a)
        tvec.append(np.round(fi - pf[a]).astype(int))
    return np.array(pmap), np.array(tvec), Msc, inv_eff, pf


def fold(fm, d_af):
    """Group the supercell pairs by primitive bond (a, b, T)."""
    af = Model(d_af)
    pmap, tvec, Msc, inv_eff, pf = site_map(fm, af)
    idx = {n: k for k, n in enumerate(af.names)}
    acc = defaultdict(list)
    for p in d_af['pairs']:
        a, b = pmap[idx[p[0]]], pmap[idx[p[1]]]
        dfrac = np.array(p[6:9]) @ inv_eff - (pf[b] - pf[a])
        T = np.round(dfrac).astype(int)
        assert np.abs(dfrac - T).max() < 0.05
        acc[(int(a), int(b), tuple(int(t) for t in T))].append(p[5])
    return acc


def unfold_K0(fm, af):
    """The FM-reference couplings placed on the supercell sites by lattice translation (real K(0) matrix)."""
    pmap, tvec, Msc, inv_eff, pf = site_map(fm, af)
    Minv = np.linalg.inv(Msc)

    def key(a, t):
        return int(a), tuple(np.round(((t @ Minv) % 1.0)*1000).astype(int) % 1000)
    lookup = {key(a, t): i for i, (a, t) in enumerate(zip(pmap, tvec))}
    K0 = np.zeros((af.n, af.n))
    for i in range(af.n):
        sel = fm.I == pmap[i]
        for b, R, V in zip(fm.J[sel], fm.R[sel].astype(int), fm.V[sel]):
            K0[i, lookup[key(b, tvec[i] + R)]] += V
    return K0


def lowest(model, k):
    return float(np.linalg.eigvalsh(model.h(k))[0, 0])


def kmin(model, n=10, starts=()):
    """k_min: minimum over k of the lowest eigenvalue of h(k) around the FM configuration (model's own cell)."""
    e_fm = np.tile([0, 0, 1.0], (model.n, 1))
    low = lambda k: float(np.linalg.eigvalsh(model.h(np.atleast_2d(k), e=e_fm))[0, 0])
    ks = grid(n)
    ev = np.linalg.eigvalsh(model.h(ks, e=e_fm))[:, 0]
    best = None
    for s in list(ks[np.argsort(ev)[:6]]) + [np.array(x, float) for x in starts]:
        r = minimize(low, s, method='Nelder-Mead', options={'xatol': 1e-7, 'fatol': 1e-10, 'maxiter': 4000})
        if best is None or r.fun < best.fun:
            best = r
    return best.x - np.round(best.x), float(best.fun)


def qabs(q):
    q = np.asarray(q, float)
    return np.abs(q - np.round(q))


def qdist(q, ref):
    """Distance in r.l.u. after reducing each component modulo 1 and taking absolute values (sign/mirror equivalence)."""
    return float(np.linalg.norm(qabs(q) - qabs(ref)))


def stabilization(fm_data, af_data, e=None):
    """(E_FM - E_conf) per magnetic ion from the FM-reference and from the converged (AFM-reference) exchange."""
    fm, af = Model(fm_data), Model(af_data)
    e = af.e if e is None else e
    ones = np.tile([0, 0, 1.0], (af.n, 1))
    Kf, Ka = unfold_K0(fm, af), af.K0()
    return energy(Kf, ones) - energy(Kf, e), energy(Ka, ones) - energy(Ka, e)
