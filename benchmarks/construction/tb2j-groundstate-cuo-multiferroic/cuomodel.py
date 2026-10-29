"""Parsing and Heisenberg/LSWT utilities for the CuO task construction (numpy/scipy only; runtime image).

Convention (TB2J, as archived): E = -sum_{i != j} J_ij e_i.e_j over ordered pairs, unit vectors, J > 0 ferromagnetic.
Wave vectors are in reciprocal-lattice units (r.l.u.) of the conventional 4-Cu cell of the FM-reference calculation.

The authors' bosonic matrix h(k) (paper Eq. 9, Toth-Lake with S = 1, J in meV, values archived in eV) is
reproduced by `Native.h`; physical magnon energies for spin length S are the positive eigenvalues of g.h / S.
"""
import json
import re
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

RAW = Path.home()/'.cache/ai-theorist/tb2j-groundstate-cuo-multiferroic/raw'
SETS = {'FM/ligand': 'ligand/FM', 'AFM/ligand': 'ligand/AFM', 'FM/no_ligand': 'noligand/FM', 'AFM/no_ligand': 'noligand/AFM'}
DISPERSION = {'FM/ligand': 'ligand/FM/dispersion_FM', 'AFM/ligand': 'ligand/AFM/dispersion_AFM',
              'FM/no_ligand': 'noligand/FM/dispersion_nol_FM', 'AFM/no_ligand': 'noligand/AFM/dispersion_nol_AFM'}
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
            columns = tokens[tokens.index('z') + 1:]   # header is 'Atom number' or 'Atom_number' then x y z ...
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


def load(key):
    return parse_exchange(RAW/SETS[key]/'exchange.out')


def dispersion(key):
    path = RAW/DISPERSION[key]
    text = path.read_text().split('\n')
    return text[0], np.loadtxt(path)


CONV = np.array(load('FM/ligand')['cell_A'])
FRAC0 = np.array([a['position_A'] for a in load('FM/ligand')['atoms'] if a['name'].startswith('Cu')]) @ np.linalg.inv(CONV)


def supercell(cell):
    return np.round(np.array(cell) @ np.linalg.inv(CONV)).astype(int)


def cu_info(d):
    """Cu names, positions, sign of the moment along z, conventional-fractional positions and plaquette type (0/1)."""
    cell = np.array(d['cell_A'])
    Mi = supercell(cell)
    inv = np.linalg.inv(cell)
    cu = [a for a in d['atoms'] if a['name'].startswith('Cu')]
    ox = np.array([a['position_A'] for a in d['atoms'] if a['name'].startswith('O')])
    pos = np.array([a['position_A'] for a in cu])
    mz = np.array([a.get('M(z)', a.get('w_magmom')) for a in cu])
    fcu, fo = pos @ inv @ Mi, ox @ inv @ Mi
    types = []
    for f in fcu:
        v = np.unique(np.round((fo - f) - np.round(fo - f), 2), axis=0)
        near = v[np.argsort(np.linalg.norm(v @ CONV, axis=1))[:4]]
        key = [np.sign(n[1]) for n in near if n[0] < 0 and n[2] > 0]
        types.append(0 if key and key[0] > 0 else 1)
    return {'names': [a['name'] for a in cu], 'pos': pos, 'spins': np.sign(mz), 'fconv': fcu,
            'types': np.array(types), 'Mi': Mi, 'cell': cell}


class Native:
    """Isotropic model exactly as archived: the calculation's own (super)cell and pair list."""

    def __init__(self, key, rmax=None, data=None):
        d = data or load(key)
        info = cu_info(d)
        self.info, self.n, self.cell = info, len(info['names']), info['cell']
        idx = {n: k for k, n in enumerate(info['names'])}
        P = [p for p in d['pairs'] if rmax is None or p[9] < rmax]
        self.I = np.array([idx[p[0]] for p in P])
        self.J = np.array([idx[p[1]] for p in P])
        self.R = np.array([p[2:5] for p in P], float)
        self.V = np.array([p[5] for p in P])
        self.spins = info['spins']

    def K(self, k):
        """K_ij(k) = sum_R J_ij(R) exp(2 pi i k.R), k in this cell's reciprocal fractional coordinates."""
        k = np.atleast_2d(np.asarray(k, float))
        ph = np.exp(2j*np.pi*(k @ self.R.T))*self.V
        out = np.zeros((len(k), self.n, self.n), complex)
        for a in range(self.n):
            for b in range(self.n):
                m = (self.I == a) & (self.J == b)
                if m.any():
                    out[:, a, b] = ph[:, m].sum(1)
        return out

    def h(self, k, spins=None):
        """Authors' bosonic matrix h(k) in meV (S = 1) for a collinear state with spins +-1 along z."""
        return bosonic(self.K(k), self.K([0, 0, 0])[0].real, self.spins if spins is None else spins)

    def energy(self, spins=None):
        s = self.spins if spins is None else np.asarray(spins, float)
        return float(-s @ self.K([0, 0, 0])[0].real @ s/self.n)


def bosonic(Kk, K0, spins):
    s = np.asarray(spins, float)
    C = np.diag(-(K0*np.outer(s, s)).sum(1))
    out = []
    for Kq in Kk:
        A = -0.5*Kq*(1 + np.outer(s, s))
        B = -0.5*Kq*(1 - np.outer(s, s))
        out.append(2*np.block([[A - C, B], [B.conj().T, A - C]]))
    return np.array(out)


def omegas(H):
    """Positive eigenvalues of g.h for each k (paraunitary problem); also the largest imaginary part."""
    n = H.shape[1]//2
    g = np.diag(np.r_[np.ones(n), -np.ones(n)])
    ev = np.linalg.eigvals(g @ H)
    return np.sort(ev.real, axis=1)[:, n:], float(np.abs(ev.imag).max())


def grid(n, shift=0.5):
    g = (np.arange(n) + shift)/n
    return np.array(np.meshgrid(g, g, g, indexing='ij')).reshape(3, -1).T


def magnon_stats(model, n=8, spin=0.5, spins=None):
    """Mean over the BZ of all branches and band maximum of the LSWT energies (meV) for spin length `spin`."""
    ks = grid(n)
    H = model.h(ks, spins=spins)
    w, imag = omegas(H)
    hev = np.linalg.eigvalsh(H)
    i = np.unravel_index(np.argmax(w), w.shape)[0]
    r = minimize(lambda k: -omegas(model.h(k, spins=spins))[0].max(), ks[i], method='Nelder-Mead',
                 options={'xatol': 1e-6, 'fatol': 1e-8})
    return {'mean': float(w.mean()/spin), 'max': float(max(w.max(), -r.fun)/spin), 'min_h_eig': float(hev.min()),
            'max_imag': imag, 'h_mean': float(hev[hev > 0].mean()/spin), 'h_max': float(hev.max()/spin),
            'max_grid_only': float(w.max()/spin)}


def min_h_eig(model, n=8):
    return float(np.linalg.eigvalsh(model.h(grid(n, 0.0)))[:, 0].min())


# ---------- symmetric (folded) two-sublattice model on the C-centred lattice ----------

def folded(key, rmax=None, data=None):
    """Average the pair list over Cu of the same plaquette type and equal ideal displacement (conventional r.l.u.)."""
    d = data or load(key)
    info = cu_info(d)
    idx = {n: k for k, n in enumerate(info['names'])}
    inv = np.linalg.inv(info['cell'])
    acc = defaultdict(list)
    for p in d['pairs']:
        if rmax is not None and p[9] >= rmax:
            continue
        i, j = idx[p[0]], idx[p[1]]
        df = np.array(p[6:9]) @ inv @ info['Mi']
        acc[(int(info['types'][i]), int(info['types'][j]), tuple(np.round(df*4).astype(int)))].append(p[5])
    nsite = [int((info['types'] == 0).sum()), int((info['types'] == 1).sum())]
    return {k: float(np.sum(v))/nsite[k[0]] for k, v in acc.items()}


class Folded:
    def __init__(self, model):
        items = list(model.items())
        self.A = np.array([k[0] for k, v in items])
        self.B = np.array([k[1] for k, v in items])
        self.D = np.array([k[2] for k, v in items])/4.0
        self.V = np.array([v for k, v in items])

    def K(self, q):
        q = np.atleast_2d(np.asarray(q, float))
        out = np.zeros((len(q), 2, 2), complex)
        for a in (0, 1):
            for b in (0, 1):
                m = (self.A == a) & (self.B == b)
                out[:, a, b] = np.exp(2j*np.pi*(q @ self.D[m].T)) @ self.V[m]
        return out

    def top(self, q):
        return float(np.linalg.eigvalsh(self.K(q))[0, -1])

    def lt(self, starts=((0.5, 0, 0.5), (0.49, 0, 0.53), (0, 0.5, 0), (0.5, 0.5, 0.5), (0.5, 0, 0), (0.1, 0.1, 0.1))):
        """Luttinger-Tisza: maximise the top eigenvalue of K(q) (equivalently minimise the FM-reference h(k))."""
        g = np.arange(16)/16
        qs = np.array(np.meshgrid(g, g, g, indexing='ij')).reshape(3, -1).T
        lam = np.linalg.eigvalsh(self.K(qs))[:, -1]
        starts = list(starts) + [tuple(qs[np.argmax(lam)])]
        best = None
        for s in starts:
            r = minimize(lambda q: -self.top(q), s, method='Nelder-Mead',
                         options={'xatol': 1e-8, 'fatol': 1e-12, 'maxiter': 8000})
            if best is None or -r.fun > best[1]:
                best = (r.x, -r.fun)
        return np.asarray(best[0], float), best[1]   # raw maximiser (the 2-site spectrum is periodic only on the C-centred reciprocal lattice)


class FoldedOnSupercell:
    """A folded model evaluated on the 16-site AF1 supercell of the AFM/ligand calculation (k in conventional r.l.u.)."""

    def __init__(self, model, sites):
        self.pos, self.types, self.spins, self.Mi = sites
        self.n = len(self.pos)
        Minv = np.linalg.inv(self.Mi)
        I, J, D, V = [], [], [], []
        for i in range(self.n):
            for (a, b, k), val in model.items():
                if a != self.types[i]:
                    continue
                d = np.array(k)/4.0
                diff = (self.pos[i] + d - self.pos) @ Minv
                ok = np.where((np.abs(diff - np.round(diff)) < 1e-6).all(1))[0]
                assert len(ok) == 1 and self.types[ok[0]] == b
                I.append(i); J.append(ok[0]); D.append(d); V.append(val)
        self.I, self.J, self.D, self.V = map(np.array, (I, J, D, V))

    def K(self, k):
        k = np.atleast_2d(np.asarray(k, float))
        ph = np.exp(2j*np.pi*(k @ self.D.T))*self.V
        out = np.zeros((len(k), self.n, self.n), complex)
        for a in range(self.n):
            for b in range(self.n):
                m = (self.I == a) & (self.J == b)
                if m.any():
                    out[:, a, b] = ph[:, m].sum(1)
        return out

    def h(self, k_sc, spins=None):
        k = np.atleast_2d(k_sc) @ np.linalg.inv(self.Mi).T  # supercell fractional -> conventional r.l.u.
        return bosonic(self.K(k), self.K([0, 0, 0])[0].real, self.spins if spins is None else spins)


def af1_sites():
    info = cu_info(load('AFM/ligand'))
    return np.round(info['fconv']*4)/4, info['types'], info['spins'], info['Mi']


def canonical(q):
    """Representative of +-q modulo integer vectors, with components in (-0.5, 0.5] after sign choice."""
    q = np.asarray(q, float)
    q = q - np.round(q)
    return q


def q_distance(q, ref):
    """Smallest Euclidean distance (r.l.u.) between q and any of +-ref + G, G integer."""
    q, ref = np.asarray(q, float), np.asarray(ref, float)
    best = np.inf
    for sgn in (1, -1):
        d = q - sgn*ref
        best = min(best, float(np.linalg.norm(d - np.round(d))))
    return best


def chain_J(key, data=None):
    """Mean TB2J J of the dominant (most negative) bond class, with its distance and spread."""
    d = data or load(key)
    P = np.array([[p[5], p[9]] for p in d['pairs']])
    jmin = P[:, 0].min()
    dist = P[np.argmin(P[:, 0]), 1]
    sel = P[np.abs(P[:, 1] - dist) < 0.05]
    sel = sel[sel[:, 0] < 0.5*jmin]
    return float(sel[:, 0].mean()), float(dist), float(sel[:, 0].min()), float(sel[:, 0].max()), len(sel)
