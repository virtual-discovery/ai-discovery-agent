"""Teacher spin model for amaranta-2d-triangular-helimagnets (numpy only).

The model is the AMaRaNTA Hamiltonian of Putatunda/Orlando et al. with the authors' archived
parameters (pool.json, from the AiiDA provenance):

    E = sum_{i<j, 1NN} e_i . J1(ij) . e_j + J2 sum_{i<j, 2NN} e_i.e_j + J3 sum_{i<j, 3NN} e_i.e_j
        + K sum_i (e_i^z)^2

with unit vectors e_i, each pair counted once, and the 1NN/2NN/3NN classes being the exact
bond lengths of the three pairs ('couples') AMaRaNTA used (it groups distances within 1 %).
J1 is the archived tensor on the archived 1NN pair; the other 1NN bonds carry the tensor
transformed by the site point group of the magnetic ion (determined from the full structure,
ligands included): J1(P d) = Q J1(d) Q^T with Q = det(P) P for spins (axial vectors), and
J1(-d) = J1(d)^T. K = 2 A_archived: AMaRaNTA's SIA formula (E+z - E+x - E-x + E-z)/4 returns
half the coefficient of (e^z)^2, and K = 2A reproduces the archived SIA energies (validate()).
"""
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
MAGNETIC = ('Ni', 'V', 'Mn')
COMPOUNDS = ['NiCl2', 'NiBr2', 'VCl2', 'VBr2', 'VI2', 'MnCl2', 'MnBr2', 'MnI2']


def load_pool():
    return json.loads((HERE/'pool.json').read_text())


def point_group(cell, sites):
    """Site point group (3x3 polar matrices) of the first magnetic site, ligands included."""
    cell = np.array(cell, float)
    inv = np.linalg.inv(cell)
    origin = np.array(next(s['position'] for s in sites if s['species'] in MAGNETIC), float)
    X = np.array([s['position'] for s in sites], float) - origin
    kinds = [s['species'] for s in sites]
    ops = []
    for n in range(6):
        th = n * np.pi / 3
        Rz = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
        for my in (1, -1):
            for mz in (1, -1):
                P = Rz @ np.diag([1, my, mz])
                ok = True
                for y, k in zip(X @ P.T, kinds):
                    hit = False
                    for x, k2 in zip(X, kinds):
                        if k2 == k:
                            f = (y - x) @ inv
                            if np.linalg.norm((f - np.round(f)) @ cell) < 0.02:
                                hit = True
                                break
                    if not hit:
                        ok = False
                        break
                if ok:
                    ops.append(P)
    return ops


class Compound:
    def __init__(self, name, rec):
        self.name = name
        st = rec['structure']
        self.cell = np.array(st['cell'], float)
        self.basis = np.array([s['position'] for s in st['sites'] if s['species'] in MAGNETIC], float)
        self.species = next(s['species'] for s in st['sites'] if s['species'] in MAGNETIC)
        J = rec['J_dict']
        self.T = np.array([[J['J_xx'], J['J_xy'], J['J_xz']], [J['J_yx'], J['J_yy'], J['J_yz']], [J['J_zx'], J['J_zy'], J['J_zz']]])
        self.J2, self.J3, self.A = J['J_scalar_2NN'], J['J_scalar_3NN'], J['A_scalar']
        self.K = 2.0 * self.A
        c = rec['couples']
        self.d = [np.subtract(c[k][1], c[k][0]) for k in ('1NN couple', '2NN couple', '3NN couple')]
        self.L = [float(np.linalg.norm(x)) for x in self.d]
        self.ops = point_group(st['cell'], st['sites'])
        self._bond_cache = {}

    # ---- bonds -----------------------------------------------------------------
    def tensor(self, d):
        for P in self.ops:
            v = P @ self.d[0]
            Q = np.linalg.det(P) * P
            if np.allclose(v[:2], d[:2], atol=0.02):
                return Q @ self.T @ Q.T
            if np.allclose(v[:2], -d[:2], atol=0.02):
                return (Q @ self.T @ Q.T).T
        raise ValueError('bond not in the 1NN orbit')

    def coupling(self, d, iso=False):
        """3x3 coupling for a bond vector d (Cartesian, A), or None."""
        L = np.linalg.norm(d)
        if abs(L - self.L[0]) < 0.01:
            return np.eye(3) * np.trace(self.T) / 3 if iso else self.tensor(d)
        if abs(L - self.L[1]) < 0.01:
            return np.eye(3) * self.J2
        if abs(L - self.L[2]) < 0.01:
            return np.eye(3) * self.J3
        return None

    def lattice_bonds(self, rmax=8.5, iso=False):
        """All bonds from basis site 0 to (basis k, lattice R): list of (k, n1, n2, d, Jmat)."""
        key = (rmax, iso)
        if key in self._bond_cache:
            return self._bond_cache[key]
        out = []
        for k, b in enumerate(self.basis):
            for n1 in range(-4, 5):
                for n2 in range(-4, 5):
                    d = b + n1 * self.cell[0] + n2 * self.cell[1] - self.basis[0]
                    if 0.1 < np.linalg.norm(d) < rmax:
                        Jm = self.coupling(d, iso)
                        if Jm is not None:
                            out.append((k, n1, n2, d, Jm))
        self._bond_cache[key] = out
        return out

    # ---- supercells ------------------------------------------------------------
    def sites(self, n1, n2):
        """Lab ordering: basis atom k, then i < n1, then j < n2; position b_k + i a1 + j a2."""
        return np.array([b + i * self.cell[0] + j * self.cell[1] for b in self.basis for i in range(n1) for j in range(n2)])

    def supercell_couplings(self, n1, n2, iso=False):
        """Sparse list (i, j, J) over ordered pairs i != j (images summed), for E = 1/2 sum e_i J e_j."""
        A1, A2 = n1 * self.cell[0], n2 * self.cell[1]
        acc = {}
        pos = self.sites(n1, n2)
        N = len(pos)
        shifts = [a * A1 + b * A2 for a in range(-3, 4) for b in range(-3, 4)]
        for i in range(N):
            for j in range(N):
                for s in shifts:
                    d = pos[j] + s - pos[i]
                    L = np.linalg.norm(d)
                    if L < 0.1 or L > 8.5:
                        continue
                    Jm = self.coupling(d, iso)
                    if Jm is not None:
                        acc[(i, j)] = acc.get((i, j), 0) + Jm
        return [(i, j, Jm) for (i, j), Jm in acc.items()]

    def energy_meV(self, spins, couplings, iso=False):
        e = np.asarray(spins, float)
        e = e / np.linalg.norm(e, axis=1, keepdims=True)
        E = 0.0
        for i, j, Jm in couplings:
            E += 0.5 * e[i] @ Jm @ e[j]
        if not iso:
            E += self.K * np.sum(e[:, 2] ** 2)
        return E

    # ---- reciprocal space ------------------------------------------------------
    def bvec(self):
        return 2 * np.pi * np.linalg.inv(self.cell[:2, :2]).T

    def Jq_iso(self, q):
        """Energy per site of a flat single-q spiral of the isotropic model (q Cartesian, 1/A)."""
        b = self.lattice_bonds(iso=True)
        D = np.array([x[3][:2] for x in b])
        Jv = np.array([x[4][0, 0] for x in b])
        return 0.5 * float(Jv @ np.cos(D @ q))


def all_compounds():
    pool = load_pool()
    return {n: Compound(n, pool[n]) for n in COMPOUNDS}


def validate(pool=None):
    """Within-set residuals (meV) of the teacher model against the archived DFT energies."""
    pool = pool or load_pool()
    out = {}
    for n in COMPOUNDS:
        c = Compound(n, pool[n])
        soc, col = [], []
        for fs in pool[n]['four_state_sets']:
            st = fs['supercell']
            cell = np.array(st['cell'])
            pos = np.array([s['position'] for s in st['sites'] if s['species'] in MAGNETIC])
            nmag = len(pos)
            shifts = [a * cell[0] + b * cell[1] for a in range(-2, 3) for b in range(-2, 3)]
            cp = {}
            for i in range(nmag):
                for j in range(nmag):
                    for s in shifts:
                        d = pos[j] + s - pos[i]
                        if np.linalg.norm(d) > 0.1:
                            Jm = c.coupling(d)
                            if Jm is not None:
                                cp[(i, j)] = cp.get((i, j), 0) + Jm
            cpl = [(i, j, J) for (i, j), J in cp.items()]
            Em, Ed = [], []
            for conf, E in zip(fs['moment_directions'], fs['energies_eV']):
                conf = np.array(conf, float)
                m = conf.reshape(-1, 3)[:nmag] if fs['soc'] else np.c_[np.zeros((nmag, 2)), conf[:nmag]]
                Em.append(sum(0.5 * m[i] @ J @ m[j] for i, j, J in cpl) + (c.K * np.sum(m[:, 2] ** 2) if fs['soc'] else 0))
                Ed.append(E * 1000)
            r = (np.array(Ed) - Ed[0]) - (np.array(Em) - Em[0])
            (soc if fs['soc'] else col).append(r[1:])
        soc = np.concatenate(soc)
        out[n] = {'soc_rms_meV': float(np.sqrt(np.mean(soc ** 2))), 'soc_max_meV': float(np.abs(soc).max()),
                  'collinear_single_flip_resid_meV': [float(x[0]) for x in col]}
    return out


if __name__ == '__main__':
    print(json.dumps(validate(), indent=1))
