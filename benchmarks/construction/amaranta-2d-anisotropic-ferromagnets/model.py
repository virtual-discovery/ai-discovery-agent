"""Teacher spin model for amaranta-2d-anisotropic-ferromagnets (numpy only).

The model is the AMaRaNTA Hamiltonian (Putatunda, Cuono, Orlando, Picozzi) with the authors'
archived parameters (pool.json, from the AiiDA provenance):

    E = sum_{i<j, 1NN} e_i . J1(ij) . e_j + J2 sum_{i<j, 2NN} e_i.e_j + J3 sum_{i<j, 3NN} e_i.e_j
        + K sum_i (e_i^z)^2

with unit vectors e_i and each pair counted once. The 1NN/2NN/3NN classes are the bonds of
exactly the lengths of the three pairs ("couples") AMaRaNTA used. J1 is the archived 3x3 tensor
on the archived 1NN couple (first atom -> second atom); every other bond of that class carries the
tensor transformed by the layer's space-group operation g = (R, t) that maps the archived couple
onto it: J(g a -> g b) = Q J Q^T with Q = det(R) R (spins are axial vectors), and a reversed bond
carries the transpose. The operations are found from the full structure (ligands included), so
sublattice-staggered antisymmetric parts (for example the DMI on the two Cr chains of CrSBr) are
kept. K = 2 A_archived: AMaRaNTA's single-ion formula (E+z - E+x - E-x + E-z)/4 returns half the
coefficient of (e^z)^2 (validate() checks both choices against the archived DFT energies).
"""
import itertools
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
MAGNETIC = ('Cr',)
COMPOUNDS = ['CrI3', 'CrBr3', 'CrCl3', 'CrGeTe3', 'CrSiTe3', 'CrSiSe3', 'CrSBr', 'CrPS4', 'CrClO', 'CrBrO']
SYM_TOL = 0.03          # A, tolerance for symmetry detection
PSEUDO_TOL = 0.15       # A, tolerance for pseudo-symmetry (only used where no exact operation exists: CrCl3)
CLASS_REL_TOL = 0.01    # AMaRaNTA groups neighbour distances within 1 %


def load_pool():
    return json.loads((HERE/'pool.json').read_text())


def unwrap_layer(cell, sites):
    """Shift atoms by multiples of c so the layer is contiguous in z (CrSiSe3 straddles the cell edge)."""
    c = np.array(cell, float)[2]
    pos = np.array([s['position'] for s in sites], float)
    z = pos[:, 2]
    # put every atom within c/2 of the first magnetic atom
    z0 = next(p[2] for p, s in zip(pos, sites) if s['species'] in MAGNETIC)
    shift = np.round((z - z0) / c[2])
    pos = pos - shift[:, None] * c[None, :]
    return pos


def lattice_automorphisms(A2, tol=2e-3):
    """In-plane orthogonal 2x2 matrices that map the 2D lattice (rows of A2) onto itself."""
    out = []
    for m in itertools.product(range(-2, 3), repeat=4):
        M = np.array(m, float).reshape(2, 2)
        if abs(abs(np.linalg.det(M)) - 1) > 1e-9:
            continue
        # new basis vectors b_i = sum_j M_ij a_j; R maps a_i -> b_i
        B = M @ A2
        R = np.linalg.solve(A2, B).T   # R a_i^T = b_i^T  ->  A2 R^T = B
        if np.allclose(R @ R.T, np.eye(2), atol=tol):
            u, _, vt = np.linalg.svd(R)
            R = u @ vt                      # nearest orthogonal matrix
            if not any(np.allclose(R, X, atol=1e-6) for X in out):
                out.append(R)
    return out


class Structure:
    def __init__(self, cell, sites):
        self.cell = np.array(cell, float)
        self.species = [s['species'] for s in sites]
        self.pos = unwrap_layer(cell, sites)
        self.A2 = self.cell[:2, :2]
        self.ops = self.space_group(SYM_TOL, 2e-3)
        self.pseudo_ops = None

    def get_pseudo_ops(self):
        if self.pseudo_ops is None:
            self.pseudo_ops = self.space_group(PSEUDO_TOL, 3e-2)
        return self.pseudo_ops

    def frac2(self, xy):
        return np.linalg.solve(self.A2.T, np.asarray(xy, float).T).T

    def same_site(self, x, k, tol=SYM_TOL):
        """Index of the atom equivalent to Cartesian position x (in-plane lattice translations), species k, or None."""
        for i, (p, s) in enumerate(zip(self.pos, self.species)):
            if s != k or abs(p[2] - x[2]) > tol:
                continue
            f = self.frac2((x - p)[:2])
            if np.linalg.norm((f - np.round(f)) @ self.A2) < tol:
                return i
        return None

    def space_group(self, tol, otol):
        ops = []
        for R2 in lattice_automorphisms(self.A2, otol):
            for sz in (1.0, -1.0):
                R = np.eye(3)
                R[:2, :2] = R2
                R[2, 2] = sz
                # map atom 0 onto every atom of the same species to propose translations
                x0 = self.pos[0]
                for p, s in zip(self.pos, self.species):
                    if s != self.species[0]:
                        continue
                    t = p - R @ x0
                    ok = all(self.same_site(R @ x + t, k, tol) is not None for x, k in zip(self.pos, self.species))
                    if ok:
                        f = self.frac2(t[:2])
                        t = t.copy()
                        t[:2] = (f - np.floor(f + 1e-9)) @ self.A2
                        if not any(np.allclose(R, R1) and np.allclose(t, t1, atol=tol) for R1, t1 in ops):
                            ops.append((R, t))
        return ops


class Compound:
    def __init__(self, name, rec, K_factor=2.0):
        self.name = name
        st = rec['structure']
        self.st = Structure(st['cell'], st['sites'])
        mag = [i for i, s in enumerate(self.st.species) if s in MAGNETIC]
        self.basis = self.st.pos[mag]
        self.cell = self.st.cell
        J = rec['J_dict']
        self.T = np.array([[J['J_xx'], J['J_xy'], J['J_xz']], [J['J_yx'], J['J_yy'], J['J_yz']], [J['J_zx'], J['J_zy'], J['J_zz']]])
        self.J2, self.J3, self.A = J['J_scalar_2NN'], J['J_scalar_3NN'], J['A_scalar']
        self.K = K_factor * self.A
        c = rec['couples']
        # couples in the unwrapped frame (fix z wrap by the cell's c)
        cz = self.cell[2, 2]
        def unwrap(p):
            p = np.array(p, float)
            p[2] -= np.round((p[2] - self.basis[0][2]) / cz) * cz
            return p
        self.couples = [(unwrap(c[k][0]), unwrap(c[k][1])) for k in ('1NN couple', '2NN couple', '3NN couple')]
        self.d = [b - a for a, b in self.couples]
        self.L = [float(np.linalg.norm(x)) for x in self.d]
        self._cache = {}

    # ---- bonds -----------------------------------------------------------------
    def site_index(self, x):
        """(basis index k, lattice translation (n1, n2)) of a magnetic position x."""
        for k, b in enumerate(self.basis):
            if abs(b[2] - x[2]) > SYM_TOL:
                continue
            f = self.st.frac2((x - b)[:2])
            if np.linalg.norm((f - np.round(f)) @ self.st.A2) < SYM_TOL:
                return k, tuple(int(v) for v in np.round(f))
        raise ValueError('not a magnetic site')

    def _equiv(self, x, y, tol=SYM_TOL):
        if abs(x[2] - y[2]) > tol:
            return False
        f = self.st.frac2((x - y)[:2])
        return np.linalg.norm((f - np.round(f)) @ self.st.A2) < tol

    def _images(self, ops, tol, xi, d):
        """Transformed archived tensors Q T Q^T (or transposes) for every op mapping the couple onto bond (xi, xi + d)."""
        a, b = self.couples[0]
        out = []
        for R, t in ops:
            Q = np.linalg.det(R) * R
            ga, gb = R @ a + t, R @ b + t
            gd = gb - ga
            if np.allclose(gd, d, atol=tol) and self._equiv(ga, xi, tol):
                out.append(('f', Q))
            if np.allclose(-gd, d, atol=tol) and self._equiv(gb, xi, tol):
                out.append(('r', Q))
        return out

    def symmetry_ops(self):
        """Exact space-group operations if they connect every 1NN-class bond to the couple, else the
        pseudo-symmetric set (CrCl3: monoclinic distortion of 0.1 %, pseudo-D3d)."""
        if 'ops' in self._cache:
            return self._cache['ops']
        chosen = (self.st.ops, SYM_TOL, 'exact')
        for k0, b0 in enumerate(self.basis):
            for k, b in enumerate(self.basis):
                for m1 in range(-3, 4):
                    for m2 in range(-3, 4):
                        d = b + m1 * self.cell[0] + m2 * self.cell[1] - b0
                        if abs(np.linalg.norm(d) / self.L[0] - 1) < CLASS_REL_TOL and not self._images(self.st.ops, SYM_TOL, b0, d):
                            chosen = (self.st.get_pseudo_ops(), PSEUDO_TOL, 'pseudo')
        self._cache['ops'] = chosen
        return chosen

    def symmetrized_T(self):
        """Archived tensor averaged over the stabilizer of the archived bond (removes the
        symmetry-forbidden part, 0.003 meV in CrI3, up to 0.08 meV in CrGeTe3)."""
        if 'Tsym' in self._cache:
            return self._cache['Tsym']
        ops, tol, _ = self.symmetry_ops()
        a, b = self.couples[0]
        acc = []
        for kind, Q in self._images(ops, tol, a, b - a):
            X = Q @ self.T @ Q.T
            acc.append(X if kind == 'f' else X.T)
        Ts = np.mean(acc, axis=0)
        self._cache['Tsym'] = Ts
        return Ts

    def tensor(self, xi, d):
        """1NN coupling matrix for the bond from magnetic position xi to xi + d (J in e_i . J . e_j)."""
        ops, tol, _ = self.symmetry_ops()
        im = self._images(ops, tol, xi, d)
        if not im:
            raise ValueError('bond %s from %s not in the 1NN orbit of %s' % (d, xi, self.name))
        kind, Q = im[0]
        X = Q @ self.symmetrized_T() @ Q.T
        return X if kind == 'f' else X.T

    def coupling(self, xi, d, iso=False):
        L = np.linalg.norm(d)
        if abs(L / self.L[0] - 1) < CLASS_REL_TOL:
            return np.eye(3) * np.trace(self.T) / 3 if iso else self.tensor(xi, d)
        if abs(L / self.L[1] - 1) < CLASS_REL_TOL:
            return np.eye(3) * self.J2
        if abs(L / self.L[2] - 1) < CLASS_REL_TOL:
            return np.eye(3) * self.J3
        return None

    def bond_table(self, iso=False, rmax=None):
        """[(k0, k, m1, m2, J 3x3)]: coupling from basis site k0 (cell 0) to basis site k in cell (m1, m2)."""
        key = ('tab', iso)
        if key in self._cache:
            return self._cache[key]
        rmax = rmax or max(self.L) * 1.05
        out = []
        for k0, b0 in enumerate(self.basis):
            for k, b in enumerate(self.basis):
                for m1 in range(-6, 7):
                    for m2 in range(-6, 7):
                        d = b + m1 * self.cell[0] + m2 * self.cell[1] - b0
                        if 0.1 < np.linalg.norm(d) < rmax:
                            Jm = self.coupling(b0, d, iso)
                            if Jm is not None:
                                out.append((k0, k, m1, m2, np.array(Jm, float)))
        self._cache[key] = out
        return out


def supercell_energy(table, K, nbasis, n1, n2, spins):
    """Energy (meV) of a diagonal n1 x n2 supercell (site order k, i, j); periodic images summed exactly."""
    e = np.asarray(spins, float)
    e = e / np.linalg.norm(e, axis=1, keepdims=True)
    S = e.reshape(nbasis, n1, n2, 3)
    E = 0.0
    for k0, k, m1, m2, J in table:
        Sj = np.roll(np.roll(S[k], -m1, axis=0), -m2, axis=1)
        E += 0.5 * np.einsum('ija,ab,ijb->', S[k0], J, Sj)
    return E + K * np.sum(e[:, 2] ** 2)


def all_compounds(K_factor=2.0):
    pool = load_pool()
    return {n: Compound(n, pool[n], K_factor) for n in COMPOUNDS}


def archived_set_energies(c, fs):
    """Model energies (meV) of the archived four-state configurations of one set (general cell)."""
    st = fs['supercell']
    cell = np.array(st['cell'], float)
    sites = st['sites']
    pos = unwrap_layer(cell, sites)
    mag = [i for i, s in enumerate(sites) if s['species'] in MAGNETIC]
    pos = pos[mag]
    # align with the compound frame (same Cartesian frame; z unwrap relative to the compound's basis)
    cz = cell[2, 2]
    pos[:, 2] -= np.round((pos[:, 2] - c.basis[0][2]) / cz) * cz
    n = len(pos)
    shifts = [a * cell[0] + b * cell[1] for a in range(-3, 4) for b in range(-3, 4)]
    cp = {}
    for i in range(n):
        for j in range(n):
            for s in shifts:
                d = pos[j] + s - pos[i]
                if np.linalg.norm(d) > 0.1:
                    Jm = c.coupling(pos[i], d)
                    if Jm is not None:
                        cp[(i, j)] = cp.get((i, j), 0) + Jm
    out = []
    for conf in fs['moment_directions']:
        conf = np.array(conf, float)
        if len(conf) == 3 * len(sites):
            m = conf.reshape(-1, 3)[mag]
            soc = True
        else:
            m = np.c_[np.zeros((n, 2)), conf[mag]]
            soc = False
        E = sum(0.5 * m[i] @ J @ m[j] for (i, j), J in cp.items())
        if soc:
            E += c.K * np.sum(m[:, 2] ** 2)
        out.append((E, soc))
    return out


def validate(pool=None, K_factor=2.0):
    """Within-set residuals (meV) of the teacher model against the archived DFT energies."""
    pool = pool or load_pool()
    res = {}
    for n in COMPOUNDS:
        c = Compound(n, pool[n], K_factor)
        soc, col = [], []
        for fs in pool[n]['four_state_sets']:
            em = archived_set_energies(c, fs)
            Em = np.array([e for e, s in em])
            Ed = np.array(fs['energies_eV']) * 1000
            r = (Ed - Ed[0]) - (Em - Em[0])
            (soc if em[0][1] else col).append(r[1:])
        soc = np.concatenate(soc)
        res[n] = {'soc_rms_meV': float(np.sqrt(np.mean(soc ** 2))), 'soc_max_meV': float(np.abs(soc).max()),
                  'collinear_resid_meV': [float(x[0]) for x in col]}
    return res


if __name__ == '__main__':
    pool = load_pool()
    for n in COMPOUNDS:
        c = Compound(n, pool[n])
        print(n, 'ops', len(c.st.ops), 'pseudo', len(c.st.get_pseudo_ops()), 'L', np.round(c.L, 4), 'nbonds', len(c.bond_table()))
    print(json.dumps({'K=2A': validate(pool), 'K=A': validate(pool, 1.0)}, indent=1))
