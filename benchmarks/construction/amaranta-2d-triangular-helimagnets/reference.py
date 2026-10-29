"""Reference analysis and wrong-path variants for amaranta-2d-triangular-helimagnets.

Builds the lab's bond tables from model.py, and computes, for the teacher model and for the
realistic wrong paths: Luttinger-Tisza wave vectors, ground-state classes, energies of the
hidden prediction configurations, and the VI2 easy axis.
"""
import numpy as np
from scipy.optimize import minimize
import model as M

TYPES = ('ferromagnet', '120-degree', 'collinear stripe', 'incommensurate spiral')


def bond_table(c, iso=False, rmax=8.5):
    """[(k0, k, m1, m2, J 3x3)]: coupling from basis site k0 (cell 0) to basis site k in cell (m1, m2)."""
    out = []
    for k0, b0 in enumerate(c.basis):
        for k, b in enumerate(c.basis):
            for m1 in range(-4, 5):
                for m2 in range(-4, 5):
                    d = b + m1 * c.cell[0] + m2 * c.cell[1] - b0
                    if 0.1 < np.linalg.norm(d) < rmax:
                        Jm = c.coupling(d, iso)
                        if Jm is not None:
                            out.append((k0, k, m1, m2, np.array(Jm, float)))
    return out


def energy(table, K, nbasis, n1, n2, spins):
    """Energy (meV) of a diagonal n1 x n2 supercell (lab site order k, i, j); images summed exactly."""
    e = np.asarray(spins, float)
    e = e / np.linalg.norm(e, axis=1, keepdims=True)
    S = e.reshape(nbasis, n1, n2, 3)
    E = 0.0
    for k0, k, m1, m2, J in table:
        Sj = np.roll(np.roll(S[k], -m1, axis=0), -m2, axis=1)
        E += 0.5 * np.einsum('ija,ab,ijb->', S[k0], J, Sj)
    return E + K * np.sum(e[:, 2] ** 2)


def iso_bonds_from_classes(c, classes):
    """Isotropic model from a list of (bond vector representative set, J): every lattice bond whose
    vector matches one in a class gets that J. classes: list of (list of 2D vectors, J)."""
    out = []
    for k0, b0 in enumerate(c.basis):
        for k, b in enumerate(c.basis):
            for m1 in range(-4, 5):
                for m2 in range(-4, 5):
                    d = (b + m1 * c.cell[0] + m2 * c.cell[1] - b0)[:2]
                    for vecs, J in classes:
                        if any(np.linalg.norm(d - v) < 0.02 for v in vecs):
                            out.append((k0, k, m1, m2, np.eye(3) * J))
    return out


def Jq(table, c, q):
    """Energy per site of a flat single-q spiral for an isotropic table (q Cartesian)."""
    e = 0.0
    for k0, k, m1, m2, J in table:
        if k0 == 0:
            d = (c.basis[k] + m1 * c.cell[0] + m2 * c.cell[1] - c.basis[0])[:2]
            e += 0.5 * J[0, 0] * np.cos(q @ d)
    return e


def lt(table, c, grid=121):
    B = c.bvec()
    best = None
    for h in np.linspace(0, 1, grid):
        for k in np.linspace(0, 1, grid):
            e = Jq(table, c, np.array([h, k]) @ B)
            if best is None or e < best[0] - 1e-12:
                best = (e, h, k)
    r = minimize(lambda x: Jq(table, c, np.array(x) @ B), best[1:], method='Nelder-Mead',
                 options={'xatol': 1e-10, 'fatol': 1e-13})
    return np.array(r.x), float(r.fun)


def primitive_G(c):
    """Integer (h, k) translations of the primitive reciprocal lattice in catalog r.l.u."""
    frac = [np.linalg.solve(c.cell[:2, :2].T, (b - c.basis[0])[:2]) for b in c.basis]
    G = []
    for h in range(-3, 4):
        for k in range(-3, 4):
            g = np.array([h, k], float)
            if all(np.isclose(np.cos(2 * np.pi * (g @ f)), 1.0) for f in frac):
                G.append(g)
    return np.array(G)


def equivalents(c, q):
    """All r.l.u. vectors equivalent to q under the site point group, sign and primitive G."""
    B = c.bvec()
    qc = np.asarray(q, float) @ B
    out = []
    for P in c.ops:
        for s in (1, -1):
            v = s * (P[:2, :2] @ qc)
            h = np.linalg.solve(B.T, v)
            for g in primitive_G(c):
                out.append(h + g)
    return np.array(out)


def q_distance(c, q1, q2):
    """Smallest Cartesian distance (1/A) between q1 and any equivalent of q2 (both r.l.u.)."""
    B = c.bvec()
    E = equivalents(c, q2)
    return float(np.min(np.linalg.norm((E - np.asarray(q1, float)) @ B, axis=1)))


def classify(c, q):
    """Ground-state class of a single-q state with wave vector q (r.l.u.)."""
    B = c.bvec()
    if q_distance(c, q, [0, 0]) < 1e-3:
        return 'ferromagnet'
    G = primitive_G(c)
    two_q = 2 * np.asarray(q, float)
    if np.min(np.linalg.norm((G - two_q) @ B, axis=1)) < 1e-3:
        return 'collinear stripe'
    three_q = 3 * np.asarray(q, float)
    if np.min(np.linalg.norm((G - three_q) @ B, axis=1)) < 1e-3:
        return '120-degree'
    return 'incommensurate spiral'


def canonical_q(c, q):
    """A canonical representative: the equivalent with the smallest |q|, then lexicographic."""
    B = c.bvec()
    E = equivalents(c, q)
    n = np.round(np.linalg.norm(E @ B, axis=1), 6)
    cand = E[n == n.min()]
    cand = cand[np.lexsort((-cand[:, 0], -cand[:, 1]))]
    return [float(round(x, 6)) + 0.0 for x in cand[0]]


def easy_axis(table, K, nbasis, n1, n2, pattern):
    """Lowest-energy direction n for a collinear pattern (signs per site): returns unit vector."""
    def f(x):
        n = np.array([np.sin(x[0]) * np.cos(x[1]), np.sin(x[0]) * np.sin(x[1]), np.cos(x[0])])
        return energy(table, K, nbasis, n1, n2, np.outer(pattern, n))
    best = min((f([t, p]), t, p) for t in np.linspace(0, np.pi, 37) for p in np.linspace(0, 2 * np.pi, 73))
    r = minimize(f, best[1:], method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-13})
    n = np.array([np.sin(r.x[0]) * np.cos(r.x[1]), np.sin(r.x[0]) * np.sin(r.x[1]), np.cos(r.x[0])])
    return n * np.sign(n[2] if abs(n[2]) > 1e-9 else n[0])


def axis_angle_deg(n1, n2):
    return float(np.degrees(np.arccos(min(1.0, abs(np.dot(n1, n2)) / np.linalg.norm(n1) / np.linalg.norm(n2)))))
