"""Reference analysis for amaranta-2d-anisotropic-ferromagnets (teacher model in model.py).

Ground states (Luttinger-Tisza on the isotropic part, checked against the full model), the
anisotropy of the ordered state and its single-ion / exchange split, Kitaev-Gamma
decompositions in the ideal cubic frame of the octahedra, DM vectors, the zone-centre spin-wave
gap and the magnetostatic dipolar anisotropy.
"""
import numpy as np
from scipy.optimize import minimize
import model as M

S_SPIN = 1.5


def table_energy(c, table, K, n1, n2, spins):
    return M.supercell_energy(table, K, len(c.basis), n1, n2, spins)


def Jq_matrix(c, table, q):
    """nb x nb isotropic-part (trace/3) Fourier matrix; E/site of a flat spiral = min eigenvalue / 2 (unit vectors)."""
    nb = len(c.basis)
    Mq = np.zeros((nb, nb), complex)
    for k0, k, m1, m2, J in table:
        d = (c.basis[k] + m1 * c.cell[0] + m2 * c.cell[1] - c.basis[k0])[:2]
        Mq[k0, k] += np.trace(J) / 3 * np.exp(1j * q @ d)
    return 0.5 * (Mq + Mq.conj().T) / 2


def bvec(c):
    return 2 * np.pi * np.linalg.inv(c.cell[:2, :2]).T


def luttinger_tisza(c, grid=61):
    tab = c.bond_table(iso=True)
    B = bvec(c)
    nb = len(c.basis)
    def lam(x):
        return float(np.linalg.eigvalsh(Jq_matrix(c, tab, np.array(x) @ B))[0])
    best = min((lam([h, k]), h, k) for h in np.linspace(0, 1, grid) for k in np.linspace(0, 1, grid))
    r = minimize(lam, best[1:], method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-13})
    q = np.array(r.x)
    w, V = np.linalg.eigh(Jq_matrix(c, tab, q @ B))
    return q, float(r.fun) / nb * nb, V[:, 0], w


def collinear_patterns(c):
    """Candidate collinear patterns on a 2x2 supercell (sign per site, lab order k, i, j)."""
    nb = len(c.basis)
    pats = {}
    import itertools
    n1 = n2 = 2
    N = nb * n1 * n2
    seen = []
    for signs in itertools.product([1, -1], repeat=N - 1):
        s = np.array((1,) + signs, float)
        pats[tuple(s)] = s
    return pats, n1, n2


def best_collinear(c, iso=True):
    """Lowest-energy collinear state among all sign patterns of a 2x2 and a 1x2/2x1 supercell (isotropic part)."""
    tab = c.bond_table(iso=iso)
    nb = len(c.basis)
    best = None
    import itertools
    for n1, n2 in ((1, 1), (2, 1), (1, 2), (2, 2)):
        N = nb * n1 * n2
        for signs in itertools.product([1, -1], repeat=N - 1):
            s = np.array((1,) + signs, float)
            spins = np.outer(s, [0, 0, 1.0])
            e = table_energy(c, tab, 0.0, n1, n2, spins) / N
            if best is None or e < best[0] - 1e-9:
                best = (e, n1, n2, s)
    return best


def anisotropy_matrix(c, n1, n2, pattern, part='all'):
    """Symmetric 3x3 M with E(n)/site = n.M.n + const for the collinear pattern with moments +-n.
    part: 'all', 'sia' (single-ion only) or 'exchange'."""
    tab = c.bond_table()
    N = len(pattern)
    K = c.K if part in ('all', 'sia') else 0.0
    if part == 'sia':
        tab = []
    E = lambda n: table_energy(c, tab, K, n1, n2, np.outer(pattern, n)) / N
    Mm = np.zeros((3, 3))
    ex = np.eye(3)
    e0 = 0.0
    for a in range(3):
        Mm[a, a] = E(ex[a])
    for a in range(3):
        for b in range(a + 1, 3):
            v = (ex[a] + ex[b]) / np.sqrt(2)
            w = (ex[a] - ex[b]) / np.sqrt(2)
            Mm[a, b] = Mm[b, a] = (E(v) - E(w)) / 2
    return Mm


def principal(Mm):
    w, V = np.linalg.eigh(Mm)
    return w, V


def gap_uniform(Mm, S=S_SPIN):
    """Zone-centre LSWT gap of the uniform precession (acoustic mode) for E(n)/site = n.M.n (unit vectors)."""
    w, V = np.linalg.eigh(Mm)
    return 2.0 / S * np.sqrt(max(w[1] - w[0], 0) * max(w[2] - w[0], 0)), w, V


# ---------------------------------------------------------------- Kitaev frame
def bridging_ligands(c, xi, xj, nlig=2):
    """The ligands closest to the midpoint of the bond xi -> xj (all images)."""
    st = c.st
    mid = (xi + xj) / 2
    cands = []
    for p, s in zip(st.pos, st.species):
        if s in M.MAGNETIC:
            continue
        for m1 in range(-2, 3):
            for m2 in range(-2, 3):
                x = p + m1 * c.cell[0] + m2 * c.cell[1]
                cands.append((np.linalg.norm(x - mid), s, x))
    cands.sort(key=lambda t: t[0])
    return [x for _, s, x in cands[:nlig]]


def kitaev_frame(c, xi, xj):
    """Ideal cubic frame (e1, e2, e3) of the bond xi -> xj: e3 perpendicular to the bond, at arccos(1/sqrt3)
    from +z, on the side of the normal of the Cr-L-Cr plane; e1 + e2 + e3 = sqrt(3) z."""
    d = xj - xi
    u = d / np.linalg.norm(d)
    z = np.array([0, 0, 1.0])
    L1, L2 = bridging_ligands(c, xi, xj)
    nrm = np.cross(u, L2 - L1)
    nrm /= np.linalg.norm(nrm)
    if nrm[2] < 0:
        nrm = -nrm
    perp = np.cross(z, u)                  # in-plane, perpendicular to the bond
    perp /= np.linalg.norm(perp)
    sgn = np.sign(nrm @ perp)
    e3 = sgn * np.sqrt(2 / 3) * perp + np.sqrt(1 / 3) * z
    # e1, e2: rotate e3 by +-120 deg about z
    def rz(v, th):
        R = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
        return R @ v
    e1, e2 = rz(e3, 2 * np.pi / 3), rz(e3, -2 * np.pi / 3)
    return e1, e2, e3, nrm


def jkgg(J, e1, e2, e3):
    """Least-squares (exact for C2h bonds) J, K, Gamma, Gamma' of the symmetric part of J in the frame e1, e2, e3."""
    Js = (J + J.T) / 2
    basis = [np.eye(3), np.outer(e3, e3), np.outer(e1, e2) + np.outer(e2, e1),
             np.outer(e1, e3) + np.outer(e3, e1) + np.outer(e2, e3) + np.outer(e3, e2)]
    X = np.array([b.ravel() for b in basis]).T
    coef, res, *_ = np.linalg.lstsq(X, Js.ravel(), rcond=None)
    fit = (X @ coef).reshape(3, 3)
    return coef, float(np.abs(fit - Js).max())


def dm_vector(J):
    """D with e_i.J_A.e_j = D.(e_i x e_j), J_A = (J - J^T)/2."""
    return np.array([(J[1, 2] - J[2, 1]) / 2, (J[2, 0] - J[0, 2]) / 2, (J[0, 1] - J[1, 0]) / 2])


# ---------------------------------------------------------------- dipolar
MU0_MUB2_OVER_4PI = 0.05368  # meV A^3 for two moments of 1 muB (mu0 muB^2 / 4 pi)


def dipolar_mae(c, moment_muB=3.0, rmax=200.0, pattern=None):
    """E(in-plane, along x) - E(along z) per magnetic site (meV) of point dipoles, collinear pattern on the
    primitive cell (signs per basis site; default ferromagnetic)."""
    nb = len(c.basis)
    pattern = np.ones(nb) if pattern is None else np.asarray(pattern, float)
    A = c.cell[:2, :2]
    nmax = int(rmax / min(np.linalg.norm(A, axis=1))) + 1
    m1, m2 = np.meshgrid(np.arange(-nmax, nmax + 1), np.arange(-nmax, nmax + 1), indexing='ij')
    R = m1.ravel()[:, None] * c.cell[0][None, :2] + m2.ravel()[:, None] * c.cell[1][None, :2]
    Ez = Ex = 0.0
    for k0 in range(nb):
        for k in range(nb):
            d = R + (c.basis[k] - c.basis[k0])[None, :2]
            r = np.linalg.norm(d, axis=1)
            ok = (r > 0.1) & (r < rmax)
            d, r = d[ok], r[ok]
            s = pattern[k0] * pattern[k]
            # E = 1/2 sum C m^2 [mi.mj - 3 (mi.r)(mj.r)/r^2] / r^3
            Ez += 0.5 * s * np.sum(1.0 / r ** 3)
            Ex += 0.5 * s * np.sum((1.0 - 3 * (d[:, 0] / r) ** 2) / r ** 3)
    C = MU0_MUB2_OVER_4PI * moment_muB ** 2
    return C * (Ex - Ez) / nb


def gap_linearized(c, table, K, pattern, n0, S=S_SPIN):
    """Zone-centre LSWT frequencies from the linearised Landau-Lifshitz equations of the primitive cell
    (images summed), de_i/dt = -(1/S) e_i x dE/de_i, about the collinear state e_k = pattern_k n0."""
    nb = len(c.basis)
    Jsum = np.zeros((nb, nb, 3, 3))
    for k0, k, m1, m2, J in table:
        Jsum[k0, k] += J
    n0 = np.asarray(n0, float) / np.linalg.norm(n0)
    u = np.cross(n0, [0.3, 0.5, 0.8]); u /= np.linalg.norm(u)
    v = np.cross(n0, u)
    def field(e):
        h = np.einsum('klab,lb->ka', Jsum, e)
        h[:, 2] += 2 * K * e[:, 2]
        return h
    def rhs(x):
        e = np.array([pattern[k] * n0 + x[2 * k] * u + x[2 * k + 1] * v for k in range(nb)])
        de = -np.cross(e, field(e)) / S
        return np.array([[d @ u, d @ v] for d in de]).ravel()
    eps = 1e-6
    A = np.zeros((2 * nb, 2 * nb))
    for i in range(2 * nb):
        dx = np.zeros(2 * nb); dx[i] = eps
        A[:, i] = (rhs(dx) - rhs(-dx)) / (2 * eps)
    w = np.linalg.eigvals(A)
    return np.sort(np.abs(w.imag))
