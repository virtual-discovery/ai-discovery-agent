"""Build Ca3Ti2O7 domain-wall supercells from the bulk Cmc2_1 (A2_1am) structure.

The bulk displacement field u (relative to ideal I4/mmm sites) is mapped onto the sites of a
supercell; the sites on one side of the wall planes take a symmetry-transformed field.
"""
import numpy as np


def ideal_sites(bulk_sym, bulk_pos, bulk_cell, ref_sym, ref_frac):
    """Reference I4/mmm sites of the 48-atom orthorhombic cell (a = polar axis along x, b along y).
    ref_frac: tetragonal-cell fractional positions (24 atoms). Returns sites (48,3) matched to atoms."""
    a, b, c = np.diag(bulk_cell)
    t1 = np.array([a / 2, b / 2, 0.0]); t2 = np.array([-a / 2, b / 2, 0.0]); t3 = np.array([0, 0, c])
    base = ref_frac[:, 0:1] * t1 + ref_frac[:, 1:2] * t2 + ref_frac[:, 2:3] * t3
    cand, csym = [], []
    for i in range(-2, 3):
        for j in range(-2, 3):
            cand.append(base + i * t1 + j * t2); csym += list(ref_sym)
    cand = np.concatenate(cand); csym = np.array(csym)
    sym = np.array(bulk_sym)
    # origin: best global shift (scan Ti candidates), then nearest-site matching per species
    best = None
    ti = np.nonzero(sym == 'Ti')[0]
    for k in np.nonzero(csym == 'Ti')[0]:
        sh = bulk_pos[ti[0]] - cand[k]
        d = mic(bulk_pos[ti][:, None] - cand[csym == 'Ti'][None] - sh, bulk_cell)
        err = np.linalg.norm(d, axis=2).min(1).sum()
        if best is None or err < best[0]:
            best = (err, sh)
    sh = best[1]
    sites = np.zeros_like(bulk_pos)
    for s in set(bulk_sym):
        m = sym == s
        st = cand[csym == s] + sh
        d = mic(bulk_pos[m][:, None] - st[None], bulk_cell)
        k = np.linalg.norm(d, axis=2).argmin(1)
        sites[m] = bulk_pos[m] - d[np.arange(m.sum()), k]
    # remove mean displacement so the origin sits on the ideal inversion centre lattice
    u = bulk_pos - sites
    sites = sites + 0.0
    # snap the site origin so that one Ca rock-salt site is exactly an inversion centre (0,0,0)
    return sites, u


def mic(d, cell):
    f = d @ np.linalg.inv(cell)
    f -= np.round(f)
    return f @ cell


class Bulk:
    def __init__(self, sym, pos, cell, ref_sym, ref_frac):
        self.sym = np.array(sym)
        self.cell = np.asarray(cell, float)
        sites, u = ideal_sites(sym, np.asarray(pos, float), self.cell, ref_sym, ref_frac)
        # put an inversion centre of the ideal structure at the origin: the Ca site nearest z=0 on the rock-salt layer
        ca = np.nonzero(self.sym == 'Ca')[0]
        fz = (sites[ca] @ np.linalg.inv(self.cell))[:, 2] % 1
        # rock-salt Ca sit at the Ca site whose z coincides with an inversion centre; choose the site with
        # fractional z nearest 0 among those with a partner at -z
        k = ca[np.argmin(np.minimum(fz, 1 - fz))]
        self.origin = sites[k].copy()
        self.sites = sites - self.origin
        self.u = u
        self.inv = np.linalg.inv(self.cell)
        self.frac = (self.sites @ self.inv) % 1.0

    def lookup(self, s):
        """Displacement and species at Cartesian site(s) s (orthorhombic frame)."""
        f = (np.atleast_2d(s) @ self.inv) % 1.0
        d = f[:, None, :] - self.frac[None]
        d -= np.round(d)
        dist = np.linalg.norm(d @ self.cell, axis=2)
        k = dist.argmin(1)
        if dist[np.arange(len(k)), k].max() > 0.05:
            raise ValueError('site not on the reference lattice')
        return self.u[k], self.sym[k]

    def variant(self, op, s):
        """Displacement at site s for a domain variant.
        op: 'A' (original), 'I' (inversion through the origin), 'IT' (inversion then tetragonal
        translation t1 = (a/2, b/2, 0)), 'T' (translation t1 only)."""
        t1 = np.array([self.cell[0, 0] / 2, self.cell[1, 1] / 2, 0.0])
        s = np.atleast_2d(s)
        if op == 'A':
            u, sp = self.lookup(s)
        elif op == 'I':
            u, sp = self.lookup(-s); u = -u
        elif op == 'T':
            u, sp = self.lookup(s - t1)
        elif op == 'IT':
            u, sp = self.lookup(t1 - s); u = -u
        else:
            raise ValueError(op)
        return u, sp


def supercell_sites(bulk, frame, n, normal_axis=0, nz=1):
    """Ideal sites of a supercell. frame 'tet': cell = n*t1' x 2*t2' x c (t' = tetragonal axes);
    'orth_a': n*a x b x c; 'orth_b': n*b x a x c rotated so the wall normal is x; returns
    (cell, sites in supercell frame, sites in the bulk orthorhombic frame)."""
    a, b, c = np.diag(bulk.cell)
    t1 = np.array([a / 2, b / 2, 0.0]); t2 = np.array([-a / 2, b / 2, 0.0]); t3 = np.array([0, 0, c])
    if frame == 'tet':
        A = np.array([n * t1, 2 * t2, nz * t3])
        reps = [(i, j, k) for i in range(n) for j in range(2) for k in range(nz)]
        unit = [t1, t2, t3]
    elif frame == 'orth_a':
        A = np.array([n * np.array([a, 0, 0]), [0, b, 0], nz * t3])
        unit = [np.array([a, 0, 0]), np.array([0, b, 0]), t3]
        reps = [(i, j, k) for i in range(n) for j in range(1) for k in range(nz)]
    elif frame == 'orth_b':
        A = np.array([[0, n * b, 0], [-a, 0, 0], nz * t3])
        unit = [np.array([0, b, 0]), np.array([-a, 0, 0]), t3]
        reps = [(i, j, k) for i in range(n) for j in range(1) for k in range(nz)]
    elif frame == 'orth_c':
        A = np.array([[a, 0, 0], [0, b, 0], n * t3])
        unit = [np.array([a, 0, 0]), np.array([0, b, 0]), t3]
        reps = [(0, 0, k) for k in range(n)]
    else:
        raise ValueError(frame)
    # basis: all bulk sites inside one unit cell spanned by `unit`
    U = np.array(unit)
    Uinv = np.linalg.inv(U)
    cand = []
    for i in range(-3, 4):
        for j in range(-3, 4):
            for k in range(-1, 2):
                cand.append(bulk.sites + i * bulk.cell[0] + j * bulk.cell[1] + k * bulk.cell[2])
    cand = np.concatenate(cand)
    fu = cand @ Uinv
    inside = np.all((fu > -1e-6) & (fu < 1 - 1e-6), axis=1)
    basis = cand[inside]
    sites = np.concatenate([basis + i * U[0] + j * U[1] + k * U[2] for i, j, k in reps])
    # rotate into a frame with the first supercell vector along x
    e1 = A[0] / np.linalg.norm(A[0]); e3 = np.array([0, 0, 1.0]); e2 = np.cross(e3, e1)
    R = np.array([e1, e2, e3])
    return A @ R.T, sites @ R.T, sites, R


def build_wall(bulk, frame, n, op_right, offset=0.0, nz=1, axis=0):
    """Two walls per cell: sites with fractional coordinate (along `axis`) in [offset, offset+1/2)
    take variant A, the others variant op_right."""
    cell, s_sc, s_bulk, R = supercell_sites(bulk, frame, n, nz=nz)
    f = (s_sc @ np.linalg.inv(cell))[:, axis]
    left = ((f - offset) % 1.0) < 0.5
    u = np.zeros_like(s_bulk); sym = np.empty(len(s_bulk), dtype=object)
    uA, spA = bulk.variant('A', s_bulk[left]); u[left] = uA; sym[left] = spA
    if op_right is None:
        uB, spB = bulk.variant('A', s_bulk[~left])
    else:
        uB, spB = bulk.variant(op_right, s_bulk[~left])
    u[~left] = uB; sym[~left] = spB
    pos = (s_bulk + u) @ R.T
    return list(sym), pos, cell
