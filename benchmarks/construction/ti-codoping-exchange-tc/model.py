"""Dilute-impurity spin model of doped Bi2Te3 built from the archived Jij table (construction utilities).

Lattice: the Bi sublattice of Bi2Te3 (archive Bi2Te3.cif: hexagonal a = 4.38296 A, c = 30.48701 A,
Bi1 at z = 0.4, Bi2 at z = 0.6 of the rhombohedral cell) in an orthorhombic supercell
(7L a) x (4L sqrt(3) a) x (L c) with 336 L^3 Bi sites.
Model: H = -sum_{i<j} J_ij e_i.e_j (the paper's H = -1/2 sum_{i,j}), unit vectors, J(meV) of the archived
shell and species pair applied to every neighbour of that shell; J = 0 beyond the 17 tabulated shells.
"""
import csv
import json
import os
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
CACHE = Path(os.environ.get('TI_CACHE', Path.home()/'.cache/ai-theorist/ti-codoping-exchange-tc'))
A_HEX, C_HEX = 4.382963164015663, 30.48701213585697
NAMES = 'Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd'.split()
ZOF = {n: z for n, z in zip(NAMES, list(range(21, 31)) + list(range(39, 49)))}
NAMEOF = {z: n for n, z in ZOF.items()}
KB = 0.08617333262  # meV/K


def load_table(path=None):
    path = Path(path) if path else CACHE/'Jij_table.csv'
    rows = list(csv.DictReader(open(path)))
    out = []
    for r in rows:
        out.append({'shell': int(r['i_shell']), 'R': float(r['R']), 'Rvec': [float(r['Rx']), float(r['Ry']), float(r['Rz'])],
                    'ilayer2': int(float(r['ilayer2'])), 'Z1': int(float(r['Zimp1'])), 'Z2': int(float(r['Zimp2'])),
                    'num_neighbors': int(float(r['num_neighbors'])), 'm1': float(r['m1']), 'm2': float(r['m2']),
                    'm01': float(r['m01']), 'm02': float(r['m02']), 'J': float(r['J']), 'uuid': r['J_data_uuid']})
    return out


TABLE = load_table() if (CACHE/'Jij_table.csv').exists() else None


def shells(table=None):
    table = table or TABLE
    sh = {}
    for r in table:
        sh.setdefault(r['shell'], (r['R'], r['ilayer2'], r['num_neighbors']))
    return [sh[k] for k in sorted(sh)]


def single_moments(table=None):
    table = table or TABLE
    m = {}
    for r in table:
        m[NAMEOF[r['Z1']]] = abs(r['m01']); m[NAMEOF[r['Z2']]] = abs(r['m02'])
    return m


def pair_J(a, b, table=None, shells_used=None):
    """J(shell) for species names a, b (None where the archive has no calculation)."""
    table = table or TABLE
    z1, z2 = sorted((ZOF[a], ZOF[b]))
    J = [None] * 17
    for r in table:
        if r['Z1'] == z1 and r['Z2'] == z2:
            J[r['shell']] = r['J']
    if shells_used is not None:
        J = [j if k in shells_used else 0.0 for k, j in enumerate(J)]
    return J


def bi_sites(L):
    """Cartesian Bi positions (n,3), sublattice labels (1/2) and orthorhombic box lengths."""
    a, c = A_HEX, C_HEX
    # orthorhombic cell a x sqrt(3) a x c: hexagonal lattice points (0,0),(1/2,1/2) in-plane, R-centring along c
    base = []
    for (fx, fy) in [(0, 0), (0.5, 0.5)]:
        for k, (dx, dy) in enumerate([(0, 0), (2/3, 1/3), (1/3, 2/3)]):
            # hexagonal fractional (dx, dy) -> cartesian with a1=(a,0), a2=(-a/2, a sqrt3/2)
            x = a * (dx - dy / 2); y = a * np.sqrt(3) / 2 * dy
            for zf, lab in [(0.4, 1), (0.6, 2)]:
                base.append((fx * a + x, fx * a * np.sqrt(3) + y, ((zf + k / 3) % 1) * c, lab))
    base = np.array(base)
    nx, ny, nz = 7 * L, 4 * L, L
    box = np.array([nx * a, ny * a * np.sqrt(3), nz * c])
    g = np.array([(i, j, k) for i in range(nx) for j in range(ny) for k in range(nz)], float)
    shifts = g * np.array([a, a * np.sqrt(3), c])
    pos = (shifts[:, None, :] + base[None, :, :3]).reshape(-1, 3) % box
    lab = np.tile(base[:, 3], len(g)).astype(int)
    return pos, lab, box


_SITE_CACHE = {}


def sample(L, comp, seed, Jfun, rcut=12.3, tol=0.02):
    """Random dilute configuration. comp: {species: concentration per Bi site}; Jfun(a, b) -> list of 17 J (meV).
    Returns species index per spin (order of comp keys), CSR neighbour list and J."""
    if L not in _SITE_CACHE:
        _SITE_CACHE[L] = bi_sites(L)
    pos, lab, box = _SITE_CACHE[L]
    rng = np.random.default_rng(seed)
    nsite = len(pos)
    counts = [int(round(c * nsite)) for c in comp.values()]
    chosen = rng.choice(nsite, size=sum(counts), replace=False)
    spec = np.repeat(np.arange(len(counts)), counts)
    p = pos[chosen]
    tree = cKDTree(p, boxsize=box)
    pairs = tree.query_pairs(rcut, output_type='ndarray')
    radii = np.array([s[0] for s in shells()])
    names = list(comp)
    Jtab = {(i, j): Jfun(names[i], names[j]) for i in range(len(names)) for j in range(len(names))}
    rows, cols, vals = [], [], []
    if len(pairs):
        d = p[pairs[:, 1]] - p[pairs[:, 0]]
        d -= box * np.round(d / box)
        r = np.linalg.norm(d, axis=1)
        k = np.abs(r[:, None] - radii[None, :]).argmin(axis=1)
        ok = np.abs(r - radii[k]) < tol
        assert ok.all(), 'unassigned distance'
        for (i, j), kk in zip(pairs, k):
            Jv = Jtab[(spec[i], spec[j])][kk]
            if Jv is None:
                raise ValueError(f'missing J for {names[spec[i]]}-{names[spec[j]]} shell {kk}')
            if Jv != 0.0:
                rows += [i, j]; cols += [j, i]; vals += [Jv, Jv]
    n = len(p)
    order = np.lexsort((cols, rows)) if rows else np.array([], int)
    rows, cols, vals = np.array(rows, int)[order], np.array(cols, int)[order], np.array(vals, float)[order]
    offs = np.zeros(n + 1, int)
    np.add.at(offs, rows + 1, 1)
    offs = np.cumsum(offs)
    return spec, offs, cols, vals, lab[chosen]


def write_graph(path, spec, offs, cols, vals):
    with open(path, 'wb') as f:
        np.array([len(spec)], '<i4').tofile(f)
        np.asarray(spec, '<i4').tofile(f)
        np.asarray(offs, '<i4').tofile(f)
        np.asarray(cols, '<i4').tofile(f)
        np.asarray(vals, '<f8').tofile(f)


def mean_field_tc(comp, Jfun, nshell=17):
    """Two-species (multi-sublattice) mean-field Tc in K: k_B T = (1/3) lambda_max(c_b J0_ab)."""
    names = list(comp)
    N = [s[2] for s in shells()]
    M = np.zeros((len(names), len(names)))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            J = Jfun(a, b)
            M[i, j] = comp[b] * sum((J[k] or 0.0) * N[k] for k in range(nshell))
    w, v = np.linalg.eig(M)
    k = np.argmax(w.real)
    return float(w[k].real / 3 / KB), v[:, k].real
