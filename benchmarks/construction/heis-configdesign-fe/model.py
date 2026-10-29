"""Heisenberg-model utilities for the Fe task construction (numpy only).

Convention (the authors'): E/atom - E_FM/atom = (1/N) sum_i sum_{j != i, all images} J_ij (1 - e_i.e_j),
i.e. H = -sum_{i != j} J_ij e_i.e_j with unit vectors and each pair counted twice.
For a conical spiral with cone angle theta this gives sin^2(theta) sum_j J_0j (1 - cos(2 pi q.n_j)).
"""
import json
from pathlib import Path
import numpy as np

POOL = json.loads((Path(__file__).with_name('pool.json')).read_text())
A = np.array(POOL['primitive_cell_A'])
_n = np.array([(i, j, k) for i in range(-7, 8) for j in range(-7, 8) for k in range(-7, 8) if (i, j, k) != (0, 0, 0)])
_d = np.round(np.linalg.norm(_n @ A, axis=1), 4)
RADII = sorted(set(_d))[:12]
COUNTS = [int((_d == r).sum()) for r in RADII]


def spiral_design(q, theta_deg, nshell):
    q = np.atleast_2d(q)
    s2 = np.sin(np.radians(np.atleast_1d(theta_deg))) ** 2
    cols = [(1 - np.cos(2 * np.pi * (q @ _n[_d == r].T))).sum(axis=1) * s2 for r in RADII[:nshell]]
    return np.array(cols).T


def lattice_constant(state):
    """bcc lattice constant of a random-state supercell from its shortest interatomic vector."""
    cell, frac = np.array(state['cell_A']), np.array(state['fractional_positions'])
    pos = frac @ cell
    shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]) @ cell
    d = np.linalg.norm(pos[None, :, :] + shifts[:, None, :] - pos[0], axis=2)
    return d[d > 1e-6].min() * 2 / np.sqrt(3)


def random_design(state, nshell):
    cell, frac = np.array(state['cell_A']), np.array(state['fractional_positions'])
    m = np.array(state['moments_muB'])
    e = m / np.linalg.norm(m, axis=1, keepdims=True)
    pos = frac @ cell
    nat = len(pos)
    a_state = lattice_constant(state)
    radii = [r * a_state / POOL['lattice_A'] for r in RADII]
    reach = radii[nshell - 1] + 1e-3
    span = int(np.ceil(reach / min(np.linalg.norm(cell, axis=1)))) + 1
    shifts = np.array([(i, j, k) for i in range(-span, span + 1) for j in range(-span, span + 1) for k in range(-span, span + 1)]) @ cell
    row = np.zeros(nshell)
    for i in range(nat):
        d = pos[None, :, :] + shifts[:, None, :] - pos[i]
        dist = np.round(np.linalg.norm(d, axis=2), 3)
        for s, r in enumerate(radii[:nshell]):
            mask = np.abs(dist - r) < 2e-3
            row[s] += (mask * (1 - e @ e[i])[None, :]).sum()
    return row / nat


def fit(X, y, intercept=False):
    M = np.c_[np.ones(len(y)), X] if intercept else X
    c, *_ = np.linalg.lstsq(M, y, rcond=None)
    return (c[0], c[1:]) if intercept else (0.0, c)


def stiffness(J):
    """D in meV A^2: E_flat(q) - E_flat(0) ~ D q^2 per atom, q in 1/A (isotropic cubic average)."""
    return sum(Js * (np.linalg.norm(_n[_d == r] @ A, axis=1) ** 2).sum() / 6 for Js, r in zip(J, RADII))


def mean_field_tc_K(J):
    """Classical mean-field Tc for unit vectors: k_B Tc = (2/3) sum_j J_0j (pair counted twice convention)."""
    return 2 / 3 * sum(Js * c for Js, c in zip(J, COUNTS)) / 0.08617333262


def direction_of(q):
    q = np.array(q)
    if np.allclose(q, 0):
        return 'gamma'
    u = q / np.abs(q[np.abs(q) > 1e-9]).min()
    return '[' + ''.join(str(int(round(x))) for x in u) + ']'


def rpa_tc_K(J, ngrid=40):
    """Classical RPA (Tyablikov) Tc: k_B Tc = (1/3) / <1/(J'(0) - J'(q))>_BZ with J' = 2J (pair counted twice)."""
    B = 2 * np.pi * np.linalg.inv(A).T
    f = (np.arange(ngrid) + 0.5) / ngrid
    qs = np.array(np.meshgrid(f, f, f, indexing='ij')).reshape(3, -1).T  # reciprocal-lattice coordinates
    Jq = np.zeros(len(qs))
    J0 = 0.0
    for Js, r in zip(J, RADII):
        nn = _n[_d == r]
        Jq += 2 * Js * np.cos(2 * np.pi * qs @ nn.T).sum(axis=1)
        J0 += 2 * Js * len(nn)
    return (1 / 3) / np.mean(1 / (J0 - Jq)) / 0.08617333262
