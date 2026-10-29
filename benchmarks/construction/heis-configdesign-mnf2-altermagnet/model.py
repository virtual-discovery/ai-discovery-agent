"""Two-sublattice Heisenberg utilities for the MnF2 task construction (numpy only).

Convention (the authors'): E_cell(q) - E_cell(0) = sum_{i in cell} sum_{j != i} J_ij (1 - e_i.e_j), unit vectors,
every ordered pair counted, energies per 6-atom primitive cell (2 Mn). Per Mn this is
(1/N_Mn) sum_i sum_j J_ij (1 - e_i.e_j), so the J values are the same in the per-Mn convention.

Released spirals are flat and 'uniform': the moment on an Mn at position r (fractional) has azimuth
2 pi q.r + const, including the basis offset (verified from the in-cell constrained moments). Hence for any
Mn-Mn bond vector d (fractional) the pair term is J(d) (1 - cos 2 pi q.d) and the cell energy is
E_cell(q) = 2 sum_shells J_s X_s(q), X_s(q) = sum_{d in shell, from one site} (1 - cos 2 pi q.d).
"""
import json
from pathlib import Path
import numpy as np

POOL = json.loads((Path(__file__).with_name('pool.json')).read_text())
CELL = np.array(POOL['cell_A'])
FRAC = np.array(POOL['fractional_positions'])
TAU = {'A': FRAC[0], 'B': FRAC[1]}     # Mn1 (CSV sublattice 1) at (1/2,1/2,1/2), Mn2 at the origin


def _bonds(span=4):
    """All Mn-Mn bond vectors (fractional) from sublattice A, labelled same/inter, sorted by length."""
    out = []
    rng = range(-span, span + 1)
    for n1 in rng:
        for n2 in rng:
            for n3 in rng:
                R = np.array([n1, n2, n3], float)
                for kind, d in (('same', R), ('inter', R + TAU['B'] - TAU['A'])):
                    if np.allclose(d, 0):
                        continue
                    out.append((round(float(np.linalg.norm(d @ CELL)), 4), kind, d))
    return sorted(out, key=lambda x: (x[0], x[1]))


BONDS = _bonds()
SHELLS = []          # list of dicts: distance, kind, vectors (fractional, from one site), count
for dist, kind, d in BONDS:
    if SHELLS and abs(SHELLS[-1]['distance'] - dist) < 2e-3 and SHELLS[-1]['kind'] == kind:
        SHELLS[-1]['vectors'].append(d)
    else:
        SHELLS.append({'distance': dist, 'kind': kind, 'vectors': [d]})
SHELLS = [dict(s, vectors=np.array(s['vectors']), count=len(s['vectors'])) for s in SHELLS if s['distance'] < 12.0]


def altermagnetic(shell):
    """True when the same-distance shell splits into two sublattice-inequivalent bond classes (rutile 4/mmm, site mmm):
    same-sublattice lattice vectors (n1, n2, n3) with n1 != 0 and n2 != 0."""
    v = shell['vectors']
    return shell['kind'] == 'same' and bool(np.any((np.abs(v[:, 0]) > 1e-9) & (np.abs(v[:, 1]) > 1e-9)))


def X(q, nshell, convention='uniform'):
    """Per-Mn design rows X_s(q) for the first nshell shells. E_cell = 2 * X @ J.

    convention: 'uniform' (released spirals), 'lattice_phase' (wrong: both Mn of a cell share the angle
    2 pi q.R, ignoring the basis offset), 'neel_reference' (wrong: q measured from the Neel state, i.e. an
    extra pi on every inter-sublattice bond, so the q=0 reference is taken to be the Neel state)."""
    q = np.atleast_2d(np.asarray(q, float))
    cols = []
    for s in SHELLS[:nshell]:
        v = s['vectors']
        if convention == 'lattice_phase' and s['kind'] == 'inter':
            v = np.round(v - (TAU['B'] - TAU['A']))      # cell translation only
        ph = 2 * np.pi * q @ v.T
        if convention == 'neel_reference' and s['kind'] == 'inter':
            cols.append((-(1 - np.cos(ph))).sum(axis=1))
        else:
            cols.append((1 - np.cos(ph)).sum(axis=1))
    return np.array(cols).T


def fit(Xm, y, intercept=False, w=None):
    M = np.c_[np.ones(len(y)), Xm] if intercept else Xm
    if w is not None:
        M, y = M * w[:, None], y * w
    c, *_ = np.linalg.lstsq(M, y, rcond=None)
    return (c[0], c[1:]) if intercept else (0.0, c)


def jq_matrix(J, q):
    """2x2 exchange Fourier matrix J_ss'(q) = sum_d J(d) exp(i 2 pi q.d) (sublattice-symmetrised model)."""
    q = np.asarray(q, float)
    aa = sum(Js * np.cos(2 * np.pi * s['vectors'] @ q).sum() for Js, s in zip(J, SHELLS) if s['kind'] == 'same')
    ab = sum(Js * np.exp(2j * np.pi * s['vectors'] @ q).sum() for Js, s in zip(J, SHELLS) if s['kind'] == 'inter')
    return aa, ab


def neel_mean_field_K(J):
    """Classical mean-field Neel temperature (unit vectors, ordered-pair convention): k_B T = (2/3) sum_j J_0j sigma_j."""
    tot = sum(Js * s['count'] * (1 if s['kind'] == 'same' else -1) for Js, s in zip(J, SHELLS))
    return 2 / 3 * tot / 0.08617333262
