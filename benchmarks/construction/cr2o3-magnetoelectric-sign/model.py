"""Lattice-mediated magnetoelectric response from the parsed archive (pool.json). Construction-side only.

alpha_xx = mu0 * dM_x/dE_x  with  M = (net spin moment per cell) / Omega,
dM/dE = sum_m (dm_x/du_m) (p_m,x / k_m) / Omega,   p_m,a = sum_{kappa,b} Z*_{kappa,ab} e_{m,kappa b}
Units: dm/du in muB/A, p in e, k in eV/A^2, E in V/A -> alpha in ps/m via PREF.
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
POOL = json.loads((HERE/'pool.json').read_text())
NAMES = ['Cr1', 'Cr2', 'Cr3', 'Cr4', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6']
MODES = ['Eu_1_x', 'Eu_2_x', 'Eu_3_x', 'Eu_4_x']
MU0, MUB = 1.25663706212e-6, 9.2740100783e-24
CELL = np.array(POOL['structures']['phonon_POSCAR']['cell_A'])
OMEGA = abs(np.linalg.det(CELL))
PREF = MU0 * MUB * 1e-10 / (OMEGA * 1e-30) * 1e12
AMU_EV = 1.66053906660e-27 / 1.602176634e-19 * 1e-20   # amu * (rad/s)^2 -> eV/A^2


def amp(tag):
    return float(tag.replace('_', '.'))


def frac0():
    return np.array(POOL['born']['0_0']['reference']['poscar']['frac'])


def born_dipole(rec):
    """Physical dipole change (e*A) of one VASP LCALCPOL run: p is printed in 'electrons Angst' (charge -|e|)."""
    return -(np.array(rec['p_ion']) + np.array(rec['p_elc']))


def born_charges(tags=('0_02', '-0_02'), displacement='poscar', transpose=False, asr=True):
    """Z[kappa, a(polarization), b(displacement)] from finite differences against the 0_0 reference."""
    ref = POOL['born']['0_0']['reference']
    d0 = born_dipole(ref)
    f0 = frac0()
    pos0 = np.array(ref['outcar_positions_cartesian_A'])
    Z = np.zeros((10, 3, 3))
    for k, n in enumerate(NAMES):
        for b, c in enumerate('xyz'):
            vals = []
            for t in tags:
                rec = POOL['born'][t][f'{n}_{c}']
                if displacement == 'poscar':
                    du = np.array(rec['poscar']['frac']) - f0
                    du -= np.round(du)
                    u = (du @ CELL)[k, b]
                else:                                 # displacement read from the OUTCAR positions
                    u = (np.array(rec['outcar_positions_cartesian_A']) - pos0)[k, b]
                vals.append((born_dipole(rec) - d0) / u)
            Z[k, :, b] = np.mean(vals, axis=0)
    if transpose:
        Z = np.transpose(Z, (0, 2, 1))
    if asr:
        Z = Z - Z.mean(axis=0)
    return Z


def fc_modes():
    fc = POOL['fc_eigen']
    return {m: (np.array(fc['vectors'][m]), fc['values'][m]) for m in MODES}


def polarity(Z, vec):
    return np.einsum('kab,kb->a', Z, vec)


def elk_slope(mode, tags, quantity='total_moment_muB'):
    a = np.array([amp(t) for t in tags])
    m = np.array([POOL['elk'][mode][t][quantity][0] for t in tags])
    return float(np.polyfit(a, m, 1)[0])


def vasp_slope(mode, tags, quantity='total'):
    a = np.array([amp(t) for t in tags])
    recs = [POOL['vasp_moments'][mode][t] for t in tags]
    if quantity == 'total':
        m = np.array([r['total_magnetization_muB'][0] for r in recs])
    else:
        m = np.array([sum(row[-1] for row in r['ion_magnetization_x']) for r in recs])
    return float(np.polyfit(a, m, 1)[0])


def alpha(Z, slopes, stiff=None, vecs=None):
    """Total alpha (ps/m) and per-mode contributions for the domain of `slopes`."""
    modes = fc_modes()
    contrib = {}
    for m in MODES:
        vec, k = modes[m]
        if vecs is not None:
            vec = vecs[m]
        if stiff is not None:
            k = stiff[m]
        contrib[m] = PREF * polarity(Z, vec)[0] * slopes[m] / k
    tot = sum(contrib.values())
    return tot, contrib


def fractions(contrib):
    tot = sum(contrib.values())
    return {m: c / tot for m, c in contrib.items()}


def cr_pairs(frac, cell=CELL):
    """Cr-Cr distances along the trigonal axis for the 4 Cr (first four atoms): returns sorted list of (d, i, j)."""
    pos = np.array(frac)[:4] @ cell
    out = []
    for i in range(4):
        for j in range(i + 1, 4):
            best = min(np.linalg.norm(pos[j] - pos[i] + np.array(s) @ cell)
                       for s in [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)])
            out.append((best, i, j))
    return sorted(out)


def domain(frac, cr_mz, cell=CELL, pair=None):
    """'out_pointing' if the moments of a Cr pair stacked along the trigonal axis point away from each other.

    pair=None uses the nearest (face-sharing, 2.62 A) pair; pair=(i, j) forces a given pair of the same column.
    """
    pos = np.array(frac)[:4] @ cell
    if pair is None:
        d, i, j = cr_pairs(frac, cell)[0]
    else:
        i, j = pair
    rij = pos[j] - pos[i]                      # same trigonal column: positions (t, t, t) in the cell
    si, sj = np.sign(cr_mz[i] * rij[2]), np.sign(cr_mz[j] * -rij[2])
    assert si == sj
    return 'out_pointing' if si < 0 else 'in_pointing'
