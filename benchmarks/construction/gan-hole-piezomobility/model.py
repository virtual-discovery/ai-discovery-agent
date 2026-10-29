"""Shared analysis utilities for gan-hole-piezomobility construction (reference, ablations, robustness).

All mobilities are the zero-field 300 K tensors of pool.json, symmetrised: mu_ij = (mu_ij + mu_ji) / 2.
"""
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
POOL = json.loads((HERE/'pool.json').read_text())
S = POOL['states']
K = np.array(POOL['kgrid_crystal'])
KEYS = {'xx': (0, 0), 'yy': (1, 1), 'zz': (2, 2), 'yz': (1, 2), 'xz': (0, 2), 'xy': (0, 1)}


def mu(label, T=300, solver='ibte'):
    m = np.array(S[label][solver][str(T)]['mu'])
    return 0.5 * (m + m.T)


def as_dict(m):
    return {k: float(m[i, j]) for k, (i, j) in KEYS.items()}


def label(sign, pct, comp):
    return 'unstrained' if pct == 0 else f"{'compress' if sign < 0 else 'expand'}_{pct}%_e{comp}"


def series(sign, comp, solver='ibte'):
    """Tensors at 0, 1, 2 % along one signed component."""
    return [mu(label(sign, p, comp), solver=solver) for p in (0, 1, 2)]


def geometric_next(m0, m1, m2):
    """Third point of a saturating series: increments shrink geometrically when they keep their sign and shrink,
    otherwise the last increment is continued."""
    d1, d2 = m1 - m0, m2 - m1
    out = m2 + d2
    same = (np.sign(d1) == np.sign(d2)) & (np.abs(d2) < np.abs(d1)) & (np.abs(d1) > 1e-9)
    out[same] = (m2 + d2 ** 2 / np.where(np.abs(d1) > 1e-9, d1, 1.0))[same]
    return out


def linear_next(m0, m1, m2):
    return m2 + (m2 - m1)


def mirror(m, comp):
    """Mirror partner of a shear state: the diagonal is unchanged, the strained off-diagonal changes sign."""
    out = m.copy()
    i, j = {4: (1, 2), 5: (0, 2), 6: (0, 1)}[comp]
    out[i, j] = out[j, i] = -m[i, j]
    return out


def swap_xy(m):
    p = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    return p @ m @ p.T


BANDS = json.loads((HERE.parents[1]/'tasks/gan-hole-piezomobility/environment/lab/bands.json').read_text())


def gamma_levels(label_):
    """Top three Kramers levels at Gamma relative to the VBM (meV, descending) and the top three at kz = 0.1."""
    e = np.array(BANDS['energies_eV'][label_])[:, :6]
    g = int(np.where(np.all(np.abs(K) < 1e-9, axis=1))[0][0])
    z = int(np.where(np.all(np.abs(K - [0, 0, 0.1]) < 1e-9, axis=1))[0][0])
    vbm = e[g].max()
    return (np.sort(e[g])[::-1][0::2] - vbm) * 1000, (np.sort(e[z])[::-1][0::2] - vbm) * 1000


def bracket_zz():
    """mu_zz(compress 3% e3): past the band crossing the 2->3 % increment lies between zero (saturation) and the
    1->2 % increment (no saturation); the reference estimate is the geometric mean of the two bounds."""
    y1, y2 = mu('compress_1%_e3')[2, 2], mu('compress_2%_e3')[2, 2]
    low, high = y2, y2 + (y2 - y1)
    return float(np.sqrt(low * high)), float(low), float(high)


def inplane_max(m):
    return float(np.linalg.eigvalsh(m[:2, :2]).max())
