"""Forward model of the HYSPEC fixed-orientation data and the IN20 constant-E scans of MnTe."""
import numpy as np
import lswt
import sym

# HYSPEC kinematics (derived from the authors' trajectory in Fig1c.mat; see README): Ei = 37 meV, omega = 32 deg,
# reduction lattice a = 4.1001 A, c = 6.8429 A. E = Ei - 2.0721 |k_i - Q|^2 with Q in that lattice.
_A = np.sqrt(6.48814446)
_KX = 26.27916711 / (2 * _A)
_KZ = -_KX * np.tan(np.radians(32.0))
_G = min(np.roots([1.4884, 2.44 * _KZ, 7.79528992]).real)    # root giving c = 6.843 A
EI = _KX ** 2 + _KZ ** 2


def hyspec_energy(u, v):
    h = np.asarray(u) + np.asarray(v)
    l = -2.23 * np.asarray(u) + 1.6 * np.asarray(v)
    return EI - (_KX - h * _A) ** 2 - (_KZ - l * _G) ** 2


def uv_to_hkl(u, v):
    return np.array([u + v, 0.0, -2.23 * u + 1.6 * v])


def mn_formfactor(qA):
    s2 = (qA / (4 * np.pi)) ** 2
    A, a, B, b, C, c, D = 0.4220, 17.684, 0.5948, 6.005, 0.0043, -0.609, -0.0219   # Mn2+ <j0>
    return A * np.exp(-a * s2) + B * np.exp(-b * s2) + C * np.exp(-c * s2) + D


ASTAR_DIR_DEG = 30.0      # a* lies 30 deg from Cartesian x (|| a) in the ab plane


def spectral(hkl, J, D, domains=((0, 1 / 3), (120, 1 / 3), (240, 1 / 3)), sign=+1):
    """Mode energies (2,), unpolarized S_perp (2,), chiral term (2,) for Neel order along a*-type axes.

    domains: (angle from a* in deg, weight); sign: antiphase domain (+1: sublattice A spin along +n).
    Chiral term = i * Qhat . (S x S^dagger) convention: M_n = -2 Im(S^{yz}_n) in the Blume-Maleyev frame,
    reported per unit net Neel sign; multiply by the imbalance to get the measured value."""
    q = lswt.qcart(hkl)
    qh = q / np.linalg.norm(q)
    E = None
    Sp = np.zeros(2)
    Mc = np.zeros(2)
    for ang, w in domains:
        fr = lswt.frames(ASTAR_DIR_DEG + ang, sign)
        e, Sab = lswt.correlations(hkl, J, D, fr)
        if E is None:
            E = e
        P = np.eye(3) - np.outer(qh, qh)
        Sp += w * np.einsum('ab,nab->n', P, Sab).real
        # chiral: Qhat . (i * eps_{abc} S^{bc}) ~ 2 Im part of antisymmetric part
        anti = Sab - np.transpose(Sab, (0, 2, 1))
        vec = np.stack([anti[:, 1, 2], anti[:, 2, 0], anti[:, 0, 1]], axis=1)   # S^{yz}-S^{zy} etc
        Mc += w * (1j * vec @ qh).real
    return E, Sp, Mc


def window_points(axis_vals, fixed_lo, fixed_hi, along='u', step=0.004):
    """Sub-sampled (u, v) points for each cut bin: along 'u' (v integrated) or along 'v' (u integrated)."""
    fixed = np.arange(fixed_lo + step / 2, fixed_hi, step)
    pts = []
    for x in axis_vals:
        sub = np.linspace(x - 0.004 + 0.001, x + 0.004 - 0.001, 4) if along == 'u' else np.linspace(x - 0.01 + 0.0025, x + 0.01 - 0.0025, 4)
        if along == 'u':
            pts.append([(a, b) for a in sub for b in fixed])
        else:
            pts.append([(b, a) for a in sub for b in fixed])
    return pts


class CutModel:
    """Precomputes mode energies, S_perp and chiral weights on the integration points of one cut."""

    def __init__(self, pts, J, D, domains=((0, 1 / 3), (120, 1 / 3), (240, 1 / 3))):
        self.pts = pts
        self.data = []
        for group in pts:
            rows = []
            for (u, v) in group:
                hkl = uv_to_hkl(u, v)
                E, Sp, Mc = spectral(hkl, J, D, domains)
                _, _, Mc1 = spectral(hkl, J, D, ((0, 1.0),))
                ff = mn_formfactor(np.linalg.norm(lswt.qcart(hkl))) ** 2
                rows.append((hyspec_energy(u, v), E, Sp * ff, Mc1 * ff))
            self.data.append(rows)

    def predict(self, sigma, scale, bg=(0.0, 0.0), eta=0.0, kind='sym'):
        out = []
        for rows in self.data:
            acc = []
            for (Ep, E, Sp, Mc) in rows:
                g = np.exp(-0.5 * ((Ep - E) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
                if kind == 'sym':
                    acc.append(scale * (Sp * g).sum() + bg[0] + bg[1] * (Ep - 31.0))
                else:
                    acc.append(scale * eta * (Mc * g).sum())
            out.append(np.mean(acc))
        return np.array(out)
