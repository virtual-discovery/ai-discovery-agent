"""Defect-thermodynamics utilities for the gan-irradiation-defect-levels construction (numpy only).

Conventions
- Charge q = NELECT(neutral) - NELECT. Energies in eV.
- VASP's 'average (electrostatic) potential at core' is the potential energy of an electron (eV); the
  electrostatic potential of the defect is therefore -(V_defect - V_perfect).
- Point-charge (PC) potential of a periodic array of unit charges with a neutralising background in an
  anisotropic dielectric (Kumagai & Oba, PRB 89, 195205 (2014)), Ewald-summed.
- eFNV correction: E_corr = E_PC + q * < (V_defect - V_perfect) + q * phi_PC >_far, with
  E_PC = -(q^2/2) * Madelung potential (positive for a cubic-like cell), and the average taken over atoms
  outside the sphere inscribed in the supercell around the defect.
"""
import itertools
import json
import math
from pathlib import Path
import numpy as np
from scipy.special import erfc

HERE = Path(__file__).resolve().parent
KE = 14.399645  # e^2/(4 pi eps0) in eV*Angstrom
NEUTRAL_NELECT = {'Va_Ga1': 381, 'Va_N1': 379, 'Ga_i1': 387, 'Ga_i2': 387, 'N_i1': 389, 'N_i2': 389,
                  'Ga_N1': 382, 'N_Ga1': 386, 'VV_GaN1': 376, 'C_N1': 383}


def load_pool():
    return json.loads((HERE/'pool.json').read_text())


def split_name(name):
    base, q = name.rsplit('_', 1)
    return base, int(q)


def charge(rec, name):
    base, q = split_name(name)
    assert NEUTRAL_NELECT[base] - rec['nelect'] == q, name
    return q


# ---------------------------------------------------------------- structure mapping
def cart(rec, key='frac_final'):
    return np.array(rec[key]) @ np.array(rec['cell_A'])


def min_image(d, cell):
    f = d @ np.linalg.inv(cell)
    f -= np.round(f)
    return f @ cell


def map_to_perfect(rec, perfect):
    """Match atoms of a defect cell to perfect-cell sites using the starting (near-ideal) geometry.

    Returns (pairs, vacancies, extras, substitutions, defect_center_cart): pairs maps defect-cell atom index ->
    perfect-cell site index; the defect centre is the mean of the vacant sites, the substituted sites and the
    final positions of the extra (interstitial) atoms."""
    cell = np.array(perfect['cell_A'])
    pp, ps = cart(perfect), np.array(perfect['symbols'])
    ds = np.array(rec['symbols'])
    start = np.array(rec['frac_start']) @ np.array(rec['cell_A'])
    d = np.linalg.norm(min_image(start[:, None, :] - pp[None, :, :], cell), axis=2)
    pairs, used, extras, subs = {}, set(), [], []
    for i in np.argsort(d.min(axis=1)):
        free = np.where(np.isin(np.arange(len(pp)), list(used)), np.inf, d[i])
        j = int(np.argmin(free))
        if free[j] < 0.7:
            pairs[int(i)] = j
            used.add(j)
            if ds[i] != ps[j]:
                subs.append((int(i), j))
        else:
            extras.append(int(i))
    vacancies = [j for j in range(len(pp)) if j not in used]
    final = cart(rec)
    # an atom displaced at the start that relaxed back onto the site it left is a host atom
    vac_c, ext_c = list(vacancies), list(extras)
    for i in extras:
        for j in vac_c:
            if ds[i] == ps[j] and np.linalg.norm(min_image(final[i] - pp[j], cell)) < 0.6:
                vac_c.remove(j)
                ext_c.remove(i)
                break
    pts = [pp[j] for j in vac_c] + [final[i] for i in ext_c] + [pp[j] for _, j in subs]
    ref = pts[0]
    pts = [ref + min_image(p - ref, cell) for p in pts]
    return pairs, vacancies, extras, subs, np.mean(pts, axis=0)


# ---------------------------------------------------------------- anisotropic Ewald
class Ewald:
    def __init__(self, cell, eps, gamma=None, nreal=4, nrec=6):
        self.cell = np.array(cell, float)
        self.eps = np.array(eps, float)
        self.vol = abs(np.linalg.det(self.cell))
        self.einv = np.linalg.inv(self.eps)
        self.sdet = math.sqrt(np.linalg.det(self.eps))
        self.g = gamma or math.sqrt(np.trace(self.eps) / 3) / self.vol ** (1 / 3)
        rng = range(-nreal, nreal + 1)
        self.R = np.array(list(itertools.product(rng, rng, rng))) @ self.cell
        rec = 2 * np.pi * np.linalg.inv(self.cell).T
        rng = range(-nrec, nrec + 1)
        G = np.array([g for g in itertools.product(rng, rng, rng) if g != (0, 0, 0)]) @ rec
        GeG = np.einsum('ij,jk,ik->i', G, self.eps, G)
        self.G, self.wG = G, 4 * np.pi / self.vol * np.exp(-GeG / (4 * self.g ** 2)) / GeG

    def _real(self, r, exclude_origin):
        d = r[None, :] - self.R
        s = np.sqrt(np.einsum('ij,jk,ik->i', d, self.einv, d))
        if exclude_origin:
            s = s[s > 1e-8]
        return float(np.sum(erfc(self.g * s) / s)) / self.sdet

    def potential(self, r):
        """Potential (in e/(4 pi eps0) per Angstrom units) at r of unit charges at the lattice points + background."""
        r = np.asarray(r, float)
        return self._real(r, False) + float(np.sum(self.wG * np.cos(self.G @ r))) - math.pi / (self.vol * self.g ** 2)

    def madelung(self):
        return (self._real(np.zeros(3), True) + float(np.sum(self.wG)) - 2 * self.g / (math.sqrt(math.pi) * self.sdet)
                - math.pi / (self.vol * self.g ** 2))


def inscribed_radius(cell):
    cell = np.array(cell)
    vol = abs(np.linalg.det(cell))
    return min(vol / np.linalg.norm(np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])) for i in range(3)) / 2


def efnv(rec, perfect, q, eps, region='outside_inscribed', alignment=True, ewald=None, potential_sign=-1):
    """Return dict with E_PC, alignment term and total correction (eV) for charge q.

    potential_sign = -1 treats the VASP core potentials as electron potential energies (correct);
    +1 treats them as electrostatic potentials (the sign error)."""
    if q == 0:
        return {'E_pc': 0.0, 'alignment': 0.0, 'total': 0.0, 'n_sites': 0}
    cell = np.array(perfect['cell_A'])
    ew = ewald or Ewald(cell, eps)
    e_pc = -0.5 * q * q * KE * ew.madelung()
    pairs, vac, extras, subs, center = map_to_perfect(rec, perfect)
    vd, vp = np.array(rec['core_potentials_V']), np.array(perfect['core_potentials_V'])
    pd = cart(rec)
    rows = []
    for i, j in pairs.items():
        if (i, j) in subs:
            continue
        r = min_image(pd[i] - center, cell)
        rows.append((np.linalg.norm(r), vd[i] - vp[j], q * KE * ew.potential(r)))
    rows = np.array(rows)
    if region == 'outside_inscribed':
        mask = rows[:, 0] > inscribed_radius(cell)
    elif region == 'all':
        mask = np.ones(len(rows), bool)
    else:
        mask = rows[:, 0] > float(region)
    align = float(np.mean(-potential_sign * rows[mask, 1] + rows[mask, 2])) if alignment else 0.0
    return {'E_pc': e_pc, 'alignment': q * align, 'total': e_pc + q * align, 'n_sites': int(mask.sum()),
            'raw_far_shift': float(np.mean(rows[mask, 1]))}


def makov_payne_pc(cell, q, eps_scalar):
    """Isotropic point-charge (Madelung) term only: q^2 * alpha / (2 eps L) via the isotropic Ewald."""
    ew = Ewald(cell, np.eye(3) * eps_scalar)
    return -0.5 * q * q * KE * ew.madelung()


# ---------------------------------------------------------------- thermodynamics
def levels_from_energies(E, vbm, emax):
    """Thermodynamic transition levels from corrected energies E[q] (eV, without the q*E_F term).

    E_f(q, E_F) = E[q] + q*(vbm + E_F) with E_F measured from the VBM. Returns (inside, all): the levels
    ((q_low_EF, q_high_EF), position above the VBM) between consecutive ground states, restricted to
    [0, emax] for `inside` and over an extended window for `all`."""
    qs = sorted(E)
    grid = np.linspace(-1.0, emax + 1.0, 60001)
    form = np.array([E[q] + q * (vbm + grid) for q in qs])
    gs = [qs[k] for k in form.argmin(axis=0)]
    levels = []
    for k in range(1, len(gs)):
        a, b = gs[k - 1], gs[k]
        if a != b:
            levels.append(((a, b), (E[b] - E[a]) / (a - b) - vbm))
    inside = [(p, x) for p, x in levels if 0 <= x <= emax]
    return inside, levels


def ground_state(E, vbm, ef):
    """(q, formation-energy part) of the lowest E[q] + q*(vbm + ef)."""
    return min(((q, E[q] + q * (vbm + ef)) for q in E), key=lambda t: t[1])
