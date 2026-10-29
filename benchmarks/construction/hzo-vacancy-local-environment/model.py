"""Local-environment and population utilities for the hzo-vacancy-local-environment construction (numpy only).

Conventions (stated in the task instruction):
  E_f = E(vacancy cell) - E(pristine cell, same composition) + mu_O, mu_O = -10.09 eV (O-poor).
  First-shell cations: Hf/Zr within 2.3 A of the oxygen in the staged (ideal) structures.
  Population: dilute neutral vacancies at 700 K in a macroscopic random solid solution (each cation site Zr with
  probability x); site fraction = <exp(-E_f/kT)> over all O sites; distribution over the number of Zr neighbours.
"""
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
POOL = json.loads((HERE/'pool.json').read_text())
MU_POOR = POOL['mu_O_eV']['poor']
KB = 8.617333262e-5
T = 700.0
KT = KB * T
CUTOFF = 2.3
XS = ['0.25', '0.50', '0.75']
CLASSES = [(4, k) for k in range(5)] + [(3, k) for k in range(4)]   # (coordination, n_Zr)


def parse_cif(text):
    lat, sp, frac = {}, [], []
    for line in text.splitlines():
        s = line.split()
        if line.startswith('_cell_length_') or line.startswith('_cell_angle_'):
            lat[s[0]] = float(s[1])
        elif len(s) == 7 and s[0] in ('Hf', 'Zr', 'O'):
            sp.append(s[0])
            frac.append([float(v) for v in s[3:6]])
    assert all(abs(lat[k] - 90) < 1e-9 for k in lat if 'angle' in k)
    cell = np.diag([lat['_cell_length_a'], lat['_cell_length_b'], lat['_cell_length_c']])
    return cell, sp, np.array(frac)


def environments(cell, sp, frac, cutoff=CUTOFF):
    """{O atom index: (coordination, n_Zr)} using cations within `cutoff` (minimum image)."""
    cat = [i for i, s in enumerate(sp) if s != 'O']
    out = {}
    for i, s in enumerate(sp):
        if s != 'O':
            continue
        d = frac[cat] - frac[i]
        d -= np.round(d)
        r = np.linalg.norm(d @ cell, axis=1)
        nb = [cat[j] for j in np.where(r <= cutoff)[0]]
        out[i] = (len(nb), sum(sp[j] == 'Zr' for j in nb))
    return out


CELLS = {x: parse_cif(POOL['cif'][x]) for x in XS}
ENV = {x: environments(*CELLS[x]) for x in XS}
EF = {x: {int(s): e - POOL['vasp']['bulk_eV'][x] + MU_POOR for s, e in POOL['vasp']['vacancy_eV'][x].items()} for x in XS}


def binom(x, cls):
    n, k = cls
    return 0.5 * math.comb(n, k) * x ** k * (1 - x) ** (n - k)


def hyper(x, cls, ncat=32):
    """Environment frequency in a 32-cation cell with exactly x*32 Zr (fixed-composition sampling)."""
    nz, (n, k) = round(x * ncat), cls
    return 0.5 * math.comb(nz, k) * math.comb(ncat - nz, n - k) / math.comb(ncat, n)


def population(x, factor, freq=binom):
    """factor: {cls: <exp(-E_f/kT)> for that class}. Returns distribution over n_Zr, mean n_Zr, site fraction."""
    w = {cls: freq(x, cls) * factor[cls] for cls in CLASSES}
    z = sum(w.values())
    dist = [sum(v for cls, v in w.items() if cls[1] == k) / z for k in range(5)]
    return {'distribution': dist, 'mean_nZr': sum(k * p for k, p in enumerate(dist)), 'site_fraction': z,
            'o3_share': sum(v for cls, v in w.items() if cls[0] == 3) / z}


def factors_from_energies(energy_of_class):
    return {cls: math.exp(-energy_of_class(cls) / KT) for cls in CLASSES}


def features(cls, x, kind):
    n, z = cls
    o3 = float(n == 3)
    if kind == 'pooled':                 # n_Zr only, sublattices pooled
        return [1.0, z]
    if kind == 'local':                  # sublattice offset + per-sublattice Zr slopes, composition-blind
        return [1.0, o3, o3 * z, (1 - o3) * z]
    if kind == 'local_x':                # ... + linear composition term
        return [1.0, o3, o3 * z, (1 - o3) * z, x]
    if kind == 'local_offsets':          # ... + one offset per composition
        return [1.0, o3, o3 * z, (1 - o3) * z, float(x == 0.25), float(x == 0.75)]
    raise KeyError(kind)


def fit(rows, kind):
    """rows: list of (x, cls, E_f). Returns least-squares coefficients."""
    A = np.array([features(cls, x, kind) for x, cls, e in rows])
    y = np.array([e for x, cls, e in rows])
    return np.linalg.lstsq(A, y, rcond=None)[0]


def predict(coef, x, cls, kind):
    return float(np.dot(coef, features(cls, x, kind)))


def truth_population(x):
    """Reference population from all 64 archived sites of the composition's cell: per class, the Boltzmann factor
    averaged over its archived sites; classes absent from the cell take the per-sublattice additive fit of that
    cell (their random-alloy weight is below 0.3 % of the population)."""
    rows = [(float(x), ENV[x][s], e) for s, e in EF[x].items()]
    coef = fit(rows, 'local')
    factor, counts = {}, {}
    for cls in CLASSES:
        es = [e for s, e in EF[x].items() if ENV[x][s] == cls]
        counts[cls] = len(es)
        factor[cls] = float(np.mean(np.exp(-np.array(es) / KT))) if es else math.exp(-predict(coef, float(x), cls, 'local') / KT)
    pop = population(float(x), factor)
    absent = sum(binom(float(x), cls) * factor[cls] for cls in CLASSES if counts[cls] == 0) / pop['site_fraction']
    return pop, absent
