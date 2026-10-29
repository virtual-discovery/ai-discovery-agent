"""Task definitions shared by the construction scripts."""
import numpy as np
import kmclib as K

COMPOSITIONS = {'A': {'Cr': 0.2, 'Nb': 0.2, 'Mo': 0.2, 'Ta': 0.2, 'W': 0.2},
                'B': {'V': 0.25, 'Cr': 0.25, 'Nb': 0.25, 'W': 0.25},
                'C': {'Nb': 0.50, 'Ta': 0.25, 'W': 0.25}}
T_ETA = 2000.0
T_ARRHENIUS = (1700.0, 2000.0, 2300.0)
N_HIDDEN = 100


def xvec(c):
    return [COMPOSITIONS[c].get(e, 0.0) for e in K.EL]


def hidden_indices():
    """100 clean alloy hops (not pure/dilute, saddle above both end states), never served."""
    rng = np.random.default_rng(2026)
    cand = np.where(~K.PATHO)[0]
    cand = cand[cand >= 35]
    return np.sort(rng.choice(cand, N_HIDDEN, replace=False))


def model(mask, n_orb=6, quad=True, lam=0.03, nsh=3, ce_mode='total'):
    K0, W, Q = K.train_kra(mask, n_orb=n_orb, quad=quad, lam=lam)
    return {'K0': K0, 'W': W, 'Q': Q, 'V': K.train_ce(mask, nsh=nsh, mode=ce_mode)}


def thermo_ideal(m):
    return dict(m, V=None)


def kra_simple(m):
    return {'K0': K.PURE_KRA.copy(), 'W': np.zeros((6, 0, 6)), 'Q': None, 'V': m['V']}


def arrhenius(Ts, Ds):
    """Least-squares slope of ln D vs 1/T -> activation energy (eV)."""
    x, y = 1.0 / np.asarray(Ts), np.log(np.asarray(Ds))
    return -np.polyfit(x, y, 1)[0] * K.KB
