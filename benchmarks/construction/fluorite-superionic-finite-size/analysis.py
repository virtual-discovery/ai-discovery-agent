"""Reference analysis for fluorite-superionic-finite-size (numpy only; mirrored in solution/solve.py).

Protocol (the graded definitions):
- Mobile-ion self-diffusion D [A^2/ps]: slope/6 of a straight-line fit to the mobile-species MSD over lag
  times 20-100 ps, where lag time = lag frames x (timestep x dump interval) read from each run's input deck.
  Frame: relative to the immobile sublattice, MSD_lab = MSD + a(2 + a) MSD_cm(mobile), a = N_mob m_mob / (N_imm m_imm)
  (the author's "laboratory" convention; frame choice changes D by < 5%).
- c_V per formula unit [k_B]: central finite difference of the mean total energy per formula unit over the
  run temperatures (mean measured temperature). T_c = maximum of a parabola through the largest c_V point and
  its two neighbours.
- Activation energies: least-squares fit of ln D vs 1/(k_B T) on the 768-atom (4x4x4) cell, below-T_c window
  [T_low, T_c - 100 K] and above-T_c window [T_c + 100 K, T_max] (temperatures on the 50 K grid).
"""
import numpy as np

KB = 8.617333262e-5  # eV/K
MASS = {'CaF2': (40.078, 18.998403), 'PbF2': (207.2, 18.998403), 'UO2': (238.028909, 15.9994),
        'AgI': (126.90447, 107.8682)}   # (type 1, type 2) masses from the input decks
MOBILE_RATIO = {'CaF2': 2.0, 'PbF2': 2.0, 'UO2': 2.0, 'AgI': 1.0}   # N(type 2) / N(type 1); type 2 is mobile
FU_ATOMS = {'CaF2': 3, 'PbF2': 3, 'UO2': 3, 'AgI': 2}
D_WINDOW = (20.0, 100.0)


def lab_msd(material, msd):
    """Mobile-species MSD relative to the immobile sublattice from msd_CM_B4 columns."""
    m1, m2 = MASS[material]
    a = MOBILE_RATIO[material] * m2 / m1
    return msd[:, 2] + a * (2 + a) * msd[:, 6]


def diffusion(t, msd, window=D_WINDOW):
    sel = (t >= window[0]) & (t <= window[1])
    return float(np.polyfit(t[sel], msd[sel], 1)[0] / 6)


def heat_capacity(T, E_fu):
    """Central finite differences of <E>/f.u. vs T, in k_B per formula unit (one-sided at the ends)."""
    T, E = np.asarray(T, float), np.asarray(E_fu, float)
    return np.gradient(E, T) / KB


def peak_temperature(T, cv):
    T, cv = np.asarray(T, float), np.asarray(cv, float)
    k = int(np.argmax(cv))
    if k == 0 or k == len(T) - 1:
        return float(T[k])
    c = np.polyfit(T[k - 1:k + 2], cv[k - 1:k + 2], 2)
    return float(np.clip(-c[1] / (2 * c[0]), T[k - 1], T[k + 1]))


def arrhenius(T, D):
    T, D = np.asarray(T, float), np.asarray(D, float)
    slope, intercept = np.polyfit(1 / (KB * T), np.log(D), 1)
    return float(-slope), float(np.exp(intercept))
