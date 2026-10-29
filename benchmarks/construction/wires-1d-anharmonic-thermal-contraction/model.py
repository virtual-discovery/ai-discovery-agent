"""Reference numerics for wires-1d-anharmonic-thermal-contraction (construction side).

All inputs are the archived texts in pool.json (dynamical matrices, SSCHA log excerpts, DFT files).
"""
import io
import json
import re
from pathlib import Path
import numpy as np
import dyn

RY_EV = 13.605693122994
KB_CM = 0.69503476            # k_B in cm^-1 / K
EV_J = 1.602176634e-19
BOHR_A = 0.529177210903


def load_pool(path):
    return json.loads(Path(path).read_text())


def dyn_from_texts(texts):
    """Parse a dict {filename: text} (ordered by q-star index) into one dynamical-matrix set."""
    import tempfile, os
    names = sorted(texts, key=lambda n: int(re.findall(r'(\d+)$', n)[0]))
    names = [n for n in names if 'Dynamical' in texts[n]]
    with tempfile.TemporaryDirectory() as tmp:
        for i, n in enumerate(names):
            Path(tmp, f'd{i + 1}').write_text(texts[n])
        return dyn.read_set(os.path.join(tmp, 'd'), len(names))


def grid_frequencies(d):
    """Frequencies (cm^-1) per commensurate q (rows, sorted ascending), with the q list (2pi/alat)."""
    return np.array([q for q, _ in d['blocks']]), np.array([np.sort(dyn.frequencies(d, phi)) for _, phi in d['blocks']])


def cv_mode(w, T):
    x = w / (KB_CM * T)
    return x ** 2 * np.exp(x) / np.expm1(x) ** 2


def zero_mask(qs, W, n_zero):
    """True for the n_zero lowest-|omega| modes at Gamma (rigid translations and, for a non-straight wire, the axial rotation)."""
    m = np.zeros(W.shape, bool)
    g = int(np.argmin(np.linalg.norm(qs, axis=1)))
    m[g, np.argsort(np.abs(W[g]))[:n_zero]] = True
    return m


def cv_ratio(freqs, T0, temps, deg, n_zero, sign=-1.0, derivative=True, per='vibrational'):
    """SSCHA heat capacity of the simulation supercell, C_V = sum c(omega)(1 - T/omega domega/dT), per vibrational mode in k_B.

    freqs: {T: (qs, W[nq, nmode])} auxiliary frequencies; mode-by-mode polynomial fit of omega(T) over `temps`.
    sign=+1 reproduces the archived notebook / paper Eq. (2) sign; derivative=False gives the harmonic c(omega(T)).
    per='vibrational' divides by (3N - n_zero); per='atoms' divides by 3N.
    """
    qs, W0 = freqs[T0]
    zm = zero_mask(qs, W0, n_zero)
    stack = np.array([freqs[T][1] for T in temps])
    total = 0.0
    for iq in range(W0.shape[0]):
        for j in range(W0.shape[1]):
            if zm[iq, j]:
                continue
            w = W0[iq, j]
            corr = 1.0
            if derivative:
                fit = np.polyfit(temps, stack[:, iq, j], deg)
                corr = 1 + sign * T0 / w * np.polyval(np.polyder(fit), T0)
            total += cv_mode(w, T0) * corr
    n = W0.size - (n_zero if per == 'vibrational' else 0)
    return total / n


def axial_stiffness(d, ncell, frac=0.01):
    """Long-wavelength LA stiffness k = M_cell v_LA^2 (eV per unit cell) from Fourier-interpolated force constants."""
    Rz, P = dyn.real_space_ifc(d, ncell)
    c = d['cell'][2, 2]
    qz = frac * 2 * np.pi / c
    w2, v = np.linalg.eigh(dyn.dynmat_at(d, Rz, P, qz))
    zc = np.array([np.sum(np.abs(v[2::3, k]) ** 2) for k in range(len(w2))])
    k = int(np.argmax(zc[:6]))
    return float(d['masses'].sum() * w2[k] / qz ** 2 * RY_EV)


def cell_volume_A3(d):
    return float(abs(np.linalg.det(d['cell'])) * BOHR_A ** 3)


def stress_zz(excerpt, which='sscha'):
    """zz component (GPa) of the SSCHA stress (or of the ab initio average stress) from a log excerpt."""
    lines = excerpt.splitlines()
    num = r'[-+]?\d+\.\d+(?:[eE][-+]?\d+)?'
    for i, l in enumerate(lines):
        if which == 'sscha' and 'STRESS TENSOR' in l:
            return float(re.findall(num, lines[i + 3])[2]), float(re.findall(num, lines[i + 3])[5])
        if which == 'average' and 'Ab initio average stress' in l:
            return float(re.findall(num, lines[i + 3])[2]), None
    raise ValueError('no stress block')


def converged(excerpt):
    return 'satisfy the convergence criteria' in excerpt


def slope(T, y, T0, deg=1):
    return float(np.polyval(np.polyder(np.polyfit(T, y, deg)), T0))


def alpha_per_K(dPdT_GPa_per_K, V_A3, k_eV):
    """alpha = V (dP_zz/dT)_L / k  (P_zz pressure-like: positive = compressive)."""
    return dPdT_GPa_per_K * 1e9 * V_A3 * 1e-30 / (k_eV * EV_J)


def softest_hessian(qs, W, n_zero):
    """Signed lowest frequency (cm^-1) of the Hessian excluding the n_zero rigid-body modes at Gamma."""
    zm = zero_mask(qs, W, n_zero)
    return float(W[~zm].min())


def crossing(T1, w1, T2, w2):
    """Temperature where the signed squared frequency crosses zero (linear interpolation in s = sign(w) w^2)."""
    s1, s2 = np.sign(w1) * w1 ** 2, np.sign(w2) * w2 ** 2
    return float(T1 + (T2 - T1) * s1 / (s1 - s2))
