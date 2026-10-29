"""Shubnikov-de Haas Fourier analysis used to build and calibrate hgte-gated-subband-densities.

The authors' chart (FFT_experiment.ipynb) is reproduced exactly by `chart(..., quantity='sxx',
window='hamming', bmax=0.8, nmax=8e11, nres=0.025e11, pad=1)` for the 70 nm sample (max relative
difference 1e-14 against 70nm_FT.csv). Densities are e*f/h (one spin-resolved Fermi pocket per peak),
in 1e11 cm^-2.
"""
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import periodogram

E_CHARGE, H_PLANCK = 1.60217663e-19, 6.62607004e-34      # values used in the authors' notebook
LBYW = 600 / 200


def sheet(vi, vxx, vxy, rpre, quantity='sxx'):
    rxx, rxy = vxx / vi * rpre, vxy / vi * rpre
    rho = rxx / LBYW
    if quantity == 'sxx':
        return rho / (rho ** 2 + rxy ** 2)
    if quantity == 'rxx':
        return rho
    raise ValueError(quantity)


def inv_b_grid(bmax, nmax, nres, bmin=None):
    dbi = 1 / (nmax * 1e4) * (E_CHARGE / H_PLANCK) / 2
    if bmin is None:
        bmin = 1 / (1 / bmax + (nmax / nres) * dbi)
    return 1 / np.arange(1 / bmax, 1 / bmin + dbi, dbi), 1 / dbi


def spectrum(B, y, bmax=0.8, nmax=8e11, nres=0.025e11, window='hamming', pad=1, bmin=None):
    ib, fs = inv_b_grid(bmax, nmax, nres, bmin)
    x = interp1d(B, y, kind='linear', bounds_error=False, fill_value=0)(ib)
    nfft = len(x) * pad if pad > 1 else None
    f, p = periodogram(x, fs=fs, window=window, return_onesided=True, detrend='linear', scaling='density', nfft=nfft)
    return (E_CHARGE / H_PLANCK * f / 1e4) * 1e-11, p


def chart(vg, B, data, rpre, quantity='sxx', **kw):
    """vg: (G,), B: (F,), data: (G, F, 3) = VI, Vxx, Vxy. Returns n (1e11 cm^-2) and power (G, N)."""
    out = []
    for g in range(len(vg)):
        y = sheet(data[g, :, 0], data[g, :, 1], data[g, :, 2], rpre, quantity)
        n, p = spectrum(B, y, **kw)
        out.append(p)
    return n, np.array(out)


def local_maxima(n, logp, lo=0.15, hi=8.0):
    """Local maxima of a log spectrum with parabolic refinement: list of (n, logp)."""
    out = []
    dn = n[1] - n[0]
    for i in range(1, len(n) - 1):
        if lo <= n[i] <= hi and logp[i] > logp[i - 1] and logp[i] >= logp[i + 1]:
            a, b, c = logp[i - 1], logp[i], logp[i + 1]
            den = a - 2 * b + c
            d = 0.5 * (a - c) / den if den != 0 else 0.0
            out.append((n[i] + d * dn, b - 0.25 * (a - c) * d))
    return out


def averaged_log(vg, power, center, half_width):
    m = np.abs(vg - center) <= half_width + 1e-9
    return np.log10(power[m]).mean(0)


def peak_near(n, logp, guess, window):
    """The strongest local maximum within +-window of guess."""
    cands = [p for p in local_maxima(n, logp) if abs(p[0] - guess) <= window]
    if not cands:
        return None
    return max(cands, key=lambda p: p[1])[0]


def hall_density(B, vi, vxy, rpre, bmax=0.2):
    """Single-carrier low-field Hall density (1e11 cm^-2, electrons positive) from the slope of Rxy vs B."""
    m = np.abs(B) <= bmax + 1e-9
    slope = np.polyfit(B[m], vxy[m] / vi[m] * rpre, 1)[0]
    return 1 / (slope * 1.602176634e-19) * 1e-15
