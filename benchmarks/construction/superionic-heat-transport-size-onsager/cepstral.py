"""Multicomponent cepstral Green-Kubo analysis (numpy only), after Ercole, Marcolongo, Baroni,
Sci. Rep. 7, 15835 (2017) and Bertossa, Grasselli, Ercole, Baroni, PRL 122, 255901 (2019); the same
algorithm as the thermocepstrum code shipped in the source record.

Inputs are flux time series of shape (n_steps, 3) per current. Current 0 is the energy flux; the
other currents (e.g. the O velocity sum) are treated as the "non-heat" currents of the
multicomponent formula: the zero-frequency value of the reduced spectrum
    S_red(f) = 1 / [S^{-1}(f)]_00
gives the Onsager-corrected conductivity, independent of the energy zero of each species.
"""
import numpy as np
from scipy.special import polygamma, digamma

KB_EV = 8.617333262e-5          # eV/K
EV_PER_PS_A_TO_W_PER_MK = 1.602176634e-19 / 1e-12 / 1e-10   # eV/(ps A K) -> W/(m K)


def spectra(currents, dt):
    """Matrix periodogram averaged over Cartesian components.

    currents: list of arrays (n, 3). Returns freqs (THz, 1/ps), S (nf, M, M) one-sided-normalised
    so that S[0,0,0] -> int_{-inf}^{inf} <J(t)J(0)> dt = 2 int_0^inf <J(t) J(0)> dt (per component).
    """
    n = min(len(c) for c in currents)
    n -= n % 2
    F = np.array([np.fft.rfft(np.asarray(c[:n], float), axis=0) for c in currents])   # no mean removal: keeps f=0                                   # (M, nf, 3)
    S = np.einsum('ifl,jfl->fij', np.conj(F), F) * dt / n / 3.0     # (nf, M, M) Hermitian, averaged over l
    freqs = np.fft.rfftfreq(n, dt)
    return freqs, S, 3


def reduced(S):
    """Schur complement 1/(S^-1)_00 for each frequency (S: (nf, M, M))."""
    if S.shape[1] == 1:
        return S[:, 0, 0].real.copy()
    return 1.0 / np.linalg.inv(S)[:, 0, 0].real


def cepstral(spec, n_comp, n_currents, f_cut_index=None, P=None):
    """Cepstral estimate of the zero-frequency value of a (reduced) periodogram.

    spec: periodogram on rfft frequencies (each value ~ S(f) chi2_nu / (2 n_comp) with nu = 2(n_comp - n_currents + 1)).
    f_cut_index: keep only the first f_cut_index frequencies (low-pass, the F* of thermocepstrum).
    P: number of cepstral coefficients (default: AIC).
    Returns (S0, sigma_S0, P*, aic).
    """
    s = np.real(spec[:f_cut_index] if f_cut_index else spec)
    nu = 2 * (n_comp - n_currents + 1)
    L = np.log(s)
    # symmetric (even) extension and real IFFT -> cepstral coefficients
    ext = np.concatenate([L, L[-2:0:-1]])
    N = len(ext)
    C = np.fft.ifft(ext).real
    var_l = polygamma(1, nu / 2.0)
    # periodogram of the component-averaged (Schur-reduced) matrix ~ S chi2_nu / (2 n_comp)  (complex Wishart)
    mean_l = digamma(nu / 2.0) - np.log(n_comp)
    Nh = len(L)
    aic = None
    if P is None:
        tail = np.cumsum((C[:Nh] ** 2)[::-1])[::-1]
        aic = np.array([N * tail[p] / var_l + 2 * p for p in range(1, Nh - 1)])
        P = int(np.argmin(aic)) + 1
    L0 = C[0] + 2 * C[1:P].sum() - mean_l
    var_L0 = var_l * (4 * P - 2) / N
    S0 = np.exp(L0)
    return S0, S0 * np.sqrt(var_L0), P, aic


def kappa_prefactor(volume_A3, T):
    """kappa = V/(kB T^2) * (1/2) S(0) with S in (flux units)^2 ps and J an intensive flux in eV A^-2 ps^-1.

    If the stored flux is the extensive J*V (eV A / ps) divide S by V^2 first (done by the caller).
    Returns the factor converting S0/2 [eV^2 A^-4 ps^-2 * ps] * V/(kB T^2) to W/(m K).
    """
    return volume_A3 / (KB_EV * T ** 2) * EV_PER_PS_A_TO_W_PER_MK


def running_integral(j, dt, tmax):
    """Direct Green-Kubo running integral of the autocorrelation (averaged over components) up to tmax.
    Returns (t, I(t)) with I = int_0^t <J(t')J(0)> dt'. j: (n, 3)."""
    j = np.asarray(j, float) - np.asarray(j, float).mean(0)
    n = len(j)
    m = int(tmax / dt)
    F = np.fft.rfft(j, n=2 * n, axis=0)
    acf = np.fft.irfft(np.abs(F) ** 2, axis=0)[:m].mean(1) / (n - np.arange(m))
    I = np.concatenate([[0.0], np.cumsum(0.5 * (acf[1:] + acf[:-1]) * dt)])
    return np.arange(m) * dt, I
