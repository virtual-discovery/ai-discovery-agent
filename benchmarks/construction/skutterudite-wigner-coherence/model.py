"""Reference Wigner (SMA) thermal-conductivity evaluator on the archived per-mode arrays.

kappa_P: populations (Peierls-Boltzmann, single-mode RTA) with degenerate subspaces summed as blocks;
kappa_C: coherences between non-degenerate bands (Simoncelli, Marzari, Mauri), Lorentzian overlap
with full widths Gamma = 2*gamma (phono3py gamma = imaginary part of the self-energy, THz).
With the archived inputs this reproduces every archived kappa_P_RTA and kappa_C to < 1e-7.
"""
import numpy as np

THZ_TO_EV = 4.135667696923859e-3   # h * 1 THz in eV
HBAR_EV_PS = 6.582119569e-4         # hbar in eV ps
EV = 1.602176634e-19
CUTOFF_THZ = 0.01                   # phono3py default cutoff frequency
CORRUPT_THZ = 1.0                   # archived gamma_isotope above this is a numerical failure (see README)


def volume(lattice):
    return abs(np.linalg.det(np.array(lattice, float)))


def kappa(fr, gam, gam_iso, cv, G, w, vol, iso=True, deg_thr=1e-4, shape='lorentzian', width_factor=2.0,
          coh_degenerate=False, pop_block=True, tau_factor=1.0, drop_q=None):
    """Return (kappa_P, kappa_C) in W/m/K. G[q,s,s'] = (1/3) sum_a |V^a_ss'(q)|^2 per irreducible q (A^2 THz^2)."""
    fr = np.asarray(fr, float); w = np.asarray(w, float).copy()
    if drop_q is not None:
        w[list(drop_q)] = 0.0
    g = np.asarray(gam, float) + (np.asarray(gam_iso, float) if iso and gam_iso is not None else 0.0)
    conv = EV * 1e4 * 1e-12 / (vol * 1e-30)        # eV/K * (A THz)^2 * ps / A^3 -> W/m/K
    ok = fr > CUTOFF_THZ
    with np.errstate(divide='ignore', invalid='ignore'):
        tau = np.where(ok, tau_factor / (4 * np.pi * g), 0.0)
    nb = fr.shape[1]
    same = np.abs(fr[:, :, None] - fr[:, None, :]) < deg_thr
    pair = ok[:, :, None] & ok[:, None, :]
    G = np.asarray(G, float)
    if G.ndim == 2:                      # only diagonal velocities available
        Gd = G
        G = None
    if pop_block and G is not None:
        pop = np.where(same & pair, G, 0.0).sum(axis=2) * cv * tau
    else:
        diag = np.einsum('qss->qs', G) if G is not None else Gd
        pop = diag * cv * tau
    kp = (pop.sum(axis=1) * w).sum() * conv / w.sum()
    if G is None:
        return kp, float('nan')
    hw = fr * THZ_TO_EV
    hg = width_factor * g * THZ_TO_EV / tau_factor
    with np.errstate(divide='ignore', invalid='ignore'):
        s = hg[:, :, None] + hg[:, None, :]
        dw = hw[:, :, None] - hw[:, None, :]
        if shape == 'lorentzian':
            lor = 0.5 * s / (dw ** 2 + 0.25 * s ** 2)
        else:                            # Gaussian with the same FWHM, same normalisation (integral pi)
            sig = s / 2.3548200450309493
            lor = np.pi * np.exp(-0.5 * (dw / sig) ** 2) / (np.sqrt(2 * np.pi) * sig)
        x = cv / hw
        pref = 0.25 * (hw[:, :, None] + hw[:, None, :]) * (x[:, :, None] + x[:, None, :])
        term = np.where(pair, G * pref * lor, 0.0)
    off = ~np.eye(nb, dtype=bool)[None]
    mask = off if coh_degenerate else (off & ~same)
    kc = (np.where(mask, term, 0.0).sum(axis=(1, 2)) * w).sum() * conv * HBAR_EV_PS / w.sum()
    return kp, kc


def repaired_isotope(gam_iso, mode='zero', fr=None):
    gi = np.array(gam_iso, float)
    bad = gi > CORRUPT_THZ
    if mode == 'zero':
        gi[bad] = 0.0
    elif mode == 'band_median':
        for q, s in zip(*np.where(bad)):
            gi[q, s] = np.median(np.delete(gi[:, s], q))
    return gi, bad


def t_half(T, kp, kc, mode='ratio'):
    """Lowest T with kC/kP = 0.5 by linear interpolation between grid temperatures (nan if not reached;
    'extrapolate' continues the last two points)."""
    T = np.asarray(T, float); r = np.asarray(kc) / np.asarray(kp)
    idx = np.where(r >= 0.5)[0]
    if len(idx) == 0 or idx[0] == 0:
        if mode == 'extrapolate':
            return float(T[-2] + (0.5 - r[-2]) * (T[-1] - T[-2]) / (r[-1] - r[-2]))
        return float('nan')
    i = idx[0]
    if mode in ('ratio', 'extrapolate'):
        return float(T[i - 1] + (0.5 - r[i - 1]) * (T[i] - T[i - 1]) / (r[i] - r[i - 1]))
    tt = np.linspace(T[i - 1], T[i], 20001)
    p = np.interp(tt, T, kp); c = np.interp(tt, T, kc)
    return float(tt[np.argmax(c >= 0.5 * p)])
