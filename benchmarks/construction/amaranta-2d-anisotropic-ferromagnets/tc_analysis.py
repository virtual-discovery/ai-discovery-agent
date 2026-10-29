"""Ordering-temperature estimators from the reference Monte Carlo (mc_reference.py output)."""
import numpy as np


def _fit(T, y, w=None, deg=3):
    return np.polyfit(T, y, deg, w=w)


def crossing(data, L1, L2, window=None):
    """Binder-cumulant crossing of sizes L1 < L2: the first temperature (going up) where U_L1 - U_L2 changes
    from negative (ordered side: the larger cell has the larger U) to positive; located on the cubic fits of
    U(T) inside the bracketing data interval (linear interpolation if the fits do not cross there)."""
    def series(L):
        r = sorted([x for x in data if x['L'] == L and (window is None or window[0] <= x['T_K'] <= window[1])],
                   key=lambda x: x['T_K'])
        return {x['T_K']: x for x in r}
    s1, s2 = series(L1), series(L2)
    T = sorted(set(s1) & set(s2))
    diff = [s1[t]['U'] - s2[t]['U'] for t in T]
    idx = next((i for i in range(len(T) - 1) if diff[i] < 0 <= diff[i + 1]), None)
    if idx is None:
        return None
    lo, hi = T[idx], T[idx + 1]
    def fit(sr):
        t = np.array(sorted(sr))
        return _fit(t, [sr[x]['U'] for x in t], [1 / max(sr[x]['U_err'], 1e-4) for x in t], deg=min(3, len(t) - 1))
    roots = np.roots(np.polysub(fit(s1), fit(s2)))
    roots = [x.real for x in roots if abs(x.imag) < 1e-9 and lo <= x.real <= hi]
    if roots:
        return float(roots[0])
    return float(lo + (hi - lo) * (-diff[idx]) / (diff[idx + 1] - diff[idx]))


def chi_peak(data, L):
    r = sorted([x for x in data if x['L'] == L], key=lambda x: x['T_K'])
    T = np.array([x['T_K'] for x in r])
    chi = np.array([x['chi'] for x in r])
    i = int(np.argmax(chi))
    lo, hi = max(0, i - 2), min(len(T), i + 3)
    p = np.polyfit(T[lo:hi], chi[lo:hi], 2)
    return float(-p[1] / (2 * p[0])) if p[0] < 0 else float(T[i])


def analyse(mc, window=None):
    data = mc['results']
    sizes = sorted({r['L'] for r in data})
    pairs = {f'{a}-{b}': crossing(data, a, b, window) for i, a in enumerate(sizes) for b in sizes[i + 1:]}
    out = {'sizes': sizes, 'crossings': pairs,
           'Tc_binder_largest_pair': pairs[f'{sizes[-2]}-{sizes[-1]}'],
           'Tc_binder_smallest_pair': pairs[f'{sizes[0]}-{sizes[1]}'],
           'Tc_chi_peak_largest_L': chi_peak(data, sizes[-1]),
           'chi_peaks': {L: chi_peak(data, L) for L in sizes}}
    return out


def mean_field_tc(J0, K):
    """Classical mean-field ordering temperature (K) for unit vectors along z: 1 = -J0_zz <e_z^2>_K / T."""
    x = np.linspace(-1, 1, 4001)
    def f(T):
        w = np.exp(-K * x ** 2 / T)
        ez2 = np.trapz(x ** 2 * w, x) / np.trapz(w, x)
        return -J0[2, 2] * ez2 / T - 1
    lo, hi = 0.1, 100.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return float(0.5 * (lo + hi) * 11.604518)
