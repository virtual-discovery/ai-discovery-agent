"""Profile-prediction recipes used for calibration, robustness variants and ablations (construction side).

The reference recipe is the one in solution/solve.py:
  * every purchased (003) profile (narrow scans and wide-scan regions, T >= 140 K) is fitted jointly with
    y_T(q) = b0 + b1 (T - 200) + s q + A_T / (1 + (q / kappa_T)^2);
  * a withheld temperature gets A and kappa interpolated linearly between the nearest purchased temperatures
    on either side, and the background from the linear law;
  * the (101) profile at 180 K is the (101) background (linear-in-q fits at 165 K and 300 K, interpolated in T)
    plus the (003) critical term at 180 K scaled by |F(101)|^2/|F(003)|^2 = 1/9 (ferromagnetic correlations,
    isotropic spin components: the polarisation sums are equal at the two reflections).
"""
import numpy as np
from scipy.optimize import least_squares

EDGES = 0.015*np.array([4, 6, 8, 10, 12, 15, 17, 20, 27])
QC = 0.5*(EDGES[1:] + EDGES[:-1])


def lorentz(q, A, k):
    return A/(1 + (q/k)**2)


def global_fit(profiles, fix_b1=None):
    """profiles: {T: (values, errors)} in units of 1e-6."""
    Ts = sorted(profiles)
    n = len(Ts)

    def unpack(p):
        b0, b1, s = p[0], (p[1] if fix_b1 is None else fix_b1), p[2]
        return b0, b1, s, p[3:3+n], p[3+n:]

    def resid(p):
        b0, b1, s, A, k = unpack(p)
        return np.concatenate([(profiles[T][0] - (b0 + b1*(T - 200) + s*QC + lorentz(QC, A[i], k[i])))/profiles[T][1]
                               for i, T in enumerate(Ts)])
    p0 = np.r_[8.2, 0.002, 0.0, np.full(n, 0.5), np.full(n, 0.1)]
    lb = np.r_[5, -1, -20, np.zeros(n), np.full(n, 0.02)]
    ub = np.r_[12, 1, 20, np.full(n, 30), np.full(n, 0.5)]
    r = least_squares(resid, p0, bounds=(lb, ub), x_scale='jac')
    b0, b1, s, A, k = unpack(r.x)
    return {'Ts': Ts, 'b0': b0, 'b1': b1, 's': s, 'A': dict(zip(Ts, A)), 'k': dict(zip(Ts, k)),
            'chi2r': float(np.sum(r.fun**2)/(len(r.fun) - len(r.x)))}


def bracket(Ts, T):
    lo = max(t for t in Ts if t < T)
    hi = min(t for t in Ts if t > T)
    return lo, hi, (T - lo)/(hi - lo)


def critical(g, T):
    if T in g['A']:
        return lorentz(QC, g['A'][T], g['k'][T])
    lo, hi, w = bracket(g['Ts'], T)
    return lorentz(QC, (1 - w)*g['A'][lo] + w*g['A'][hi], (1 - w)*g['k'][lo] + w*g['k'][hi])


def background(g, T):
    return g['b0'] + g['b1']*(T - 200) + g['s']*QC


def predict(g, T):
    return background(g, T) + critical(g, T)


def predict_mirror(g, T, tc=177.5):
    """Wrong path: critical scattering symmetric in |T - Tc|, background at the true temperature."""
    return background(g, T) + critical(g, 2*tc - T)


def linear_q(values, errors):
    c = np.polyfit(QC, values, 1, w=1/errors)
    return np.polyval(c, QC)


def predict_101(p101, crit180, ratio=1/9, smooth=True, temps=(165, 300)):
    """p101: {T: (values, errors)} at (101); crit180: (003) critical term at 180 K."""
    lo, hi = temps
    w = (180 - lo)/(hi - lo)
    f = linear_q if smooth else (lambda v, e: v)
    bg = (1 - w)*f(*p101[lo]) + w*f(*p101[hi])
    return bg + ratio*crit180


def interp_bins(profiles, T):
    """Robustness variant: per-annulus linear interpolation between the nearest purchased temperatures."""
    lo, hi, w = bracket(sorted(profiles), T)
    return (1 - w)*profiles[lo][0] + w*profiles[hi][0]


def chi2r(pred, truth, err):
    return float(np.mean(((np.asarray(pred) - truth)/err)**2))


def independent_fits(profiles):
    """Robustness recipe: fit each temperature on its own (constant + slope + Lorentzian)."""
    par = {}
    for T, (y, e) in profiles.items():
        r = least_squares(lambda p: (y - (p[0] + p[1]*QC + lorentz(QC, p[2], p[3])))/e, [8.2, 0, 0.5, 0.1],
                          bounds=([5, -20, 0, 0.02], [12, 20, 30, 0.5]))
        par[T] = r.x
    return par


def predict_independent(par, T):
    lo, hi, w = bracket(sorted(par), T)
    p = (1 - w)*par[lo] + w*par[hi]
    return p[0] + p[1]*QC + lorentz(QC, p[2], p[3])
