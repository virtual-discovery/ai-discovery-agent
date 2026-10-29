"""Relaxations with the numpy DeePMD evaluator: fixed-cell ionic relaxation and full (cell + ions) enthalpy
minimisation at a target pressure. scipy L-BFGS-B with analytic gradients."""
import numpy as np
from scipy.optimize import minimize

EV_A3_PER_GPA = 1 / 160.21766208


def relax_ions(dp, cell, pos, fmax=1e-3, maxiter=500):
    cell = np.asarray(cell, float)
    pos = np.asarray(pos, float)
    n = len(pos)

    def f(x):
        r = dp.compute(cell, x.reshape(n, 3))
        return r['energy'], -r['forces'].ravel()
    res = minimize(f, pos.ravel(), jac=True, method='L-BFGS-B',
                   options={'gtol': fmax, 'maxiter': maxiter, 'ftol': 1e-14, 'maxcor': 30})
    p = res.x.reshape(n, 3)
    r = dp.compute(cell, p)
    return {'positions': p, 'energy': r['energy'], 'max_force': float(np.abs(r['forces']).max()),
            'stress_eV_A3': r['stress_eV_A3'], 'pressure_GPa': r['pressure_GPa'], 'iterations': int(res.nit)}


def relax_full(dp, cell, pos, pressure_GPa, fmax=1e-3, maxiter=1000, hydrostatic_only=False):
    """Minimise H = E + P V over fractional coordinates and a symmetric cell strain (lower-triangular free)."""
    cell0 = np.asarray(cell, float)
    frac0 = np.asarray(pos, float) @ np.linalg.inv(cell0)
    n = len(frac0)
    P = pressure_GPa * EV_A3_PER_GPA
    scale = 10.0   # strain variables are scaled for conditioning

    def unpack(x):
        e = x[3 * n:] / scale
        eps = np.array([[e[0], e[3], e[4]], [e[3], e[1], e[5]], [e[4], e[5], e[2]]])
        if hydrostatic_only:
            eps = np.eye(3) * e[0]
        M = np.eye(3) + eps
        cell = cell0 @ M.T
        return cell, x[:3 * n].reshape(n, 3) @ cell, M

    def f(x):
        cell, p, M = unpack(x)
        r = dp.compute(cell, p)
        V = abs(np.linalg.det(cell))
        H = r['energy'] + P * V
        # gradient wrt fractional coordinates: dE/ds = -F @ cell^T... positions p = s @ cell
        gs = (-r['forces']) @ cell.T
        # gradient wrt strain (relative to current cell, mapped back through M)
        sig = r['dE_dstrain'] + P * V * np.eye(3)       # dH/d(eps_current)
        gE = sig @ np.linalg.inv(M).T                    # chain to eps relative to cell0 (small-strain form)
        gE = 0.5 * (gE + gE.T)
        if hydrostatic_only:
            ge = np.array([np.trace(gE), 0, 0, 0, 0, 0])
        else:
            ge = np.array([gE[0, 0], gE[1, 1], gE[2, 2], 2 * gE[0, 1], 2 * gE[0, 2], 2 * gE[1, 2]])
        return H, np.concatenate([gs.ravel(), ge / scale])
    x0 = np.concatenate([frac0.ravel(), np.zeros(6)])
    res = minimize(f, x0, jac=True, method='L-BFGS-B',
                   options={'gtol': fmax, 'maxiter': maxiter, 'ftol': 1e-15, 'maxcor': 50})
    cell, p, M = unpack(res.x)
    r = dp.compute(cell, p)
    V = abs(np.linalg.det(cell))
    return {'cell': cell, 'positions': p, 'energy': r['energy'], 'volume': V, 'enthalpy': r['energy'] + P * V,
            'pressure_GPa': r['pressure_GPa'], 'stress_eV_A3': r['stress_eV_A3'],
            'max_force': float(np.abs(r['forces']).max()), 'iterations': int(res.nit)}
