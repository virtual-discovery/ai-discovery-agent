"""Structure relaxation with the numpy ACE evaluator (scipy L-BFGS-B on positions and optional cell strain).

mode: 'positions' (fixed cell), 'volume' (positions + isotropic cell scaling, zero pressure),
      'cell' (positions + full symmetric cell strain, zero stress).
fixed: indices of atoms held fixed (positions mode and with cell modes in fractional coordinates).
"""
import numpy as np
from scipy.optimize import minimize


def relax(ace, cell, symbols, positions, mode='positions', fixed=(), fmax=1e-3, max_steps=400):
    cell0 = np.asarray(cell, float)
    pos0 = np.asarray(positions, float)
    N = len(pos0)
    free = np.ones(N, bool)
    free[list(fixed)] = False
    frac0 = pos0 @ np.linalg.inv(cell0)
    vol0 = abs(np.linalg.det(cell0))
    scale = vol0 ** (1 / 3)
    nstrain = {'positions': 0, 'volume': 1, 'cell': 6}[mode]
    calls = {'n': 0}
    last = {}

    def unpack(x):
        u = x[:3 * free.sum()].reshape(-1, 3)
        e = x[3 * free.sum():] / scale          # strain variables are carried in length units
        if mode == 'volume':
            F = np.eye(3) * (1 + e[0])
        elif mode == 'cell':
            F = np.eye(3) + np.array([[e[0], e[5] / 2, e[4] / 2], [e[5] / 2, e[1], e[3] / 2], [e[4] / 2, e[3] / 2, e[2]]])
        else:
            F = np.eye(3)
        cell = cell0 @ F.T
        frac = frac0.copy()
        frac[free] += u @ np.linalg.inv(cell0)    # displacements in reference Cartesian coordinates
        return cell, frac @ cell, F

    def fun(x):
        cell, pos, F = unpack(x)
        out = ace.compute(cell, symbols, pos)
        calls['n'] += 1
        last.update(out=out, cell=cell, pos=pos)
        # gradient w.r.t. reference displacements: dE/du = -F_force @ F (chain rule through cell deformation)
        g_u = -(out['forces'][free] @ F)
        g = [g_u.ravel()]
        if nstrain:
            vol = abs(np.linalg.det(cell))
            dEdF = vol * out['stress'] @ np.linalg.inv(F).T
            if mode == 'volume':
                g.append([np.trace(dEdF)])
            else:
                g.append([dEdF[0, 0], dEdF[1, 1], dEdF[2, 2], (dEdF[1, 2] + dEdF[2, 1]) / 2,
                          (dEdF[0, 2] + dEdF[2, 0]) / 2, (dEdF[0, 1] + dEdF[1, 0]) / 2])
        grad = np.concatenate([np.ravel(a) for a in g])
        grad[3 * free.sum():] /= scale
        return out['energy'], grad

    x0 = np.zeros(3 * free.sum() + nstrain)
    bounds = [(None, None)] * (3 * free.sum()) + [(-0.15 * scale, 0.15 * scale)] * nstrain
    res = minimize(fun, x0, jac=True, method='L-BFGS-B', bounds=bounds,
                   options={'maxiter': max_steps, 'gtol': fmax * 0.1, 'ftol': 1e-14, 'maxcor': 30})
    cell, pos, F = unpack(res.x)
    out = ace.compute(cell, symbols, pos)
    f = out['forces'].copy()
    f[~free] = 0
    return {'energy': out['energy'], 'cell': cell, 'positions': pos, 'forces': out['forces'], 'stress': out['stress'],
            'max_force': float(np.abs(f).max()) if N else 0.0, 'steps': int(res.nit), 'evaluations': calls['n'],
            'converged': bool(np.abs(f).max() < fmax * 5)}
