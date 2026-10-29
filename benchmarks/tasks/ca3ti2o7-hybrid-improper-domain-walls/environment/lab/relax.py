"""Structure relaxation with a NEP model: L-BFGS on positions, optionally the cell length along one axis."""
import numpy as np
from scipy.optimize import minimize


class Relaxer:
    def __init__(self, model, types, pos, cell, relax_axis=None, skin=1.0):
        self.m, self.types = model, np.asarray(types, int)
        self.pos0, self.cell0 = np.asarray(pos, float).copy(), np.asarray(cell, float).copy()
        self.axis = relax_axis
        self.skin = skin
        self.calls = 0
        self.pairs, self.ref = None, None
        if relax_axis is not None:
            a, b = self.cell0[(relax_axis + 1) % 3], self.cell0[(relax_axis + 2) % 3]
            n = np.cross(a, b)
            self.n = n / np.linalg.norm(n)
            self.L = abs(np.dot(self.cell0[relax_axis], self.n))

    def geometry(self, x):
        N = len(self.types)
        pos = x[:3 * N].reshape(N, 3)
        cell = self.cell0
        eps = 0.0
        if self.axis is not None:
            eps = x[-1] / self.L
            F = np.eye(3) + eps * np.outer(self.n, self.n)
            pos = pos @ F.T
            cell = self.cell0 @ F.T
        return pos, cell, eps

    def evaluate(self, x, want_bec=False):
        pos, cell, eps = self.geometry(x)
        if self.pairs is None or np.abs(pos - self.ref[0]).max() > 0.4 * self.skin or abs(eps - self.ref[1]) > 0.02:
            self.pairs = self.m.neighbours(pos, cell, self.skin)
            self.ref = (pos.copy(), eps)
        self.calls += 1
        return self.m.compute(self.types, pos, cell, want_bec=want_bec, pairs=self.pairs), pos, cell, eps

    def fun(self, x):
        o, pos, cell, eps = self.evaluate(x)
        N = len(self.types)
        g = -o['forces']
        if self.axis is not None:
            F = np.eye(3) + eps * np.outer(self.n, self.n)
            g = g @ F                     # gradient wrt reference coordinates
            dEdeps = -self.n @ o['virial'] @ self.n / (1 + eps)
            return o['energy'], np.r_[g.ravel(), dEdeps / self.L]
        return o['energy'], g.ravel()

    def run(self, fmax=0.02, max_steps=500):
        x0 = self.pos0.ravel()
        if self.axis is not None:
            x0 = np.r_[x0, 0.0]
        res = minimize(self.fun, x0, jac=True, method='L-BFGS-B',
                       options={'maxiter': max_steps, 'maxfun': int(max_steps * 1.5), 'gtol': fmax,
                                'ftol': 1e-14, 'maxcor': 30})
        o, pos, cell, eps = self.evaluate(res.x)
        f = np.abs(o['forces']).max() if len(o['forces']) else 0.0
        fnorm = np.sqrt((o['forces'] ** 2).sum(1)).max()
        return {'energy': o['energy'], 'positions': pos, 'cell': cell, 'strain': eps,
                'max_force': float(fnorm), 'virial': o['virial'], 'converged': bool(fnorm <= fmax * 1.75),
                'iterations': int(res.nit), 'evaluations': self.calls, 'message': str(res.message)}
