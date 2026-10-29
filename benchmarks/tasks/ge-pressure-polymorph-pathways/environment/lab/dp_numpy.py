"""Pure-numpy evaluator of the Fantasia et al. (2024) DeePMD-kit Ge model (Fantasy-frz.pb).

Model (from the frozen graph's training script): hybrid descriptor se_e3 (rcut 3.3 A, rcut_smth 2.0, sel 16,
embedding [8, 16, 32]) + se_e2_a (rcut 6.6 A, rcut_smth 3.3, sel 96, embedding [16, 32, 64], axis 32), both with
zero mean; fitting net [160, 120, 80] tanh; atom_ener = -4.52 eV, so E_i = NN(D_i) - NN(0) - 4.52.
Energies, analytic forces and stress (virial) are returned. Weights are read from potential/dp_weights.npz and
potential/dp_meta.json (extracted from the frozen graph by a minimal protobuf reader; no TensorFlow).
"""
import json
from pathlib import Path
import numpy as np


class DP:
    def __init__(self, folder):
        folder = Path(folder)
        with np.load(folder/'dp_weights.npz') as z:
            self.W = {k: np.array(z[k], dtype=float) for k in z.files}
        m = json.loads((folder/'dp_meta.json').read_text())
        self.d0, self.d1, self.bias = m['d0'], m['d1'], m['bias_atom_e']
        self.std0 = self.W['descrpt_attr_0/t_std'].reshape(-1, 4)[0]
        self.std1 = self.W['descrpt_attr_1/t_std'].reshape(-1, 4)[0]
        self.fit0 = self._fit_forward(np.zeros((1, self.W['layer_0_type_0/matrix'].shape[0])))[0][0]

    # ---------------------------------------------------------------- neighbour lists
    @staticmethod
    def neighbours(cell, pos, rcut, sel):
        """Per atom: indices j and displacement vectors d = r_j + shift - r_i (r < rcut), sorted, <= sel."""
        cell = np.asarray(cell, float)
        pos = np.asarray(pos, float)
        rec = np.linalg.inv(cell).T
        nrep = np.ceil(rcut * np.linalg.norm(rec, axis=1)).astype(int)
        grid = np.array([(a, b, c) for a in range(-nrep[0], nrep[0] + 1) for b in range(-nrep[1], nrep[1] + 1)
                         for c in range(-nrep[2], nrep[2] + 1)])
        shifts = grid @ cell
        n = len(pos)
        J = np.full((n, sel), -1, dtype=int)
        D = np.zeros((n, sel, 3))
        for i in range(n):
            d = (pos[None, :, :] + shifts[:, None, :] - pos[i]).reshape(-1, 3)
            r = np.linalg.norm(d, axis=1)
            idx = np.where((r > 1e-8) & (r < rcut))[0]
            idx = idx[np.argsort(r[idx], kind='stable')][:sel]
            J[i, :len(idx)] = idx % n
            D[i, :len(idx)] = d[idx]
        return J, D

    @staticmethod
    def _sw(r, rmin, rmax):
        u = np.clip((r - rmin) / (rmax - rmin), 0.0, 1.0)
        s = np.where(r < rmin, 1.0, np.where(r < rmax, u ** 3 * (-6 * u * u + 15 * u - 10) + 1, 0.0))
        ds = np.where((r >= rmin) & (r < rmax), (-30 * u ** 4 + 60 * u ** 3 - 30 * u ** 2) / (rmax - rmin), 0.0)
        return s, ds

    def _env(self, D, rcut, rsmth, std):
        """Normalized environment matrix env (n, sel, 4) and its derivative wrt D: denv (n, sel, 4, 3)."""
        r = np.linalg.norm(D, axis=2)
        mask = r > 0
        rs = np.where(mask, r, 1.0)
        s, ds = self._sw(rs, rsmth, rcut)
        q = s / rs                          # env0 = q
        dq = (ds * rs - s) / rs ** 2
        h = q / rs                          # env_a = h * d_a
        dh = (dq * rs - q) / rs ** 2
        env = np.zeros(D.shape[:2] + (4,))
        env[..., 0] = q
        env[..., 1:] = h[..., None] * D
        u = D / rs[..., None]
        denv = np.zeros(D.shape[:2] + (4, 3))
        denv[..., 0, :] = dq[..., None] * u
        denv[..., 1:, :] = h[..., None, None] * np.eye(3) + dh[..., None, None] * D[..., :, None] * u[..., None, :]
        env *= mask[..., None]
        denv *= mask[..., None, None]
        return env / std, denv / std[None, None, :, None]

    # ---------------------------------------------------------------- networks
    def _embed(self, x, prefix, suffix, sizes):
        """Embedding net on scalar inputs x (...,): returns g (..., M) and dg/dx (..., M)."""
        W = self.W
        h = None
        for k in range(len(sizes)):
            Wk, bk = W[f'{prefix}/matrix_{k + 1}{suffix}'], W[f'{prefix}/bias_{k + 1}{suffix}']
            if k == 0:
                y = np.tanh(x[..., None] * Wk[0] + bk)
                dy = (1 - y * y) * Wk[0]
                h, dh = y, dy
                continue
            y = np.tanh(h @ Wk + bk)
            dy = (1 - y * y) * (dh @ Wk)
            if sizes[k] == sizes[k - 1]:
                h, dh = h + y, dh + dy
            elif sizes[k] == 2 * sizes[k - 1]:
                h, dh = np.concatenate([h, h], -1) + y, np.concatenate([dh, dh], -1) + dy
            else:
                h, dh = y, dy
        return h, dh

    def _fit_forward(self, x):
        W = self.W
        acts = [x]
        for k in range(3):
            x = np.tanh(x @ W[f'layer_{k}_type_0/matrix'] + W[f'layer_{k}_type_0/bias'])
            acts.append(x)
        return (x @ W['final_layer_type_0/matrix'] + W['final_layer_type_0/bias'])[:, 0], acts

    def _fit_backward(self, acts):
        W = self.W
        g = np.repeat(W['final_layer_type_0/matrix'][:, 0][None, :], acts[0].shape[0], axis=0)
        for k in (2, 1, 0):
            y = acts[k + 1]
            g = (g * (1 - y * y)) @ W[f'layer_{k}_type_0/matrix'].T
        return g

    # ---------------------------------------------------------------- main
    def compute(self, cell, pos, forces=True):
        cell = np.asarray(cell, float)
        pos = np.asarray(pos, float)
        n = len(pos)
        d0, d1 = self.d0, self.d1
        J0, D0v = self.neighbours(cell, pos, d0['rcut'], d0['sel'])
        J1, D1v = self.neighbours(cell, pos, d1['rcut'], d1['sel'])
        e0, de0 = self._env(D0v, d0['rcut'], d0['rcut_smth'], self.std0)
        e1, de1 = self._env(D1v, d1['rcut'], d1['rcut_smth'], self.std1)
        # se_e3
        v = e0[..., 1:]
        eij = np.einsum('ijm,ikm->ijk', v, v)
        g3, dg3 = self._embed(eij, 'filter_type_all_0', '_0_0', d0['neuron'])
        s0 = 1.0 / d0['sel'] ** 2
        Dsc3 = np.einsum('ijk,ijkm->im', eij, g3) * s0
        # se_e2_a
        G, dG = self._embed(e1[..., 0], 'filter_type_0_1', '_0', d1['neuron'])
        T = np.einsum('ijc,ijm->icm', e1, G) / d1['sel']
        ax = d1['axis']
        Dsa = np.einsum('icm,icl->iml', T, T[:, :, :ax])
        desc = np.concatenate([Dsc3, Dsa.reshape(n, -1)], axis=1)
        out, acts = self._fit_forward(desc)
        atomic = out - self.fit0 + self.bias
        E = float(atomic.sum())
        if not forces:
            return {'energy': E, 'atomic': atomic}
        gd = self._fit_backward(acts)                      # dE/ddesc (n, 2080)
        m3 = Dsc3.shape[1]
        gD3, gDa = gd[:, :m3], gd[:, m3:].reshape(Dsa.shape)
        # se_e3 backward
        de_jk = (np.einsum('im,ijkm->ijk', gD3, g3) + eij * np.einsum('im,ijkm->ijk', gD3, dg3)) * s0
        gv = np.einsum('ijk,ika->ija', de_jk + de_jk.transpose(0, 2, 1), v)
        genv0 = np.zeros_like(e0)
        genv0[..., 1:] = gv
        # se_e2_a backward
        gT = np.einsum('icl,iml->icm', T[:, :, :ax], gDa)
        gT[:, :, :ax] += np.einsum('icm,iml->icl', T, gDa)
        genv1 = np.einsum('icm,ijm->ijc', gT, G) / d1['sel']
        gG = np.einsum('ijc,icm->ijm', e1, gT) / d1['sel']
        genv1[..., 0] += np.einsum('ijm,ijm->ij', gG, dG)
        # chain to displacement vectors
        gd0 = np.einsum('ijc,ijcb->ijb', genv0, de0)
        gd1 = np.einsum('ijc,ijcb->ijb', genv1, de1)
        F = np.zeros((n, 3))
        vir = np.zeros((3, 3))
        for J, Dv, gdv in ((J0, D0v, gd0), (J1, D1v, gd1)):
            valid = J >= 0
            # dE/dr_j += gdv ; dE/dr_i -= gdv
            np.add.at(F, J[valid], -gdv[valid])
            F += gdv.sum(axis=1)
            vir += np.einsum('ija,ijb->ab', Dv, gdv)       # dE/d(strain) = sum d_a dE/dd_b
        vol = abs(np.linalg.det(cell))
        stress = vir / vol                                    # dE/deps / V  (eV/A^3), tensile positive
        return {'energy': E, 'atomic': atomic, 'forces': F, 'dE_dstrain': vir,
                'stress_eV_A3': 0.5 * (stress + stress.T), 'pressure_GPa': float(-np.trace(stress) / 3 * 160.21766208)}

    def energy(self, cell, pos):
        return self.compute(cell, pos, forces=False)['energy']
