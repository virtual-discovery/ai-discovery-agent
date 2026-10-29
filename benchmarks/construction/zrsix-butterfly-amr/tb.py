"""Tight-binding utilities for the archived ZrSiX Wannier models (construction and lab)."""
import numpy as np

HBAR_EV_S = 6.582119569e-16
LATT = {'ZrSiS': (3.5496073769164020, 8.1636955490019929),   # POSCAR / wannier90.win of each compound
        'ZrSiSe': (3.623, 8.365),
        'ZrSiTe': (3.6974, 9.508)}


class Model:
    def __init__(self, path, name):
        z = np.load(path)
        self.R = z['R'].astype(float)
        self.H = z['H'] / z['deg'][:, None, None]
        self.nw = self.H.shape[1]
        a, c = LATT[name]
        self.A = np.diag([a, a, c])
        self.B = 2 * np.pi * np.linalg.inv(self.A).T
        self.Rc = self.R @ self.A
        self.Hf = self.H.reshape(len(self.R), -1)

    def hk(self, kfrac):
        ph = np.exp(2j * np.pi * kfrac @ self.R.T)
        return (ph @ self.Hf).reshape(-1, self.nw, self.nw)

    def bands(self, kfrac, chunk=4000):
        out = []
        for s in range(0, len(kfrac), chunk):
            out.append(np.linalg.eigvalsh(self.hk(kfrac[s:s + chunk])))
        return np.concatenate(out)

    def bands_vel(self, kfrac, bands=None, chunk=2000):
        """Energies (eV) and Hellmann-Feynman velocities dE/dk (eV*Angstrom) of selected bands."""
        E, V = [], []
        for s in range(0, len(kfrac), chunk):
            k = kfrac[s:s + chunk]
            ph = np.exp(2j * np.pi * k @ self.R.T)
            h = (ph @ self.Hf).reshape(-1, self.nw, self.nw)
            e, u = np.linalg.eigh(h)
            if bands is not None:
                e, u = e[:, bands], u[:, :, bands]
            v = []
            for d in range(3):
                dh = ((ph * (1j * self.Rc[:, d])) @ self.Hf).reshape(-1, self.nw, self.nw)
                v.append(np.einsum('kin,kij,kjn->kn', u.conj(), dh, u).real)
            E.append(e); V.append(np.stack(v, -1))
        return np.concatenate(E), np.concatenate(V)
