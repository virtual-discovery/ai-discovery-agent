"""Band energies of an archived model on a uniform Gamma-centred grid (exact at the grid points).

For each k1 = i/n1 the hoppings are summed along R1 with the phase exp(2 pi i k1 R1); the remaining (R2, R3)
sum is an inverse FFT on the n2 x n3 grid (folding R mod n is exact at grid points)."""
import numpy as np


def eig_grid(model, n1, n2, n3):
    nw = model.nw
    R = model.R.astype(int)
    E = np.empty((n1, n2, n3, nw))
    for i in range(n1):
        ph = np.exp(2j * np.pi * (i / n1) * R[:, 0])
        HR = np.zeros((n2, n3, nw, nw), complex)
        np.add.at(HR, (R[:, 1] % n2, R[:, 2] % n3), model.H * ph[:, None, None])
        Hk = np.fft.ifftn(HR, axes=(0, 1)) * (n2 * n3)
        E[i] = np.linalg.eigvalsh(Hk)
    return E
