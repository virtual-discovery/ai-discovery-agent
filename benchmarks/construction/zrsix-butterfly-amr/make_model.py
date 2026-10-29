"""Convert the three archived wannier90_hr files to compact npz models in /data/model.

    bash run_in_image.sh make_model.py

Each model stores the integer lattice vectors R, the Wigner-Seitz degeneracies and H(R) in eV (complex128),
plus the lattice and orbital positions read from the matching WannierTools input.dat / POSCAR.
"""
import numpy as np
from pathlib import Path

SRC = Path('/data/ZrSiX data/input_MR')
OUT = Path('/data/model')
FILES = {'ZrSiS': 'ZrSiS/wannier90_hr.dat_nsymm16', 'ZrSiSe': 'ZrSiSe/wannier90_hr.dat_nsymm16',
         'ZrSiTe': 'ZrSiTe/wannier90_hr.dat'}


def read_hr(path):
    with open(path) as f:
        f.readline()
        nw = int(f.readline()); nr = int(f.readline())
        deg = []
        while len(deg) < nr:
            deg += [int(x) for x in f.readline().split()]
        data = np.loadtxt(f)
    data = data.reshape(nr, nw * nw, 7)
    R = data[:, 0, :3].astype(int)
    i = data[0, :, 3].astype(int) - 1; j = data[0, :, 4].astype(int) - 1
    H = np.zeros((nr, nw, nw), complex)
    H[:, i, j] = data[:, :, 5] + 1j * data[:, :, 6]
    return R, np.array(deg), H


def main():
    OUT.mkdir(exist_ok=True)
    for name, rel in FILES.items():
        R, deg, H = read_hr(SRC/rel)
        herm = max(np.abs(H[k] - H[np.where((R == -R[k]).all(1))[0][0]].conj().T).max() for k in range(len(R)))
        np.savez_compressed(OUT/f'{name}.npz', R=R, deg=deg, H=H)
        print(name, 'nR', len(R), 'deg sum', deg.sum(), 'unique deg', np.unique(deg), 'hermiticity err', herm,
              'onsite diag', np.round(np.diag(H[(R == 0).all(1)][0]).real, 3))


if __name__ == '__main__':
    main()
