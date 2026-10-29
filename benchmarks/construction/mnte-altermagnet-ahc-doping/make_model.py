"""Convert the archived wannier90_hr.dat to npz (R, deg, H as complex128) for the evaluator.

    bash run_in_image.sh make_model.py

Reads /data/dft_data/2_wannier90/wannier90_hr.dat (fetch.py) and writes /data/model/archived.npz.
The conversion is exact (the hr values are read as printed, 6 decimals).
"""
from pathlib import Path
import numpy as np

SRC = Path('/data/dft_data/2_wannier90/wannier90_hr.dat')
OUT = Path('/data/model/archived.npz')


def read_hr(path):
    with open(path) as f:
        f.readline()
        nw = int(f.readline()); nr = int(f.readline())
        deg = []
        while len(deg) < nr:
            deg += [int(x) for x in f.readline().split()]
        data = np.loadtxt(f)
    R = data[:, :3].astype(int).reshape(nr, nw * nw, 3)[:, 0, :]
    m = data[:, 3].astype(int).reshape(nr, nw, nw)
    n = data[:, 4].astype(int).reshape(nr, nw, nw)
    assert (m[:, 0, :] == np.arange(1, nw + 1)).all() and (n[:, :, 0] == np.arange(1, nw + 1)).all()
    H = np.transpose((data[:, 5] + 1j * data[:, 6]).reshape(nr, nw, nw), (0, 2, 1))   # H[R, m, n] = <m0|H|nR>
    return R, np.array(deg), H


if __name__ == '__main__':
    R, deg, H = read_hr(SRC)
    index = {tuple(r): i for i, r in enumerate(R)}
    assert max(np.abs(H[index[tuple(-r)]] - H[i].conj().T).max() for i, r in enumerate(R)) == 0
    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez(OUT, R=R, deg=deg, H=H)
    print(OUT, len(R), 'R vectors', H.shape[1], 'Wannier functions; sum 1/deg =', (1 / deg).sum())
