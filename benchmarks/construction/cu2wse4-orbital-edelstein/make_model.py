"""Convert the two archived OpenMX Wannier Hamiltonians (*.HWR) to npz.

    bash run_in_image.sh make_model.py

Reads /data/x/Cu2WSe4_{wosoc,soc}.HWR (fetch.py) and writes /data/model/{wosoc,soc}.npz with
R (integer lattice vectors), deg (Wigner-Seitz degeneracy), H (complex, Hartree, as in the file: H[r, m, n] is the
line 'm n re im' of block R), lattice (Bohr), fermi_Ha and the header lines.
"""
import re
from pathlib import Path
import numpy as np

SRC, OUT = Path('/data/x'), Path('/data/model')


def parse(path):
    lines = path.read_text().splitlines()
    nw = int(lines[1].split()[-1]); nr = int(lines[2].split()[-1])
    lat = np.array([[float(x) for x in lines[4 + i].split()] for i in range(3)])
    header = lines[:9]
    fermi = float(lines[8].split()[-1])
    R, deg, H = np.zeros((nr, 3), int), np.zeros(nr, int), np.zeros((nr, nw, nw), complex)
    i = 9
    for r in range(nr):
        m = re.match(r'R\s*\(\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\s*\)\s+(\d+)', lines[i])
        R[r] = [int(m.group(k)) for k in (1, 2, 3)]; deg[r] = int(m.group(4))
        block = np.array(' '.join(lines[i + 1:i + 1 + nw * nw]).split(), float).reshape(nw * nw, 4)
        assert (block[:, 0].astype(int) == np.repeat(np.arange(1, nw + 1), nw)).all()
        assert (block[:, 1].astype(int) == np.tile(np.arange(1, nw + 1), nw)).all()
        H[r] = (block[:, 2] + 1j * block[:, 3]).reshape(nw, nw)
        i += 1 + nw * nw
    return dict(R=R, deg=deg, H=H, lattice_bohr=lat, fermi_Ha=np.array(fermi), header=np.array(header))


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    for name in ('wosoc', 'soc'):
        d = parse(SRC/f'Cu2WSe4_{name}.HWR')
        np.savez(OUT/f'{name}.npz', **d)
        R, H, deg = d['R'], d['H'], d['deg']
        idx = {tuple(r): k for k, r in enumerate(R)}
        herm = max(np.abs(H[idx[tuple(-r)]] - H[k].conj().T).max() for k, r in enumerate(R))
        print(name, 'nw', H.shape[1], 'nr', len(R), 'sum 1/deg', (1 / deg).sum(), 'fermi Ha', d['fermi_Ha'],
              'max |H(-R) - H(R)^+|', herm, 'max|H|', np.abs(H).max(), 'R range', R.min(0), R.max(0))
