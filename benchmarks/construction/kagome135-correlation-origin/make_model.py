"""Convert the archived Wannier90 hr files into compact npz models (run in the runtime image).

    bash run_in_image.sh make_model.py

Writes /data/models/<key>.npz with R (n_R,3) int, degen (n_R,), H (n_R,33,33) complex in eV (as archived,
absolute energies), cell (3,3) Angstrom, frac positions of the 3 transition-metal sites, Fermi level (eV)
from Fermi_level.dat, plus the lattice read from the .wout.
"""
import re
from pathlib import Path
import numpy as np

X = Path('/data/x')
OUT = Path('/data/models')
DIRS = {'CsTi3Sb5': ('TiSb', 'CsTiSb'), 'CsV3Sb5': ('VSb', 'CsVSb'), 'CsCr3Sb5': ('CrSb', 'CsCrSb'),
        'CsTi3Bi5': ('TiBi', 'CsTiBi'), 'CsCr3Bi5': ('CrBi', 'CsCrBi'), 'CsCr3Bi5_9.2': ('CrBi_9.2', 'CsCrBi')}


def read_hr(path):
    lines = path.read_text().splitlines()
    nw, nr = int(lines[1]), int(lines[2])
    nline = (nr + 14) // 15
    degen = np.array(' '.join(lines[3:3 + nline]).split(), int)
    data = np.loadtxt(lines[3 + nline:])
    R = data[::nw * nw, :3].astype(int)
    H = (data[:, 5] + 1j * data[:, 6]).reshape(nr, nw, nw).transpose(0, 2, 1)   # rows: m (col 4), n (col 3)
    # Wannier90 hr lines are "R1 R2 R3 m n Re Im" for H_mn(R) with m fastest: reshape gives [R][n][m]; transpose.
    return R, degen, H


def read_wout(path):
    text = path.read_text()
    cell = np.array([[float(v) for v in re.search(r'a_%d\s+(\S+)\s+(\S+)\s+(\S+)' % i, text).groups()] for i in (1, 2, 3)])
    return cell


def main():
    OUT.mkdir(exist_ok=True)
    for key, (folder, stem) in DIRS.items():
        R, degen, H = read_hr(X/folder/f'turn_{stem}_hr.dat')
        cell = read_wout(X/folder/f'turn_{stem}.wout')
        ef = float((X/folder/'Fermi_level.dat').read_text().split()[-2])
        herm = np.abs(H - H[[np.where((R == -r).all(1))[0][0] for r in R]].conj().transpose(0, 2, 1)).max()
        np.savez_compressed(OUT/f'{key}.npz', R=R, degen=degen, H=H, cell=cell, efermi=ef)
        print(key, H.shape, 'sum(1/deg)=%.1f' % (1 / degen).sum(), 'herm %.1e' % herm, 'EF', ef)


if __name__ == '__main__':
    main()
