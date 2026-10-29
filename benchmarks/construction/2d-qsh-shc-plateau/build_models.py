"""Build the lab's Wannier models from the archived Wannier90 files (runs in the runtime image).

    python3 -B build_models.py SRC_ROOT OUT_DIR folder [folder ...]

For each archive folder (below Data/Wannier_Hamiltonians, e.g. insulators/Bi2_P-3m1) it reads aiida_hr.dat, the cell
(aiida.in CELL_PARAMETERS), mp_grid (aiida.win) and the final Wannier centres (cut_aiida.wout), applies Wannier90's
use_ws_distance image selection (the archived runs set use_ws_distance = .true.), writes OUT_DIR/<name>.npz
(R, deg, H complex64) and checks that the model reproduces the archived Wannier bands (aiida_band.dat).
"""
import json
import re
import sys
from pathlib import Path
import numpy as np
import convert_hr
import qsh_tb


def cell_of(folder):
    t = (folder/'aiida.in').read_text()
    return np.array([[float(x) for x in l.split()] for l in t.split('CELL_PARAMETERS')[1].split('\n')[1:4]])


def band_path(win, B):
    m = re.search(r'Begin explicit_kpath\n(.*?)End explicit_kpath\n', win, re.S)
    if m:
        return np.array([[float(x) for x in l.split()] for l in m.group(1).strip().split('\n')])
    m = re.search(r'Begin kpoint_path(.*?)End kpoint_path', win, re.S | re.I)
    segs = []
    for l in m.group(1).strip().split('\n'):
        p = l.split()
        segs.append((np.array(p[1:4], float), np.array(p[5:8], float)))
    nb = re.search(r'\n\s*bands_num_points\s*=\s*(\d+)', win)
    nb = int(nb.group(1)) if nb else 100
    L = [np.linalg.norm((b - a) @ B) for a, b in segs]
    pts = []
    for (a, b), l in zip(segs, L):
        n = max(int(round(nb * l / L[0])), 1)
        pts += [a + (b - a) * j / n for j in range(n)]
    pts.append(segs[-1][1])
    return np.array(pts)


def build(src, out, folder):
    f = src/folder
    win = (f/'aiida.win').read_text()
    A = cell_of(f)
    mp = [int(x) for x in re.search(r'mp_grid\s*=\s*(\d+)\s+(\d+)\s+(\d+)', win).groups()]
    R, deg, H = convert_hr.read_hr(f/'aiida_hr.dat')
    cen = convert_hr.wannier_centres(f/'cut_aiida.wout', H.shape[1])
    R2, deg2, H2 = convert_hr.ws_expand(R, deg, H, A, cen, mp)
    assert (deg2 == 1).all()
    # 2D: k_z = 0, so images that differ only by R_z (centres reported outside the home cell along z) are merged
    keys = sorted({(int(a), int(b)) for a, b in R2[:, :2]})
    pos = {k: i for i, k in enumerate(keys)}
    Hm = np.zeros((len(keys),) + H2.shape[1:], complex)
    for r, h in zip(R2, H2):
        Hm[pos[(int(r[0]), int(r[1]))]] += h
    R2, H2 = np.array([[a, b, 0] for a, b in keys]), Hm
    name = folder.split('/')[-1].replace('|', '_')
    np.savez_compressed(out/f'{name}.npz', R=R2.astype(np.int16), H=H2.astype(np.complex64), centres=cen)
    m = qsh_tb.Model(out/f'{name}.npz', A)
    bd = np.loadtxt(f/'aiida_band.dat')
    nk = len(bd) // m.nw
    bd = bd.reshape(m.nw, nk, 2)
    pts = band_path(win, m.B)
    err = float(abs(m.bands(pts) - bd[:, :, 1].T).max()) if len(pts) == nk else None
    return {'name': name, 'folder': folder, 'num_wann': m.nw, 'nR_archived': int(len(R)), 'nR_ws': int(len(R2)),
            'mp_grid': mp, 'band_dat_max_err_eV': err}


if __name__ == '__main__':
    src, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    for folder in sys.argv[3:]:
        print(json.dumps(build(src, out, folder)), flush=True)
