"""Readers for archived Wannier90 files and the use_ws_distance image selection (used by build_models.py)."""
import numpy as np


def read_hr(path):
    with open(path) as f:
        f.readline()
        nw = int(f.readline())
        nr = int(f.readline())
        deg = []
        while len(deg) < nr:
            deg += [int(x) for x in f.readline().split()]
        data = np.loadtxt(f)
    data = data.reshape(nr, nw * nw, 7)
    R = data[:, 0, :3].astype(int)
    m = data[0, :, 3].astype(int) - 1
    n = data[0, :, 4].astype(int) - 1
    H = np.zeros((nr, nw, nw), complex)
    H[:, m, n] = data[:, :, 5] + 1j * data[:, :, 6]
    return R, np.array(deg), H


def wannier_centres(wout_path, nw):
    """Final-state WF centres (Angstrom) from the archived (head/tail) .wout."""
    lines = open(wout_path).read().split('\n')
    start = max(i for i, l in enumerate(lines) if 'Final State' in l)
    cen = []
    for l in lines[start + 1:start + 1 + nw]:
        inside = l.split('(')[1].split(')')[0]
        cen.append([float(x) for x in inside.split(',')])
    return np.array(cen)


def ws_expand(R, deg, H, cell_A, centres, mp_grid, search=2, tol=1e-5):
    """Wannier90 use_ws_distance: for every (R, i, j) choose the supercell images R + T (T in mp_grid supercell
    lattice) that minimise |R + T + tau_j - tau_i|, splitting degenerate images evenly (ws_distance.F90).
    Returns an equivalent hr on an enlarged R set with deg = 1."""
    A = np.array(cell_A)
    nw = H.shape[1]
    sc = np.array(mp_grid)
    shifts = np.array([(i, j, k) for i in range(-search, search + 1) for j in range(-search, search + 1)
                       for k in range(-search, search + 1)]) * sc
    out = {}
    d0 = centres[None, :, :] - centres[:, None, :]          # tau_j - tau_i  [i, j, 3]
    for r, (Rv, dg) in enumerate(zip(R, deg)):
        vec = (Rv[None, :] + shifts) @ A                     # [s, 3]
        dist = np.linalg.norm(vec[:, None, None, :] + d0[None], axis=-1)   # [s, i, j]
        dmin = dist.min(0)
        sel = np.abs(dist - dmin[None]) < tol               # degenerate images
        nd = sel.sum(0)
        for s in np.where(sel.any((1, 2)))[0]:
            key = tuple(Rv + shifts[s])
            w = np.where(sel[s], H[r] / dg / nd, 0)
            out[key] = out.get(key, 0) + w
    keys = sorted(out)
    return np.array(keys), np.ones(len(keys), int), np.array([out[k] for k in keys])
