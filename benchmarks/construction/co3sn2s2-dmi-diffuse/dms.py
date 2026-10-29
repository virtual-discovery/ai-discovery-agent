"""Diffuse-scattering utilities for co3sn2s2-dmi-diffuse (construction side).

The archived DMC arrays (DMCpy Cut2D) hold four channels per pixel: counts, monitor, normalization (vanadium)
and normcount; readcombined.py defines the normalized intensity counts*normcount/(monitor*normalization).
Bin edges are a regular 0.015 1/A grid in the (H0L) plane: Qx along a*, Qy along c* (1/A).
"""
import os
from pathlib import Path
import numpy as np

CACHE = Path(os.environ.get('CO3SN2S2_CACHE', str(Path.home()/'.cache/ai-theorist/co3sn2s2-dmi-diffuse')))
A_LAT, C_LAT = 5.357, 13.127                      # lattice constants used in the authors' SpinW script
ASTAR, CSTAR = 4*np.pi/(np.sqrt(3)*A_LAT), 2*np.pi/C_LAT
STEP = 0.015
NARROW = {140: 3896, 170: 3885, 175: 3879, 178: 3882, 185: 3883, 192: 3887, 200: 3884, 220: 3903, 250: 3905}
WIDE = {1.5: 3869, 165: 3871, 180: 3874, 300: 3875}   # 3871 combines 3870-3871, 3874 combines 3872-3874 (README)
EDGES_PIX = [4, 6, 8, 10, 12, 15, 17, 20, 27]      # annulus edges in pixel steps (x 0.015 1/A)
EDGES = [STEP*e for e in EDGES_PIX]
UNIT = 1e6                                          # profiles are reported in 1e-6 of the normalized intensity


def nominal(h, l):
    return (h*ASTAR, l*CSTAR)


def load(kind, run):
    folder = CACHE/'x'/('narrowscans' if kind == 'narrow' else 'widescans')
    d = np.load(folder/f'{run}_Cut2D_Data.npy')
    b = np.load(folder/f'{run}_Cut2D_BinEdges.npy')
    qx, qy = b[0][:, 0], b[1][0, :]
    assert np.allclose(b[0], qx[:, None]) and np.allclose(b[1], qy[None, :])   # separable regular grid
    assert np.allclose(np.diff(qx), STEP) and np.allclose(np.diff(qy), STEP)
    return {'counts': d[0], 'monitor': d[1], 'normalization': d[2], 'normcount': d[3],
            'qx_edges': qx, 'qy_edges': qy}


def crop(scan, center, half):
    """Pixels whose centres lie within +-half (1/A) of center in both directions."""
    cx = scan['qx_edges'][:-1] + STEP/2
    cy = scan['qy_edges'][:-1] + STEP/2
    ix = np.where(np.abs(cx - center[0]) <= half)[0]
    iy = np.where(np.abs(cy - center[1]) <= half)[0]
    out = {k: scan[k][ix[0]:ix[-1]+1, iy[0]:iy[-1]+1] for k in ('counts', 'monitor', 'normalization', 'normcount')}
    out['qx_edges'] = scan['qx_edges'][ix[0]:ix[-1]+2]
    out['qy_edges'] = scan['qy_edges'][iy[0]:iy[-1]+2]
    return out


def profile(scan, center, edges_pix=EDGES_PIX):
    """Graded annular profile (see instruction.md): pooled counts / pooled exposure, Poisson errors.

    Annulus membership uses integer pixel offsets (lo^2 <= di^2 + dj^2 < hi^2), so it is exact."""
    c, m, n, nc = (np.asarray(scan[k], float) for k in ('counts', 'monitor', 'normalization', 'normcount'))
    cx = (np.asarray(scan['qx_edges'])[:-1] + STEP/2)[:, None] * np.ones(c.shape)
    cy = (np.asarray(scan['qy_edges'])[:-1] + STEP/2)[None, :] * np.ones(c.shape)
    valid = (nc > 0) & (m > 0) & np.isfinite(n) & (n > 0)
    with np.errstate(all='ignore'):
        inten = np.where(valid, c*nc/(m*n), -np.inf)
        expo = np.where(valid, m*n/nc, 0.0)
    near = valid & (np.hypot(cx - center[0], cy - center[1]) < 0.05)
    k = np.argmax(np.where(near, inten, -np.inf))
    i, j = np.unravel_index(k, c.shape)
    ii, jj = np.indices(c.shape)
    di, dj = ii - i, jj - j
    r2 = di**2 + dj**2
    keep = valid & ~((np.abs(di) <= 3) & (np.abs(dj) <= 3))
    vals, errs, npix = [], [], []
    for lo, hi in zip(edges_pix[:-1], edges_pix[1:]):
        s = keep & (r2 >= lo*lo) & (r2 < hi*hi)
        C, X = c[s].sum(), expo[s].sum()
        vals.append(UNIT*C/X)
        errs.append(UNIT*np.sqrt(C)/X)
        npix.append(int(s.sum()))
    return np.array(vals), np.array(errs), np.array(npix), (float(cx[i, j]), float(cy[i, j]))
