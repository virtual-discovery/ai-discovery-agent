"""Extract the archived UppASD inputs and outputs used by the task (run after fetch.py).

    python3 -B extract.py [CACHE]     (default ~/.cache/ai-theorist/mnwo4-halfmoon-frustration)

Reads CACHE/raw/x/... and writes CACHE/extracted.npz plus CACHE/model.json:
  - model.json: the authors' Hamiltonian exactly as archived (jfile rows in mRy, kfile, posfile, momfile, cell),
    identical in every directory used (asserted, apart from the reduced-model jfiles of Fig. 5);
  - spectra: archived 1.5 K S(q,w) of paper Fig. 1b,c (101 path points; Sxx, Syy, Szz) for E < 5 meV;
  - maps_T{14,20,30,60}: archived static S(q) of paper Fig. 3 on its (k, l) grid, q = (-l/2, k, l);
  - ko_{5b..5f}: archived 20 K static S(q) of the reduced models of paper Fig. 5 on their (k, l) grid;
  - bz_{36_24_36,48_24_24}: archived 11 K static S(q) on the full commensurate grid (paper SI Figs. 7/9).
No value is altered except the cast to float64 arrays; |S| is the archived fourth column.
"""
import json
import sys
from pathlib import Path
import numpy as np

CACHE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/mnwo4-halfmoon-frustration'
X = CACHE/'raw/x'
CELL = None


def read_cell(inpsd):
    lines = inpsd.read_text().splitlines()
    i = [k for k, l in enumerate(lines) if l.startswith('cell')][0]
    rows = [lines[i].split()[1:4], lines[i + 1].split()[:3], lines[i + 2].split()[:3]]
    return np.array(rows, float)


def read_q(path):
    rows = [l.split() for l in path.read_text().splitlines() if l.strip()]
    pts = [r for r in rows if len(r) == 3]
    return np.array(pts, float)


def read_sq(path):
    return np.loadtxt(path, comments='#', skiprows=1)


def check_q(qfrac, sq, cell):
    cart = qfrac @ np.linalg.inv(cell).T
    assert np.abs(cart - sq[:, 1:4]).max() < 2e-6, 'q grid does not match the archived output'


def plane_grid(qfrac, sq):
    k = np.round(qfrac[:, 1], 6)
    l = np.round(qfrac[:, 2], 6)
    assert np.allclose(qfrac[:, 0], -qfrac[:, 2] / 2, atol=1e-6), 'map points must lie on q = (-l/2, k, l)'
    ku, lu = np.unique(k), np.unique(l)
    img = np.full((len(ku), len(lu), 4), np.nan)
    img[np.searchsorted(ku, k), np.searchsorted(lu, l)] = sq[:, 5:9]
    assert not np.isnan(img).any()
    return ku, lu, img


def main():
    out = {}
    base = X/'Fig_1_b_c/calculation'
    cell = read_cell(base/'inpsd.dat')
    model = {'cell_A': cell.tolist(),
             'jfile': [[float(x) for x in l.split()] for l in (base/'jfile').read_text().splitlines() if l.strip()],
             'kfile': [[float(x) for x in l.split()] for l in (base/'kfile').read_text().splitlines() if l.strip()],
             'posfile': [[float(x) for x in l.split()] for l in (base/'posfile').read_text().splitlines() if l.strip()],
             'momfile': [[float(x) for x in l.split()] for l in (base/'momfile').read_text().splitlines() if l.strip()],
             'units': 'UppASD: jfile and kfile in mRy; maptype 2 (integer cell offsets); posfile direct coordinates'}
    for d in ['Fig_3/T14/Base', 'Fig_3/T20/Base', 'Fig_3/T30/Base', 'Fig_3/T60/Base', 'Fig_5/Fig_5b',
              'S_Fig_7_and_Fig_9/36_24_36', 'S_Fig_7_and_Fig_9/48_24_24']:
        assert (X/d/'jfile').read_text().split() == (base/'jfile').read_text().split(), d
        assert np.allclose(read_cell(X/d/'inpsd.dat'), cell), d
    # 1.5 K S(q,w) (48x48x48, AF1 start)
    sqw = np.loadtxt(base/'sqw.loopTemp.out', comments='#', skiprows=1)
    nq = 101
    E = np.loadtxt(base/'swdos.loopTemp.out', comments='#', skiprows=1)[:, 0]
    ne = len(E)
    assert sqw.shape[0] == nq * ne
    qpath = read_q(base/'qfile')[:nq]
    check_q(qpath, sqw[::ne], cell)
    keep = E < 5.0
    out['spectra_q'] = qpath
    out['spectra_E'] = E[keep]
    out['spectra_S'] = sqw[:, 5:8].reshape(nq, ne, 3)[:, keep, :]
    out['spectra_absS'] = sqw[:, 8].reshape(nq, ne)[:, keep]
    tot = np.loadtxt(base/'totenergy.loopTemp.out', comments='#')
    out['thermo_1p5K_mRy'] = tot[-2000:, 1:4].mean(axis=0)
    # Fig. 3 maps (36x24x36)
    for T in (14, 20, 30, 60):
        d = X/f'Fig_3/T{T}'
        qf = read_q(d/'qfile')
        sq = read_sq(d/'sq.loopTemp.out')
        check_q(qf, sq, cell)
        ku, lu, img = plane_grid(qf, sq)
        out[f'map_T{T}'] = img
        out['map_k'], out['map_l'] = ku, lu
    # Fig. 5 reduced models at 20 K (36x24x36)
    jsets = {}
    for f in 'bcdef':
        d = X/f'Fig_5/Fig_5{f}'
        qf = read_q(d/'qfile')
        sq = read_sq(d/'sq.loopTemp.out')
        check_q(qf, sq, cell)
        ku, lu, img = plane_grid(qf, sq)
        out[f'ko_5{f}'] = img
        out['ko_k'], out['ko_l'] = ku, lu
        jsets[f'5{f}'] = {'jfile': [[float(x) for x in l.split()] for l in (d/'jfile').read_text().splitlines() if l.strip()],
                          'has_kfile': (d/'kfile').exists()}
    model['reduced_models_fig5'] = jsets
    # 11 K full commensurate grids
    for box in ('36_24_36', '48_24_24'):
        d = X/f'S_Fig_7_and_Fig_9/{box}'
        sq = read_sq(d/'sq.loopTemp.out')
        qd = np.loadtxt(d/'qpointsdir.out')
        assert np.abs(qd[:, 1:4] - sq[:, 1:4]).max() < 2e-6
        out[f'bz_{box}_q'] = qd[:, 4:7]
        out[f'bz_{box}_S'] = sq[:, 5:9]
    np.savez_compressed(CACHE/'extracted.npz', **out)
    (CACHE/'model.json').write_text(json.dumps(model, indent=1))
    print({k: np.shape(v) for k, v in out.items()})


if __name__ == '__main__':
    main()
