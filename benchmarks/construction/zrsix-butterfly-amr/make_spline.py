"""Exact band-group eigenvalues on a 128 x 128 x 64 grid -> spline coefficient files /data/spline/<name>.npz."""
import numpy as np, sys, time
from pathlib import Path
sys.path.insert(0, '.')
import tb, grid, spline
GROUP = {'ZrSiS': [11, 12], 'ZrSiSe': [11, 12], 'ZrSiTe': [11, 12, 13]}
EREF = {'ZrSiS': 9.70, 'ZrSiSe': 8.92, 'ZrSiTe': 8.86}
NG = (128, 128, 64)
Path('/data/spline').mkdir(exist_ok=True)
for name in sys.argv[1:] or GROUP:
    t = time.time()
    m = tb.Model(f'/data/model/{name}.npz', name)
    E = grid.eig_grid(m, *NG)
    coef = spline.prefilter(spline.coeffs_from_eigs(E[..., GROUP[name]] - EREF[name]))
    np.savez(f'/data/spline/{name}.npz', coef=coef, eref=EREF[name], group=np.array(GROUP[name]) + 1, B=m.B,
             E_all=E[..., 9:16].astype(np.float32))
    print(name, round(time.time() - t), 's')
