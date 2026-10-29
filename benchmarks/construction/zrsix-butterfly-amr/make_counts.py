"""Band-energy histograms of every model on the Gamma-centred N x N x N/2 meshes (exact FFT eigenvalues) for the
lab's electron-count requests: /data/counts/<name>.npz with hist[mesh, band, bin] (0.2 meV bins) and nlow[mesh, band]
(states below the histogram)."""
import numpy as np, sys
from pathlib import Path
sys.path.insert(0, '.')
import tb, grid
from common import COUNT_RANGE, COUNT_BIN, COUNT_MESHES

Path('/data/counts').mkdir(exist_ok=True)
for name in sys.argv[1:] or COUNT_RANGE:
    m = tb.Model(f'/data/model/{name}.npz', name)
    lo, hi = COUNT_RANGE[name][0] - 0.3, COUNT_RANGE[name][1] + 0.3
    nb = int(round((hi - lo) / COUNT_BIN))
    hist = np.zeros((len(COUNT_MESHES), m.nw, nb), np.int32); nlow = np.zeros((len(COUNT_MESHES), m.nw), np.int64)
    for i, n in enumerate(COUNT_MESHES):
        E = grid.eig_grid(m, n, n, n // 2).reshape(-1, m.nw)
        for b in range(m.nw):
            idx = np.floor((E[:, b] - lo) / COUNT_BIN).astype(np.int64)
            nlow[i, b] = (idx < 0).sum()
            ok = (idx >= 0) & (idx < nb)
            hist[i, b] = np.bincount(idx[ok], minlength=nb)
        print(name, n, 'states', E.size, flush=True)
    np.savez_compressed(f'/data/counts/{name}.npz', hist=hist, nlow=nlow, lo=lo, ebin=COUNT_BIN, meshes=np.array(COUNT_MESHES))
