"""Compute the uniform-mesh Berry-curvature histogram caches served by the lab's `ahc` and `electron_count` requests.

    bash run_in_image.sh make_cache.py N [N ...]        # meshes N x N x N/2, e.g. 48 72 100 150 200 300 400

Each cache holds, for the archived Wannier model (Neel vector +y), the Omega-weighted energy histograms of all
32 bands on the Gamma-centred N x N x N/2 mesh (bins of 0.1 meV between 4.80 and 8.80 eV, plus the summed
curvature of all states below the window) and the state-count histogram. sigma(mu) and the electron count for
any chemical potential in the window and any Fermi smearing follow exactly (mnte_tb.sigma_from_cache). The Neel-
reversed model is the time reverse of the archived one; on these inversion-symmetric meshes its caches are the
same histograms with the curvature negated, so they are not stored separately.
Writes CACHE/mesh_N.npz (CACHE = /data/cache inside the image, i.e. ~/.cache/ai-theorist/<id>/cache).
The 400 mesh takes about 3 h on 4 cores.
"""
import sys
import time
from pathlib import Path
import numpy as np
import mnte_tb as T

CACHE = Path('/data/cache')
MODEL = Path('/data/model/archived.npz')

if __name__ == '__main__':
    CACHE.mkdir(parents=True, exist_ok=True)
    for n in map(int, sys.argv[1:]):
        out = CACHE/f'mesh_{n}.npz'
        if out.exists():
            print('exists', out); continue
        t = time.time()
        c = T.mesh_cache(str(MODEL), (n, n, n // 2), workers=4)
        c['seconds'] = np.array(time.time() - t)
        np.savez_compressed(out, **c)
        s = T.sigma_from_cache(c, 6.13013 + np.array([-0.8, -0.5, -0.3]), 0.001)
        print(out, 'seconds', round(time.time() - t), 'sigma_xy', s[2].round(2), 'yz', s[0].round(3), 'zx', s[1].round(3), flush=True)
