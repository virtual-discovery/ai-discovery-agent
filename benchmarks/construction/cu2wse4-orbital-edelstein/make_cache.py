"""Compute the uniform-mesh histogram caches served by the lab's Fermi-surface requests.

    bash run_in_image.sh make_cache.py MODEL N [N ...]      # MODEL wosoc|soc; mesh N x N x N/2, e.g. 24 32 48 64 96

Each cache holds, for one archived Wannier model on the Gamma-centred N x N x N/2 mesh of the conventional cell,
the energy histograms (1 meV bins, absolute eV of the HWR scale, EMIN..EMAX) of the state count and of
m^S_i v_j, m^O_i v_j, m^L_i v_j and v_a Omega_b for all i, j (see cws_tb.py). Writes /data/cache/MODEL_N.npz.
"""
import sys
import time
from pathlib import Path
import numpy as np
import cws_tb as T

CACHE = Path(__import__('os').environ.get('CWS_CACHE_DIR', '/data/cache'))
EMIN, EMAX = -7.6, -2.2

if __name__ == '__main__':
    CACHE.mkdir(parents=True, exist_ok=True)
    model = sys.argv[1]
    workers = int(__import__('os').environ.get('WORKERS', '4'))
    for n in map(int, sys.argv[2:]):
        out = CACHE/f'{model}_{n}.npz'
        if out.exists():
            print('exists', out); continue
        t = time.time()
        c = T.mesh_cache(f'/data/model/{model}.npz', (n, n, n // 2), EMIN, EMAX, workers=workers)
        c['seconds'] = np.array(time.time() - t)
        np.savez_compressed(out, **c)
        print(out, 'seconds', round(time.time() - t), flush=True)
