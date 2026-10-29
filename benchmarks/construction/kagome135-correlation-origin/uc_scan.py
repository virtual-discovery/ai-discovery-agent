"""Lab-path metallic/non-metallic outcome on a U grid (construction check of the Mott boundary)."""
import json, sys
import numpy as np
sys.path.insert(0, '.')
import kag_ssmf as K
key, N, JU, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
grid = [float(x) for x in sys.argv[5].split(',')]
m = K.load_model(f'/data/models/{key}.npz')
s = K.SSMF(m, mesh=K.LAB_MESH, kT=K.LAB_KT)
res = {}
for U in grid:
    r = K.lab_ssmf(m, N, U, JU, solver=s)
    res[U] = r
    print(key, N, JU, U, r['metallic'], np.round(r.get('Z', []), 4), flush=True)
json.dump({str(k): v for k, v in res.items()}, open(out, 'w'))
