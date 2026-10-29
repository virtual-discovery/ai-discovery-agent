"""Validate the SSMF evaluator against the numbers quoted in arXiv:2512.22576 (both filling conventions)."""
import json, sys, time
import numpy as np
sys.path.insert(0, '.')
import kag_ssmf as K
jobs = [a.split(':') for a in sys.argv[2:]]
out = {}
for key, N, U, JU in jobs:
    m = K.load_model(f'/data/models/{key}.npz')
    s = K.SSMF(m)
    t = time.time()
    r = s.solve(float(N), float(U), float(JU) * float(U), maxiter=400)
    out[f'{key}:{N}:{U}:{JU}'] = {'Z': r['Z'].round(4).tolist(), 'n': r['n'].round(4).tolist(), 'conv': r['converged'],
                                  'it': r['iterations'], 's': round(time.time() - t, 1)}
    print(f'{key} N={N} U={U} J/U={JU}', np.round(r['Z'], 3), 'n', np.round(r['n'], 3), r['converged'], r['iterations'], flush=True)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
