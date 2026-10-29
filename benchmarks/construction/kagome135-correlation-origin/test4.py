import sys, time, numpy as np
sys.path.insert(0, '.')
import kag_ssmf as K
key, N, JU = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
m = K.load_model(f'/data/models/{key}.npz')
s = K.SSMF(m, mesh=(9, 9, 6))
t = time.time()
r = s.scan(N, JU, [2, 3, 4, 5, 5.5, 6, 6.5, 7])
for u, x in sorted(r['results'].items()):
    print(u, np.round(x['Z'], 4), np.round(x['n'], 4), x['iterations'], '%.1e' % x['residual'], flush=True)
print(key, N, JU, 'stopped', r['stopped_at'], 'path', [(round(a, 3), b) for a, b in r['path']], 'time %.0f' % (time.time() - t))
