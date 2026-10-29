import sys, time, numpy as np
sys.path.insert(0, '.')
import kag_ssmf as K
key, N, U, JU = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
mesh = tuple(int(x) for x in sys.argv[5].split(',')) if len(sys.argv) > 5 else (12, 12, 8)
m = K.load_model(f'/data/models/{key}.npz')
s = K.SSMF(m, mesh=mesh)
t = time.time()
r = s.solve(N, U, JU * U)
print(r['ramp'], flush=True)
print(key, N, U, JU, mesh, 'Z', np.round(r['Z'], 4), 'n', np.round(r['n'], 4), 'conv', r['converged'], r['iterations'], 'res %.1e' % r['residual'], 'time %.1f' % (time.time() - t))
