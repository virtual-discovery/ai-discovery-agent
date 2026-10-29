import sys, time, numpy as np
sys.path.insert(0, '.')
import kag_ssmf as K
key, N, U, JU, mix = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
m = K.load_model(f'/data/models/{key}.npz')
s = K.SSMF(m)
t = time.time()
r = s.run(N, U, JU * U, mix=mix, verbose=True, maxiter=300)
print(key, N, U, 'Z', np.round(r['Z'], 4), 'n', np.round(r['n'], 4), 'conv', r['converged'], r['iterations'], 'time %.1f' % (time.time() - t))
