import sys, time, numpy as np
sys.path.insert(0, '.')
import kag_ssmf as K
key, N, U = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
JU = float(sys.argv[4]) if len(sys.argv) > 4 else 0.18
m = K.load_model(f'/data/models/{key}.npz')
t = time.time()
s = K.SSMF(m, mesh=(12, 12, 8))
mu, occ = s.noninteracting(N)
print('noninteracting n_d', np.round(occ[:5], 3), 'mu-EF', round(mu - m['efermi'], 3), 'setup %.1fs' % (time.time() - t))
r = s.run(N, U, JU * U, verbose='-v' in sys.argv)
print(key, N, U, 'Z', np.round(r['Z'], 4), 'n', np.round(r['n'], 4), 'conv', r['converged'], r['iterations'], 'time %.1f' % (time.time() - t))
