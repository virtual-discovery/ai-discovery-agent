import sys, numpy as np
sys.path.insert(0, '.')
import kag_ssmf as K
np.set_printoptions(precision=3, suppress=True, linewidth=150)
for key in ['CsTi3Sb5','CsV3Sb5','CsCr3Sb5','CsTi3Bi5','CsCr3Bi5','CsCr3Bi5_9.2']:
    m = K.load_model(f'/data/models/{key}.npz')
    sig = K.sigma_exact(m)
    out = [key, 'sigma_d', sig[:5]]
    for n in [(12,12,8), (24,24,16)]:
        ks = K.kmesh(*n); e = np.linalg.eigvalsh(K.hk(m, ks))
        out.append('N(EF)@%s=%.3f/%.3f' % (n[0], K.electron_count(e, m['efermi'], 0.01), K.electron_count(e, m['efermi'], 0.002)))
    print(*out)
