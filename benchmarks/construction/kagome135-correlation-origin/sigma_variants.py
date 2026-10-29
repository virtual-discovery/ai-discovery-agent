import sys, numpy as np
sys.path.insert(0, '.')
import kag_ssmf as K
np.set_printoptions(precision=3, suppress=True, linewidth=150)
for key in ['CsTi3Sb5','CsV3Sb5','CsCr3Sb5','CsTi3Bi5','CsCr3Bi5','CsCr3Bi5_9.2']:
    m = K.load_model(f'/data/models/{key}.npz')
    full = K.sigma_exact(m)[:5]
    H = m['H'].copy(); H[:, :15, 15:] = 0; H[:, 15:, :15] = 0
    donly = K.sigma_exact(dict(m, H=H))[:5]
    ks, w = K.kmesh_tr(18, 18, 12)
    e, v = np.linalg.eigh(K.hk(m, ks)); ef = m['efermi']
    wgt = (np.abs(v[:, :5, :]) ** 2) * w[:, None, None]
    E = e[:, None, :] - ef
    win = (np.abs(E) < 4)
    mom1 = (wgt * E * win).sum((0, 2)) / (wgt * win).sum((0, 2))
    mom2 = (wgt * (E - mom1[None, :, None]) ** 2 * win).sum((0, 2)) / (wgt * win).sum((0, 2))
    m1 = (wgt * E).sum((0, 2)) / wgt.sum((0, 2)); m2 = (wgt * (E - m1[None, :, None]) ** 2).sum((0, 2)) / wgt.sum((0, 2))
    print(key, 'full', full, 'mesh18', np.sqrt(m2), 'd-only', donly, 'win4', np.sqrt(mom2), 'mean-EF', m1)
