import numpy as np
from pathlib import Path
np.set_printoptions(precision=3, suppress=True, linewidth=150)
for key in ['CsTi3Sb5','CsV3Sb5','CsCr3Sb5','CsTi3Bi5','CsCr3Bi5','CsCr3Bi5_9.2']:
    d = np.load(f'/data/models/{key}.npz'); R, deg, H, ef = d['R'], d['degen'], d['H'], float(d['efermi'])
    i0 = np.where((R == 0).all(1))[0][0]
    H0 = H[i0]
    print('=====', key, 'EF', ef)
    print('onsite diag', np.real(np.diag(H0)) - ef)
    for s in range(3):
        blk = H0[5*s:5*s+5, 5*s:5*s+5].real - ef*np.eye(5)
        print('site', s, 'block\n', blk)
    n = 12
    ks = np.array([(i/n, j/n, l/8) for i in range(n) for j in range(n) for l in range(8)])
    ph = np.exp(2j*np.pi*ks@R.T)/deg
    Hk = np.einsum('kr,rmn->kmn', ph, H)
    e, v = np.linalg.eigh(Hk)
    for N in [27, 30, 33]:
        es = np.sort(e.ravel()); mu = 0.5*(es[len(ks)*N//2-1] + es[len(ks)*N//2])
        print('N', N, 'mu-EF %.3f' % (mu-ef))
    ne = 2*(e < ef).sum()/len(ks)
    print('electrons below archive EF: %.3f' % ne)
