import sys, time, numpy as np
sys.path.insert(0, '.')
import kag_ssmf as K
m = K.load_model('/data/models/CsCr3Sb5.npz')
s = K.SSMF(m)
H = s.Hoff.copy()
t=time.time(); e,v=np.linalg.eigh(H); print('eigh %.3f'%(time.time()-t))
t=time.time(); mu=K.find_mu(e,34,0.01); print('mu %.3f'%(time.time()-t))
t=time.time(); f=K.fermi((e-mu)/0.01); nf = np.einsum('kan,kn->a', np.abs(v) ** 2, f); print('nf %.3f'%(time.time()-t))
t=time.time(); AV=np.matmul(H[:, :15, :], v); print('mm %.3f'%(time.time()-t))
t=time.time(); h = np.real(np.einsum('kan,kn,kan->a', AV, f, v[:, :15, :].conj(), optimize=True)); print('h %.3f'%(time.time()-t))
