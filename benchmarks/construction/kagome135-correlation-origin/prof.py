import sys, time, numpy as np
sys.path.insert(0, '.')
import kag_ssmf as K
m = K.load_model('/data/models/CsCr3Sb5.npz')
s = K.SSMF(m)
t=time.time(); r = s.fermion_step(np.full(5,0.5), np.zeros(5), 34); print('fermion %.2f'%(time.time()-t))
site = K.SlaveSpinSite(6, 1.08)
t=time.time()
for i in range(20): site.ground(np.full(5,-1.0), np.zeros(5)+0.1*i)
print('ground x20 %.2f'%(time.time()-t))
