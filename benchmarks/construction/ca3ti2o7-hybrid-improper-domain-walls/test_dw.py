import numpy as np, time, glob
from nep_numpy import NEP
from xyzio import read_xyz
R='/raw/'
pots={p:NEP(R+'Potentials/Potential-%s.txt'%p) for p in ['NEP','qNEP','qNEP-w-BECs']}
def ev(m,s,bec=False):
    ty=[m.elements.index(x) for x in s['symbols']]; return m.compute(ty,s['pos'],s['lattice'],want_bec=bec)
bulk={}
for fn in ['Cmc2_1-qNEP-w-BECs','Cmc2_1-DFT','I4mmm-qNEP-w-BECs','I4mmm-DFT']:
    s=read_xyz(R+'Bulk-structures/%s.xyz'%fn)[0]
    for p,m in pots.items():
        o=ev(m,s); n=len(s['symbols'])
        print(fn,p,'E/fu %.5f'%(o['energy']/n*12),'maxF %.3f'%abs(o['forces']).max(),'virial diag/N',np.round(np.diag(o['virial'])/n,4))
        bulk[(fn,p)]=o['energy']/n*12
for kind in ['qNEP-w-BECs','DFT']:
  for fn in sorted(glob.glob(R+'Domain-walls/%s/*.xyz'%kind)):
    s=read_xyz(fn)[0]; L=s['lattice']; n=len(s['symbols'])
    lens=np.linalg.norm(L,axis=1); k=int(np.argmax(lens)); A=np.linalg.norm(np.cross(L[(k+1)%3],L[(k+2)%3]))
    out=[]
    for p,m in pots.items():
        t0=time.time(); o=ev(m,s); dt=time.time()-t0
        e=(o['energy']-n/12*bulk[('Cmc2_1-qNEP-w-BECs',p)])/(2*A)*16021.77
        out.append('%s %.1f mJ/m2 (maxF %.3f, %.1fs)'%(p,e,abs(o['forces']).max(),dt))
    print(kind,fn.split('/')[-1],n,'A=%.1f'%A,' | '.join(out))
