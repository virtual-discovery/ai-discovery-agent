import numpy as np, time, sys
from nep_numpy import NEP
from xyzio import read_xyz
R='/raw/'
tr=read_xyz(R+'Training-sets/Training-set-qNEP-w-BECs.xyz')
for pot in ['Potential-NEP.txt','Potential-qNEP.txt','Potential-qNEP-w-BECs.txt']:
    m=NEP(R+'Potentials/'+pot)
    de=[];df=[];db=[]
    sel=list(range(0,len(tr),7))
    t0=time.time()
    for j in sel:
        c=tr[j]; ty=[m.elements.index(s) for s in c['symbols']]
        o=m.compute(ty,c['pos'],c['lattice'],want_bec=('bec' in c and m.charge))
        de.append((o['energy']-c['energy'])/len(ty)); df.append(np.sqrt(((o['forces']-c['forces'])**2).mean()))
        if 'bec' in o: db.append(np.sqrt(((o['bec'].reshape(-1,9)-c['bec'])**2).mean()))
    de=np.array(de)
    print(pot,'n',len(sel),'E mean off %.4f rmse(after mean) %.4f eV/at'%(de.mean(),np.std(de)),'F rmse %.3f'%np.mean(df),'BEC rmse',np.mean(db) if db else None,'time',time.time()-t0)
