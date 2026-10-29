import numpy as np, time
from nep_numpy import NEP
for d in ['test_qnep','test_nep']:
    base='nepcpu/NEP_CPU-main/'+d+'/'
    try: m=NEP(base+'nep.txt')
    except ValueError as e: print(d,e); continue
    t=open(base+'xyz.in').read().split()
    N=int(t[0]); cell=np.array(t[1:10],float).reshape(3,3)
    rows=np.array(t[10:10+4*N]).reshape(N,4)
    types=[m.elements.index(s) for s in rows[:,0]]; pos=rows[:,1:].astype(float)
    t0=time.time(); o=m.compute(types,pos,cell,want_bec=True); print('time',time.time()-t0)
    fref=np.loadtxt(base+'force_analytical_ref.out'); vref=np.loadtxt(base+'virial_ref.out').sum(0).reshape(3,3)
    print(d,'max|dF|',abs(o['forces']-fref).max(),'max|F|',abs(fref).max())
    print('virial mine\n',o['virial'],'\nref\n',vref)
    dref=np.loadtxt(base+'descriptor_ref.out'); print('desc shape',dref.shape)
