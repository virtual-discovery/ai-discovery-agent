import numpy as np, json, sys, time
from nep_numpy import NEP
from xyzio import read_xyz
from builder import Bulk, build_wall
from relax import Relaxer
args=json.loads(sys.argv[1]); tag=args['tag']
m=NEP('/raw/Potentials/Potential-%s.txt'%args.get('pot','qNEP-w-BECs'))
B=read_xyz('/raw/Bulk-structures/Cmc2_1-qNEP-w-BECs.xyz')[0]
I4=read_xyz('/raw/Bulk-structures/I4mmm-qNEP-w-BECs.xyz')[0]
cell=np.diag(np.diag(B['lattice']))
bulk=Bulk(B['symbols'],B['pos'],cell,I4['symbols'],I4['pos']/np.diag(I4['lattice']))
if 'file' in args:
    s=read_xyz(args['file'])[0]; sym,pos,c=s['symbols'],s['pos'],s['lattice']
else:
    sym,pos,c=build_wall(bulk,args['frame'],args['n'],args['op'],offset=args.get('offset',0.0),nz=args.get('nz',1),axis=args.get('axis',0))
ty=[m.elements.index(x) for x in sym]
eb=m.compute([m.elements.index(x) for x in B['symbols']],B['pos'],cell)['energy']/48
t0=time.time()
r=Relaxer(m,ty,pos,c,relax_axis=args.get('relax_axis'))
res=r.run(fmax=args.get('fmax',0.02),max_steps=args.get('steps',600))
c2=res['cell']; ax=args.get('axis',0)
A=np.linalg.norm(np.cross(c2[(ax+1)%3],c2[(ax+2)%3]))
e0=r.evaluate(np.r_[pos.ravel(),[0.0]*(args.get('relax_axis') is not None)])[0]['energy']
out={'tag':tag,'args':args,'N':len(sym),'E':res['energy'],'E_initial':e0,'Ebulk_per_atom':eb,'A':A,
     'E_DW_mJm2':(res['energy']-len(sym)*eb)/(2*A)*16021.77,'E_DW_initial':(e0-len(sym)*eb)/(2*A)*16021.77,
     'maxF':res['max_force'],'converged':res['converged'],'iters':res['iterations'],'evals':res['evaluations'],
     'strain':res['strain'],'seconds':time.time()-t0,'msg':res['message']}
print(json.dumps(out))
np.savez('runs/%s.npz'%tag,sym=np.array(sym),pos=res['positions'],cell=c2)
json.dump(out,open('runs/%s.json'%tag,'w'))
