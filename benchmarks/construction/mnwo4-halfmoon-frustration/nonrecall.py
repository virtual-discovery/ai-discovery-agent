"""Non-recall fitting route (no literature values): joint differential evolution over J1..J12, D and the SCGA
temperature scale (10 spectra + block-averaged 20/30 K maps + 11 K ordering-peak constraint), joint polish,
J3 - J11 null-direction scan, predictions. Uses the lab-rounded data (4 significant digits).

    python3 -B nonrecall.py CACHE SEED      (~13 min single-threaded; set OMP_NUM_THREADS=1)

Writes nonrecall_seedSEED.json (parameters after the null scan, tau, DE objective) and prints the graded scores.
"""
import numpy as np, sys, json, time, os
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
from model import *
from fitlib import Joint, null_scan, predictions
seed=int(sys.argv[2]); CACHE=sys.argv[1]
shells,J0,D0,_=truth(CACHE+'/model.json'); M=Model(shells)
X=np.load(CACHE+'/extracted.npz')
r4=np.vectorize(lambda x: float('%.4g'%x))
IDX=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,88,94,100]
spectra=[(X['spectra_q'][i],X['spectra_E'][:np.shape(X['spectra_S'])[1]],r4(X['spectra_S'][i])) for i in IDX]
maps={T:half_to_full(norm3(r4(X['map_T%d'%T][186:,:,:3]))) for T in (20,30)}
t0=time.time()
J=Joint(M,spectra,maps,[-8/36,0.5,16/36])
p,f=J.search(seed)
q=null_scan(M,p[:13],maps); pr=predictions(M,q)
ref14=map14_blocks(X['map_T14'][:,:,3]); red={m:ko_blocks(X['ko_'+m][:,:,3]) for m in ('5d','5e')}
def pear(a,b): a=np.ravel(a)-np.mean(a); b=np.ravel(b)-np.mean(b); return float(a@b/np.sqrt((a@a)*(b@b)))
kref,_=M.lt(J0); kd=min(np.abs((np.array(pr['k'])-s*kref)-np.round(np.array(pr['k'])-s*kref)).max() for s in (1,-1))
json.dump({'p':q.tolist(),'tau':float(p[13]),'f':f},open(str(__import__('pathlib').Path(__file__).resolve().parent/('nonrecall_seed%d.json'%seed)),'w'))
print('seed %d f %.3f [%.0fs] J1 %.3f J6 %.3f J9 %.3f D %.3f | kdist %.4f dE %.4f r14 %.3f r5d %.3f r5e %.3f'%(seed,f,time.time()-t0,q[0],q[5],q[8],q[12],kd,pr['dE'],pear(pr['map14'],ref14),pear(pr['reduced']['5d'],red['5d']),pear(pr['reduced']['5e'],red['5e'])),flush=True)
