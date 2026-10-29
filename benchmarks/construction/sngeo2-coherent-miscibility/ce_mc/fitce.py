import ce, json, numpy as np, sys
from pos import parse, X
aG,cG,aS,cS=4.4827,2.9075,4.8327,3.2422
def load(fn='PBE', cut={1:0.1,2:7.0,3:5.0,4:4.0}):
    d=json.load(open(X+f'data_CE_{fn}.json'))
    orbs=ce.enumerate_orbits((aG+aS)/2,(cG+cS)/2,cut)
    Phi=[];E=[];xs=[];Ts=[]
    eG=eS=None
    recs=[]
    for r in d:
        A,sp,fr=parse(r['structure']); n=sp.count('Ge')+sp.count('Sn'); x=sp.count('Ge')/n
        T,sites=ce.struct_sites(A,sp,fr,x*aG+(1-x)*aS,x*cG+(1-x)*cS)
        assert abs(round(abs(np.linalg.det(T)))*2-n)<1e-9,(T,n)
        recs.append((x,r['energy']/n,ce.correlations(T,sites,orbs),T))
    eG=[r[1] for r in recs if r[0]==1][0]; eS=[r[1] for r in recs if r[0]==0][0]
    xs=np.array([r[0] for r in recs]); ef=np.array([r[1]-r[0]*eG-(1-r[0])*eS for r in recs])
    Phi=np.array([r[2] for r in recs])
    return xs,ef,Phi,orbs,recs
if __name__=='__main__':
    xs,ef,Phi,orbs,recs=load()
    print(Phi.shape); np.save('phi.npy',Phi); np.save('ef.npy',ef); np.save('xs.npy',xs)
    # check: point correlation = 2x-1
    print('point check',np.abs(Phi[:,1]-(2*xs-1)).max())
    # duplicates of correlation vectors?
    from itertools import combinations
    dup=[(i,j) for i,j in combinations(range(len(xs)),2) if np.allclose(Phi[i],Phi[j])]
    print('duplicate corr vectors',dup)
    print('rank',np.linalg.matrix_rank(Phi))
