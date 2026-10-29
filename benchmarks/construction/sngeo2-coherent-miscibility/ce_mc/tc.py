import numpy as np, pickle, fittools as F, mc, sys, json
KB=mc.KB
def binder(xs,N,T,mu0):
    """reweight in mu to the max-variance point; return U4, mu*, x*"""
    beta=1/(KB*T); best=None
    for dmu in np.linspace(-6,6,241):
        lw=beta*dmu*N*xs; w=np.exp(lw-lw.max()); w/=w.sum()
        m=(w*xs).sum(); v=(w*(xs-m)**2).sum()
        if best is None or v>best[0]:
            m4=(w*(xs-m)**4).sum(); best=(v,1-m4/(3*v*v),mu0+dmu,m)
    return best
if __name__=='__main__':
    lam=float(sys.argv[1]); L=int(sys.argv[2]); Lz=int(sys.argv[3]); Ts=[float(t) for t in sys.argv[4].split(',')]
    nsw=int(sys.argv[5]); mu0=float(sys.argv[6]); tag=sys.argv[7]
    xs_,ef,Phi,orbs=pickle.load(open('big.pkl','rb'))
    w=F.lasso(Phi,ef*1000,lam)
    P,J,groups=mc.build(orbs,w,L,Lz); N=P.shape[0]
    rng=np.random.default_rng(L*1000+int(Ts[0]))
    out=[]
    s=rng.choice([-1.,1.],N)
    for T in Ts:
        mu=mu0
        # short pre-run to tune mu toward symmetric point
        for it in range(3):
            x=mc.sgc(P,J,groups,T,mu,s,150,rng)
            v,u,mus,m=binder(x,N,T,mu); mu=mus
        x=mc.sgc(P,J,groups,T,mu,s,nsw,rng,nequil=200)
        v,u,mus,m=binder(x,N,T,mu)
        out.append(dict(T=T,L=L,mu=mus,U4=u,var=v*N,x=m,xmin=float(x.min()),xmax=float(x.max())))
        print(json.dumps(out[-1]),flush=True)
    json.dump(out,open(f'tc_{tag}.json','w'))
