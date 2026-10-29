import numpy as np, ce
KB=8.617333e-5*1000  # meV/K
def build(orbs, eci, L, Lz):
    """returns partner array P (N,K,3), couplings J (N,K) for local fields; eci per cation (meV), eci[0]=V0, eci[1..]=orbits"""
    ncell=L*L*Lz; N=2*ncell
    def idx(cell,sub):
        i,j,k=cell; return 2*(((i%L)*L+(j%L))*Lz+(k%Lz))+sub
    inst=[[] for _ in range(N)]
    lin=np.zeros(N)
    for o,V in zip(orbs,eci[1:]):
        if abs(V)<1e-9: continue
        m=len(o['clusters']); Jc=2*V/m
        for cl in o['clusters']:
            pts=[ce.site_index(np.array(p)) for p in cl]
            for i in range(L):
                for j in range(L):
                    for k in range(Lz):
                        ids=[idx((c[0]+i,c[1]+j,c[2]+k),s) for c,s in pts]
                        for a in range(len(ids)):
                            others=ids[:a]+ids[a+1:]
                            inst[ids[a]].append((others,Jc))
    K=max(len(v) for v in inst)
    P=np.full((N,K,3),N,dtype=np.int64); J=np.zeros((N,K))
    for n,v in enumerate(inst):
        for q,(oth,Jc) in enumerate(v):
            P[n,q,:len(oth)]=oth; J[n,q]=Jc
    # coloring
    adj=[set() for _ in range(N)]
    for n,v in enumerate(inst):
        for oth,_ in v: adj[n].update(oth)
    color=-np.ones(N,int)
    for n in range(N):
        used={color[m] for m in adj[n] if color[m]>=0}
        c=0
        while c in used: c+=1
        color[n]=c
    groups=[np.where(color==c)[0] for c in range(color.max()+1)]
    return P,J,groups
def energy(s,P,J):
    se=np.append(s,1)
    h=(J*se[P].prod(axis=2)).sum(1)
    return 0  # unused
def sgc(P,J,groups,T,mu,s,nsweep,rng,nequil=0):
    """mu in meV per Ge (exchange chem pot): H - mu*N_Ge"""
    N=len(s); se=np.append(s,1).astype(float); beta=1/(KB*T); xs=[]
    for sw in range(nsweep):
        for g in groups:
            h=(J[g]*se[P[g]].prod(axis=2)).sum(1)
            # flip s->-s: dE = -2 s h  (cluster energies), chem: -mu*(dN_Ge) with dN_Ge = -s
            sg=se[g]; dE=-2*sg*h + mu*sg
            acc=rng.random(len(g))<np.exp(-beta*np.clip(dE,-50/beta,None))
            se[g[acc]]*=-1
        if sw>=nequil: xs.append((se[:N].mean()+1)/2)
    s[:]=se[:N]
    return np.array(xs)
