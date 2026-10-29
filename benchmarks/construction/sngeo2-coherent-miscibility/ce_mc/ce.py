"""Rutile (Ge,Sn)O2 cation cluster expansion utilities (numpy only)."""
import numpy as np, itertools, json
from pos import parse, X
U=0.3066
SITES=np.array([[0,0,0],[.5,.5,.5]])
OX=np.array([[U,U,0],[-U,-U,0],[.5+U,.5-U,.5],[.5-U,.5+U,.5]])%1
def symops():
    ops=[]
    for perm in [(0,1),(1,0)]:
        for sx in (1,-1):
            for sy in (1,-1):
                for sz in (1,-1):
                    R=np.zeros((3,3)); R[0,perm[0]]=sx; R[1,perm[1]]=sy; R[2,2]=sz
                    for t in itertools.product([0,.5],repeat=3):
                        t=np.array(t); new=(OX@R.T+t)%1
                        ok=all(np.min(np.linalg.norm(((new[i]-OX+0.5)%1)-0.5,axis=1))<1e-6 for i in range(4))
                        ok2=all(np.min(np.linalg.norm(((((SITES@R.T+t)%1)[i]-SITES+0.5)%1)-0.5,axis=1))<1e-6 for i in range(2))
                        if ok and ok2: ops.append((R,t))
    return ops
OPS=symops()
def site_index(f):
    """prim-frac coordinate -> (cell int vector, sublattice)"""
    for s in (0,1):
        d=f-SITES[s]; r=np.round(d)
        if np.allclose(d,r,atol=1e-6): return tuple(r.astype(int)), s
    raise ValueError(f)
def canon(cl):
    """cluster = list of prim-frac coords; canonical translation-invariant key"""
    pts=sorted([tuple(np.round(p,6)) for p in cl])
    base=np.array(pts[0])-SITES[site_index(np.array(pts[0]))[1]]
    pts=sorted([tuple(np.round(np.array(p)-base,6)) for p in pts])
    return tuple(pts)
def orbit(cl):
    out=set()
    for R,t in OPS:
        out.add(canon([R@np.array(p)+t for p in cl]))
    return out
def enumerate_orbits(a,c,cut):
    """cut: dict n->cutoff (Angstrom, max pair distance in cluster)"""
    P=np.diag([a,a,c])
    rng=range(-3,4)
    pts=[np.array(i)+SITES[s] for i in itertools.product(rng,rng,rng) for s in (0,1)]
    pts=[p for p in pts if np.linalg.norm(p@P)<=max(cut.values())+1e-6]
    orbs=[]; seen=set()
    for n,rc in sorted(cut.items()):
        for s0 in (0,1):
            p0=SITES[s0]
            near=[p for p in pts if 1e-6<np.linalg.norm((p-p0)@P)<=rc+1e-6]
            for comb in itertools.combinations(range(len(near)),n-1):
                cl=[p0]+[near[i] for i in comb]
                dm=max([np.linalg.norm((x-y)@P) for x,y in itertools.combinations(cl,2)]+[0.0])
                if dm>rc+1e-6: continue
                k=canon(cl)
                if k in seen: continue
                ob=orbit(cl); seen|=ob
                orbs.append({'n':n,'diam':dm,'clusters':sorted(ob)})
    orbs.sort(key=lambda o:(o['n'],round(o['diam'],3)))
    return orbs
def struct_sites(A,sp,fr,a,c):
    """map a relaxed POSCAR onto prim lattice; return supercell matrix T and list of (cell,sub,spin)"""
    P=np.diag([a,a,c])
    T=np.round(A@np.linalg.inv(P)).astype(int)
    out=[]
    for f,s in zip(fr,sp):
        if s=='O': continue
        pf=f@T
        best=None
        for sub in (0,1):
            d=pf-SITES[sub]; r=np.round(d); e=np.linalg.norm((d-r)@P)
            if best is None or e<best[0]: best=(e,tuple(r.astype(int)),sub)
        assert best[0]<0.6, best
        out.append((best[1],best[2],1 if s=='Ge' else -1))
    return T,out
def correlations(T,sites,orbs):
    """average cluster functions per orbit (point term included as orbit n=1)"""
    T=np.array(T); Ti=np.linalg.inv(T)
    occ={}
    for cell,sub,spin in sites:
        f=np.array(cell)@Ti   # supercell-frac of the cell origin
        key=(tuple(np.round((f%1)*1e4).astype(int)%10000),sub)
        occ[key]=spin
    def spin_at(p):
        cell,sub=site_index(p)
        f=(np.array(cell)@Ti)%1
        return occ[(tuple(np.round(f*1e4).astype(int)%10000),sub)]
    N=len(sites)
    phi=[1.0]
    for o in orbs:
        tot=0; cnt=0
        for cell,sub,_ in sites:
            base=np.array(cell,float)
            for cl in o['clusters']:
                # cluster canonical form starts at a site of sublattice of its first point; translate so that first point sits on this site if sublattice matches
                p0=np.array(cl[0]); s0=site_index(p0)[1]
                if s0!=sub: continue
                shift=base+SITES[sub]-p0
                pr=1
                for p in cl: pr*=spin_at(np.array(p)+shift)
                tot+=pr; cnt+=1
        phi.append(tot/cnt)
    return np.array(phi)
