import re, numpy as np
def read_xyz(fn):
    lines=open(fn).read().split('\n'); i=0; out=[]
    while i<len(lines) and lines[i].strip():
        n=int(lines[i]); h=lines[i+1]
        lat=np.array(re.search(r'Lattice="([^"]+)"',h).group(1).split(),float).reshape(3,3)
        props=re.search(r'Properties=(\S+)',h).group(1).split(':')
        cols=[];k=0
        for j in range(0,len(props),3):
            name,typ,cnt=props[j],props[j+1],int(props[j+2]); cols.append((name,k,cnt)); k+=cnt
        rows=[l.split() for l in lines[i+2:i+2+n]]
        d={'lattice':lat,'header':h}
        for name,k0,cnt in cols:
            if name=='species': d['symbols']=[r[k0] for r in rows]
            else: d[name]=np.array([[float(x) for x in r[k0:k0+cnt]] for r in rows])
        for key in ('energy','free_energy'):
            m=re.search(r'(?<![\w_])'+key+r'=([-\d.eE+]+)',h)
            if m: d[key]=float(m.group(1))
        m=re.search(r'stress="([^"]+)"',h)
        if m: d['stress']=np.array(m.group(1).split(),float).reshape(3,3)
        m=re.search(r'virial="([^"]+)"',h)
        if m: d['virial']=np.array(m.group(1).split(),float).reshape(3,3)
        out.append(d); i+=n+2
    return out
