import numpy as np, json
import os
X=os.path.expanduser('~/.cache/ai-theorist/sngeo2-coherent-miscibility/x/')
def parse(s):
    L=s.strip('\n').split('\n')
    sc=float(L[1]); A=np.array([[float(v) for v in L[i].split()] for i in (2,3,4)])*sc
    el=L[5].split(); n=[int(v) for v in L[6].split()]
    N=sum(n); fr=np.array([[float(v) for v in L[8+i].split()[:3]] for i in range(N)])
    sp=[e for e,k in zip(el,n) for _ in range(k)]
    return A,sp,fr
