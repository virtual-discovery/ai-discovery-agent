import numpy as np, fittools as F, fitce, pickle, sys
cut={1:0.1,2:9.0,3:6.6,4:5.6}
xs,ef,Phi,orbs,recs=fitce.load(cut=cut)
pickle.dump((xs,ef,Phi,orbs),open('big.pkl','wb'))
print(Phi.shape, 'rank',np.linalg.matrix_rank(Phi))
for o in orbs: print(o['n'],round(o['diam'],3),len(o['clusters']))
