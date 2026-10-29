import numpy as np 
import matplotlib.pyplot as plt 

fnumber = 3884 #only edit me please :)


data = np.load('{}_Cut2D_Data.npy'.format(fnumber))
bins = np.load('{}_Cut2D_BinEdges.npy'.format(fnumber))
norms = (data[0]*data[3])/(data[1]*data[2])
errs = (np.sqrt(data[0])*data[3])/(data[1]*data[2])

Qxs = bins[0, :-1, :-1] + 0.015
Qys = bins[1, :-1, :-1] + 0.015
Is = norms
errs = errs

print(np.shape(Qxs))
np.savetxt('{}_Qxs.csv'.format(fnumber), Qxs, delimiter=',')
np.savetxt('{}_Qys.csv'.format(fnumber), Qys, delimiter=',')
np.savetxt('{}_Is.csv'.format(fnumber), Is, delimiter=',')
np.savetxt('{}_errs.csv'.format(fnumber), errs, delimiter=',')


#errs = (np.sqrt(data[0])*data[3])/(data[1]*data[2])
plt.imshow(errs, vmax=1e-5)
plt.figure()
plt.imshow(norms, vmax=1e-5)
plt.close()
#print(np.shape(data))