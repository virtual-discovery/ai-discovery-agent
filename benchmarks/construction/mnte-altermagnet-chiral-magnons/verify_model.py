"""Checks behind two build decisions (run in the runtime image after fetch.py):

1. The LSWT of lswt.py with the authors' Hamiltonian (J1 3.99, J2 -0.120, J3 0.472, J10 0.0681, J11 -0.0221, D 0.0482)
   reproduces the archived dispersion Fig1c.mat along (h, 0, -1.22).
2. forward.hyspec_energy reproduces the archived HYSPEC trajectory Etraj(h) of Fig1c.mat along (h, 0, -1.22).
"""
import numpy as np
import scipy.io as sio
import forward as F
import lswt
import reduce as R

d = sio.loadmat(R.X + '/Fig1c.mat')
J = {1: 3.99, 2: -0.120, 3: 0.472, 10: 0.0681, 11: -0.0221}
fr = lswt.frames(30.0)
q, Et = d['q'][0], d['Et']
dev = [np.abs(lswt.modes([q[i], 0, -1.22], J, 0.0482, fr)[0] - np.sort(Et[:2, i])).max() for i in range(0, 4001, 100)]
print('max |LSWT - Fig1c.mat| along (h,0,-1.22): %.4f meV' % max(dev))
h, E = d['h'][0], d['Etraj'][0]
v = (-1.22 + 2.23 * h) / 3.83
print('max |E_kinematic - Etraj|: %.2e meV' % np.abs(F.hyspec_energy(h - v, v) - E).max())
