import numpy as np, sys
sys.path.insert(0, '.')
import chambers as C
HB2_2M = 3.80998    # hbar^2/(2 m_e) in eV A^2
class Para:
    def __init__(self, mass, a=10.0):
        self.mass = mass; self.Binv = np.linalg.inv(2*np.pi*np.eye(3)/a); self.B = 2*np.pi*np.eye(3)/a
    def eig_vel(self, kf):
        kc = (kf - np.round(kf)) @ self.B
        e = HB2_2M/self.mass * (kc**2).sum(1)
        v = 2*HB2_2M/self.mass * kc
        return e[:, None], v[:, None, :]
sp = Para(0.2)
k0 = np.array([[0.01, 0.0, 0.0], [0.0, 0.012, 0.003]])
X = np.array([0, 0.5, 1.0, 2.0, 5.0])
tr = C.integrate(sp, k0, np.zeros(2, int), [0, 0, 1], 200, dl=0.0005)
avg, per = C.orbit_average(tr, X)
m = 0.2 * 9.1093837e-31
S_exact = 2*np.pi*m/1.602176634e-19 * 1e12
print('period', per, 'exact', S_exact, 'closed', tr['closed'], tr['nsteps'])
wct = 2*np.pi*X/S_exact
for i in range(2):
    r = avg[i, :, 0, 0] / avg[i, 0, 0, 0] if avg[i,0,0,0] else avg[i,:,1,1]/avg[i,0,1,1]
    print('sigma_xx ratio', r.round(4), 'expected', (1/(1+wct**2)).round(4))
    print('xy/xx0', (avg[i, :, 0, 1] / (avg[i, 0, 0, 0]+avg[i,0,1,1])*2).round(4), 'expected +-', (wct/(1+wct**2)).round(4))
