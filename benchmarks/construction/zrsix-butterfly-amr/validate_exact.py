"""Validate the spline-based orbit integrator against exact diagonalization of the archived Hamiltonian.

For random Fermi-surface start points (ZrSiSe, mu = neutral, field along c and tilted), integrate the same orbits
with velocities from exact Hellmann-Feynman diagonalization of H(k) and compare periods and orbit averages."""
import numpy as np, sys
sys.path.insert(0, '.')
import tb, spline, chambers as C

class Exact:
    def __init__(self, m, group):
        self.m, self.g = m, group
        self.Binv = np.linalg.inv(m.B)
    def eig_vel(self, kf):
        return self.m.bands_vel(kf, self.g)

name, mu = sys.argv[1], float(sys.argv[2])
GROUP = {'ZrSiS': [11, 12], 'ZrSiSe': [11, 12], 'ZrSiTe': [11, 12, 13]}[name]
z = np.load(f'/data/spline/{name}.npz')
sp = spline.GroupSpline(z['coef'], z['B'], float(z['eref']))
m = tb.Model(f'/data/model/{name}.npz', name)
ex = Exact(m, GROUP)
rng = np.random.default_rng(3)
K = rng.random((30000, 3))
e, _ = sp.eig_vel(K)
ks, bs = np.nonzero(np.abs(e - mu) < 0.01)
pick = rng.choice(len(ks), 24, replace=False)
X = np.array([0, 1, 5, 10])
for th in (0, 45, 90):
    b = [0, np.sin(np.radians(th)), np.cos(np.radians(th))]
    ta = C.integrate(sp, K[ks[pick]], bs[pick], b, 120)
    tb_ = C.integrate(ex, K[ks[pick]], bs[pick], b, 120)
    aa, pa = C.orbit_average(ta, X); ab, pb = C.orbit_average(tb_, X)
    both = ta['closed'] & tb_['closed']
    rel = np.abs(pa - pb)[both] / pb[both]
    tr_a = np.trace(aa, axis1=2, axis2=3); tr_b = np.trace(ab, axis1=2, axis2=3)
    d = np.abs(tr_a - tr_b) / np.maximum(np.abs(tr_b[:, :1]), 1e-9)
    print(f'theta {th}: closed spline {ta["closed"].sum()} exact {tb_["closed"].sum()} same-flag {(ta["closed"]==tb_["closed"]).mean():.2f};'
          f' period rel err median {np.median(rel):.4f} max {rel.max():.4f}; trace<v vbar> rel err (to X=0) median {np.median(d[:,1:]):.4f} 90% {np.percentile(d[:,1:],90):.4f}')
