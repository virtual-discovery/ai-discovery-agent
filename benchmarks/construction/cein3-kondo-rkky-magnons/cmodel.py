"""CeIn3 magnon forward model (numpy/scipy only).

LSWT of a Heisenberg antiferromagnet on the simple-cubic Bravais lattice ordered at R=(1/2,1/2,1/2):
  H = sum_<ij> J_ij S_i.S_j (each bond once, J>0 antiferromagnetic), J(q) = sum_n J(n) exp(2 pi i q.n)
  omega(q) = S sqrt[(J(q)-J(R)) (J(q+R)-J(R))],  chi''(q,w) ~ sqrt[(J(q+R)-J(R))/(J(q)-J(R))] delta(w-omega)
Exchange = short-range shells SJ_s (S*J, meV) plus an optional long-range term -A*exp(-|q-R'|^2/2kappa^2)
(R' the nearest R-equivalent point), i.e. a staggered, Gaussian-decaying real-space tail of range ~1/(2 pi kappa).
Optional gap: omega -> sqrt(omega^2 + gap^2). Each data bin is averaged over its momentum box.
"""
import numpy as np
from scipy.signal import fftconvolve

R = np.array([0.5, 0.5, 0.5])
def _shells(nmax=30):
    N = np.array([(i, j, k) for i in range(-6, 7) for j in range(-6, 7) for k in range(-6, 7) if (i, j, k) != (0, 0, 0)])
    out = []
    for r2 in sorted(set((N ** 2).sum(1))):
        vecs = N[(N ** 2).sum(1) == r2]
        groups = {}
        for v in vecs:
            groups.setdefault(tuple(sorted(np.abs(v))), []).append(v)
        for key in sorted(groups):
            out.append((int(r2), key, np.array(groups[key])))
    return out[:nmax]
SHELLS = _shells()
NSH = len(SHELLS)
EF = np.arange(-1.0, 12.0, 0.01)

def design(P):
    return np.stack([np.cos(2 * np.pi * P @ v.T).sum(-1) for (_, _, v) in SHELLS], -1)
CR = design(R[None])[0]
def d2R(P):
    Pr = P - np.floor(P)
    return ((Pr - 0.5) ** 2).sum(1)
def unit(v):
    v = np.array(v, float); return v / np.linalg.norm(v)
def perp(u):
    a = np.array([1., 0, 0]) if abs(u[0]) < 0.9 else np.array([0, 1., 0])
    e1 = unit(np.cross(u, a)); return e1, unit(np.cross(u, e1))

class Bins:
    """Momentum bins; each is a box (half_along along u, +-half_perp along the two perpendicular axes)."""
    def __init__(self, centers, dirs, half_along, half_perp, n_perp=14, n_along=3):
        pts = []
        for c, u in zip(centers, dirs):
            e1, e2 = perp(u)
            a = ((np.arange(n_perp) + 0.5) / n_perp * 2 - 1) * half_perp if half_perp > 0 else np.zeros(1)
            t = ((np.arange(n_along) + 0.5) / n_along * 2 - 1) * half_along
            A_, B_, T_ = np.meshgrid(a, a, t, indexing='ij')
            pts.append(c + T_.reshape(-1, 1) * u + A_.reshape(-1, 1) * e1 + B_.reshape(-1, 1) * e2)
        self.nb, self.npt = len(pts), len(pts[0])
        P = np.concatenate(pts)
        self.Ca, self.Cb = design(P) - CR, design(P + R) - CR
        self.da, self.db = d2R(P), d2R(P + R)
    def disp(self, p):
        J = np.r_[p['SJ'], np.zeros(NSH - len(p['SJ']))]
        a, b = self.Ca @ J, self.Cb @ J
        if p.get('A', 0):
            k2 = 2 * p['kappa'] ** 2
            a = a - p['A'] * (np.exp(-self.da / k2) - 1); b = b - p['A'] * (np.exp(-self.db / k2) - 1)
        a = np.clip(a, 1e-9, None); b = np.clip(b, 0, None)
        w = np.sqrt(a * b + p.get('gap', 0.0) ** 2)
        return w, np.sqrt(b / a)
    def spectra(self, p):
        w, A = self.disp(p)
        de = EF[1] - EF[0]
        x = np.clip((w - EF[0]) / de, 0, len(EF) - 1.001); k = x.astype(int); f = x - k
        base = np.repeat(np.arange(self.nb), self.npt) * len(EF)
        n = self.nb * len(EF)
        H = np.bincount(base + k, A * (1 - f), n) + np.bincount(base + k + 1, A * f, n)
        return H.reshape(self.nb, len(EF)) / (self.npt * de)

def kernel(sigma, gamma):
    half = min(4.0, 6 * sigma + 15 * gamma + 0.05)
    x = np.arange(-half, half + 1e-9, EF[1] - EF[0])
    g = np.exp(-0.5 * (x / sigma) ** 2); g /= g.sum()
    if gamma > 1e-4:
        L = gamma / np.pi / (x ** 2 + gamma ** 2) * (EF[1] - EF[0])
        return fftconvolve(L, g, mode='same')
    return g
def broaden(H, sigma, gamma):
    return fftconvolve(H, kernel(sigma, gamma)[None, :], mode='same', axes=1)
def energy_bins(Hc, Es, ebin):
    c = np.concatenate([np.zeros((Hc.shape[0], 1)), np.cumsum(Hc, 1)], 1)
    lo = np.searchsorted(EF, Es - ebin / 2 - 1e-9); hi = np.searchsorted(EF, Es + ebin / 2 - 1e-9)
    return (c[:, hi] - c[:, lo]) / np.maximum(hi - lo, 1)

# ---------- geometry of the archived data sets
HR_DIR = {'RG': unit([-1, -1, -1]), 'RX': unit([0, -1, -1]), 'RM': unit([0, 0, -1])}
HR_SIGMA = 0.106 / 2.3548
def he_center(path, s):
    r3 = np.sqrt(0.75)
    if path == 'RG': return R * (1 - s / r3), unit([-1, -1, -1])
    if path == 'GX': return np.array([s - r3, 0, 0]), unit([1, 0, 0])
    if path == 'XM': return np.array([0.5, s - r3 - 0.5, 0]), unit([0, 1, 0])
    f = (s - r3 - 1.0) / np.sqrt(0.5); return np.array([0.5, 0.5, 0]) * (1 - f), unit([-1, -1, 0])
HE_PERP = {'RG': 0.09, 'MG': 0.09, 'GX': 0.17, 'XM': 0.17}
def path_of(s):
    r3 = np.sqrt(0.75)
    return 'RG' if s <= r3 + 1e-9 else 'GX' if s <= r3 + 0.5 + 1e-9 else 'XM' if s <= r3 + 1.0 + 1e-9 else 'MG'

class HRRows:
    """High-resolution constant-q energy spectra along RG/RX/RM (q = signed distance from R, r.l.u.)."""
    def __init__(self, path, qs, Es, half_perp=0.08, dq=0.02, n_perp=14):
        self.B = Bins([R + t * HR_DIR[path] for t in qs], [HR_DIR[path]] * len(qs), dq / 2, half_perp, n_perp)
        self.Es = np.array(Es); self.ebin = self.Es[1] - self.Es[0]
    def model(self, p):
        return energy_bins(broaden(self.B.spectra(p), HR_SIGMA, p.get('gamma_hr', 0)), self.Es, self.ebin)
class HRCut:
    """High-resolution constant-energy cut (energy window E0 +- 0.1 meV) along RG/RX/RM."""
    def __init__(self, path, qs, E0, half_perp=0.08, dq=0.01, ewin=0.1, n_perp=14):
        self.B = Bins([R + t * HR_DIR[path] for t in qs], [HR_DIR[path]] * len(qs), dq / 2, half_perp, n_perp)
        self.E0, self.ewin = E0, ewin
    def model(self, p):
        Hc = broaden(self.B.spectra(p), HR_SIGMA, p.get('gamma_hr', 0))
        m = (EF >= self.E0 - self.ewin - 1e-9) & (EF < self.E0 + self.ewin - 1e-9)
        return Hc[:, m].mean(1)
class HERows:
    """12-meV constant-Q energy spectra at path coordinates s (R=0, Gamma=0.866, X=1.366, M=1.866, Gamma=2.573)."""
    def __init__(self, ss, Es, ds=0.025, n_perp=12, perp_scale=1.0):
        cs = [he_center(path_of(s), s) for s in ss]
        hp = [HE_PERP[path_of(s)] * perp_scale for s in ss]
        self.groups = []
        for h in sorted(set(hp)):
            idx = [i for i in range(len(ss)) if hp[i] == h]
            self.groups.append((idx, Bins([cs[i][0] for i in idx], [cs[i][1] for i in idx], ds / 2, h, n_perp)))
        self.n = len(ss); self.Es = np.array(Es); self.ebin = self.Es[2] - self.Es[1]
    def model(self, p):
        out = np.zeros((self.n, len(self.Es)))
        for idx, B in self.groups:
            out[idx] = energy_bins(broaden(B.spectra(p), p.get('sigma_he', 0.2), p.get('gamma_he', 0)), self.Es, self.ebin)
        return out
class Integrated:
    """Energy-integrated chi''(Q) at path coordinates s (12-meV binning)."""
    def __init__(self, ss, n_perp=16, perp_scale=1.0):
        self.Bs = [Bins([he_center(path_of(s), s)[0]], [he_center(path_of(s), s)[1]], 0.0125, HE_PERP[path_of(s)] * perp_scale, n_perp) for s in ss]
    def model(self, p):
        return np.array([B.spectra(p)[0].sum() * (EF[1] - EF[0]) for B in self.Bs])

def props(p):
    """velocity at R (meV per r.l.u., averaged over [111],[110],[001]), bandwidth W (max omega), stability."""
    J = np.r_[p['SJ'], np.zeros(NSH - len(p['SJ']))]
    def ab(P):
        a, b = (design(P) - CR) @ J, (design(P + R) - CR) @ J
        if p.get('A', 0):
            k2 = 2 * p['kappa'] ** 2
            a = a - p['A'] * (np.exp(-d2R(P) / k2) - 1); b = b - p['A'] * (np.exp(-d2R(P + R) / k2) - 1)
        return a, b
    vs = []
    for u in [unit([1, 1, 1]), unit([0, 1, 1]), unit([0, 0, 1])]:
        a, b = ab(np.array([R + 2e-3 * u])); vs.append(float(np.sqrt(max(a[0] * b[0], 0)) / 2e-3))
    f = np.linspace(0, 1, 24, endpoint=False)
    G = np.array(np.meshgrid(f, f, f, indexing='ij')).reshape(3, -1).T
    a, b = ab(G)
    W = float(np.sqrt(np.clip(a * b, 0, None)).max())
    return {'v': float(np.mean(vs)), 'W': W, 'eta': float(np.mean(vs)) / (2 * np.pi * W),
            'stable': bool((a > -1e-6).all() and (b > -1e-6).all())}
