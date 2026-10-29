"""Evaluators of the archived ZrSiX Wannier models used by the lab (copied unchanged into environment/lab).

- Model: exact Fourier evaluation of the archived H(R) (24 orbitals, no spin-orbit coupling, spinless):
  band energies and Hellmann-Feynman velocities dE/dk (eV Angstrom) at arbitrary k.
- orbits: semiclassical k-space orbits on one band for a field direction (RK4 in arc length, exact velocities).
- counts: Fermi-Dirac electron counts from precomputed band-energy histograms of the uniform meshes.
- transport: Chambers magnetoconductivity from the precomputed orbit tables (see make_cache.py in the
  construction folder): sigma/tau = (2 e^2/(V N_k)) sum (-df/dE) <v_i vbar_j>_orbit.
"""
import numpy as np

HBAR = 6.582119569e-16
E = 1.602176634e-19
KB = 8.617333262e-5
ME = 9.1093837015e-31
S_PER_L = HBAR ** 2 * 1e32


class Model:
    def __init__(self, path):
        z = np.load(path)
        self.R = z['R'].astype(float)
        self.H = (z['H_re'].astype(np.complex128) + 1j * z['H_im'].astype(np.complex128)).reshape(len(self.R), -1)
        self.nw = int(z['nw'])
        self.A = z['A'].astype(float)
        self.B = 2 * np.pi * np.linalg.inv(self.A).T
        self.Binv = np.linalg.inv(self.B)
        self.Rc = self.R @ self.A
        self.volume = abs(np.linalg.det(self.A))

    def bands(self, kfrac, chunk=500):
        out = []
        for s in range(0, len(kfrac), chunk):
            ph = np.exp(2j * np.pi * kfrac[s:s + chunk] @ self.R.T)
            out.append(np.linalg.eigvalsh((ph @ self.H).reshape(-1, self.nw, self.nw)))
        return np.concatenate(out)

    def bands_vel(self, kfrac, bands=None, chunk=500):
        Es, Vs = [], []
        for s in range(0, len(kfrac), chunk):
            ph = np.exp(2j * np.pi * kfrac[s:s + chunk] @ self.R.T)
            e, u = np.linalg.eigh((ph @ self.H).reshape(-1, self.nw, self.nw))
            if bands is not None:
                e, u = e[:, bands], u[:, :, bands]
            v = [np.einsum('kin,kij,kjn->kn', u.conj(), ((ph * (1j * self.Rc[:, d])) @ self.H).reshape(-1, self.nw, self.nw), u).real
                 for d in range(3)]
            Es.append(e); Vs.append(np.stack(v, -1))
        return np.concatenate(Es), np.concatenate(Vs)


class OneBand:
    """Adapter: energies and velocities of a single band (0-based index) for the orbit integrator."""
    def __init__(self, model, band):
        self.m, self.band, self.Binv = model, band, model.Binv

    def eig_vel(self, kf):
        e, v = self.m.bands_vel(kf, [self.band])
        return e, v


def trace_orbits(model, band, k0, bvec, max_length, dl=0.004):
    """Orbits through start points k0 (frac) on `band` for field direction bvec. Arc length <= max_length (1/A).
    Returns per orbit: closed flag, G (reciprocal-lattice displacement per period), period in s = B*t (T ps),
    cyclotron mass |m*|/m_e = e S / (2 pi m_e), arc length, energy, and the path (frac)."""
    b = np.asarray(bvec, float); b = b / np.linalg.norm(b)
    ad = OneBand(model, band)
    M = len(k0)
    k = np.array(k0, float)
    e0, v = ad.eig_vel(k); e0 = e0[:, 0]; v = v[:, 0]
    s = np.zeros(M); active = np.ones(M, bool); closed = np.zeros(M, bool); G = np.zeros((M, 3), int)
    paths = [[k[i].copy()] for i in range(M)]
    length = np.zeros(M)

    def f(vv):
        w = np.cross(vv, b)
        n = np.maximum(np.linalg.norm(w, axis=1), 1e-3)
        return (w / n[:, None]) @ model.Binv, n
    step = 0
    nmax = int(max_length / dl)
    while active.any() and step < nmax:
        a = np.where(active)[0]
        ka, va = k[a], v[a]
        d1, n1 = f(va)
        _, v2 = ad.eig_vel(ka + 0.5 * dl * d1); d2, n2 = f(v2[:, 0])
        _, v3 = ad.eig_vel(ka + 0.5 * dl * d2); d3, n3 = f(v3[:, 0])
        _, v4 = ad.eig_vel(ka + dl * d3); d4, n4 = f(v4[:, 0])
        kn = ka + dl * (d1 + 2 * d2 + 2 * d3 + d4) / 6
        en, vn = ad.eig_vel(kn); en, vn = en[:, 0], vn[:, 0]
        kn = kn - (((en - e0[a]) / np.maximum((vn ** 2).sum(1), 1e-6))[:, None] * vn) @ model.Binv
        _, vn = ad.eig_vel(kn); vn = vn[:, 0]
        s[a] += S_PER_L * dl * (1 / n1 + 2 / n2 + 2 / n3 + 1 / n4) / 6
        k[a], v[a] = kn, vn
        length[a] += dl
        step += 1
        for t, i in enumerate(a):
            paths[i].append(kn[t].copy())
        delta = kn - np.array(k0)[a]
        Gn = np.round(delta)
        dist = np.linalg.norm((delta - Gn) @ model.B, axis=1)
        done = (dist < 0.75 * dl) & (step > 8)
        closed[a[done]] = True; G[a[done]] = Gn[done].astype(int); active[a[done]] = False
    period = np.where(closed, s, np.nan)
    mass = period * 1e-12 * E / (2 * np.pi * ME)
    return dict(closed=closed, G=G, period=period, mass=mass, length=length, energy=e0, paths=paths, s_total=s)


def fermi_weights(e, mu, kT):
    """-df/dE (1/eV) and occupations f for energies e (array) and chemical potential mu."""
    if kT <= 0:
        return None, (e < mu).astype(float)
    x = np.clip((e - mu) / kT, -60, 60)
    f = 1.0 / (np.exp(x) + 1.0)
    return f * (1 - f) / kT, f


def rho_from_sigma(sig):
    return np.linalg.inv(sig)
