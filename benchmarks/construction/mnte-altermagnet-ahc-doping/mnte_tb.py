"""Numpy evaluator of the archived MnTe Wannier Hamiltonian (copied unchanged into the task lab).

Model: wannier90_hr.dat of Zenodo record 18863304 (dft_data/2_wannier90), 32 spinor Wannier functions
(Mn d x 2 atoms, Te p x 2 atoms; spin-major order: functions 1-16 spin up along z, 17-32 spin down), stored as
npz (R in lattice coordinates, Wigner-Seitz degeneracies, H(R) in eV; the lab copy stores H(R) as integer micro-eV,
which is exact for the 6-decimal hr.dat).

Conventions (those of WannierTools 2.6.1, sigma_AHC, used by the authors):
  H(k) = sum_R exp(i k.R) H(R) / deg(R), k Cartesian in 1/Angstrom (2 pi included);
  velocity in the atomic gauge: Wannier centres at the atomic positions of the projections (wt.in);
  Omega^c_n = -2 Im sum_{m != n} V^a_nm V^b_mn / (E_n - E_m)^2 (Omega = curl_k i<u|grad_k u>), in Angstrom^2,
  c = x (a, b = y, z), y (z, x), z (x, y); pairs with |E_n - E_m| <= 1e-12 eV are skipped;
  sigma_ab = +(e^2/hbar) (1/(N_k V)) sum_k sum_n f(E_n - mu) Omega^c_n  in S/cm (sign of the archived output),
  Gamma-centred N1 x N2 x N3 mesh of fractional coordinates, Fermi-Dirac occupations f with temperature kT.
Neel reversal: the time-reversed Hamiltonian H'(R) = (i s_y) H(R)* (i s_y)^dagger (moments -> -moments).
"""
import numpy as np

CELL = np.array([[4.14, 0.0, 0.0], [-2.07, 3.585345171667576, 0.0], [0.0, 0.0, 6.71]])
B = 2 * np.pi * np.linalg.inv(CELL).T
VOL = abs(np.linalg.det(CELL))
POS = np.array([[0, 0, 0], [0, 0, 0.5], [1/3, 2/3, 0.25], [2/3, 1/3, 0.75]])
SITE = np.array([0] * 5 + [1] * 5 + [2] * 3 + [3] * 3)
SITE = np.r_[SITE, SITE]
TAU = POS[SITE] @ CELL
E2HBAR = 1.602176634e-19 ** 2 / 1.054571817e-34
EMIN, EMAX, DE = 4.80, 8.80, 1e-4          # histogram window (eV, absolute scale of the hr) and bin width
NBIN = int(round((EMAX - EMIN) / DE))
PAIRS = {0: (1, 2), 1: (2, 0), 2: (0, 1)}   # c = 0: sigma_yz (Omega_x), 1: sigma_zx (Omega_y), 2: sigma_xy (Omega_z)
NW = 32
SY = np.array([[0, -1j], [1j, 0]])
TR = np.kron(1j * SY, np.eye(16))           # time reversal on the spin-major spinor basis (with complex conjugation)
SPIN = [np.kron(s, np.eye(16)) for s in (np.array([[0, 1], [1, 0]], complex), SY, np.array([[1, 0], [0, -1]], complex))]


def load(path, neel=+1):
    z = np.load(path)
    H = z['H'] if 'H' in z else (z['H_re_ueV'] + 1j * z['H_im_ueV']) * 1e-6   # compact lab copy: integer micro-eV
    R, H = z['R'].astype(int), H / z['deg'][:, None, None]
    if neel < 0:
        H = np.einsum('ij,rjk,lk->ril', TR, H.conj(), TR.conj())
    return R, H


class Model:
    """Direct Fourier sums at arbitrary k-points (for k-lists)."""

    def __init__(self, path, neel=+1):
        R, H = load(path, neel)
        self.Rc = R @ CELL
        self.Hflat = H.reshape(len(R), -1)

    def hd(self, kc):
        kc = np.atleast_2d(np.asarray(kc, float))
        ph = np.exp(1j * kc @ self.Rc.T)
        out = np.empty((len(kc), 4, NW, NW), complex)
        out[:, 0] = (ph @ self.Hflat).reshape(-1, NW, NW)
        for a in range(3):
            out[:, 1 + a] = ((1j * self.Rc[:, a] * ph) @ self.Hflat).reshape(-1, NW, NW)
        return out

    def bands(self, kc, chunk=1000):
        kc = np.atleast_2d(np.asarray(kc, float))
        out = []
        for s in range(0, len(kc), chunk):
            ph = np.exp(1j * kc[s:s + chunk] @ self.Rc.T)
            out.append(np.linalg.eigvalsh((ph @ self.Hflat).reshape(-1, NW, NW)))
        return np.concatenate(out)

    def spin(self, kc):
        """Band energies and spin expectation values <s_x>, <s_y>, <s_z> (Pauli matrices, range -1..1)."""
        kc = np.atleast_2d(np.asarray(kc, float))
        ph = np.exp(1j * kc @ self.Rc.T)
        E, U = np.linalg.eigh((ph @ self.Hflat).reshape(-1, NW, NW))
        s = np.array([np.real(np.einsum('kmn,mp,kpn->kn', U.conj(), S, U)) for S in SPIN])
        return E, s

    def curvature(self, kc, chunk=500):
        kc = np.atleast_2d(np.asarray(kc, float))
        Es, Os = [], []
        for s in range(0, len(kc), chunk):
            E, om = curvature_from(self.hd(kc[s:s + chunk]))
            Es.append(E); Os.append(om)
        return np.concatenate(Es), np.concatenate(Os, axis=1)


def curvature_from(Hd, gauge='atomic', comps=(0, 1, 2)):
    """Hd: (nk, 4, nw, nw) = H, dH/dkx, dH/dky, dH/dkz in the lattice gauge. Returns E (nk, nw), Omega (nc, nk, nw)."""
    E, U = np.linalg.eigh(Hd[:, 0])
    Ud = np.conj(np.transpose(U, (0, 2, 1)))
    dE = E[:, :, None] - E[:, None, :]          # E_n - E_m
    need = sorted({a for c in comps for a in PAIRS[c]})
    V = {}
    for a in need:
        Va = Ud @ Hd[:, 1 + a] @ U
        if gauge == 'atomic':                   # V^a_nm = U^dag (dH + i [H, tau_a]) U = V^l_nm + i (E_n - E_m) tau_nm
            Va = Va + 1j * dE * (Ud @ (TAU[:, a][None, :, None] * U))
        V[a] = Va
    den = dE ** 2
    den[np.abs(dE) <= 1e-12] = np.inf
    om = [-2 * np.imag(np.sum(V[PAIRS[c][0]] * np.transpose(V[PAIRS[c][1]], (0, 2, 1)) / den, axis=2)) for c in comps]
    return E, np.array(om)


def fermi(x, kT):
    if kT <= 0:
        return (x < 0).astype(float)
    return 0.5 * (1 - np.tanh(np.clip(x / (2 * kT), -350, 350)))


def bin_centres():
    return EMIN + (np.arange(NBIN) + 0.5) * DE


def sigma_from_cache(c, mu, kT):
    """sigma_yz, sigma_zx, sigma_xy (S/cm) at chemical potentials mu (absolute eV) from a mesh histogram cache."""
    mu = np.atleast_1d(np.asarray(mu, float))
    f = fermi(bin_centres()[None, :] - mu[:, None], kT)
    tot = c['below'][:, None] + c['hist'] @ f.T
    return tot / float(c['nk']) / VOL * E2HBAR * 1e8


def electrons_from_cache(c, mu, kT):
    """Electrons per cell in the 32 Wannier bands at chemical potentials mu."""
    mu = np.atleast_1d(np.asarray(mu, float))
    f = fermi(bin_centres()[None, :] - mu[:, None], kT)
    return (float(c['nlow']) + f @ c['count']) / float(c['nk'])


# ---------------------------------------------------------------- mesh engine (construction only)
class GridModel:
    """Separable Fourier transforms for Gamma-centred uniform meshes."""

    def __init__(self, path, neel=+1):
        R, H = load(path, neel)
        self.lo = R.min(0)
        n = R.max(0) - self.lo + 1
        Rc = R @ CELL
        T = np.zeros((n[0], n[1], n[2], 4, NW * NW), complex)
        for r, rc, h in zip(R, Rc, H):
            i = tuple(r - self.lo)
            T[i][0] = h.ravel()
            for a in range(3):
                T[i][1 + a] = 1j * rc[a] * h.ravel()
        self.T = T.reshape(n[0], n[1], n[2], -1)
        self.r = [np.arange(self.lo[d], self.lo[d] + n[d]) for d in range(3)]

    def plane(self, k3):
        return np.tensordot(self.T, np.exp(2j * np.pi * k3 * self.r[2]), axes=([2], [0]))

    def row(self, P, k1, k2s):
        Tr = np.tensordot(np.exp(2j * np.pi * k1 * self.r[0]), P, axes=([0], [0]))
        ph = np.exp(2j * np.pi * np.outer(k2s, self.r[1]))
        return (ph @ Tr).reshape(len(k2s), 4, NW, NW)


def _accumulate(E, om, acc):
    idx = np.floor((E - EMIN) / DE).astype(int)
    low = idx < 0
    keep = (~low) & (idx < NBIN)
    for c in range(om.shape[0]):
        acc['below'][c] += om[c][low].sum()
        acc['hist'][c] += np.bincount(idx[keep], weights=om[c][keep], minlength=NBIN)
    acc['count'] += np.bincount(idx[keep], minlength=NBIN)
    acc['nlow'] += int(low.sum())


def _mesh_worker(args):
    path, neel, N, planes = args
    m = GridModel(path, neel)
    acc = {'hist': np.zeros((3, NBIN)), 'below': np.zeros(3), 'count': np.zeros(NBIN), 'nlow': 0}
    k2s = np.arange(N[1]) / N[1]
    for l in planes:
        P = m.plane(l / N[2])
        for i in range(N[0]):
            E, om = curvature_from(m.row(P, i / N[0], k2s))
            _accumulate(E, om, acc)
    return acc


def mesh_cache(path, N, neel=+1, workers=4):
    from multiprocessing import Pool
    jobs = [(path, neel, N, list(range(w, N[2], workers))) for w in range(workers)]
    jobs = [j for j in jobs if j[3]]
    with Pool(len(jobs)) as pool:
        res = pool.map(_mesh_worker, jobs)
    return {'N': np.array(N), 'nk': np.array(int(np.prod(N))), 'hist': sum(r['hist'] for r in res),
            'below': sum(r['below'] for r in res), 'count': sum(r['count'] for r in res),
            'nlow': np.array(sum(r['nlow'] for r in res))}
