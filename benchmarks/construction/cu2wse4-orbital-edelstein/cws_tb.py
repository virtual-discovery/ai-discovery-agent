"""Numpy evaluator of the archived Cu2WSe4 OpenMX Wannier Hamiltonians (copied unchanged into the task lab).

Models (Zenodo 20131063, CWS_PRM.zip): Cu2WSe4_wosoc.HWR (54 Wannier functions, no spin-orbit coupling) and
Cu2WSe4_soc.HWR (108 spinor Wannier functions, fully relativistic). Basis order: Cu0-Cu3 d (5 each), W4, W5 d
(5 each), Se6-Se13 p (3 each) in the atom order of Fig1a.cif; OpenMX real harmonics d = (z2, x2-y2, xy, xz, yz),
p = (x, y, z). The SOC model holds spin up (functions 1-54) then spin down (55-108), quantised along z.
Stored as npz (R in lattice coordinates of the conventional cell, Wigner-Seitz degeneracy, H(R) in Hartree as in
the HWR file; the lab copy stores H(R) as integers in units of 1e-9 Hartree, exact to the file's last digit
for |H| > 1e-9 Hartree).

Conventions:
  H(k) = sum_R exp(i k.(R + tau_n - tau_m)) H_mn(R) / deg(R) (atomic gauge: Wannier centres at the atomic
  positions), k Cartesian in 1/Angstrom with 2 pi; energies in eV (1 Hartree = 27.211386245988 eV).
  V^a = dH/dk_a in that gauge (eV Angstrom); band velocity v_n = V^a_nn / hbar.
  Berry curvature Omega^c_n = -2 Im sum_{m != n} V^a_nm V^b_mn / (E_n - E_m)^2 (Angstrom^2), (a, b, c) cyclic,
    i.e. Omega = curl_k i<u|grad_k u>.
  Orbital moment (modern theory, self-rotation): m^c_n = (|e| / 2 hbar) X^c_n with
    X^c_n = Im eps_abc <d_a n|(H - E_n)|d_b n> = -2 Im sum_{m != n} V^a_nm V^b_mn / (E_n - E_m) (eV Angstrom^2),
    in Bohr magnetons m/mu_B = X / (hbar^2 / m_e) = X / 7.61996 eV A^2.
  Spin moment m^S = -(g mu_B / 2) <sigma>, g = 2.
  Atomic (intra-atomic) orbital moment m^L = -mu_B <L> with L the angular momentum of the d and p orbitals
    about their own atom.
  Pairs with |E_n - E_m| < DEGEN (5 meV) are treated as degenerate and skipped in the band sums (the off-grid
  interpolation splits symmetry-enforced degeneracies by up to ~3 meV).
"""
import numpy as np

HARTREE = 27.211386245988
BOHR = 0.529177210903
HB2M = 7.619964231        # hbar^2 / m_e in eV Angstrom^2
DEGEN = float(__import__("os").environ.get("CWS_DEGEN", "5e-3"))   # eV
A_LAT, C_LAT = 5.56010564, 11.21352958
CELL = np.diag([A_LAT, A_LAT, C_LAT])
B = 2 * np.pi * np.linalg.inv(CELL).T
VOL = float(np.linalg.det(CELL))
FRAC = np.array([[0.5, 0, 0], [0, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5],                     # Cu0-Cu3
                 [0, 0, 0], [0.5, 0.5, 0.5],                                                  # W4, W5
                 [0.24223226, 0.75776774, 0.12510258], [0.75776774, 0.24223226, 0.12510258],  # Se6-Se13
                 [0.25776773, 0.25776773, 0.37489743], [0.74223227, 0.74223227, 0.37489743],
                 [0.74223227, 0.25776773, 0.62510257], [0.25776773, 0.74223227, 0.62510257],
                 [0.75776774, 0.75776774, 0.87489743], [0.24223226, 0.24223226, 0.87489743]])
SPECIES = ['Cu'] * 4 + ['W'] * 2 + ['Se'] * 8
NORB = {'Cu': 5, 'W': 5, 'Se': 3}
ATOM = np.concatenate([[i] * NORB[s] for i, s in enumerate(SPECIES)])
NW0 = len(ATOM)                                     # 54 orbitals per spin
TAU0 = FRAC[ATOM] @ CELL


def _lmat_p():
    eps = np.zeros((3, 3, 3))
    for a, b, c in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        eps[a, b, c], eps[a, c, b] = 1, -1
    return np.array([-1j * eps[k] for k in range(3)])      # (L_k)_ij = -i eps_kij in the (x, y, z) basis


def _lmat_d():
    """L in the OpenMX real d basis (z2, x2-y2, xy, xz, yz) from the action of -i r x grad on quadratic forms."""
    s3 = 2 * np.sqrt(3)
    Q = [np.diag([-1, -1, 2]) / s3, np.diag([0.5, -0.5, 0]),
         np.array([[0, .5, 0], [.5, 0, 0], [0, 0, 0]]), np.array([[0, 0, .5], [0, 0, 0], [.5, 0, 0]]),
         np.array([[0, 0, 0], [0, 0, .5], [0, .5, 0]])]
    eps = np.zeros((3, 3, 3))
    for a, b, c in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        eps[a, b, c], eps[a, c, b] = 1, -1
    L = np.zeros((3, 5, 5), complex)
    for k in range(3):
        for j, q in enumerate(Q):
            # f = r.Q.r; (r x grad f)_k = eps_kab r_a 2 (Q r)_b = r.(2 eps_k Q).r (symmetrised)
            m = 2 * eps[k] @ q
            m = -1j * (m + m.T) / 2
            for i, p in enumerate(Q):
                L[k, i, j] = np.trace(p @ m) / np.trace(p @ p)
    return L


LP, LD = _lmat_p(), _lmat_d()


def lmat0():
    """Atomic L_x, L_y, L_z on the 54 orbitals of one spin (block diagonal per atom)."""
    L = np.zeros((3, NW0, NW0), complex)
    i = 0
    for s in SPECIES:
        n = NORB[s]
        L[:, i:i + n, i:i + n] = LD if n == 5 else LP
        i += n
    return L


def load(path):
    z = np.load(path)
    H = z['H'] if 'H' in z else (z['H_re'] + 1j * z['H_im']) * 1e-9
    return z['R'].astype(int), H * HARTREE / z['deg'][:, None, None]


class Model:
    """Bloch Hamiltonian, its k-derivatives and band-resolved properties at arbitrary k-points.

    H and dH/dk are built in the lattice gauge (one matrix product over the 147 lattice vectors); the
    atomic-gauge velocity follows as V^a = U^+ (dH^l/dk_a + i [H^l, tau_a]) U (properties_from). Energies,
    spin and atomic L expectation values do not depend on the gauge.
    """

    def __init__(self, path):
        R, H = load(path)
        self.nw = H.shape[1]
        self.soc = self.nw == 2 * NW0
        ns = 2 if self.soc else 1
        self.tau = np.tile(TAU0, (ns, 1))
        self.Rc = R @ CELL
        self.Hflat = H.reshape(len(R), -1)
        L0 = lmat0()
        if self.soc:
            self.L = np.array([np.kron(np.eye(2), l) for l in L0])
            sx, sy, sz = np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])
            self.S = np.array([np.kron(s, np.eye(NW0)) for s in (sx, sy, sz)])
        else:
            self.L, self.S = L0, None

    def hd(self, kc):
        """Lattice-gauge H(k) and dH/dk_a (eV, eV Angstrom) for Cartesian k-points (nk, 3)."""
        kc = np.atleast_2d(np.asarray(kc, float))
        nw = self.nw
        ph = np.exp(1j * kc @ self.Rc.T)
        out = np.empty((len(kc), 4, nw, nw), complex)
        out[:, 0] = (ph @ self.Hflat).reshape(-1, nw, nw)
        for a in range(3):
            out[:, 1 + a] = ((1j * self.Rc[:, a] * ph) @ self.Hflat).reshape(-1, nw, nw)
        return out

    def bands(self, kc, chunk=200):
        kc = np.atleast_2d(np.asarray(kc, float))
        out = []
        for s in range(0, len(kc), chunk):
            ph = np.exp(1j * kc[s:s + chunk] @ self.Rc.T)
            out.append(np.linalg.eigvalsh((ph @ self.Hflat).reshape(-1, self.nw, self.nw)))
        return np.concatenate(out)

    def properties(self, kc, chunk=100):
        """Band energies and band-resolved quantities at Cartesian k-points; see properties_from."""
        kc = np.atleast_2d(np.asarray(kc, float))
        parts = [properties_from(self.hd(kc[s:s + chunk]), self.tau, self.L, self.S) for s in range(0, len(kc), chunk)]
        out = {}
        for key in parts[0]:
            if parts[0][key] is None:
                out[key] = None
            else:
                out[key] = np.concatenate([p[key] for p in parts], axis=0 if key == 'E' else 1)
        return out


def properties_from(Hd, tau, L, S):
    """Hd: (nk, 4, nw, nw) lattice-gauge H and dH/dk; tau: (nw, 3) Wannier centres. Returns dict of E (nk, nw), v (3, nk, nw) eV A, omega (3, nk, nw) A^2,
    m_orb (3, nk, nw) mu_B, m_atomic (3, nk, nw) mu_B, m_spin (3, nk, nw) mu_B (None without SOC)."""
    E, U = np.linalg.eigh(Hd[:, 0])
    Ud = np.conj(np.transpose(U, (0, 2, 1)))
    dE = E[:, :, None] - E[:, None, :]                      # E_n - E_m
    V = [Ud @ Hd[:, 1 + a] @ U + 1j * dE * (Ud @ (tau[:, a][None, :, None] * U)) for a in range(3)]
    inv = np.where(np.abs(dE) > DEGEN, 1.0 / np.where(np.abs(dE) > DEGEN, dE, 1.0), 0.0)
    om, X = [], []
    for c, (a, b) in enumerate([(1, 2), (2, 0), (0, 1)]):
        P = V[a] * np.transpose(V[b], (0, 2, 1))           # V^a_nm V^b_mn
        om.append(-2 * np.imag(np.sum(P * inv ** 2, axis=2)))
        X.append(-2 * np.imag(np.sum(P * inv, axis=2)))
    out = {'E': E, 'v': np.array([np.real(np.diagonal(V[a], axis1=1, axis2=2)) for a in range(3)]),
           'omega': np.array(om), 'm_orb': np.array(X) / HB2M,
           'm_atomic': -np.array([np.real(np.sum(U.conj() * (l @ U), axis=1)) for l in L])}
    out['m_spin'] = (-np.array([np.real(np.sum(U.conj() * (s @ U), axis=1)) for s in S]) if S is not None else None)
    return out


# ---------------------------------------------------------------- Fermi-surface integrals from histograms
# Every mesh cache holds, per energy bin (absolute eV on the scale of the HWR file), the sum over states of
#   count, m^S_i v_j, m^O_i v_j, m^L_i v_j (mu_B eV Angstrom) and v_a Omega_b (eV Angstrom^3), i, j, a, b = x, y, z,
# with v = dE/dk (eV Angstrom). A response at chemical potential mu and temperature kT is the bin sum weighted by
# -df/dE, divided by N_k V_cell.
QUANTS = ['spin', 'orb', 'atomic', 'bcd']
EBIN = 1e-3                                  # eV


def fermi_window(centres, mu, kT):
    """-df/dE (1/eV) at bin centres for chemical potentials mu (nmu,)."""
    x = (centres[None, :] - np.asarray(mu, float)[:, None]) / kT
    x = np.clip(x, -300, 300)
    return 1.0 / (kT * (2 * np.cosh(x / 2)) ** 2)


def surface_integral(cache, quant, mu, kT):
    """(nmu, 3, 3) array: (1/(N_k V)) sum_nk (-df/dE) q_ij, q = m_i v_j (mu_B / Angstrom^2) or v_a Omega_b (Angstrom^0)."""
    w = fermi_window(cache['centres'], mu, kT)            # (nmu, nbin)
    h = cache[quant].reshape(len(cache['centres']), 9)      # (nbin, 9)
    return (w @ h).reshape(-1, 3, 3) / (float(cache["nk"]) * VOL)


def density_of_states(cache, mu, kT):
    """States per eV per Angstrom^3 (Fermi-window weighted) and electrons per cell above the window bottom."""
    w = fermi_window(cache['centres'], mu, kT)
    return (w @ cache['count']) / (float(cache['nk']) * VOL)


# ---------------------------------------------------------------- mesh engine (construction only)
class GridModel:
    """Separable Fourier sums on Gamma-centred uniform meshes (the lattice vectors fill a 7 x 7 x 3 box)."""

    def __init__(self, path):
        R, H = load(path)
        self.base = Model(path)
        self.lo = R.min(0)
        n = R.max(0) - self.lo + 1
        nw = H.shape[1]
        Rc = R @ CELL
        T4 = np.zeros((n[0], n[1], n[2], 4, nw * nw), complex)
        for r, rc, h in zip(R, Rc, H):
            i = tuple(r - self.lo)
            T4[i][0] = h.ravel()
            for a in range(3):
                T4[i][1 + a] = 1j * rc[a] * h.ravel()
        self.T = T4.reshape(n[0], n[1], n[2], -1)
        self.r = [np.arange(self.lo[d], self.lo[d] + n[d]) for d in range(3)]
        self.nw = nw

    def plane(self, k3):
        return np.tensordot(self.T, np.exp(2j * np.pi * k3 * self.r[2]), axes=([2], [0]))

    def row(self, P, k1, k2s):
        Tr = np.tensordot(np.exp(2j * np.pi * k1 * self.r[0]), P, axes=([0], [0]))
        ph = np.exp(2j * np.pi * np.outer(k2s, self.r[1]))
        return (ph @ Tr).reshape(len(k2s), 4, self.nw, self.nw)


def _accumulate(p, weight, acc, emin, nbin):
    E = p['E']
    idx = np.floor((E - emin) / EBIN).astype(int)
    keep = (idx >= 0) & (idx < nbin)
    ii = idx[keep]
    acc['count'] += weight * np.bincount(ii, minlength=nbin)
    acc['nlow'] += weight * float((idx < 0).sum())
    v = p['v']
    for q, key in (('spin', 'm_spin'), ('orb', 'm_orb'), ('atomic', 'm_atomic'), ('bcd', 'omega')):
        if p[key] is None:
            continue
        for i in range(3):
            for j in range(3):
                val = (v[i] * p[key][j]) if q == 'bcd' else (p[key][i] * v[j])
                acc[q][:, 3 * i + j] += weight * np.bincount(ii, weights=val[keep], minlength=nbin)


def _mesh_worker(args):
    path, N, planes, emin, nbin = args
    g = GridModel(path)
    m = g.base
    acc = {'count': np.zeros(nbin), 'nlow': 0.0, **{q: np.zeros((nbin, 9)) for q in QUANTS}}
    k2s = np.arange(N[1]) / N[1]
    for l, weight in planes:
        P = g.plane(l / N[2])
        for i in range(N[0]):
            p = properties_from(g.row(P, i / N[0], k2s), m.tau, m.L, m.S)
            _accumulate(p, weight, acc, emin, nbin)
    return acc


def mesh_cache(path, N, emin, emax, workers=4):
    """Histogram cache on the Gamma-centred N1 x N2 x N3 mesh. Planes l and N3 - l are related by k -> -k, under
    which every accumulated quantity is even (time reversal), so only l = 0 .. N3/2 are computed."""
    from multiprocessing import Pool
    nbin = int(round((emax - emin) / EBIN))
    half = [(l, 1.0 if (l == 0 or 2 * l == N[2]) else 2.0) for l in range(N[2] // 2 + 1)]
    # plane l = 0 (and N3/2) maps onto itself under k -> -k; counted once with weight 1 (exact)
    jobs = [(path, N, half[w::workers], emin, nbin) for w in range(workers)]
    jobs = [j for j in jobs if j[2]]
    with Pool(len(jobs)) as pool:
        res = pool.map(_mesh_worker, jobs)
    out = {'N': np.array(N), 'nk': np.array(int(np.prod(N))), 'emin': np.array(emin),
           'centres': emin + (np.arange(nbin) + 0.5) * EBIN}
    for key in ['count', 'nlow'] + QUANTS:
        out[key] = np.array(sum(r[key] for r in res))
    return out
