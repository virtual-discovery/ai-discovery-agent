"""Numpy phonon Boltzmann transport for theta-TaN (2-atom hexagonal cell), ShengBTE/FourPhonon conventions.

Shared by the construction scripts and the task lab (copied verbatim into environment/lab/).

Conventions (as in ShengBTE, Li et al. CPC 185, 1747 (2014), and FourPhonon, Han et al. CPC 270, 108179 (2022)):
- Gamma-centred n1 x n2 x n3 grid; q in reduced reciprocal coordinates; phases exp(i q.R) on cell vectors.
- Second-order IFCs in phonopy full format (eV/A^2) for a 5x5x5 supercell, minimum-image (Wigner-Seitz) weights.
- Third-order IFCs in ShengBTE format (eV/A^3), fourth-order in FourPhonon format (eV/A^4).
- Adaptive Gaussian broadening sigma = scalebroad * sqrt(sum_nu (b_nu . dv / n_nu)^2 / 6), processes kept
  for |d omega| <= 2 sigma.
- 1/tau_3 = (1/N) [sum_+ G+ + 1/2 sum_- G-],  G = (hbar pi / 4) pop / (w w' w'') |V3|^2 delta.
- 1/tau_4 = (1/N^2) sum [1/2 G++ + 1/2 G+- + 1/6 G--],  G = (hbar^2 pi / 8) pop / (w w' w'' w''') |V4|^2 delta.
- Iterative BTE for three-phonon scattering (ShengBTE Delta term); four-phonon, isotope and extra
  (phonon-electron) rates enter the diagonal only, as in FourPhonon.
Angular frequencies are rad/ps in all public outputs; rates are ps^-1.
"""
import itertools
import numpy as np

HBAR = 1.054571817e-34
KB = 1.380649e-23
EV = 1.602176634e-19
AMU = 1.66053906660e-27
ANG = 1e-10
MASSES = {'Ta': 0.00012 * 179.9474648 + 0.99988 * 180.9479958, 'N': 0.99636 * 14.0030740048 + 0.00364 * 15.0001088982}
ISO_G = {'Ta': 0.00012 * (1 - 179.9474648 / MASSES['Ta']) ** 2 + 0.99988 * (1 - 180.9479958 / MASSES['Ta']) ** 2,
         'N': 0.99636 * (1 - 14.0030740048 / MASSES['N']) ** 2 + 0.00364 * (1 - 15.0001088982 / MASSES['N']) ** 2}
W_SI = np.sqrt(EV / ANG ** 2 / AMU)          # sqrt(eV/A^2/amu) in rad/s


# ----------------------------------------------------------------------------------------------- input
def read_poscar(path):
    lines = open(path).read().split('\n')
    scale = float(lines[1])
    A = np.array([[float(x) for x in lines[i].split()] for i in (2, 3, 4)]) * scale
    names = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    species = [n for n, c in zip(names, counts) for _ in range(c)]
    pos = np.array([[float(x) for x in lines[8 + i].split()[:3]] for i in range(sum(counts))])
    return {'A': A, 'species': species, 'frac': pos}


def read_fc2(path, nat, scell):
    tok = open(path).read().split()
    if tok[1] != '1':           # "N N" header (newer phonopy) vs "N" header (older)
        ntot, pos = int(tok[0]), 2
    else:
        ntot, pos = int(tok[0]), 1
    assert ntot == nat * int(np.prod(scell))
    fc = np.zeros((nat, ntot, 3, 3))
    for _ in range(ntot * ntot):
        i, j = int(tok[pos]) - 1, int(tok[pos + 1]) - 1
        block = np.array(tok[pos + 2:pos + 11], float).reshape(3, 3)
        pos += 11
        ncell = ntot // nat
        if i % ncell == 0:
            fc[i // ncell, j] = block
    return fc


def read_fcn(path, order):
    """ShengBTE FORCE_CONSTANTS_3RD (order 3) or FourPhonon FORCE_CONSTANTS_4TH (order 4)."""
    tok = open(path).read().split()
    n = int(tok[0])
    pos = 1
    R, atoms, vals = [], [], []
    ncomp = 3 ** order
    for _ in range(n):
        pos += 1
        r = np.array(tok[pos:pos + 3 * (order - 1)], float).reshape(order - 1, 3)
        pos += 3 * (order - 1)
        at = [int(x) - 1 for x in tok[pos:pos + order]]
        pos += order
        block = np.zeros((3,) * order)
        for _ in range(ncomp):
            idx = tuple(int(x) - 1 for x in tok[pos:pos + order])
            block[idx] = float(tok[pos + order])
            pos += order + 1
        R.append(r)
        atoms.append(at)
        vals.append(block)
    return {'R': np.array(R), 'atoms': np.array(atoms), 'fc': np.array(vals)}


# ----------------------------------------------------------------------------------------------- harmonic
class Crystal:
    def __init__(self, poscar, fc2, scell=(5, 5, 5), fc3=None, fc4=None):
        self.A = poscar['A']                       # rows: lattice vectors (A)
        self.species = poscar['species']
        self.frac = poscar['frac']
        self.nat = len(self.species)
        self.nb = 3 * self.nat
        self.mass = np.array([MASSES[s] for s in self.species])
        self.B = 2 * np.pi * np.linalg.inv(self.A).T          # rows: reciprocal vectors (1/A, with 2 pi)
        self.volume = abs(np.linalg.det(self.A))
        self.scell = np.array(scell)
        self._build_fc2(fc2)
        self.fc3 = self._reduce(fc3, 3) if fc3 is not None else None
        self.fc4 = self._reduce(fc4, 4) if fc4 is not None else None

    def _build_fc2(self, fc2):
        n1, n2, n3 = self.scell
        ncell = n1 * n2 * n3
        cells = np.array([(x, y, z) for z in range(n3) for y in range(n2) for x in range(n1)])   # phonopy: x fastest
        shifts = np.array(list(itertools.product(range(-2, 3), repeat=3))) * self.scell
        terms_R, terms_ij, terms_fc = [], [], []
        for i in range(self.nat):
            for s in range(self.nat * ncell):
                j, c = divmod(s, ncell)
                if not np.any(fc2[i, s]):
                    continue
                cand = cells[c] + shifts                           # cell vector of atom s and its images
                d = (cand + self.frac[j] - self.frac[i]) @ self.A
                dist = np.linalg.norm(d, axis=1)
                ok = np.abs(dist - dist.min()) < 1e-4 * max(1.0, dist.min())
                for R in cand[ok]:
                    terms_R.append(R)
                    terms_ij.append((i, j))
                    terms_fc.append(fc2[i, s] / ok.sum())
        self.fc2_R = np.array(terms_R, float)
        self.fc2_ij = np.array(terms_ij)
        mm = np.sqrt(self.mass[self.fc2_ij[:, 0]] * self.mass[self.fc2_ij[:, 1]])
        self.fc2_val = np.array(terms_fc) / mm[:, None, None]
        self.fc2_cart = self.fc2_R @ self.A

    def _reduce(self, fcn, order):
        Ainv = np.linalg.inv(self.A)
        n = np.rint(fcn['R'] @ Ainv).astype(int)               # (nblock, order-1, 3) integer cell vectors
        assert np.allclose(n @ self.A, fcn['R'], atol=1e-4)
        m = np.prod(np.sqrt(self.mass[fcn['atoms']]), axis=1)
        return {'n': n, 'atoms': fcn['atoms'], 'fc': fcn['fc'] / m.reshape((-1,) + (1,) * order)}

    def dynmat(self, q, deriv=False):
        """q: (nq,3) reduced. Returns D (nq,nb,nb) and optionally dD/dq_cart (nq,3,nb,nb) in (eV/A^2/amu)*A."""
        q = np.atleast_2d(q)
        ph = np.exp(2j * np.pi * q @ self.fc2_R.T)               # (nq, nterm)
        nq, nat = len(q), self.nat
        D = np.zeros((nq, nat, 3, nat, 3), complex)
        dD = np.zeros((nq, 3, nat, 3, nat, 3), complex) if deriv else None
        for i in range(nat):
            for j in range(nat):
                sel = (self.fc2_ij[:, 0] == i) & (self.fc2_ij[:, 1] == j)
                D[:, i, :, j, :] = np.einsum('qt,tab->qab', ph[:, sel], self.fc2_val[sel])
                if deriv:
                    dD[:, :, i, :, j, :] = np.einsum('qt,tc,tab->qcab', ph[:, sel], 1j * self.fc2_cart[sel],
                                                     self.fc2_val[sel])
        D = D.reshape(nq, self.nb, self.nb)
        D = 0.5 * (D + np.conj(np.transpose(D, (0, 2, 1))))
        if deriv:
            return D, dD.reshape(nq, 3, self.nb, self.nb)
        return D

    def phonons(self, q):
        """Angular frequencies (rad/ps), eigenvectors e[q, iα, j] and group velocities (km/s = nm/ps)."""
        D, dD = self.dynmat(q, deriv=True)
        w2, e = np.linalg.eigh(D)
        w = np.sign(w2) * np.sqrt(np.abs(w2)) * W_SI * 1e-12         # rad/ps
        vmat = np.einsum('qaj,qcab,qbj->qjc', np.conj(e), dD, e).real    # (eV/A^2/amu) * A
        with np.errstate(divide='ignore', invalid='ignore'):
            v = vmat * W_SI ** 2 * ANG / (2 * w[..., None] * 1e12) * 1e-3   # km/s
        v[~np.isfinite(v)] = 0.0
        v[np.abs(w) < 1e-6] = 0.0
        return w, e, v


# ----------------------------------------------------------------------------------------------- grid & symmetry
def grid_points(mesh):
    n1, n2, n3 = mesh
    idx = np.array([(i, j, k) for k in range(n3) for j in range(n2) for i in range(n1)])   # i fastest
    return idx, idx / np.array(mesh, float)


def grid_index(ijk, mesh):
    ijk = np.mod(ijk, mesh)
    return ijk[..., 0] + mesh[0] * (ijk[..., 1] + mesh[1] * ijk[..., 2])


def point_group(cr, tol=1e-5):
    """Integer rotations W (lattice coords, x' = W x) of the crystal; the cell has one atom per species."""
    G = cr.A @ cr.A.T
    ops = []
    for flat in itertools.product((-1, 0, 1), repeat=9):
        W = np.array(flat).reshape(3, 3)
        if abs(abs(np.linalg.det(W)) - 1) > 1e-9:
            continue
        # x_cart = x @ A; metric preserved: W^T... use row convention x' = x @ W^T
        if not np.allclose(W.T @ G @ W, G, atol=1e-4 * np.abs(G).max()):
            continue
        ok = True
        for i, x in enumerate(cr.frac):
            y = W @ x
            j = [k for k in range(cr.nat) if cr.species[k] == cr.species[i]
                 and np.allclose(np.mod(y - cr.frac[k] + 0.5, 1) - 0.5, 0, atol=tol)]
            if not j:
                ok = False
                break
        if ok:
            ops.append(W)
    return ops


def irreducible(cr, mesh):
    """Irreducible q of the grid under the point group plus time reversal.

    Returns (irr_index, star_of, rot_of): irr_index (nirr,) grid indices; star_of (N,) which irreducible point
    each grid point maps to; rot_of (N,3,3) Cartesian rotation R with v(q_full) = R v(q_irr)."""
    mesh = np.array(mesh)
    idx, _ = grid_points(mesh)
    ops = point_group(cr)
    # reciprocal action on reduced q: q' = q @ W^{-1}  (since x' = W x)
    cart = [cr.A.T @ W @ np.linalg.inv(cr.A.T) for W in ops]
    recops = [np.rint(np.linalg.inv(W)).astype(int) for W in ops]
    N = len(idx)
    star_of = -np.ones(N, int)
    rot_of = np.zeros((N, 3, 3))
    irr = []
    for n in range(N):
        if star_of[n] >= 0:
            continue
        k = len(irr)
        irr.append(n)
        for Wr, Rc in zip(recops, cart):
            for sgn in (1, -1):
                img = grid_index(sgn * (idx[n] @ Wr), mesh)
                if star_of[img] < 0:
                    star_of[img] = k
                    rot_of[img] = sgn * Rc
    return np.array(irr), star_of, rot_of


# ----------------------------------------------------------------------------------------------- transport
def bose(w_radps, T):
    x = HBAR * w_radps * 1e12 / (KB * T)
    with np.errstate(over='ignore', divide='ignore'):
        return 1.0 / np.expm1(x)


def base_sigma(dv, cr, mesh):
    """dv (..., 3) km/s -> sigma in rad/ps (ShengBTE base_sigma)."""
    Bnm = cr.B * 10.0                                  # 1/nm
    s = sum((dv @ Bnm[nu] / mesh[nu]) ** 2 for nu in range(3))
    return np.sqrt(s / 6.0)


class Grid:
    def __init__(self, cr, mesh):
        self.cr, self.mesh = cr, np.array(mesh)
        self.idx, self.q = grid_points(self.mesh)
        self.N = len(self.q)
        self.w, self.e, self.v = cr.phonons(self.q)
        self.irr, self.star_of, self.rot_of = irreducible(cr, self.mesh)
        self.mult = np.bincount(self.star_of, minlength=len(self.irr))
        self.valid = self.w > 1e-3                      # excludes the three acoustic modes at Gamma

    def _gaussian(self, d, s):
        return np.exp(-(d / s) ** 2) / (np.sqrt(np.pi) * s)

    def isotope_rates(self, scalebroad):
        """Tamura isotope scattering (ps^-1) on irreducible modes."""
        cr = self.cr
        g = np.array([ISO_G[s] for s in cr.species])
        out = np.zeros((len(self.irr), cr.nb))
        e = self.e.reshape(self.N, cr.nat, 3, cr.nb)
        for k, iq in enumerate(self.irr):
            # overlap |sum_a e*_j(q) e_j'(q')|^2 per atom
            prod = np.einsum('iaj,niaJ->nijJ', np.conj(e[iq]), e)             # (N, nat, nb, nb)
            weight = np.einsum('i,nijJ->njJ', g, np.abs(prod) ** 2)
            sig = scalebroad * base_sigma(self.v[:, None, :, :], cr, self.mesh)
            sig = np.broadcast_to(sig, (self.N, cr.nb, cr.nb))
            d = self.w[iq][None, :, None] - self.w[:, None, :]
            ok = (np.abs(d) <= 2 * sig) & self.valid[:, None, :] & (sig > 0)
            delta = np.zeros_like(d)
            delta[ok] = self._gaussian(d[ok], sig[ok])
            rate = np.pi / 2 * self.w[iq] ** 2 * np.sum(weight * delta, axis=(0, 2)) / self.N
            out[k] = np.where(self.valid[iq], rate, 0.0)
        return out

    def three_phonon(self, scalebroad, T, keep=True):
        """Three-phonon rates on irreducible modes (ps^-1) and, if keep, the process lists for the iteration."""
        cr, N, mesh = self.cr, self.N, self.mesh
        nb = cr.nb
        fc3 = cr.fc3
        n2, n3 = fc3['n'][:, 0], fc3['n'][:, 1]
        E = np.exp(2j * np.pi * self.q @ (n2 - n3).T)                        # (N, nblock), q-independent
        groups = {}
        for b, at in enumerate(fc3['atoms']):
            groups.setdefault(tuple(at), []).append(b)
        nbose = bose(np.where(self.valid, self.w, 1.0), T)
        nbose[~self.valid] = 0.0
        # units: |V|^2 in (eV/A^3/amu^1.5)^2, omegas rad/ps -> 1e36 rad^3/s^3, delta in ps -> 1e-12 s.
        pref = HBAR * np.pi / 4 * (EV / (ANG ** 3 * AMU ** 1.5)) ** 2 / 1e36 * 1e-12 * 1e-12   # -> ps^-1
        rates_p = np.zeros((len(self.irr), nb))
        rates_m = np.zeros((len(self.irr), nb))
        procs = []
        e_c = np.conj(self.e)
        validw = self.valid
        for k, iq in enumerate(self.irr):
            qi = self.idx[iq]
            c = np.exp(-2j * np.pi * self.q[iq] @ n3.T)                       # (nblock,)
            w0 = self.w[iq][None, :, None, None]
            w1 = self.w[:, None, :, None]
            for sign in (1, -1):
                i2 = grid_index(qi + sign * self.idx, mesh)    # q'' = q + q' (+) or q - q' (-)
                w2 = self.w[i2][:, None, None, :]
                d = w0 + sign * w1 - w2
                sig = scalebroad * base_sigma(self.v[:, None, :, None, :] - self.v[i2][:, None, None, :, :], cr, mesh)
                sig = np.broadcast_to(sig, d.shape)
                ok = (np.abs(d) <= 2 * sig) & (sig > 0) & validw[iq][None, :, None, None] \
                    & validw[:, None, :, None] & validw[i2][:, None, None, :]
                rows = np.nonzero(ok.any(axis=(1, 2, 3)))[0]
                if len(rows) == 0:
                    continue
                Er = E[rows] if sign == 1 else np.conj(E[rows])
                Tt = np.zeros((len(rows), cr.nat, 3, cr.nat, 3, cr.nat, 3), complex)
                for (a1, a2, a3), bl in groups.items():
                    Wg = (c[bl][:, None] * fc3['fc'][bl].reshape(len(bl), 27))
                    Tt[:, a1, :, a2, :, a3, :] += (Er[:, bl] @ Wg).reshape(len(rows), 3, 3, 3)
                Tt = Tt.reshape(len(rows), nb, nb * nb)
                E1 = self.e[rows] if sign == 1 else e_c[rows]
                X = np.matmul(self.e[iq].T[None], Tt).reshape(len(rows), nb, nb, nb)       # (r, j, b, c)
                X = np.matmul(X.transpose(0, 1, 3, 2), E1[:, None])                         # (r, j, c, J)
                V = np.matmul(X.transpose(0, 1, 3, 2), e_c[i2[rows]][:, None])              # (r, j, J, K)
                okr = ok[rows]
                r_, j_, J_, K_ = np.nonzero(okr)
                n_ = rows[r_]
                n1 = nbose[n_, J_]
                n2b = nbose[i2[n_], K_]
                pop = (n1 - n2b) if sign == 1 else (n1 + n2b + 1)
                ww = self.w[iq][j_] * self.w[n_, J_] * self.w[i2[n_], K_]
                G = pref * pop / ww * np.abs(V[r_, j_, J_, K_]) ** 2 * self._gaussian(d[rows][okr], sig[rows][okr])
                if sign == 1:
                    np.add.at(rates_p[k], j_, G / N)
                else:
                    np.add.at(rates_m[k], j_, 0.5 * G / N)
                if keep:
                    procs.append((np.full(len(G), k, np.int32), np.full(len(G), sign, np.int8), j_.astype(np.int8),
                                  n_.astype(np.int32), J_.astype(np.int8), i2[n_].astype(np.int32),
                                  K_.astype(np.int8), G.astype(np.float64)))
        return rates_p, rates_m, procs

    def kappa(self, T, rate_diag, procs=None, iterate=False, max_iter=300, tol=1e-6):
        """Thermal conductivity tensor (W/mK) and solver info. rate_diag: (nirr, nb) total diagonal rate (ps^-1).

        Modes with zero total rate are dropped (reported as zero_rate_modes); ShengBTE would give them infinite
        lifetimes."""
        cr, N = self.cr, self.N
        nb = cr.nb
        wi = self.w[self.irr]
        vi = self.v[self.irr]
        validi = self.valid[self.irr]
        pos = validi & (rate_diag > 0)
        tau = np.where(pos, 1.0 / np.where(pos, rate_diag, 1.0), 0.0)                 # ps
        n0 = bose(np.where(validi, wi, 1.0), T)
        cv = np.where(validi, KB * (HBAR * wi * 1e12 / (KB * T)) ** 2 * n0 * (n0 + 1), 0.0)            # J/K per mode
        fac = 1e3 * 1e-9 / (self.volume_m3() * N)
        R, star = self.rot_of, self.star_of
        cvfull = cv[star]                                                                # (N, nb)
        vfull = self.v

        def kap(F):     # F (nirr, nb, 3) in nm
            Ff = np.einsum('nab,njb->nja', R, F[star])
            return np.einsum('nj,nja,njb->ab', cvfull, vfull, Ff) * fac

        F0 = tau[..., None] * vi            # km/s * ps = nm
        k_rta = kap(F0)
        info = {'zero_rate_modes': int(np.sum(self.mult[:, None] * (validi & ~pos))), 'solver': 'RTA'}
        if not iterate:
            return k_rta, info
        P = [np.concatenate([p[i] for p in procs]) for i in range(8)]
        kidx, sgn, j_, n_, J_, n2_, K_, G = P
        xi1 = self.w[n_, J_] / wi[kidx, j_]
        xi2 = self.w[n2_, K_] / wi[kidx, j_]
        plus = (sgn == 1)[:, None]
        GN = (G / N)[:, None]
        F = F0.copy()
        history = [k_rta]
        converged = False
        for it in range(max_iter):
            Ffull = np.einsum('nab,njb->nja', R, F[star])
            F1 = Ffull[n_, J_]
            F2 = Ffull[n2_, K_]
            contrib = np.where(plus, xi2[:, None] * F2 - xi1[:, None] * F1,
                               0.5 * (xi2[:, None] * F2 + xi1[:, None] * F1)) * GN
            Delta = np.zeros((len(self.irr), nb, 3))
            np.add.at(Delta, (kidx, j_), contrib)
            F = tau[..., None] * (vi + Delta)
            kn = kap(F)
            history.append(kn)
            change = np.max(np.abs(kn - history[-2])) / max(np.max(np.abs(kn)), 1e-30)
            if change < tol:
                converged = True
                break
            if not np.all(np.isfinite(kn)) or np.max(np.abs(kn)) > 1e3 * max(np.max(np.abs(k_rta)), 1.0):
                break
        info.update(solver='iterative', iterations=len(history) - 1, converged=converged,
                    kappa_rta=k_rta, history_diag=[np.diag(h) for h in history])
        return history[-1], info

    def volume_m3(self):
        return self.cr.volume * 1e-30


# ----------------------------------------------------------------------------------------------- four-phonon (Monte Carlo)
def four_phonon_mc(g, scalebroad, T, nsample, seed, irr_subset=None, batch=1024):
    """Four-phonon rates (ps^-1) on irreducible modes by uniform Monte Carlo sampling of (q', q'') pairs.

    For each irreducible q, nsample pairs (q', q'') are drawn uniformly from the N^2 grid pairs; q''' follows from
    momentum conservation for each process type. The estimator of (1/N^2) sum_{q'q''} is the sample mean, so it
    is unbiased; all branch combinations of each sampled pair are evaluated exactly.
    Returns rates (nirr, nb) for ++, +-, -- separately (with the 1/2, 1/2, 1/6 factors applied)."""
    cr, N, mesh = g.cr, g.N, g.mesh
    nb, nat = cr.nb, cr.nat
    fc4 = cr.fc4
    n2, n3, n4 = fc4['n'][:, 0], fc4['n'][:, 1], fc4['n'][:, 2]
    groups = {}
    for b, at in enumerate(fc4['atoms']):
        groups.setdefault(tuple(at), []).append(b)
    gl = [(key, np.array(bl), fc4['fc'][bl].reshape(len(bl), 81)) for key, bl in groups.items()]
    nbose = bose(np.where(g.valid, g.w, 1.0), T)
    nbose[~g.valid] = 0.0
    pref = np.pi * HBAR ** 2 / 8 * (EV / (ANG ** 4 * AMU ** 2)) ** 2 / 1e48 * 1e-12 * 1e-12
    e, ec = g.e, np.conj(g.e)
    rng = np.random.default_rng(seed)
    irr_list = range(len(g.irr)) if irr_subset is None else irr_subset
    out = {t: np.zeros((len(g.irr), nb)) for t in ('pp', 'pm', 'mm')}
    hits = np.zeros(len(g.irr), int)
    types = {'pp': ((1, 1, -1), 0.5), 'pm': ((1, -1, -1), 0.5), 'mm': ((-1, -1, -1), 1.0 / 6.0)}
    for k in irr_list:
        iq = g.irr[k]
        qi = g.idx[iq]
        q0 = g.q[iq]
        done = 0
        while done < nsample:
            m = min(batch, nsample - done)
            done += m
            i1 = rng.integers(0, N, m)
            i2 = rng.integers(0, N, m)
            for t, ((s1, s2, s3), fac) in types.items():
                # momentum: q + s1 q' + s2 q'' = -s3 q'''  ->  q''' = -(q + s1 q' + s2 q'')/s3
                i3 = grid_index(-s3 * (qi + s1 * g.idx[i1] + s2 * g.idx[i2]), mesh)
                ph = np.exp(2j * np.pi * (s1 * g.q[i1] @ n2.T + s2 * g.q[i2] @ n3.T + s3 * g.q[i3] @ n4.T))
                Tt = np.zeros((m, nat, 3, nat, 3, nat, 3, nat, 3), complex)
                for (a1, a2, a3, a4), bl, F in gl:
                    Tt[:, a1, :, a2, :, a3, :, a4, :] += (ph[:, bl] @ F).reshape(m, 3, 3, 3, 3)
                Tt = Tt.reshape(m, nb, nb * nb * nb)
                E1 = e[i1] if s1 == 1 else ec[i1]
                E2 = e[i2] if s2 == 1 else ec[i2]
                E3 = e[i3] if s3 == 1 else ec[i3]
                X = np.matmul(e[iq].T[None], Tt).reshape(m, nb, nb, nb * nb)            # (m, j, b, cd)
                X = np.matmul(X.transpose(0, 1, 3, 2), E1[:, None]).reshape(m, nb, nb, nb, nb)  # (m, j, c, d, J)
                X = np.matmul(X.transpose(0, 1, 4, 3, 2), E2[:, None, None])           # (m, j, J, d, L)
                X = np.matmul(X.transpose(0, 1, 2, 4, 3), E3[:, None, None])           # (m, j, J, L, M)
                V2 = np.abs(X) ** 2
                w0 = g.w[iq][None, :, None, None, None]
                w1 = g.w[i1][:, None, :, None, None]
                w2 = g.w[i2][:, None, None, :, None]
                w3 = g.w[i3][:, None, None, None, :]
                d = w0 + s1 * w1 + s2 * w2 + s3 * w3
                sig = scalebroad * base_sigma(g.v[i2][:, :, None, :] - g.v[i3][:, None, :, :], cr, mesh)   # (m, L, M)
                sig = np.broadcast_to(sig[:, None, None, :, :], d.shape)
                ok = (np.abs(d) <= 2 * sig) & (sig > 0) & g.valid[iq][None, :, None, None, None] \
                    & g.valid[i1][:, None, :, None, None] & g.valid[i2][:, None, None, :, None] \
                    & g.valid[i3][:, None, None, None, :]
                kk, j_, J_, L_, M_ = np.nonzero(ok)
                if len(kk) == 0:
                    continue
                n1, n2b, n3b = nbose[i1[kk], J_], nbose[i2[kk], L_], nbose[i3[kk], M_]
                if t == 'pp':
                    pop = n1 * n2b * (1 + n3b) - (1 + n1) * (1 + n2b) * n3b
                elif t == 'pm':
                    pop = n1 * (1 + n2b) * (1 + n3b) - (1 + n1) * n2b * n3b
                else:
                    pop = (1 + n1) * (1 + n2b) * (1 + n3b) - n1 * n2b * n3b
                ww = g.w[iq][j_] * g.w[i1[kk], J_] * g.w[i2[kk], L_] * g.w[i3[kk], M_]
                dd = d[ok]
                ss = sig[ok]
                G = pref * pop / ww * V2[ok] * np.exp(-(dd / ss) ** 2) / (np.sqrt(np.pi) * ss)
                np.add.at(out[t][k], j_, fac * G / nsample)
                hits[k] += len(kk)
    return out, hits


def shengbte_order(g):
    """Irreducible points in ShengBTE's listing order (third grid index fastest, first occurrence).

    Returns perm with perm[r] = index into g.irr of the r-th ShengBTE irreducible point."""
    n1, n2, n3 = g.mesh
    perm, seen = [], set()
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                s = g.star_of[grid_index(np.array([i, j, k]), g.mesh)]
                if s not in seen:
                    seen.add(s)
                    perm.append(s)
    return np.array(perm)


def read_phel(path, g):
    """BTE.sr_phel_metal (branch-major rows over ShengBTE irreducible points; omega rad/ps, rate ps^-1)
    mapped onto g.irr ordering. Returns (omega_archive, rate) with shape (nirr, nb)."""
    a = np.loadtxt(path)
    nirr = len(g.irr)
    assert a.shape[0] == nirr * g.cr.nb
    perm = shengbte_order(g)
    w = np.zeros((nirr, g.cr.nb))
    r = np.zeros((nirr, g.cr.nb))
    w[perm] = a[:, 0].reshape(g.cr.nb, nirr).T
    r[perm] = a[:, 1].reshape(g.cr.nb, nirr).T
    return w, r
