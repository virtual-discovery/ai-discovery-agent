"""Reference solution for amaranta-2d-anisotropic-ferromagnets; runs as the Harbor oracle and in check_task.

Per compound: bond classes by length (1 % grouping) up to the fourth shell; an 18-site supercell
in which no coupling reaches a site's own image. Energies of random noncollinear configurations
plus ordered states along x, y, z (ferromagnet and the partner pattern that reverses every
nearest-neighbour product) are fitted by least squares with: a full 3x3 nearest-neighbour tensor
(honeycomb: one tensor carried to the other bonds by the C3 rotation about Cr, verified on the
ligand positions; CrSBr: an independent tensor for each of the two Cr chains), scalar couplings
for shells 2-4, an on-site K (e^z)^2 term and a constant. Ground state from Luttinger-Tisza;
anisotropy from the fitted model (single-ion part = the on-site term); CrI3 bond exchange in the
ideal octahedral frame; CrSBr DM vector and zone-centre gap; CrI3 ordering temperature from
classical Monte Carlo (Binder-cumulant crossings, two sizes).
"""
import json
import math
import os
import subprocess
import time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

APP = Path(os.environ.get('APP_DIR', '/app'))
CAT = json.loads((APP/'catalog.json').read_text())
MAXC = 6
S_SPIN = 1.5
K_PER_MEV = 11.604518
evidence, notes = [], []


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    for _ in range(5):
        out = subprocess.run(cmd, capture_output=True, text=True).stdout
        try:
            result = json.loads(out)
            break
        except json.JSONDecodeError:
            time.sleep(2)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def energies(compound, n1, n2, configs):
    out = []
    for s in range(0, len(configs), MAXC):
        r = lab('evaluate', {'kind': 'energies', 'compound': compound, 'supercell': [n1, n2],
                             'spins': [np.round(c, 6).tolist() for c in configs[s:s + MAXC]]})
        while r['status'] == 'running':
            time.sleep(1)
            r = lab('job/%d' % r['id'])
        if r['status'] != 'complete':
            raise RuntimeError(r)
        evidence.append(r['id'])
        out += [e * 1000.0 for e in r['result']['total_energies_eV']]   # meV
    return np.array(out)


def rot_z(th):
    return np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1.0]])


# ---------------------------------------------------------------- geometry
class Lattice:
    def __init__(self, name):
        comp = CAT['compounds'][name]
        self.name = name
        self.cell = np.array(comp['cell_A'], float)
        cz = self.cell[2, 2]
        mag = np.array(comp['magnetic_sites'], float)
        z0 = mag[0, 2]
        mag[:, 2] -= np.round((mag[:, 2] - z0) / cz) * cz          # the layer may straddle the cell edge in z
        self.basis = mag
        self.nb = len(mag)
        self.lig = []
        for s in comp['sites']:
            if s['species'] != 'Cr':
                p = np.array(s['position_A'], float)
                p[2] -= np.round((p[2] - z0) / cz) * cz
                self.lig.append((s['species'], p))
        self.A2 = self.cell[:2, :2]
        self.B = 2 * np.pi * np.linalg.inv(self.A2).T
        vecs = self.vectors(0, 11.0)
        lengths = []
        for k, m1, m2, d in sorted(vecs, key=lambda v: np.linalg.norm(v[3])):
            L = np.linalg.norm(d)
            if not any(abs(L / x - 1) < 0.01 for x in lengths):
                lengths.append(L)
        for k0 in range(1, self.nb):
            for k, m1, m2, d in self.vectors(k0, 11.0):
                L = np.linalg.norm(d)
                if not any(abs(L / x - 1) < 0.01 for x in lengths):
                    lengths.append(L)
        self.lengths = sorted(lengths)[:4]
        self.rcut = self.lengths[-1] * 1.01

    def vectors(self, k0, rmax):
        out = []
        for k, b in enumerate(self.basis):
            for m1 in range(-6, 7):
                for m2 in range(-6, 7):
                    d = b + m1 * self.cell[0] + m2 * self.cell[1] - self.basis[k0]
                    if 0.1 < np.linalg.norm(d) <= rmax:
                        out.append((k, m1, m2, d))
        return out

    def cls(self, d):
        L = np.linalg.norm(d)
        for i, x in enumerate(self.lengths):
            if abs(L / x - 1) < 0.01:
                return i
        return None

    def bonds(self):
        """Every pair once: (k0, k, m1, m2, d, class) with d oriented canonically (first nonzero of x, y, z positive)."""
        out = []
        for k0 in range(self.nb):
            for k, m1, m2, d in self.vectors(k0, self.rcut):
                c = self.cls(d)
                if c is None:
                    continue
                key = np.round(d, 3)
                canon = next((v > 0 for v in key if abs(v) > 1e-3), True)
                if canon:
                    out.append((k0, k, m1, m2, d, c))
        return out

    def c3_ok(self):
        """Is the 120-degree rotation about Cr site 0 a symmetry of the layer (ligands included)?"""
        R = rot_z(2 * np.pi / 3)
        x0 = self.basis[0]
        pts = [(s, p) for s, p in self.lig] + [('Cr', b) for b in self.basis]
        for s, p in pts:
            q = R @ (p - x0) + x0
            if not any(s2 == s and self.same(q, p2) for s2, p2 in pts):
                return False
        return True

    def same(self, x, y, tol=0.05):
        if abs(x[2] - y[2]) > tol:
            return False
        f = np.linalg.solve(self.A2.T, (x - y)[:2])
        return np.linalg.norm((f - np.round(f)) @ self.A2) < tol


# ---------------------------------------------------------------- model and regression
class Model:
    """E = sum_pairs e_i.J.e_j + K sum (e^z)^2 + const; class 0 tensor(s), classes 1-3 scalar."""

    def __init__(self, lat):
        self.lat = lat
        self.bl = lat.bonds()
        self.honeycomb = lat.nb == 2 and lat.c3_ok()
        nn = [b for b in self.bl if b[5] == 0]
        if self.honeycomb:
            # one tensor on the bond from site 0 along d0, carried by C3 about site 0 to the other two
            d0 = next(b[4] for b in nn if b[0] == 0)
            self.d0 = d0
            self.rot = {}
            for b in nn:
                d = b[4] if b[0] == 0 else -b[4]     # bond seen from site 0
                for m in range(3):
                    R = np.linalg.matrix_power(rot_z(2 * np.pi / 3), m)
                    if np.allclose(R @ d0, d, atol=0.05):
                        self.rot[id(b)] = (R, b[0] != 0)
            self.ntensor = 1
        else:
            # independent tensor per (start site, canonical direction) bond family
            fam = []
            self.fam_of = {}
            for b in nn:
                key = (b[0], tuple(np.round(b[4], 2)))
                if key not in fam:
                    fam.append(key)
                self.fam_of[id(b)] = fam.index(key)
            self.fams = fam
            self.ntensor = len(fam)
        self.npar = 9 * self.ntensor + (len(lat.lengths) - 1) + 2

    def features(self, n1, n2, spins):
        lat = self.lat
        e = np.asarray(spins, float)
        e = e / np.linalg.norm(e, axis=1, keepdims=True)
        S = e.reshape(lat.nb, n1, n2, 3)
        f = np.zeros(self.npar)
        for b in self.bl:
            k0, k, m1, m2, d, c = b
            Si = S[k0]
            Sj = np.roll(np.roll(S[k], -m1, axis=0), -m2, axis=1)
            if c == 0:
                if self.honeycomb:
                    R, rev = self.rot[id(b)]
                    a, bb = (Sj, Si) if rev else (Si, Sj)          # oriented from the site-0 end
                    X = np.einsum('ija,ijb->ab', a @ R, bb @ R)     # (R^T e)_a (R^T e)_b summed
                    f[:9] += X.ravel()
                else:
                    t = self.fam_of[id(b)]
                    f[9 * t:9 * t + 9] += np.einsum('ija,ijb->ab', Si, Sj).ravel()
            else:
                f[9 * self.ntensor + c - 1] += np.sum(Si * Sj)
        f[-2] = np.sum(e[:, 2] ** 2)
        f[-1] = len(e)
        return f

    def set(self, p):
        self.p = np.array(p)
        self.T = [self.p[9 * t:9 * t + 9].reshape(3, 3) for t in range(self.ntensor)]
        self.J = self.p[9 * self.ntensor:9 * self.ntensor + len(self.lat.lengths) - 1]
        self.K = self.p[-2]

    def energy(self, n1, n2, spins):
        return float(self.features(n1, n2, spins)[:-1] @ self.p[:-1])

    def bond_tensor(self, b):
        if b[5] != 0:
            return np.eye(3) * self.J[b[5] - 1]
        if self.honeycomb:
            R, rev = self.rot[id(b)]
            X = R @ self.T[0] @ R.T
            return X.T if rev else X
        return self.T[self.fam_of[id(b)]]

    def table(self):
        """Directed table (k0, k, m1, m2, J) with both orientations, for Monte Carlo and spin waves."""
        out = []
        for b in self.bl:
            k0, k, m1, m2, d, c = b
            J = self.bond_tensor(b)
            out.append((k0, k, m1, m2, J))
            out.append((k, k0, -m1, -m2, J.T))
        return out


def supercell(lat):
    """Cheapest diagonal cell with >= 12 sites in which no coupling reaches a site's own image and every
    tensor family has distinct partner sites (no J + J^T aliasing)."""
    best = None
    for n1 in range(1, 13):
        for n2 in range(1, 13):
            N = lat.nb * n1 * n2
            if N > 36 or N < 13:
                continue
            ok = True
            for k0 in range(lat.nb):
                for k, m1, m2, d in lat.vectors(k0, lat.rcut):
                    if lat.cls(d) is None:
                        continue
                    if k == k0 and m1 % n1 == 0 and m2 % n2 == 0:
                        ok = False
                    if k == k0 and (2 * m1) % n1 == 0 and (2 * m2) % n2 == 0:
                        ok = False      # +d and -d reach the same site
            if ok:
                key = (math.ceil(N / 12), N, abs(n1 - n2))
                if best is None or key < best[0]:
                    best = (key, n1, n2)
    return best[1], best[2]


def partner_pattern(lat, n1, n2, gs_signs):
    """Collinear pattern with every nearest-neighbour product reversed relative to the ground-state pattern."""
    if lat.nb == 2 and lat.c3_ok():
        s = np.array([gs_signs[0], -gs_signs[1]])
        return np.repeat(s, n1 * n2)
    return np.array([gs_signs[k] * (-1) ** i for k in range(lat.nb) for i in range(n1) for j in range(n2)], float)


def fit_compound(name, nrandom):
    lat = Lattice(name)
    mod = Model(lat)
    n1, n2 = supercell(lat)
    N = lat.nb * n1 * n2
    rng = np.random.default_rng(sum(ord(ch) * 31 ** i for i, ch in enumerate(name)) % (2 ** 32))
    confs = []
    for _ in range(nrandom):
        v = rng.normal(size=(N, 3))
        confs.append(v / np.linalg.norm(v, axis=1, keepdims=True))
    fm = np.ones(N)
    alt = partner_pattern(lat, n1, n2, [1.0] * lat.nb)
    named = {'FM_x': np.outer(fm, [1, 0, 0]), 'FM_y': np.outer(fm, [0, 1, 0]), 'FM_z': np.outer(fm, [0, 0, 1]),
             'ALT_x': np.outer(alt, [1, 0, 0]), 'ALT_z': np.outer(alt, [0, 0, 1])}
    if not mod.honeycomb:
        # low symmetry: the principal axes of the anisotropy are not fixed by C3, so probe off-axis directions too
        r2 = 1 / np.sqrt(2)
        for tag, v in (('xz+', [r2, 0, r2]), ('xz-', [r2, 0, -r2]), ('yz+', [0, r2, r2]), ('yz-', [0, r2, -r2]),
                       ('xy+', [r2, r2, 0]), ('xy-', [r2, -r2, 0])):
            named['FM_' + tag] = np.outer(fm, v)
    keys = list(named)
    confs += [named[k] for k in keys]
    E = energies(name, n1, n2, confs)
    X = np.array([mod.features(n1, n2, c) for c in confs])
    p, *_ = np.linalg.lstsq(X, E, rcond=None)
    mod.set(p)
    info = {'supercell': [n1, n2], 'n_configs': len(E), 'rank': int(np.linalg.matrix_rank(X)), 'npar': mod.npar}
    if not mod.honeycomb and mod.ntensor > 1:
        # compare the independently fitted tensors of the other bond families with the first one
        T0 = mod.T[0]
        ties = []
        for t in range(1, mod.ntensor):
            same, trans = np.linalg.norm(mod.T[t] - T0), np.linalg.norm(mod.T[t] - T0.T)
            ties.append('transpose' if trans < same else 'same')
            info[f'family{t}_vs_family0'] = {'|T - T0|': float(same), '|T - T0^T|': float(trans)}
        # refit with the tie the free fit revealed (symmetry-related chains)
        A = X[:, :9].copy()
        for t, tie in enumerate(ties, start=1):
            f = X[:, 9 * t:9 * t + 9].reshape(-1, 3, 3)
            A += (f.transpose(0, 2, 1) if tie == 'transpose' else f).reshape(-1, 9)
        Xc = np.c_[A, X[:, 9 * mod.ntensor:]]
        pc, *_ = np.linalg.lstsq(Xc, E, rcond=None)
        full = [pc[:9]]
        for tie in ties:
            Tt = pc[:9].reshape(3, 3)
            full.append((Tt.T if tie == 'transpose' else Tt).ravel())
        mod.set(np.concatenate(full + [pc[9:]]))
        info['ties'] = ties
        X = Xc
        p = pc
    resid = E - X @ p
    dof = len(E) - np.linalg.matrix_rank(X)
    info['rms_resid_meV'] = float(np.sqrt(np.sum(resid ** 2) / max(dof, 1)))
    mod.fit_info = info
    mod.collective = {k: float(E[nrandom + i]) for i, k in enumerate(keys)}
    mod.Ncell = N
    return lat, mod


# ---------------------------------------------------------------- analysis
def jq_min(lat, mod):
    iso = [(b, np.trace(mod.bond_tensor(b)) / 3) for b in mod.bl]
    def lam(x):
        q = np.array(x) @ lat.B
        Mq = np.zeros((lat.nb, lat.nb), complex)
        for b, J in iso:
            k0, k, m1, m2, d, c = b
            Mq[k0, k] += J * np.exp(1j * q @ d[:2])
            Mq[k, k0] += J * np.exp(-1j * q @ d[:2])
        return np.linalg.eigh(Mq)
    best = None
    for h in np.linspace(0, 1, 41):
        for k in np.linspace(0, 1, 41):
            w, V = lam([h, k])
            if best is None or w[0] < best[0] - 1e-9:
                best = (w[0], h, k)
    r = minimize(lambda x: lam(x)[0][0], best[1:], method='Nelder-Mead', options={'xatol': 1e-9, 'fatol': 1e-12})
    q = np.array(r.x)
    q = q - np.round(q)
    w, V = lam(q)
    return q, V[:, 0], w[0]


def classify(q, vec):
    if np.linalg.norm(q) > 1e-3:
        return 'noncollinear', None
    s = np.sign(np.real(vec / vec[np.argmax(np.abs(vec))]))
    if np.all(s > 0):
        return 'ferromagnet', s
    return 'Neel antiferromagnet', s


def aniso_matrix(lat, mod, signs, onsite_only=False):
    def E(n):
        spins = np.outer(signs, n)
        if onsite_only:
            return mod.K * np.sum(spins[:, 2] ** 2) / len(signs)
        return mod.energy(1, 1, spins) / len(signs)
    ex = np.eye(3)
    Mm = np.zeros((3, 3))
    for a in range(3):
        Mm[a, a] = E(ex[a])
    for a in range(3):
        for b in range(a + 1, 3):
            Mm[a, b] = Mm[b, a] = (E((ex[a] + ex[b]) / np.sqrt(2)) - E((ex[a] - ex[b]) / np.sqrt(2))) / 2
    return Mm


def predict(lat, mod, p):
    n1, n2 = p['supercell']
    S = np.array([[math.sin(math.radians(t)) * math.cos(math.radians(f)), math.sin(math.radians(t)) * math.sin(math.radians(f)),
                   math.cos(math.radians(t))] for t, f in p['spin_angles_deg']])
    fm = np.tile([0.0, 0.0, 1.0], (len(S), 1))
    return (mod.energy(n1, n2, S) - mod.energy(n1, n2, fm)) / len(S)


def cubic_frame(lat, d, bond_dirs):
    """Ideal octahedral axes (polar angle arccos(1/sqrt3)): of the two 120-degree sets of azimuths that are exactly
    perpendicular to the three bonds of Cr site 0, the one nearest the Cr->I directions to the I above Cr site 0;
    e3 is the axis perpendicular to bond d."""
    x0 = lat.basis[0]
    az = []
    for s, p in lat.lig:
        for m1 in range(-1, 2):
            for m2 in range(-1, 2):
                v = p + m1 * lat.cell[0] + m2 * lat.cell[1] - x0
                if v[2] > 0.5 and np.linalg.norm(v) < 3.5:
                    az.append(math.atan2(v[1], v[0]))
    best = None
    for side in (1, -1):
        cand = [math.atan2(b[1], b[0]) + side * math.pi / 2 for b in bond_dirs]
        mis = sum(min(abs(math.remainder(a - c, 2 * math.pi)) for c in cand) for a in az)
        if best is None or mis < best[0]:
            best = (mis, cand)
    axes = [np.array([math.sqrt(2 / 3) * math.cos(c), math.sqrt(2 / 3) * math.sin(c), math.sqrt(1 / 3)]) for c in best[1]]
    u = d / np.linalg.norm(d)
    i3 = int(np.argmin([abs(a @ u) for a in axes]))
    e3 = axes[i3]
    e1, e2 = [axes[i] for i in range(3) if i != i3]
    return e1, e2, e3, [math.degrees(a) for a in az]


def jkgg(J, e1, e2, e3):
    Js = (J + J.T) / 2
    basis = [np.eye(3), np.outer(e3, e3), np.outer(e1, e2) + np.outer(e2, e1),
             np.outer(e1, e3) + np.outer(e3, e1) + np.outer(e2, e3) + np.outer(e3, e2)]
    X = np.array([b.ravel() for b in basis]).T
    return np.linalg.lstsq(X, Js.ravel(), rcond=None)[0]


# ---------------------------------------------------------------- Monte Carlo
class MCLattice:
    def __init__(self, nb, table, K, L):
        self.K = K
        self.N = N = nb * L * L
        idx = lambda k, i, j: (k * L + i % L) * L + j % L
        nbr = [[] for _ in range(N)]
        Jm = [[] for _ in range(N)]
        for k0, k, m1, m2, J in table:
            for i in range(L):
                for j in range(L):
                    a, b = idx(k0, i, j), idx(k, i + m1, j + m2)
                    nbr[a].append(b)
                    Jm[a].append(J)
        z = max(len(x) for x in nbr)
        self.nbr = np.array([x + [x[0]] * (z - len(x)) for x in nbr])
        self.J = np.array([y + [np.zeros((3, 3))] * (z - len(y)) for y in Jm])
        color = -np.ones(N, int)
        for a in range(N):
            used = {color[b] for b in nbr[a] if color[b] >= 0}
            cc = 0
            while cc in used:
                cc += 1
            color[a] = cc
        self.classes = [np.where(color == c)[0] for c in range(color.max() + 1)]

    def field(self, S, sites):
        return np.einsum('nzab,nzb->na', self.J[sites], S[self.nbr[sites]])


def mc_point(args):
    nb, table, K, L, T_K, sweeps, seed = args
    lat = MCLattice(nb, [(k0, k, m1, m2, np.array(J)) for k0, k, m1, m2, J in table], K, L)
    rng = np.random.default_rng(seed)
    T = T_K / K_PER_MEV
    S = np.tile([0, 0, 1.0], (lat.N, 1))
    width = 0.5
    therm = sweeps // 4
    mz = []
    for sw in range(therm + sweeps):
        for _ in range(2):              # over-relaxation about the exchange field, Metropolis on the on-site change
            for sites in lat.classes:
                old = S[sites]
                h = lat.field(S, sites)
                hn = h / np.maximum(np.linalg.norm(h, axis=1, keepdims=True), 1e-12)
                prop = 2 * np.sum(old * hn, axis=1, keepdims=True) * hn - old
                dE = K * (prop[:, 2] ** 2 - old[:, 2] ** 2)
                ok = (dE <= 0) | (rng.random(len(sites)) < np.exp(-np.clip(dE, 0, None) / T))
                S[sites[ok]] = prop[ok]
        acc = 0
        for sites in lat.classes:
            old = S[sites]
            v = rng.normal(size=old.shape)
            prop = old + width * v / np.linalg.norm(v, axis=1, keepdims=True)
            prop /= np.linalg.norm(prop, axis=1, keepdims=True)
            h = lat.field(S, sites)
            dE = np.sum((prop - old) * h, axis=1) + K * (prop[:, 2] ** 2 - old[:, 2] ** 2)
            ok = (dE <= 0) | (rng.random(len(sites)) < np.exp(-np.clip(dE, 0, None) / T))
            S[sites[ok]] = prop[ok]
            acc += ok.sum()
        if sw < therm:
            width = min(max(width * (1.05 if acc / lat.N > 0.5 else 0.95), 0.02), 2.0)
        else:
            mz.append(S[:, 2].mean())
    mz = np.array(mz)
    return {'L': L, 'T': T_K, 'U': float(1 - np.mean(mz ** 4) / (3 * np.mean(mz ** 2) ** 2)), 'mz': float(np.mean(np.abs(mz)))}


def ordering_temperature(mod, nb):
    table = [(k0, k, m1, m2, J.tolist()) for k0, k, m1, m2, J in mod.table()]
    procs = max(1, min(4, os.cpu_count() or 1))
    with Pool(procs) as pool:
        coarse = pool.map(mc_point, [(nb, table, mod.K, 12, T, 1500, 7 + i) for i, T in enumerate(range(10, 81, 5))])
        # coarse location: where U of the 12x12 cell drops through 0.4
        Ts = [r['T'] for r in coarse]
        Us = [r['U'] for r in coarse]
        i = next(j for j in range(1, len(Us)) if Us[j] < 0.4)
        t0 = Ts[i - 1] + (Ts[i] - Ts[i - 1]) * (Us[i - 1] - 0.4) / (Us[i - 1] - Us[i])
        grid = [round(t0 + dt, 2) for dt in np.linspace(-0.12 * t0, 0.06 * t0, 7)]
        jobs = [(nb, table, mod.K, L, T, sw, 100 * L + j) for L, sw in ((12, 8000), (18, 6000)) for j, T in enumerate(grid)]
        fine = pool.map(mc_point, jobs)
    ser = {L: {x['T']: x['U'] for x in fine if x['L'] == L} for L in (12, 18)}
    diff = [ser[12][t] - ser[18][t] for t in grid]
    # crossing: U_12 - U_18 turns from negative (ordered side) to positive
    i = next((j for j in range(len(grid) - 1) if diff[j] < 0 <= diff[j + 1]), None)
    if i is None:
        Tc = float(t0)
    else:
        fits = {L: np.polyfit(grid, [ser[L][t] for t in grid], 3) for L in (12, 18)}
        roots = [x.real for x in np.roots(np.polysub(fits[12], fits[18])) if abs(x.imag) < 1e-9 and grid[i] <= x.real <= grid[i + 1]]
        Tc = float(roots[0]) if roots else float(grid[i] + (grid[i + 1] - grid[i]) * (-diff[i]) / (diff[i + 1] - diff[i]))
    return Tc, {'coarse_L12': [(r['T'], round(r['U'], 3)) for r in coarse], 'fine': [(x['L'], x['T'], round(x['U'], 3)) for x in fine]}


# ---------------------------------------------------------------- main
def main():
    names = list(CAT['compounds'])
    nrand = {'CrI3': 14, 'CrGeTe3': 14, 'CrSBr': 16}   # 10 for the others
    out = {'ground_states': {}, 'easy_axis_polar_deg': {}, 'anisotropy_meV_per_Cr': {}, 'energy_predictions_meV_per_site': {}}
    models = {}
    for name in names:
        lat, mod = fit_compound(name, nrand.get(name, 10))
        models[name] = (lat, mod)
        q, vec, e = jq_min(lat, mod)
        typ, signs = classify(q, vec)
        out['ground_states'][name] = typ
        Mm = aniso_matrix(lat, mod, signs)
        w, V = np.linalg.eigh(Mm)
        out['easy_axis_polar_deg'][name] = float(np.degrees(np.arccos(min(1.0, abs(V[2, 0])))))
        Ms = aniso_matrix(lat, mod, signs, onsite_only=True)
        out['anisotropy_meV_per_Cr'][name] = {'total': float(Mm[0, 0] - Mm[2, 2]), 'single_ion': float(Ms[0, 0] - Ms[2, 2])}
        mod.signs, mod.Mm = signs, Mm
        for pid, p in CAT['to_predict'].items():
            if p['compound'] == name:
                out['energy_predictions_meV_per_site'][pid] = predict(lat, mod, p)
        Jiso = [np.trace(mod.T[0]) / 3] + list(mod.J)
        notes.append(f"- **{name}** (cell {mod.fit_info['supercell']}, {mod.fit_info['n_configs']} configurations, fit rms "
                     f"{mod.fit_info['rms_resid_meV']:.3f} meV): shells " + ', '.join(f'{L:.3f}' for L in lat.lengths)
                     + ' A; isotropic J = ' + ', '.join(f'{j:+.3f}' for j in Jiso) + f' meV; K = {mod.K:+.4f} meV; {typ}; '
                     f"E(x) - E(z) = {Mm[0, 0] - Mm[2, 2]:+.4f} meV/Cr, of which on-site {-mod.K:+.4f}; eigenvalues "
                     + ', '.join(f'{x - w[0]:.4f}' for x in w) + ' meV')
    # CrI3 bond exchange in the ideal octahedral frame
    lat, mod = models['CrI3']
    dirs = [b[4] if b[0] == 0 else -b[4] for b in mod.bl if b[5] == 0]
    e1, e2, e3, az = cubic_frame(lat, mod.d0, dirs)
    bJ, bK, bG, bGp = jkgg(mod.T[0], e1, e2, e3)
    out['CrI3_bond_exchange_meV'] = {'J': float(bJ), 'K': float(bK), 'Gamma': float(bG), 'Gamma_prime': float(bGp)}
    # CrSBr: DM vector of the bond site 0 -> site 0 + a1, staggering check, gap
    latb, modb = models['CrSBr']
    Ds = {}
    for b in modb.bl:
        if b[5] == 0 and b[0] == b[1] and (b[2], b[3]) == (1, 0):
            J = modb.bond_tensor(b)
            Ds[b[0]] = np.array([(J[1, 2] - J[2, 1]) / 2, (J[2, 0] - J[0, 2]) / 2, (J[0, 1] - J[1, 0]) / 2])
    out['CrSBr_dm_vector_meV'] = [float(x) for x in Ds[0]]
    cE, Nc = modb.collective, modb.Ncell
    Md = np.diag([cE['FM_x'], cE['FM_y'], cE['FM_z']]) / Nc        # anisotropy of the ferromagnet measured directly
    for (a, b), tag in (((0, 2), 'xz'), ((1, 2), 'yz'), ((0, 1), 'xy')):
        Md[a, b] = Md[b, a] = (cE['FM_' + tag + '+'] - cE['FM_' + tag + '-']) / (2 * Nc)
    w, V = np.linalg.eigh(Md)
    gap = 2.0 / S_SPIN * math.sqrt(max(w[1] - w[0], 0) * max(w[2] - w[0], 0))
    notes.append(f"- CrSBr ferromagnet anisotropy measured directly (18-site cell): eigenvalues {np.round(w - w[0], 4).tolist()} meV/Cr, "
                 f"easy axis {np.round(V[:, 0], 3).tolist()}, intermediate {np.round(V[:, 1], 3).tolist()}; tensor fit: {modb.fit_info}")
    out['CrSBr_spin_wave_gap_meV'] = gap
    # CrI3 ordering temperature
    Tc, mcinfo = ordering_temperature(mod, lat.nb)
    out['CrI3_Tc_K'] = Tc
    out['evidence_ids'] = evidence
    (APP/'submission.json').write_text(json.dumps(out, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# Cr-based 2D magnets (reference solution)

Convention: E = sum over pairs (each once) of e_i.J.e_j for unit spins plus on-site terms; x along a1, z normal.
Each compound was fitted by least squares to random noncollinear configurations and ordered states along x, y, z
(ferromagnet and the partner pattern that reverses every nearest-neighbour product) in an 18-site cell in which no
coupling reaches a site's own image and +d and -d bonds never land on the same site (a 2a-periodic cell would
measure J + J^T and lose the antisymmetric part). Model: full nearest-neighbour tensor (honeycomb: C3 about Cr,
checked on the ligands, carries it to the other bonds; CrSBr: independent tensors for the two Cr chains),
scalar shells 2-4, on-site K (e^z)^2. Fourth-shell couplings come out at noise level. Fit residuals are at the
noise level (0.02 meV) except where noted.

{chr(10).join(notes)}

Ground states (Luttinger-Tisza of the isotropic part): CrSiSe3 has antiferromagnetic J1 and J3, so it is a Neel
antiferromagnet, not a ferromagnet; its anisotropy must be evaluated for the Neel pattern, where the
nearest-neighbour exchange anisotropy changes sign and adds to the easy-plane single-ion term (moments in plane).
Single-ion versus exchange: ferromagnet-only energies cannot separate them (both are quadratic forms of the
moment direction); the partner pattern reverses the nearest-neighbour contribution and leaves the on-site one.
Light ligands (Cl) are single-ion dominated; Br and I have comparable exchange anisotropy; in CrGeTe3 the single-ion
term is easy-plane and the exchange anisotropy (Te SOC) overrides it.

CrI3 bond (ideal cubic frame from the Cr->I directions above Cr, azimuths {np.round(az, 1).tolist()} deg):
J = {bJ:.3f}, K = {bK:.3f}, Gamma = {bG:.3f}, Gamma' = {bGp:.3f} meV. The tensors of the three bonds differ by C3;
copying one bond's tensor to all three would put the easy axis in the plane.

CrSBr: DM vector of the bond site 0 -> site 0 + a1: {np.round(Ds[0], 3).tolist()} meV; on the other chain
{np.round(Ds[1], 3).tolist()} meV (staggered: inversion relates the two chains), so it cancels in the ferromagnet and
in co-rotating spirals. Triaxial anisotropy (easy c, intermediate a, hard b); zone-centre LSWT gap for S = 3/2:
{gap:.4f} meV = (2/S) sqrt(D1 D2), not a uniaxial value.

CrI3 ordering temperature: classical Monte Carlo of the fitted model (Metropolis + over-relaxation), Binder-cumulant
crossing of 12x12 and 18x18 cells: Tc = {Tc:.1f} K. Isotropic exchange alone gives no order at T > 0
(Mermin-Wagner); the anisotropy gap stabilizes Ising-like order and sets Tc logarithmically; mean field
overestimates it strongly. {json.dumps(mcinfo)}

Limits: no magnetic dipole-dipole interaction in these energies; for CrCl3 the shape anisotropy of a 3 muB
honeycomb layer (about 0.05 meV/Cr, easy plane) exceeds the magnetocrystalline 0.03 meV, so the real monolayer is
likely easy-plane. Tc of CrCl3/CrBr3 (anisotropy below the source precision) is not supported; quantum
fluctuations and interlayer coupling are not addressed.
""")
    lab('submit', out)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
