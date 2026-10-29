"""Pure-numpy evaluator for the released Pd-H atomic cluster expansion (ACE) potential.

Port of the C-tilde evaluator of the ACE library (ICAMS ace-evaluator: ace_evaluator.cpp, ace_radial.cpp,
ace_spherical_cart.cpp, ace_abstract_basis.cpp) for the features used by PdH_PBE_v03a.yaml:
simplified-Bessel radial basis, complex spherical harmonics with Y00 = 1, Finnis-Sinclair shifted-scaled
embedding and the ZBL inner core. Energies in eV, forces in eV/A, stress in eV/A^3 (ASE sign convention:
stress = (1/V) dE/dstrain, so pressure = -trace/3).

Validated against python-ace 0.4.0rc1 (see benchmarks/construction/pdh-miscibility-isotope-ace/validate_ace.py).
"""
import json
import math
import numpy as np

K_COUL = 14.399645351950543
ZBL_C = np.array([0.18175, 0.50986, 0.28022, 0.02817])
ZBL_E = np.array([-3.19980, -0.94229, -0.40290, -0.20162])
ZNUM = {'H': 1, 'Pd': 46}
LMAX = 4
NLM = (LMAX + 1) ** 2
L_OF = np.array([l for l in range(LMAX + 1) for m in range(-l, l + 1)])
M_OF = np.array([m for l in range(LMAX + 1) for m in range(-l, l + 1)])


def poly_cut(r, r_in, delta):
    """cutoff_func_poly: 1 below r_in - delta, 0 above r_in, smooth in between. Returns (fc, dfc)."""
    r = np.asarray(r, float)
    fc, dfc = np.ones_like(r), np.zeros_like(r)
    fc[r >= r_in] = 0.0
    mid = (r > r_in - delta) & (r < r_in)
    x = 1 - 2 * (1 + (r[mid] - r_in) / delta)
    fc[mid] = 0.5 + 7.5 / 2. * (x / 4. - x ** 3 / 6. + x ** 5 / 20.)
    dfc[mid] = -7.5 / delta * (0.25 - x * x / 2.0 + x ** 4 / 4.)
    return fc, dfc


def sinc_and_deriv(x):
    s = np.where(np.abs(x) > 1e-12, np.sin(x) / np.where(x == 0, 1, x), 1.0)
    ds = np.where(np.abs(x) > 1e-12, (np.cos(x) * x - np.sin(x)) / np.where(x == 0, 1, x) ** 2, 0.0)
    return s, ds


def sbessel(r, rc, nrad):
    """Simplified Bessel radial basis g_n(r), n = 0..nrad-1, and derivatives (zero beyond rc)."""
    r = np.asarray(r, float)
    g, dg = np.zeros((len(r), nrad)), np.zeros((len(r), nrad))
    inside = r < rc
    x = r[inside]

    def aux(n):
        pre = (-1) ** n * math.sqrt(2) * math.pi / rc ** 1.5 * (n + 1) * (n + 2) / math.sqrt((n + 1) ** 2 + (n + 2) ** 2)
        a1, da1 = sinc_and_deriv(x * (n + 1) * math.pi / rc)
        a2, da2 = sinc_and_deriv(x * (n + 2) * math.pi / rc)
        return pre * (a1 + a2), pre * (da1 * (n + 1) * math.pi / rc + da2 * (n + 2) * math.pi / rc)

    gi, dgi = np.zeros((len(x), nrad)), np.zeros((len(x), nrad))
    gi[:, 0], dgi[:, 0] = aux(0)
    d_prev = 1.0
    for n in range(1, nrad):
        en = n ** 2 * (n + 2) ** 2 / (4 * (n + 1) ** 4 + 1)
        dn = 1 - en / d_prev
        f, df = aux(n)
        gi[:, n] = (f + math.sqrt(en / d_prev) * gi[:, n - 1]) / math.sqrt(dn)
        dgi[:, n] = (df + math.sqrt(en / d_prev) * dgi[:, n - 1]) / math.sqrt(dn)
        d_prev = dn
    g[inside], dg[inside] = gi, dgi
    return g, dg


def zbl_core(r, zi, zj, cut_out, cut_in, pre):
    """ZBL core repulsion as tabulated by ACE (shifted to vanish at the outer cutoff)."""
    a = 0.46850 / (zi ** 0.23 + zj ** 0.23)

    def e_and_d(rr):
        x = rr / a
        ex = np.exp(np.multiply.outer(x, ZBL_E))
        phi = (ZBL_C * ex).sum(-1)
        dphi = (ZBL_C * ZBL_E * ex).sum(-1)
        d2phi = (ZBL_C * ZBL_E ** 2 * ex).sum(-1)
        E = phi / rr
        dE = -phi / rr ** 2 + dphi / a / rr
        d2E = 2 * phi / rr ** 3 - 2 * dphi / a / rr ** 2 + d2phi / a ** 2 / rr
        return E, dE, d2E

    r = np.asarray(r, float)
    E, dE, _ = e_and_d(r)
    Ec, dEc, d2Ec = e_and_d(np.array(cut_out))
    dr = cut_out - cut_in
    A = (-3 * dEc + dr * d2Ec) / dr ** 2
    B = (2 * dEc - dr * d2Ec) / dr ** 3
    C = -Ec + 0.5 * dr * dEc - dr ** 2 * d2Ec / 12.
    S = np.where(r <= cut_in, C, np.where(r >= cut_out, 0.0, A / 3 * (r - cut_in) ** 3 + B / 4 * (r - cut_in) ** 4 + C))
    dS = np.where((r > cut_in) & (r < cut_out), A * (r - cut_in) ** 2 + B * (r - cut_in) ** 3, 0.0)
    pref = pre * K_COUL / 2 * zi * zj
    cr, dcr = pref * (E + S), pref * (dE + dS)
    fc, dfc = poly_cut(r, cut_out, cut_out - cut_in)
    cr, dcr = cr * fc, dcr * fc + cr * dfc
    out = r >= cut_out
    cr[out], dcr[out] = 0.0, 0.0
    return cr, dcr


def ylm(rhat):
    """Complex spherical harmonics (ACE convention, Y00 = 1) for unit vectors, all m = -l..l.

    Returns Y (P, NLM) and dY (P, NLM, 3), the gradient with respect to the unit vector projected
    perpendicular to it (so that dY/dr_vec = dY / |r|).
    """
    P = len(rhat)
    rx, ry, rz = rhat[:, 0], rhat[:, 1], rhat[:, 2]
    plm = {(0, 0): np.ones(P)}
    dplm = {(0, 0): np.zeros(P)}
    plm[(1, 0)] = math.sqrt(3) * rz
    dplm[(1, 0)] = np.full(P, math.sqrt(3))
    plm[(1, 1)] = np.full(P, -math.sqrt(1.5))
    dplm[(1, 1)] = np.zeros(P)
    for l in range(2, LMAX + 1):
        lsq, ld = l * l, 2 * l
        l1, l2 = 4 * lsq - 1, lsq - ld + 1
        for m in range(0, l - 1):
            a = math.sqrt(l1 / (lsq - m * m))
            b = -math.sqrt((l2 - m * m) / (4 * l2 - 1))
            plm[(l, m)] = a * (rz * plm[(l - 1, m)] + b * plm[(l - 2, m)])
            dplm[(l, m)] = a * (plm[(l - 1, m)] + rz * dplm[(l - 1, m)] + b * dplm[(l - 2, m)])
        t = math.sqrt(2 * (l - 1) + 3) * plm[(l - 1, l - 1)]
        plm[(l, l - 1)] = t * rz
        dplm[(l, l - 1)] = t
        plm[(l, l)] = -math.sqrt(1.0 + 0.5 / l) * plm[(l - 1, l - 1)]
        dplm[(l, l)] = np.zeros(P)
    Y = np.zeros((P, NLM), complex)
    dY = np.zeros((P, NLM, 3), complex)
    phase = rx + 1j * ry
    for l in range(LMAX + 1):
        for m in range(0, l + 1):
            phm = phase ** m
            mph = m * phase ** (m - 1) if m > 0 else np.zeros(P, complex)
            y = phm * plm[(l, m)]
            dyx = mph * plm[(l, m)]
            dyy = 1j * mph * plm[(l, m)]
            dyz = phm * dplm[(l, m)]
            rdy = rx * dyx + ry * dyy + rz * dyz
            d = np.stack([dyx - rdy * rx, dyy - rdy * ry, dyz - rdy * rz], axis=-1)
            Y[:, l * l + l + m] = y
            dY[:, l * l + l + m] = d
            if m > 0:
                f = (-1) ** m
                Y[:, l * l + l - m] = f * np.conj(y)
                dY[:, l * l + l - m] = f * np.conj(d)
    return Y, dY


def fs_embed(rho, mexp):
    if abs(mexp - 1.0) < 1e-10:
        return rho.copy(), np.ones_like(rho)
    a = np.abs(rho)
    ex = np.exp(-a)
    nx = 1. / mexp
    xoff = nx ** (nx / (1.0 - nx)) * ex
    yoff = nx ** (1 / (1.0 - nx)) * ex
    sign = np.where(np.signbit(rho), -1.0, 1.0)
    F = sign * ((xoff + a) ** mexp - yoff)
    DF = yoff + mexp * (-xoff + 1.0) * (xoff + a) ** (mexp - 1.)
    return F, DF


class ACE:
    def __init__(self, path):
        P = json.load(open(path))
        self.elements = P['elements']
        self.nel = len(self.elements)
        self.cut = np.array(P['cut'])
        self.dcut = np.array(P['dcut'])
        self.cut_in = np.array(P['cut_in'])
        self.dcut_in = np.array(P['dcut_in'])
        self.prehc = np.array(P['prehc'])
        self.nradbase, self.nradmax = P['nradbase'], P['nradmax']
        assert P['lmax'] == LMAX and P['inner_cutoff_type'] == 'zbl'
        self.crad = np.array(P['crad'])            # (nel, nel, nradmax, lmax+1, nradbase)
        self.rcmax = float(self.cut.max())
        self.emb = {int(k): v for k, v in P['embeddings'].items()}
        self.nA = self.nel * self.nradmax * NLM
        self.spec = {}
        for mu in range(self.nel):
            nd = self.emb[mu]['ndensity']
            r1 = P['rank1'][str(mu)]
            r1_idx = np.array([f['mus'][0] * self.nradbase + f['ns'][0] - 1 for f in r1], int)
            r1_c = np.array([f['ctildes'][0] if isinstance(f['ctildes'][0], list) else f['ctildes'] for f in r1], float)
            by_rank = {}
            for f in P['functions'][str(mu)]:
                rank = f['rank']
                ms = np.array(f['ms_combs']).reshape(-1, rank)
                ct = np.array(f['ctildes']).reshape(len(ms), nd)
                for comb, c in zip(ms, ct):
                    idx = [(f['mus'][t] * self.nradmax + f['ns'][t] - 1) * NLM + f['ls'][t] ** 2 + f['ls'][t] + comb[t]
                           for t in range(rank)]
                    by_rank.setdefault(rank, ([], []))
                    by_rank[rank][0].append(idx)
                    by_rank[rank][1].append(c)
            self.spec[mu] = {'nd': nd, 'r1_idx': r1_idx, 'r1_c': r1_c.reshape(len(r1), nd),
                             'ranks': {k: (np.array(v[0], int), np.array(v[1], float)) for k, v in by_rank.items()}}

    def radial(self, mi, mj, r):
        g, dg = sbessel(r, self.cut[mi, mj], self.nradbase)
        c = self.crad[mi, mj]                                   # (nradmax, lmax+1, nradbase)
        R = np.einsum('nlk,pk->pnl', c, g)
        dR = np.einsum('nlk,pk->pnl', c, dg)
        zi, zj = ZNUM[self.elements[mi]], ZNUM[self.elements[mj]]
        cr, dcr = zbl_core(r, zi, zj, self.cut[mi, mj], self.cut[mi, mj] - self.dcut[mi, mj], self.prehc[mi, mj])
        return g, dg, R, dR, cr, dcr

    def neighbours(self, cell, pos):
        cell = np.asarray(cell, float)
        inv = np.linalg.inv(cell)
        heights = 1.0 / np.linalg.norm(inv, axis=0)
        nrep = np.ceil(self.rcmax / heights).astype(int)
        frac = pos @ inv
        frac -= np.floor(frac)
        pos = frac @ cell
        shifts = np.array([(a, b, c) for a in range(-nrep[0], nrep[0] + 1) for b in range(-nrep[1], nrep[1] + 1)
                           for c in range(-nrep[2], nrep[2] + 1)]) @ cell
        I, J, V = [], [], []
        N = len(pos)
        for s in shifts:
            d = pos[None, :, :] + s[None, None, :] - pos[:, None, :]
            r2 = (d ** 2).sum(-1)
            mask = (r2 < self.rcmax ** 2) & (r2 > 1e-10)
            i, j = np.nonzero(mask)
            I.append(i)
            J.append(j)
            V.append(d[i, j])
        return pos, np.concatenate(I), np.concatenate(J), np.concatenate(V)

    def compute(self, cell, symbols, positions, forces=True):
        """Energy (eV), forces (eV/A), stress (eV/A^3, ASE convention) and per-atom energies."""
        sp = np.array([self.elements.index(s) for s in symbols], int)
        pos, I, J, V = self.neighbours(cell, np.asarray(positions, float))
        N = len(sp)
        r = np.linalg.norm(V, axis=1)
        mi, mj = sp[I], sp[J]
        keep = r < self.cut[mi, mj]
        I, J, V, r, mi, mj = I[keep], J[keep], V[keep], r[keep], mi[keep], mj[keep]
        P = len(r)
        rhat = V / r[:, None]
        g = np.zeros((P, self.nradbase))
        dg = np.zeros_like(g)
        R = np.zeros((P, self.nradmax, LMAX + 1))
        dR = np.zeros_like(R)
        cr, dcr = np.zeros(P), np.zeros(P)
        for a in range(self.nel):
            for b in range(self.nel):
                sel = (mi == a) & (mj == b)
                if sel.any():
                    g[sel], dg[sel], R[sel], dR[sel], cr[sel], dcr[sel] = self.radial(a, b, r[sel])
        Y, dY = ylm(rhat)
        RL = R[:, :, L_OF]                                      # (P, nradmax, NLM)
        dRL = dR[:, :, L_OF]
        A = np.zeros((N, self.nel, self.nradmax, NLM), complex)
        np.add.at(A, (I, mj), RL * Y[:, None, :])
        A1 = np.zeros((N, self.nel, self.nradbase))
        np.add.at(A1, (I, mj), g)
        Af = A.reshape(N, -1)
        A1f = A1.reshape(N, -1)
        rho_core = np.zeros(N)
        np.add.at(rho_core, I, cr)
        # ZBL transition per atom: nearest neighbour measured from its pair-specific inner cutoff
        d = r - (self.cut_in[mi, mj] - self.dcut_in[mi, mj])
        dmin = np.full(N, self.rcmax)
        np.minimum.at(dmin, I, d)
        energies = np.zeros(N)
        G = np.zeros((N, self.nA), complex)
        G1 = np.zeros((N, self.nel * self.nradbase))
        fcut = np.ones(N)
        dfcut = np.zeros(N)
        jmin_pair = np.full(N, -1)
        if P:
            order = np.lexsort((d, I))
            first = np.r_[True, I[order][1:] != I[order][:-1]]
            jmin_pair[I[order][first]] = order[first]
        for mu in range(self.nel):
            atoms = np.nonzero(sp == mu)[0]
            if not len(atoms):
                continue
            s = self.spec[mu]
            nd = s['nd']
            rho = A1f[atoms][:, s['r1_idx']] @ s['r1_c']      # (n, nd)
            grads = []
            for rank, (idx, ct) in s['ranks'].items():
                Ag = Af[atoms][:, idx]                           # (n, ncomb, rank)
                B = Ag.prod(axis=2)
                rho += np.real(B) @ ct
                grads.append((rank, idx, ct, Ag))
            fp = self.emb[mu]['FS_parameters']
            ev = np.zeros(len(atoms))
            dF = np.zeros((len(atoms), nd))
            for p in range(nd):
                F, DF = fs_embed(rho[:, p], fp[2 * p + 1])
                ev += fp[2 * p] * F
                dF[:, p] = fp[2 * p] * DF
            # ZBL switching
            fc, dfc = np.ones(len(atoms)), np.zeros(len(atoms))
            for k, at in enumerate(atoms):
                pj = jmin_pair[at]
                if pj >= 0:
                    dci = self.dcut_in[mu, mj[pj]]
                    tc = dci - dmin[at]
                    f1, d1 = poly_cut(np.array([tc]), dci, dci)
                    fc[k], dfc[k] = f1[0], -d1[0]
            energies[atoms] = ev * fc + rho_core[atoms] * (1 - fc)
            fcut[atoms], dfcut[atoms] = fc, (ev - rho_core[atoms]) * dfc
            if forces:
                dF *= fc[:, None]
                G1[atoms[:, None], s['r1_idx'][None, :]] += dF @ s['r1_c'].T
                for rank, idx, ct, Ag in grads:
                    w = dF @ ct.T                                # (n, ncomb)
                    for t in range(rank):
                        others = np.ones(Ag.shape[:2], complex)
                        for u in range(rank):
                            if u != t:
                                others = others * Ag[:, :, u]
                        np.add.at(G, (atoms[:, None], idx[None, :, t]), w * others)
        E = float(energies.sum())
        if not forces:
            return {'energy': E, 'energies': energies}
        Gp = G.reshape(N, self.nel, self.nradmax, NLM)[I, mj]   # (P, nradmax, NLM)
        H = (Gp * RL).sum(1)
        Hd = (Gp * dRL).sum(1)
        g1 = G1.reshape(N, self.nel, self.nradbase)[I, mj]
        scal = np.real((Hd * Y).sum(1)) + (g1 * dg).sum(1) + (1 - fcut[I]) * dcr
        vec = np.real(np.einsum('pk,pkx->px', H, dY)) / r[:, None]
        gp = scal[:, None] * rhat + vec                          # dE/d r_ij (vector from i to j)
        # derivative of the switching function through the nearest-neighbour distance
        for at in np.nonzero(dfcut != 0)[0]:
            pj = jmin_pair[at]
            gp[pj] += dfcut[at] * rhat[pj]    # dfcut already holds (E_ace - rho_core) * dfc/dr_min
        F = np.zeros((N, 3))
        np.add.at(F, I, gp)
        np.add.at(F, J, -gp)
        vol = abs(np.linalg.det(np.asarray(cell, float)))
        stress = np.einsum('px,py->xy', V, gp) / vol
        stress = 0.5 * (stress + stress.T)
        return {'energy': E, 'energies': energies, 'forces': F, 'stress': stress}


def _embed_rows(ace, mu, Af, A1f):
    """Per-atom ACE energy F(rho) for rows of flattened A (complex) and rank-1 A (real); no ZBL switching."""
    s = ace.spec[mu]
    rho = A1f[:, s['r1_idx']] @ s['r1_c']
    for rank, (idx, ct) in s['ranks'].items():
        rho += np.real(Af[:, idx].prod(axis=2)) @ ct
    fp = ace.emb[mu]['FS_parameters']
    ev = np.zeros(len(Af))
    for p in range(s['nd']):
        ev += fp[2 * p] * fs_embed(rho[:, p], fp[2 * p + 1])[0]
    return ev


def scan_atom(ace, cell, symbols, positions, index, points, chunk=2000):
    """Total energy when atom `index` is placed at each of `points` (Cartesian, A) with all other atoms frozen.

    Exact for the ACE energy; points must stay outside the inner (ZBL) cutoff of every other atom.
    """
    cell = np.asarray(cell, float)
    positions = np.asarray(positions, float)
    points = np.atleast_2d(np.asarray(points, float))
    sp = np.array([ace.elements.index(s) for s in symbols], int)
    mk = sp[index]
    keep = np.arange(len(sp)) != index
    rest_sym = [s for i, s in enumerate(symbols) if i != index]
    base = ace.compute(cell, rest_sym, positions[keep], forces=False)
    E_rest = base['energy']
    rest_idx = np.nonzero(keep)[0]
    # A of the remaining atoms
    pos, I, J, V = ace.neighbours(cell, positions[keep])
    r = np.linalg.norm(V, axis=1)
    spr = sp[keep]
    mi, mj = spr[I], spr[J]
    ok = r < ace.cut[mi, mj]
    I, J, V, r, mi, mj = I[ok], J[ok], V[ok], r[ok], mi[ok], mj[ok]
    N = len(spr)
    A = np.zeros((N, ace.nel, ace.nradmax, NLM), complex)
    A1 = np.zeros((N, ace.nel, ace.nradbase))
    for a in range(ace.nel):
        for b in range(ace.nel):
            sel = (mi == a) & (mj == b)
            if sel.any():
                g, _, R, _, _, _ = ace.radial(a, b, r[sel])
                Y, _ = ylm(V[sel] / r[sel][:, None])
                np.add.at(A, (I[sel], b), R[:, :, L_OF] * Y[:, None, :])
                np.add.at(A1, (I[sel], b), g)
    # periodic images of the remaining atoms that can see the scan region
    inv = np.linalg.inv(cell)
    heights = 1.0 / np.linalg.norm(inv, axis=0)
    span = np.ptp(points, axis=0).max() if len(points) > 1 else 0.0
    nrep = np.ceil((ace.rcmax + span) / heights).astype(int) + 1
    shifts = np.array([(a, b, c) for a in range(-nrep[0], nrep[0] + 1) for b in range(-nrep[1], nrep[1] + 1)
                       for c in range(-nrep[2], nrep[2] + 1)]) @ cell
    centre = points.mean(axis=0)
    img = (pos[None, :, :] + shifts[:, None, :]).reshape(-1, 3)
    owner = np.tile(np.arange(N), len(shifts))
    near = np.linalg.norm(img - centre, axis=1) < ace.rcmax + span + 1e-6
    img, owner = img[near], owner[near]
    cand = np.unique(owner)
    E0 = {mu: None for mu in range(ace.nel)}
    Af0 = A.reshape(N, -1)
    A1f0 = A1.reshape(N, -1)
    e_base = np.zeros(N)
    for mu in range(ace.nel):
        rows = cand[spr[cand] == mu]
        if len(rows):
            e_base[rows] = _embed_rows(ace, mu, Af0[rows], A1f0[rows])
    out = np.zeros(len(points))
    for c0 in range(0, len(points), chunk):
        pts = points[c0:c0 + chunk]
        Pn = len(pts)
        d = img[None, :, :] - pts[:, None, :]                 # from moving atom to image (Pn, M, 3)
        rr = np.linalg.norm(d, axis=2)
        if (rr < ace.cut_in[mk][spr[owner]][None, :]).any():
            raise ValueError('scan point inside the inner (ZBL) cutoff of another atom')
        dE = np.zeros(Pn)
        # energy of the moving atom
        Ak = np.zeros((Pn, ace.nel, ace.nradmax, NLM), complex)
        A1k = np.zeros((Pn, ace.nel, ace.nradbase))
        # change of the neighbours' A
        dA = {}
        for b in range(ace.nel):
            pp, mm = np.nonzero((rr < ace.cut[mk, b]) & (spr[owner][None, :] == b))
            if len(pp):
                rv = rr[pp, mm]
                dv = d[pp, mm]
                g, _, R, _, _, _ = ace.radial(mk, b, rv)
                Y, _ = ylm(dv / rv[:, None])
                np.add.at(Ak, (pp, b), R[:, :, L_OF] * Y[:, None, :])
                np.add.at(A1k, (pp, b), g)
                # neighbour i sees the moving atom along -d (radial functions mu_i -> mu_k)
                g2, _, R2, _, _, _ = ace.radial(b, mk, rv)
                Y2, _ = ylm(-dv / rv[:, None])
                dA[b] = (pp, owner[mm], R2[:, :, L_OF] * Y2[:, None, :], g2)
        dE += _embed_rows(ace, mk, Ak.reshape(Pn, -1), A1k.reshape(Pn, -1))
        for b, (pp, ii, contrib, g2) in dA.items():
            key = pp * N + ii
            uk, inv_k = np.unique(key, return_inverse=True)
            Au = np.zeros((len(uk), ace.nel, ace.nradmax, NLM), complex)
            A1u = np.zeros((len(uk), ace.nel, ace.nradbase))
            np.add.at(Au, (inv_k, mk), contrib)
            np.add.at(A1u, (inv_k, mk), g2)
            atoms = uk % N
            pidx = uk // N
            Af = Af0[atoms] + Au.reshape(len(uk), -1)
            A1f = A1f0[atoms] + A1u.reshape(len(uk), -1)
            enew = _embed_rows(ace, b, Af, A1f)
            np.add.at(dE, pidx, enew - e_base[atoms])
        out[c0:c0 + Pn] = E_rest + dE
    return out
