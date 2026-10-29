"""Pure-numpy evaluator for NEP4 and NEP4-charge (qNEP) potentials.

Written for this task from the published NEP / qNEP formulas (Fan et al., PRB 104, 104309 (2021);
GPUMD documentation) so that the lab needs nothing beyond numpy. It supports the subset used by
the released Ca3Ti2O7 potentials: nep4 / nep4_charge1, one global radial and angular cutoff,
l_max 3-body up to 4 with the 4-body (q_222) and 5-body (q_1111) terms, one hidden layer.

E(positions, cell) -> energy (eV), forces (eV/A), virial (eV), and for charge models the
per-atom charges (e) and Born effective charges (e, 3x3 per atom, row = dipole component,
column = displacement component, following the NEP_CPU convention).
"""
import math
import numpy as np

K_C = 14.399645          # e^2 / (4 pi eps0) in eV A
C3B = np.array([
    0.238732414637843, 0.119366207318922, 0.119366207318922, 0.099471839432435, 0.596831036594608,
    0.596831036594608, 0.149207759148652, 0.149207759148652, 0.139260575205408, 0.104445431404056,
    0.104445431404056, 1.044454314040563, 1.044454314040563, 0.174075719006761, 0.174075719006761,
    0.011190581936149, 0.223811638722978, 0.223811638722978, 0.111905819361489, 0.111905819361489,
    1.566681471060845, 1.566681471060845, 0.195835183882606, 0.195835183882606])
C4B = np.array([-0.007499480826664, -0.134990654879954, 0.067495327439977, 0.404971964639861,
                -0.809943929279723])
C5B = np.array([0.026596810706114, 0.053193621412227, 0.026596810706114])
ZC = {1: [[0.0, 1.0], [1.0, 0.0]],
      2: [[-1.0, 0.0, 3.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]],
      3: [[0.0, -3.0, 0.0, 5.0], [-1.0, 0.0, 5.0, 0.0], [0.0, 1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0]],
      4: [[3.0, 0.0, -30.0, 0.0, 35.0], [0.0, -3.0, 0.0, 7.0, 0.0], [-1.0, 0.0, 7.0, 0.0, 0.0],
          [0.0, 1.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0]]}


def _scatter(idx, vals, n):
    """Sum rows of vals into n bins given by idx (fast replacement for np.add.at)."""
    vals = np.asarray(vals, float)
    if vals.ndim == 1:
        return np.bincount(idx, vals, minlength=n)
    flat = vals.reshape(len(idx), -1)
    out = np.stack([np.bincount(idx, flat[:, k], minlength=n) for k in range(flat.shape[1])], axis=1)
    return out.reshape((n,) + vals.shape[1:])


def _erfc(x):
    return np.array([math.erfc(v) for v in np.ravel(x)]).reshape(np.shape(x))


class NEP:
    def __init__(self, filename):
        tok = [l.split() for l in open(filename) if l.strip()]
        head = tok[0]
        if head[0] not in ('nep4', 'nep4_charge1'):
            raise ValueError('unsupported NEP model ' + head[0])
        self.charge = head[0] == 'nep4_charge1'
        self.ntypes = int(head[1])
        self.elements = head[2:2 + self.ntypes]
        self.rc_r, self.rc_a = float(tok[1][1]), float(tok[1][2])
        self.nr, self.na = int(tok[2][1]) + 1, int(tok[2][2]) + 1          # n_max + 1
        self.kr, self.ka = int(tok[3][1]) + 1, int(tok[3][2]) + 1          # basis_size + 1
        lt = [int(x) for x in tok[4][1:]]
        self.lmax = lt[0]
        if lt[1] > 1:                       # old format: 2 = has q_222; next = 1 if q_1111
            self.q222, self.q1111 = int(lt[1] == 2), int(lt[2] == 1)
        else:
            self.q222, self.q1111 = lt[1], lt[2]
        self.numL = self.lmax + self.q222 + self.q1111
        self.nneu = int(tok[5][1])
        self.dim = self.nr + self.na * self.numL
        vals = np.array([float(t[0]) for t in tok[6:]])
        T, H, D = self.ntypes, self.nneu, self.dim
        p = 0
        self.w0, self.b0, self.w1, self.w1q = [], [], [], []
        for t in range(T):
            self.w0.append(vals[p:p + H * D].reshape(H, D)); p += H * D
            self.b0.append(vals[p:p + H]); p += H
            self.w1.append(vals[p:p + H]); p += H
            if self.charge:
                self.w1q.append(vals[p:p + H]); p += H
        self.sqrt_eps_inf = 1.0
        if self.charge:
            self.sqrt_eps_inf = vals[p]; p += 1
        self.b1 = vals[p]; p += 1
        ncr = T * T * self.nr * self.kr
        nca = T * T * self.na * self.ka
        # c[n, k, t1, t2]
        self.c_r = vals[p:p + ncr].reshape(self.nr, self.kr, T, T); p += ncr
        self.c_a = vals[p:p + nca].reshape(self.na, self.ka, T, T); p += nca
        self.q_scaler = vals[p:p + D]; p += D
        if p != len(vals):
            raise ValueError('parameter count mismatch: %d read, %d in file' % (p, len(vals)))
        if self.charge:
            self.alpha = math.pi / self.rc_r
        self.nabc = sum(2 * L + 1 for L in range(1, self.lmax + 1))
        self._nl = None

    # ------------------------------------------------------------------ basis functions
    @staticmethod
    def _fn(d, rc, nk):
        x = d / rc
        fc = np.where(d < rc, 0.5 * np.cos(np.pi * x) + 0.5, 0.0)
        fcp = np.where(d < rc, -0.5 * np.pi * np.sin(np.pi * x) / rc, 0.0)
        t = 2.0 * (x - 1.0) ** 2 - 1.0
        dt = 4.0 * (x - 1.0) / rc
        T = [np.ones_like(d), t]
        Tp = [np.zeros_like(d), np.ones_like(d)]
        for m in range(2, nk):
            Tp.append(2 * T[-1] + 2 * t * Tp[-1] - Tp[-2])
            T.append(2 * t * T[-1] - T[-2])
        T, Tp = np.array(T[:nk]).T, np.array(Tp[:nk]).T
        f = 0.5 * (T + 1.0)
        fp = 0.5 * Tp * dt[:, None]
        return f * fc[:, None], fp * fc[:, None] + f * fcp[:, None]

    def _ylm(self, u):
        """Real spherical-harmonic polynomials of NEP (unit vector u) and their u-gradients."""
        x, y, z = u[:, 0], u[:, 1], u[:, 2]
        P = len(u)
        Y = np.zeros((P, self.nabc)); dY = np.zeros((P, self.nabc, 3))
        w = x + 1j * y
        idx = 0
        for L in range(1, self.lmax + 1):
            for n1 in range(L + 1):
                zp = np.zeros(P); dzp = np.zeros(P)
                for n2 in range(L - n1 + 1):
                    if (L + n1 + n2) % 2:
                        continue
                    c = ZC[L][n1][n2]
                    zp += c * z ** n2
                    if n2 > 0:
                        dzp += c * n2 * z ** (n2 - 1)
                if n1 == 0:
                    Y[:, idx] = zp; dY[:, idx, 2] = dzp; idx += 1
                else:
                    wm = w ** n1
                    wm1 = n1 * w ** (n1 - 1)
                    Y[:, idx] = zp * wm.real
                    dY[:, idx, 0] = zp * wm1.real
                    dY[:, idx, 1] = -zp * wm1.imag
                    dY[:, idx, 2] = dzp * wm.real
                    idx += 1
                    Y[:, idx] = zp * wm.imag
                    dY[:, idx, 0] = zp * wm1.imag
                    dY[:, idx, 1] = zp * wm1.real
                    dY[:, idx, 2] = dzp * wm.imag
                    idx += 1
        return Y, dY

    def _q_angular(self, s):
        """s: (N, na, nabc) -> q (N, na, numL) and dq/ds (N, na, numL, nabc)."""
        N = s.shape[0]
        q = np.zeros((N, self.na, self.numL)); dq = np.zeros((N, self.na, self.numL, self.nabc))
        for L in range(1, self.lmax + 1):
            st, nt = L * L - 1, 2 * L + 1
            w = C3B[st:st + nt] * np.r_[1.0, 2.0 * np.ones(nt - 1)]
            q[:, :, L - 1] = (w * s[:, :, st:st + nt] ** 2).sum(-1)
            dq[:, :, L - 1, st:st + nt] = 2 * w * s[:, :, st:st + nt]
        li = self.lmax
        if self.q222:
            s3, s4, s5, s6, s7 = (s[:, :, k] for k in range(3, 8))
            q[:, :, li] = (C4B[0] * s3 ** 3 + C4B[1] * s3 * (s4 ** 2 + s5 ** 2) + C4B[2] * s3 * (s6 ** 2 + s7 ** 2)
                           + C4B[3] * s6 * (s5 ** 2 - s4 ** 2) + C4B[4] * s4 * s5 * s7)
            dq[:, :, li, 3] = 3 * C4B[0] * s3 ** 2 + C4B[1] * (s4 ** 2 + s5 ** 2) + C4B[2] * (s6 ** 2 + s7 ** 2)
            dq[:, :, li, 4] = 2 * C4B[1] * s3 * s4 - 2 * C4B[3] * s6 * s4 + C4B[4] * s5 * s7
            dq[:, :, li, 5] = 2 * C4B[1] * s3 * s5 + 2 * C4B[3] * s6 * s5 + C4B[4] * s4 * s7
            dq[:, :, li, 6] = 2 * C4B[2] * s3 * s6 + C4B[3] * (s5 ** 2 - s4 ** 2)
            dq[:, :, li, 7] = 2 * C4B[2] * s3 * s7 + C4B[4] * s4 * s5
            li += 1
        if self.q1111:
            s0, s1, s2 = s[:, :, 0], s[:, :, 1], s[:, :, 2]
            a, b = s0 ** 2, s1 ** 2 + s2 ** 2
            q[:, :, li] = C5B[0] * a * a + C5B[1] * a * b + C5B[2] * b * b
            dq[:, :, li, 0] = 4 * C5B[0] * a * s0 + 2 * C5B[1] * b * s0
            dq[:, :, li, 1] = 2 * C5B[1] * a * s1 + 4 * C5B[2] * b * s1
            dq[:, :, li, 2] = 2 * C5B[1] * a * s2 + 4 * C5B[2] * b * s2
        return q, dq

    # ------------------------------------------------------------------ neighbours
    @staticmethod
    def _pairs(pos, cell, rc):
        """All ordered pairs (i, j, image shift) with |r_j + shift - r_i| < rc, i != j or shift != 0."""
        inv = np.linalg.inv(cell)
        frac = pos @ inv
        frac -= np.floor(frac)
        pos = frac @ cell
        vol = abs(np.linalg.det(cell))
        heights = [vol / np.linalg.norm(np.cross(cell[(k + 1) % 3], cell[(k + 2) % 3])) for k in range(3)]
        nmax = [int(math.ceil(rc / h)) for h in heights]
        I, J, S = [], [], []
        N = len(pos)
        for a in range(-nmax[0], nmax[0] + 1):
            for b in range(-nmax[1], nmax[1] + 1):
                for c in range(-nmax[2], nmax[2] + 1):
                    sh = np.array([a, b, c], float) @ cell
                    d = pos[None, :, :] + sh - pos[:, None, :]
                    r2 = (d * d).sum(-1)
                    m = r2 < rc * rc
                    if a == 0 and b == 0 and c == 0:
                        m &= ~np.eye(N, dtype=bool)
                    i, j = np.nonzero(m)
                    I.append(i); J.append(j); S.append(np.repeat([[a, b, c]], len(i), axis=0))
        return pos, np.concatenate(I), np.concatenate(J), np.concatenate(S).astype(float)

    # ------------------------------------------------------------------ main evaluation
    def neighbours(self, pos, cell, skin=0.0):
        """Pair list (I, J, image shifts) for cutoff rc_radial + skin; positions are not wrapped."""
        cell = np.asarray(cell, float)
        pos = np.asarray(pos, float)
        wrapped, I, J, S = self._pairs(pos, cell, self.rc_r + skin)
        # express shifts relative to the unwrapped positions
        shift = np.round((wrapped - pos) @ np.linalg.inv(cell))
        S = S + shift[J] - shift[I]
        return I, J, S

    def compute(self, types, pos, cell, want_bec=False, pairs=None):
        types = np.asarray(types, int)
        cell = np.asarray(cell, float)
        pos = np.asarray(pos, float)
        I, J, S = pairs if pairs is not None else self.neighbours(pos, cell)
        N = len(types)
        r12 = pos[J] + S @ cell - pos[I]
        d = np.sqrt((r12 * r12).sum(1))
        keep = d < self.rc_r
        if not keep.all():
            I, J, S, r12, d = I[keep], J[keep], S[keep], r12[keep], d[keep]
        u = r12 / d[:, None]
        ti, tj = types[I], types[J]
        # radial descriptors
        fr, fpr = self._fn(d, self.rc_r, self.kr)                  # (P, k)
        cr = self.c_r[:, :, ti, tj]                               # (n, k, P)
        g_r = np.einsum('pk,nkp->pn', fr, cr)
        gp_r = np.einsum('pk,nkp->pn', fpr, cr)
        qr = _scatter(I, g_r, N)
        # angular descriptors
        A = d < self.rc_a
        Ia, Ja, da, ua = I[A], J[A], d[A], u[A]
        fa, fpa = self._fn(da, self.rc_a, self.ka)
        ca = self.c_a[:, :, ti[A], tj[A]]
        g_a = np.einsum('pk,nkp->pn', fa, ca)
        gp_a = np.einsum('pk,nkp->pn', fpa, ca)
        Y, dYu = self._ylm(ua)
        # gradient of Y(u) wrt r: (I - u u^T) dY/du / d
        dYr = (dYu - ua[:, None, :] * (dYu * ua[:, None, :]).sum(-1, keepdims=True)) / da[:, None, None]
        s = _scatter(Ia, g_a[:, :, None] * Y[:, None, :], N)
        qa, dqa = self._q_angular(s)
        q = np.concatenate([qr, qa.transpose(0, 2, 1).reshape(N, -1)], axis=1)   # angular index l*na + n
        qs = q * self.q_scaler
        # neural network
        E_at = np.zeros(N); Fp = np.zeros((N, self.dim)); chg = np.zeros(N); dchg = np.zeros((N, self.dim))
        for t in range(self.ntypes):
            m = types == t
            if not m.any():
                continue
            x = np.tanh(qs[m] @ self.w0[t].T - self.b0[t])
            der = 1 - x * x
            E_at[m] = x @ self.w1[t] - self.b1
            Fp[m] = ((der * self.w1[t]) @ self.w0[t]) * self.q_scaler
            if self.charge:
                chg[m] = x @ self.w1q[t]
                dchg[m] = ((der * self.w1q[t]) @ self.w0[t]) * self.q_scaler
        energy = E_at.sum()
        force = np.zeros((N, 3)); virial = np.zeros((3, 3))
        out = {}
        if self.charge:
            chg = chg - chg.mean()
            ec, fc_, vc, D = self._coulomb(chg, pos, cell, I, J, r12, d)
            D = D - D.mean()
            energy += ec; force += fc_; virial += vc
            Fp_tot = Fp + D[:, None] * dchg
            out['charges'] = chg
        else:
            Fp_tot = Fp

        def pair_forces(G):
            """f12 = dE_i/dr12 for radial and angular pairs given per-atom dE/dq_raw G (N, dim)."""
            f_r = (G[I, :self.nr] * gp_r).sum(1)[:, None] * u
            Ga = G[:, self.nr:].reshape(N, self.numL, self.na)            # (N, l, n)
            Gs = np.einsum('iln,inlc->inc', Ga, dqa)                       # dE/ds (N, n, abc)
            Gp = Gs[Ia]                                                   # (Pa, n, abc)
            Asum = np.einsum('pnc,pc->pn', Gp, Y)
            Bsum = np.einsum('pnc,pcx->pnx', Gp, dYr)
            f_a = (gp_a * Asum).sum(1)[:, None] * ua + np.einsum('pn,pnx->px', g_a, Bsum)
            return f_r, f_a

        f_r, f_a = pair_forces(Fp_tot)
        for ii, jj, f, rr in ((I, J, f_r, r12), (Ia, Ja, f_a, r12[A])):
            force += _scatter(ii, f, N) - _scatter(jj, f, N)
            virial -= rr.T @ f
        out.update(energy=float(energy), forces=force, virial=virial, descriptor=qs)
        if self.charge and want_bec:
            b_r, b_a = pair_forces(dchg)
            bec = np.zeros((N, 3, 3))
            bec[:, [0, 1, 2], [0, 1, 2]] = chg[:, None]
            for ii, jj, f, rr in ((I, J, b_r, r12), (Ia, Ja, b_a, r12[A])):
                t = 0.5 * rr[:, :, None] * f[:, None, :]
                bec += _scatter(ii, t, N) - _scatter(jj, t, N)
            out['bec'] = bec * self.sqrt_eps_inf
        return out

    def _coulomb(self, q, pos, cell, I, J, r12, d):
        N = len(q)
        al = self.alpha
        # real space (ordered pairs within the radial cutoff) + self term
        erfc_r = _erfc(al * d) / d
        D = -2 * al / math.sqrt(math.pi) * q
        D += _scatter(I, q[J] * erfc_r, N)
        e = 0.5 * (q[I] * q[J] * erfc_r).sum() - al / math.sqrt(math.pi) * (q * q).sum()
        f2 = -(erfc_r + 2 * al / math.sqrt(math.pi) * np.exp(-al * al * d * d)) * q[I] * q[J] / (d * d)
        fpair = 0.5 * f2[:, None] * r12                 # force on i from the i-j term counted per ordered pair
        force = _scatter(I, fpair, N) - _scatter(J, fpair, N)
        vir = -0.5 * (r12.T @ (f2[:, None] * r12))   # W = -dE/d(strain)
        # reciprocal space
        two_pi = 6.2831853
        a1, a2, a3 = cell
        det = np.dot(a1, np.cross(a2, a3))
        b1, b2, b3 = (np.cross(a2, a3), np.cross(a3, a1), np.cross(a1, a2))
        b1, b2, b3 = (b * two_pi / det for b in (b1, b2, b3))
        vk = two_pi ** 3 / abs(det)
        area = lambda x, y: np.linalg.norm(np.cross(x, y))
        n1m = int(al * two_pi * area(b2, b3) / vk)
        n2m = int(al * two_pi * area(b3, b1) / vk)
        n3m = int(al * two_pi * area(b1, b2) / vk)
        ksq_max = two_pi * two_pi * al * al
        n = np.array([(a, b, c) for a in range(0, n1m + 1) for b in range(-n2m, n2m + 1)
                      for c in range(-n3m, n3m + 1)
                      if not (a == 0 and b == 0 and c == 0) and not (a == 0 and b < 0)
                      and not (a == 0 and b == 0 and c < 0)], float)
        k = n @ np.array([b1, b2, b3])
        ksq = (k * k).sum(1)
        k, ksq = k[ksq < ksq_max], ksq[ksq < ksq_max]
        G = 2.0 * abs(two_pi / det) / ksq * np.exp(-ksq * 0.25 / (al * al))
        kr = pos @ k.T                                  # (N, K)
        cs, sn = np.cos(kr), np.sin(kr)
        Sr = q @ cs; Si = -(q @ sn)
        GSE = G * (Sr * cs - Si * sn)                   # (N, K)
        e += (q[:, None] * GSE).sum()
        D += 2.0 * GSE.sum(1)
        imag = G * (Sr * sn + Si * cs)
        force += 2.0 * q[:, None] * (imag @ k)
        Etk = G * (Sr * Sr + Si * Si)                   # per-k energy
        fac = 2.0 * 0.25 / (al * al) + 2.0 / ksq
        vir_k = np.eye(3) * Etk.sum() - np.einsum('k,k,kx,ky->xy', Etk, fac, k, k)
        return K_C * e, K_C * force, K_C * vir + K_C * vir_k, K_C * D
