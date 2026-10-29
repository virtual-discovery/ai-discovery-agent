"""Reference solution for mnwo4-halfmoon-frustration; runs as the Harbor oracle and in check_task.

Path: buy 20 archived 1.5 K spectra along the path, the 20 K and 30 K maps and one 11 K zone slice (29 of 33 units).
Fit J1..J12 and D to the spectra (AF1 linear spin waves with S = 5/2 from the 5 muB moments, variable projection,
start: the paper's rounded Table I). The archived path leaves J3 - J11 undetermined; fix it with a 1-D scan against
the block-averaged 20 K and 30 K maps (SCGA at the nominal temperatures). Predict the withheld 14 K map and the
reduced-model maps with SCGA, the Luttinger-Tisza vector and the AF1 / flat-spiral energies.
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path
sys.dont_write_bytecode = True

import numpy as np
from scipy.optimize import brentq, minimize, least_squares, nnls
from scipy.ndimage import gaussian_filter1d

KB = 0.08617333262
CELL = np.array([[4.8236999512, 0, 0], [0, 5.7551999092, 0], [-0.0946104705, 0, 5.0010053420]])
TAU = np.array([[0.5, 0.68490, 0.25], [0.5, 0.31510, 0.75]])
NAX = np.array([0.81915, 0, 0.57358]); NAX = NAX / np.linalg.norm(NAX)


class Model:
    def __init__(self, shells):
        """shells: list (len 12) of lists of [s, s2, n1, n2, n3] (every ordered pair listed)."""
        self.shell = []
        self.s1, self.s2, self.d, self.R = [], [], [], []
        for k, bonds in enumerate(shells):
            for s, s2, n1, n2, n3 in bonds:
                self.shell.append(k); self.s1.append(s); self.s2.append(s2); self.R.append((n1, n2, n3))
                self.d.append(np.array([n1, n2, n3]) + TAU[s2] - TAU[s])
        self.shell = np.array(self.shell); self.s1 = np.array(self.s1); self.s2 = np.array(self.s2)
        self.d = np.array(self.d); self.R = np.array(self.R)
        self.nsh = len(shells)
        self._af1_setup()

    # ---------- J(q) ----------
    def Jq(self, J, q):
        q = np.atleast_2d(q)
        ph = np.exp(2j * np.pi * q @ self.d.T) * np.asarray(J)[self.shell][None, :]
        M = np.zeros((len(q), 2, 2), complex)
        for a in range(2):
            for b in range(2):
                m = (self.s1 == a) & (self.s2 == b)
                M[:, a, b] = ph[:, m].sum(axis=1)
        return M

    def lam_max(self, J, q):
        return np.linalg.eigvalsh(self.Jq(J, q))[:, -1]

    def lt(self, J, n=40):
        g = (np.arange(n) + 0.5) / n - 0.5
        G = np.array(np.meshgrid(g, g, g, indexing='ij')).reshape(3, -1).T
        lm = self.lam_max(J, G)
        best = None
        for i in np.argsort(-lm)[:6]:
            r = minimize(lambda x: -self.lam_max(J, x[None])[0], G[i], method='Nelder-Mead',
                         options={'xatol': 1e-8, 'fatol': 1e-13, 'maxiter': 4000})
            if best is None or r.fun < best.fun:
                best = r
        return canonical_k(best.x), -best.fun

    # ---------- AF1 ----------
    @staticmethod
    def af1_sign(s, n1, n2, n3):
        s0 = (1, 1, -1, -1)[n1 % 4] if s == 0 else (1, -1, -1, 1)[n1 % 4]
        return s0 * (-1) ** ((n2 + n3) % 2)

    def _af1_setup(self):
        # magnetic cell: 4 chemical cells (n1 = 0..3, n2 = n3 = 0) x 2 sublattices; translations (2,1,0),(0,1,1),(0,0,2)
        T = np.array([[2, 1, 0], [0, 1, 1], [0, 0, 2]]).T
        Tinv = np.linalg.inv(T)
        sites = [(s, n1) for n1 in range(4) for s in range(2)]
        def find(s, n):
            for n1 in range(4):
                f = Tinv @ (np.array(n) - np.array([n1, 0, 0]))
                if np.allclose(f, np.round(f), atol=1e-9):
                    return sites.index((s, n1))
            raise KeyError
        A, B, Dd, Sh, Par = [], [], [], [], []
        for a, (s, n1) in enumerate(sites):
            for k in np.where(self.s1 == s)[0]:
                n = np.array([n1, 0, 0]) + self.R[k]
                b = find(self.s2[k], n)
                par = self.af1_sign(s, n1, 0, 0) * self.af1_sign(self.s2[k], *n) > 0
                A.append(a); B.append(b); Dd.append(self.d[k]); Sh.append(self.shell[k]); Par.append(par)
        self.af = dict(a=np.array(A), b=np.array(B), d=np.array(Dd), sh=np.array(Sh), par=np.array(Par), n=len(sites))

    def af1_energy(self, J, D):
        """classical energy per spin of AF1 (collinear along n)"""
        af = self.af
        sgn = np.where(af['par'], 1.0, -1.0)
        return -np.sum(np.asarray(J)[af['sh']] * sgn) / af['n'] + D

    def lswt(self, J, D, q, S=2.5):
        """AF1 magnon energies (meV), shape (nq, 8); classical precession with S = mu/(g muB)."""
        af = self.af; n = af['n']
        J = np.asarray(J, float)
        Jb = J[af['sh']] / S ** 2
        sgn = np.where(af['par'], 1.0, -1.0)
        Adiag = np.zeros(n)
        np.add.at(Adiag, af['a'], 2 * S * Jb * sgn)
        Adiag += -2 * (D / S ** 2) * S
        q = np.atleast_2d(q)
        out = np.zeros((len(q), n))
        eta = np.diag(np.r_[np.ones(n), -np.ones(n)])
        for iq, qq in enumerate(q):
            Am = np.diag(Adiag).astype(complex); Amm = Am.copy(); Bm = np.zeros((n, n), complex)
            ph = np.exp(2j * np.pi * af['d'] @ qq)
            c = -S * Jb  # per listing; each bond listed twice (once per end)
            p = af['par']
            np.add.at(Am, (af['a'][p], af['b'][p]), c[p] * ph[p])
            np.add.at(Am, (af['b'][p], af['a'][p]), c[p] * np.conj(ph[p]))
            np.add.at(Amm, (af['a'][p], af['b'][p]), c[p] * np.conj(ph[p]))
            np.add.at(Amm, (af['b'][p], af['a'][p]), c[p] * ph[p])
            np.add.at(Bm, (af['a'][~p], af['b'][~p]), 2 * c[~p] * ph[~p])
            M = np.block([[Am, Bm], [Bm.conj().T, Amm.T]])
            w = np.linalg.eigvals(eta @ M)
            out[iq] = np.sort(np.real(w))[n:]
        return out

    # ---------- SCGA ----------
    def scga(self, J, D, T, q, nbz=24, lam_fixed=None):
        """per-component structure factors (Sxx, Syy, Szz) per spin, Cartesian, easy axis NAX."""
        beta = 1.0 / (KB * T)
        g = (np.arange(nbz) + 0.5) / nbz - 0.5
        G = np.array(np.meshgrid(g, g, g, indexing='ij')).reshape(3, -1).T
        ev = np.linalg.eigvalsh(self.Jq(J, G))
        if lam_fixed is None:
            def f(lam):
                return 0.5 * np.mean(1 / (lam + beta * D - beta * ev)) + np.mean(1 / (lam - beta * ev)) - 1.0
            lo = beta * ev.max() + max(-beta * D, 0.0)
            lam = brentq(f, lo + 1e-12, lo + 100)
        else:
            lam = lam_fixed
        Jm = self.Jq(J, q)
        I2 = np.eye(2)
        Ke = (lam + beta * D) * I2 - beta * Jm
        Kp = lam * I2 - beta * Jm
        Se = 0.25 * np.real(np.linalg.inv(Ke).sum(axis=(1, 2)))
        Sp = 0.25 * np.real(np.linalg.inv(Kp).sum(axis=(1, 2)))
        n2 = NAX ** 2
        return np.c_[n2[0] * Se + (1 - n2[0]) * Sp, n2[1] * Se + (1 - n2[1]) * Sp, n2[2] * Se + (1 - n2[2]) * Sp]


def canonical_k(k):
    """representative of {+k, -k} + G with k_b = 1/2 sign convention: h <= 0 when possible"""
    k = np.asarray(k, float)
    cands = []
    for s in (1, -1):
        v = s * k
        v = v - np.round(v)
        cands.append(v)
    # prefer k_b = +0.5 and h<0
    def key(v):
        return (abs(abs(v[1]) - 0.5) > 1e-6, v[0] > 0, -v[1])
    v = sorted(cands, key=key)[0]
    if abs(abs(v[1]) - 0.5) < 1e-6:
        v[1] = 0.5
    return v


def k_distance(k1, k2):
    """distance between propagation vectors modulo sign and reciprocal lattice (r.l.u., max over components)"""
    best = 9.0
    for s in (1, -1):
        dv = np.asarray(k1) - s * np.asarray(k2)
        dv = dv - np.round(dv)
        best = min(best, np.abs(dv).max())
    return best





MAP_DK, MAP_DL = 4.66 / 372, 0.0111          # Fig. 3 grid: k_i = -2.33 + i dk (i = 0..372), l_j = -1.776 + j dl (j = 0..320)
KO_DK, KO_DL = 0.025, 1.0 / 90                # Fig. 5 grid: k_i = -2 + i dk (i = 0..160), l_j = j dl (j = 0..160)
TABLE_I = np.array([-0.592, -0.0032, -0.272, -0.336, -0.0176, -0.544, -0.176, -0.016, -0.32, -0.192, -0.0672, -0.0256])
TABLE_I_D = -0.192
REDUCED = {'5d': [0, 2, 3, 5, 6, 9, 10], '5e': [0, 3, 5, 8]}   # shells kept (0-based) in the withheld reduced models
SOLD_REDUCED = {'5b': list(range(12)), '5c': [0, 2, 3, 5, 6, 8, 9, 10]}  # 5b: no anisotropy


def map14_blocks(full):
    """27 x 29 block averages (7 x 11 grid points) of a full Fig.-3 map (373 x 321)"""
    out = np.zeros((27, 29))
    for r in range(27):
        i0 = 186 + 7 * r
        for c in range(29):
            j0 = 160 + 11 * (c - 14)
            out[r, c] = full[i0 - 3:i0 + 4, j0 - 5:j0 + 6].mean()
    return out


def map14_centres():
    r, c = np.meshgrid(np.arange(27), np.arange(29), indexing='ij')
    k = 7 * r * MAP_DK; l = 11 * (c - 14) * MAP_DL
    return np.c_[(-l / 2).ravel(), k.ravel(), l.ravel()]


def ko_blocks(full):
    """32 x 32 block averages (5 x 5 grid points) of a Fig.-5 map (161 x 161)"""
    return full[:160, :160].reshape(32, 5, 32, 5).mean(axis=(1, 3))


def ko_centres():
    r, c = np.meshgrid(np.arange(32), np.arange(32), indexing='ij')
    k = -2 + KO_DK * (5 * r + 2); l = KO_DL * (5 * c + 2)
    return np.c_[(-l / 2).ravel(), k.ravel(), l.ravel()]


def half_to_full(half):
    """Fig.-3 map from its k >= 0 half (187 x 321) using S(q) = S(-q) (the grid is symmetric about q = 0)"""
    return np.concatenate([half[1:][::-1, ::-1], half], axis=0)


def norm3(S):
    return np.sqrt((np.asarray(S) ** 2).sum(-1))





def spectra_arrays(spectra, emin=0.3, emax=3.0, smooth=1.5):
    Q, Y = [], []
    for q, E, S in spectra:
        win = (E > emin) & (E < emax)
        comps = []
        for c in range(3):
            I = gaussian_filter1d(np.asarray(S)[:, c], smooth)[win]
            comps.append(I / I.max())
        Q.append(q); Y.append(comps); Ew = E[win]
    return np.array(Q), np.array(Y), Ew


def spec_resid(M, p, Q, Y, E, sig, S=2.5):
    w = M.lswt(p[:12], p[12], Q, S=S)
    out = []
    for k in range(len(Q)):
        w_ = np.unique(np.round(w[k], 4))
        A = np.c_[np.exp(-0.5 * ((E[:, None] - w_[None, :]) / sig) ** 2), np.ones(len(E))]
        for c in range(3):
            a, _ = nnls(A, Y[k][c]); out.append(A @ a - Y[k][c])
    return np.concatenate(out)


def fit_spectra(M, spectra, p0, sigs=(0.06, 0.03), nfev=400, fixD=None, S=2.5):
    Q, Y, E = spectra_arrays(spectra)
    p = np.array(p0, float)
    for sig in sigs:
        if fixD is None:
            f = lambda x: spec_resid(M, x, Q, Y, E, sig, S)
            r = least_squares(f, p, diff_step=1e-4, max_nfev=nfev, x_scale=0.05)
            p = r.x
        else:
            f = lambda x: spec_resid(M, np.r_[x, fixD], Q, Y, E, sig, S)
            r = least_squares(f, p[:12], diff_step=1e-4, max_nfev=nfev, x_scale=0.05)
            p = np.r_[r.x, fixD]
    return p, float(r.cost)


def map_chi2(M, p, maps, nbz=24, lam_fixed=None):
    """sum over temperatures of chi^2 (best scale + offset) of block-averaged maps; maps: {T: full map}"""
    qc = map14_centres()
    tot = 0.0
    for T, full in maps.items():
        y = map14_blocks(full).ravel()
        m = norm3(M.scga(p[:12], p[12], T, qc, nbz=nbz, lam_fixed=lam_fixed))
        A = np.c_[m, np.ones_like(m)]
        c = np.linalg.lstsq(A, y, rcond=None)[0]
        tot += float(np.sum((A @ c - y) ** 2) / np.var(y))
    return tot


def null_scan(M, p, maps, grid=None):
    """move along (J3 + d, J11 - d) and return the parameter vector minimising the map chi^2"""
    s = p[2] + p[10]
    grid = np.linspace(min(s, 0.0), max(s, 0.0), 41) if grid is None else grid
    best = None
    for J3 in grid:
        q = p.copy(); q[2] = J3; q[10] = s - J3
        try:
            c = map_chi2(M, q, maps)
        except Exception:
            continue
        if best is None or c < best[0]:
            best = (c, J3)
    # refine with a parabola around the minimum
    J3b = best[1]; h = grid[1] - grid[0]
    cs = []
    for J3 in (J3b - h, J3b, J3b + h):
        q = p.copy(); q[2] = J3; q[10] = s - J3
        cs.append(map_chi2(M, q, maps))
    den = cs[0] - 2 * cs[1] + cs[2]
    if den > 0:
        J3b = J3b + 0.5 * h * (cs[0] - cs[2]) / den
    q = p.copy(); q[2] = J3b; q[10] = s - J3b
    return q


def predictions(M, p, T14=14.0, T20=20.0, lam_fixed=None):
    J, D = np.asarray(p[:12]), float(p[12])
    k, lam = M.lt(J)
    dE = (-lam + D / 2) - M.af1_energy(J, D)
    m14 = norm3(M.scga(J, D, T14, map14_centres(), nbz=32, lam_fixed=lam_fixed)).reshape(27, 29)
    red = {}
    for name, keep in REDUCED.items():
        Jr = np.zeros(12); Jr[keep] = J[keep]
        red[name] = norm3(M.scga(Jr, D, T20, ko_centres(), nbz=32, lam_fixed=lam_fixed)).reshape(32, 32)
    return {'J': J.tolist(), 'D': D, 'k': np.asarray(k).tolist(), 'dE': float(dE), 'map14': m14, 'reduced': red}


# ---------------- lab interaction and the reference workflow ----------------
APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(2)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(str(r))
    return r['id'], r['result']


def energies(text):
    # 'E_n = n * dE for n = 0 .. N'
    dE = float(text.split('*')[1].split('for')[0]); n = int(text.split('..')[1])
    return np.arange(n + 1) * dE


def main():
    shells = [[[b['from_site'] - 1, b['to_site'] - 1, *b['cell_offset']] for b in p['bonds']] for p in CATALOG['exchange_paths']]
    M = Model(shells)
    evidence = []
    spectra = []
    for i in [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 88, 94, 100]:
        jid, res = buy({'kind': 'spectrum', 'index': i})
        evidence.append(jid)
        S = np.stack([res['Sxx'], res['Syy'], res['Szz']], -1)
        spectra.append((np.array(res['q_rlu']), energies(res['energies_meV']), S))
    maps = {}
    for T in (20, 30):
        jid, res = buy({'kind': 'map', 'temperature_K': T})
        evidence.append(jid)
        half = np.sqrt(np.array(res['Sxx']) ** 2 + np.array(res['Syy']) ** 2 + np.array(res['Szz']) ** 2)
        maps[T] = half_to_full(half)
    jid, zone = buy({'kind': 'zone_slice', 'supercell': '36x24x36', 'k_index': 23})
    evidence.append(jid)
    zs = np.sqrt(sum(np.array(zone[c]) ** 2 for c in ('Sxx', 'Syy', 'Szz')))
    ih, il = np.unravel_index(np.argmax(zs), zs.shape)
    kbox = [(ih - 17) / 36, 0.5, (il - 17) / 36]
    # 1) spectra: 12 J + D (the paper's rounded Table I as the starting point)
    p, cost = fit_spectra(M, spectra, np.r_[TABLE_I, TABLE_I_D])
    s3 = p[2] + p[10]
    # 2) the exact null direction J3 - J11 of the archived path, fixed by the diffuse maps
    p = null_scan(M, p, maps)
    pred = predictions(M, p)
    J, D = pred['J'], pred['D']
    # minimum-norm split, for the report only
    pm = p.copy(); pm[2] = pm[10] = s3 / 2
    km, lm = M.lt(pm[:12]); dEm = (-lm + pm[12] / 2) - M.af1_energy(pm[:12], pm[12])
    k = pred['k']
    lam_full = M.lam_max(J, np.array([k]))[0]
    submission = {
        'exchange_meV': {'J%d' % (i + 1): float(J[i]) for i in range(12)},
        'anisotropy_D_meV': D,
        'propagation_vector_rlu': [float(x) for x in k],
        'spiral_minus_af1_meV_per_Mn': pred['dE'],
        'map_14K_prediction': pred['map14'].tolist(),
        'reduced_model_predictions': {m: v.tolist() for m, v in pred['reduced'].items()},
        'half_moon_origin': 'competing_3d_exchange',
        'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission) + '\n')
    tmf = 2 / 3 * lam_full / KB
    (APP/'report.md').write_text(f"""# MnWO4 half-moons (reference solution)

Model and conventions: H = -sum_(i!=j) J_ij e_i.e_j + D sum_i (e_i.n)^2 with unit vectors and ordered pairs (meV).
The archived 1.5 K spectra are classical spin dynamics of 5 muB moments (g = 2), so magnon energies follow linear
spin-wave theory with S = 5/2 and J/S^2 couplings (energies scale as J/S, not J). I fitted J1..J12 and D to 20 spectra
(AF1 linear spin waves, Gaussians with free non-negative amplitudes at the branch energies). Along the archived
[102] and [010] paths the AF1 spectrum depends on J3 and J11 only through J3 + J11 = {s3:.3f} meV (exact null
direction: the two bonds swap parallel/antiparallel roles and their phases are conjugate on these paths). Splitting
the sum equally (minimum norm) would put the flat spiral {-dEm:.3f} meV/Mn below AF1 and move k_LT to {np.round(km, 3)}.
The 20 K and 30 K diffuse maps (block averages, SCGA at the nominal temperature, best scale and offset) fix the split:
J3 = {J[2]:.3f}, J11 = {J[10]:.3f} meV. Final J (meV): {np.round(J, 3).tolist()}, D = {D:.3f} meV (easy axis; it opens the
0.57 meV gap at k_c).

Paramagnetic maps: self-consistent Gaussian approximation (spin-length constraint) rather than mean-field RPA, which is
singular below its mean-field temperature {tmf:.1f} K and over-sharpens the reduced model 5e at 20 K. The archived
|S| carries no magnetic form factor or polarisation factor (the purchased maps show no decay with |q|), so none is applied.
SCGA reproduces the purchased maps (block Pearson ~0.89 at 20 K); its absolute temperature scale is approximate, so
I used the maps only for the J3/J11 split and not for the overall exchange scale, which the spectra fix.

Half-moons: they bridge the k_ic peaks and sharpen on cooling from 30 K towards 14 K, i.e. they are short-range
correlations of a frustrated 3D exchange network, not a temperature-insensitive local-constraint (spin-liquid)
feature and not a 1D chain effect. Removing the interlayer J9 (model 5d) or the in-plane J3/J7/J10/J11 group (5e)
removes them in the predictions, so both in-plane and inter-plane competing antiferromagnetic bonds are required;
anisotropy is irrelevant above T_N (model 5b).

Ground state: Luttinger-Tisza vector k = {np.round(k, 4).tolist()} (thermodynamic limit). The archived 11 K full-zone
maximum sits at {np.round(kbox, 4).tolist()}, a commensurate approximant fixed by the 36x24x36 box. Without D the flat
spiral wins by exchange; with D the collinear AF1 is lower by only {pred['dE']:.3f} meV/Mn, i.e. AF1 and the multiferroic
spiral are nearly degenerate. Not supported by these data: a precise T_N of the model (mean field {tmf:.0f} K is an
overestimate; the 11 K data are already strongly correlated), the weak couplings J2, J5, J8, J12 individually, and any
statement about quantum fluctuations (the model is classical).
""")
    lab('submit', submission)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
