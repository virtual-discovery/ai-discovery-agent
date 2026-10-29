"""Reference solution for tmd-moire-chern-twist; runs as the Harbor oracle and in check_task."""
import cmath
import json
import math
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
AR = APP/'archive'
CATALOG = json.loads((APP/'catalog.json').read_text())
HB2M = 3809.982  # hbar^2 / 2 m_e (meV A^2)


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    while result.get('status') == 'running':
        time.sleep(2)
        result = lab(f"job/{result['id']}")
    return result


# ---------------------------------------------------------------- continuum model (Eq. 1)
def fourier(path, n):
    d = np.loadtxt(path)
    z = (d[:, 0] + 1j * d[:, 1]).reshape(n, n)[:-1, :-1]       # grids repeat their first row/column
    F = np.fft.fft2(z) / (n - 1) ** 2
    out = {}
    for i, j in zip(*np.where(abs(F) > 1e-3)):
        out[(int(i) if i <= (n - 1) // 2 else int(i) - (n - 1), int(j) if j <= (n - 1) // 2 else int(j) - (n - 1))] = F[i, j]
    return out


class Continuum:
    """K-valley spin-up plane-wave model. The map harmonics (1,0),(0,1),(1,1) form the first star, so the grid basis
    is 60 deg apart; Delta_T = w(1 + e^{-i g2.r} + e^{-i g3.r}) with right-handed axes fixes b1 = g1, b2 = g3."""

    def __init__(self, theta, a0, mstar, cb, ct, cT, nshell=5):
        aM = a0 / (2 * math.sin(math.radians(theta) / 2))
        g = 4 * math.pi / (math.sqrt(3) * aM)
        kk = 4 * math.pi / (3 * aM)
        self.kp, self.km = kk * np.array([-math.sqrt(3) / 2, -0.5]), kk * np.array([-math.sqrt(3) / 2, 0.5])
        self.b = np.array([[g, 0.0], [g * math.cos(2 * math.pi / 3), g * math.sin(2 * math.pi / 3)]])
        self.t = HB2M / mstar
        G = [(i, j) for i in range(-nshell, nshell + 1) for j in range(-nshell, nshell + 1) if abs(i - j) <= nshell]
        idx = {x: a for a, x in enumerate(G)}
        self.G, self.n = np.array(G), len(G)
        self.Gc = self.G @ self.b
        self.M = []
        for coef in (cb, ct, cT):
            M = np.zeros((self.n, self.n), complex)
            for a, ga in enumerate(G):
                for dg, c in coef.items():
                    j = idx.get((ga[0] - dg[0], ga[1] - dg[1]))
                    if j is not None:
                        M[a, j] += c
            self.M.append(M)
        self.shift = [np.array([idx.get((ga[0] + d[0], ga[1] + d[1]), -1) for ga in G]) for d in ((1, 0), (0, 1))]

    def bands(self, k, nb=3):
        n = self.n
        kb, kt = k + self.Gc - self.kp, k + self.Gc - self.km
        H = np.zeros((2 * n, 2 * n), complex)
        H[:n, :n] = self.M[0] - np.diag(self.t * (kb ** 2).sum(1))
        H[n:, n:] = self.M[1] - np.diag(self.t * (kt ** 2).sum(1))
        H[:n, n:] = self.M[2]
        H[n:, :n] = self.M[2].conj().T
        e, v = np.linalg.eigh(H)
        return e[::-1][:nb], v[:, ::-1][:, :nb]

    def shifted(self, u, d):
        out = np.zeros_like(u)
        p = self.shift[d]
        ok = p >= 0
        for L in range(2):
            out[L * self.n:(L + 1) * self.n][ok] = u[L * self.n:(L + 1) * self.n][p[ok]]
        return out

    def chern(self, nk=18, nb=2):
        U = np.empty((nk + 1, nk + 1, 2 * self.n, nb), complex)
        for i in range(nk):
            for j in range(nk):
                U[i, j] = self.bands(i / nk * self.b[0] + j / nk * self.b[1], nb)[1]
        for j in range(nk):
            U[nk, j] = self.shifted(U[0, j], 0)
        for i in range(nk + 1):
            U[i, nk] = self.shifted(U[i, 0], 1)
        link = lambda a, b: np.einsum('ija,ija->ij', a.conj(), b)
        out = []
        for m in range(nb):
            u = U[..., m]
            F = np.angle(link(u[:-1, :-1], u[1:, :-1]) * link(u[1:, :-1], u[1:, 1:]) *
                         link(u[1:, 1:], u[:-1, 1:]) * link(u[:-1, 1:], u[:-1, :-1]))
            out.append(int(round(-F.sum() / (2 * math.pi))))   # arg<u|u+du> = -A.dk
        return out


def model(s, theta):
    p = CATALOG['continuum_parameter_sets'][s]
    d = AR/f'fig1/set_{s}'
    return Continuum(theta, p['monolayer_lattice_constant_A'], p['effective_mass_me'], fourier(d/'Delta-bot.dat', 100),
                     fourier(d/'Delta-top.dat', 100), fourier(d/'Delta-T.dat', 10))


def es_band(path):
    blocks = path.read_text().strip().split('\n\n')
    return np.array([[float(l.split()[1]) for l in b.strip().split('\n')] for b in blocks])[::-1] * 1000


# ---------------------------------------------------------------- symmetry indicator
def chern_from_c3(states):
    """exp(2 pi i C / 3) = -xi_gamma xi_kappa xi_kappa' (spinful C3z), K-valley spin-up states only."""
    out = []
    for n in range(2):
        z = -1
        for p in ('gamma', 'kappa', 'kappa_prime'):
            up = [s for s in states[p] if s['spin'] == 'up'][n]
            z *= complex(*up['C3z_eigenvalue'])
        out.append(int(round(cmath.phase(z) / (2 * math.pi / 3))))
    return out


# ---------------------------------------------------------------- Hartree potential
def xm_frac_from_structure(poscar_text):
    """Return which of (1/3,1/3), (2/3,2/3) is XM (top-layer chalcogen above bottom-layer metal)."""
    L = poscar_text.split('\n')
    A = np.array([[float(x) for x in L[i].split()] for i in (2, 3, 4)])
    sp, cnt = L[5].split(), [int(x) for x in L[6].split()]
    P = np.array([[float(x) for x in L[8 + i].split()[:3]] for i in range(sum(cnt))])
    lab_ = np.array(sum([[s] * c for s, c in zip(sp, cnt)], []))
    metal = [s for s in sp if s not in ('Se', 'Te', 'S')][0]
    z = P[:, 2] * A[2, 2]
    zm = np.median(z[lab_ == metal])
    zt, zb = z[(lab_ == metal) & (z > zm)].mean(), z[(lab_ == metal) & (z < zm)].mean()
    top = abs(z - zt) < abs(z - zb)

    def below(r, sp_top):
        m = (lab_ == sp_top) & top
        d = P[m, :2] - r
        d -= np.round(d)
        i = np.argmin(np.linalg.norm(d @ A[:2, :2], axis=1))
        rs = P[m][i, :2]
        res = {}
        for s in (metal, 'X'):
            mb = (~top) & ((lab_ == metal) if s == metal else (lab_ != metal))
            d = P[mb, :2] - rs
            d -= np.round(d)
            res[s] = np.linalg.norm(d @ A[:2, :2], axis=1).min()
        return res
    chalc = [s for s in sp if s != metal][0]
    for f in ((1 / 3, 1 / 3), (2 / 3, 2 / 3)):
        if below(np.array(f), chalc)[metal] < 0.5:
            return f
    raise RuntimeError('XM not found')


def periodic_fft(d):
    N = d.shape[0] - 1
    return np.fft.fft2(d[:-1, :-1]) / N ** 2


def value_at(F, r, star=None):
    N = F.shape[0]
    k = np.fft.fftfreq(N, 1 / N)
    ph = np.exp(2j * np.pi * (k[:, None] * r[0] + k[None, :] * r[1]))
    if star is not None:
        mask = np.zeros_like(F, bool)
        for g in star:
            mask[g] = True
        return float((F * ph)[mask].sum().real)
    return float((F * ph).sum().real)


STAR1 = [(1, 0), (0, 1), (1, 1), (-1, 0), (0, -1), (-1, -1)]


def main():
    evidence = []
    # 1. continuum sets at 3.89 deg, and which reproduces the staged DFT bands
    dft389 = es_band(AR/'fig2/tMoTe2_3.89/Es-band.dat')[[0, 2]].T
    cont, rms = {}, {}
    for s in ('P', 'Q'):
        M = model(s, 3.89)
        cont[s] = M.chern()
        pts = [np.zeros(2), (M.kp + M.km) / 2, M.kp, np.zeros(2)]
        ks = [pts[i] + (pts[i + 1] - pts[i]) * t for i in range(3) for t in np.linspace(0, 1, 21)[:(20 if i < 2 else 21)]]
        b = np.array([M.bands(k, 2)[0] for k in ks])
        b -= b[:, 0].max()
        rms[s] = float(np.sqrt(np.mean((b - dft389) ** 2)))
    match = min(rms, key=rms.get)
    fixed_214 = model(match, 2.14).chern()
    # 2. C3z eigenvalues of the DFT doublets (gamma and kappa everywhere; kappa' checked once)
    c3 = {}
    for mat, info in CATALOG['systems'].items():
        for angle in info['angles_deg']:
            st = {}
            for p in ('gamma', 'kappa'):
                r = lab('evaluate', {'kind': 'c3', 'material': mat, 'angle': angle, 'point': p})
                evidence.append(r['id'])
                st[p] = r['result']['states']
            c3[(mat, angle)] = st
    r = lab('evaluate', {'kind': 'c3', 'material': 'tMoTe2', 'angle': '3.89', 'point': 'kappa_prime'})
    evidence.append(r['id'])
    same = [s['C3z_eigenvalue'] for s in r['result']['states'] if s['spin'] == 'up'] == \
           [s['C3z_eigenvalue'] for s in c3[('tMoTe2', '3.89')]['kappa'] if s['spin'] == 'up']
    dft = {}
    for (mat, angle), st in c3.items():
        st['kappa_prime'] = st['kappa'] if same else None
        dft[f'{mat}_{angle}'] = chern_from_c3(st)
    # 3. Hartree potential drop at XM; XM located from the relaxed structures
    xm_315 = xm_frac_from_structure((AR/'fig2/tWSe2_relaxation_3.15/CONTCAR').read_text())
    xm_125 = xm_frac_from_structure((AR/'fig2/tWSe2_relaxation_1.25/CONTCAR').read_text())
    assert xm_315 == xm_125
    xm = xm_315
    n = 320
    d315 = (np.loadtxt(AR/'fig3/tWSe2_3.15/pot_2d_top_sm.dat') - np.loadtxt(AR/'fig3/tWSe2_3.15/pos_2d_bot_sm.dat')).reshape(n, n) * 1000
    hx = {'3.15': value_at(periodic_fft(d315), xm)}
    first = {'3.15': value_at(periodic_fft(d315), xm, STAR1)}
    for angle in ('1.70', '1.25'):
        r = lab('evaluate', {'kind': 'hartree_map', 'angle': angle})
        evidence.append(r['id'])
        d = (np.array(r['result']['v_top_surface_eV']) - np.array(r['result']['v_bottom_surface_eV'])) * 1000
        F = periodic_fft(d)
        hx[angle], first[angle] = value_at(F, xm), value_at(F, xm, STAR1)
    r = lab('evaluate', {'kind': 'hartree_linecut', 'angle': '1.54'})
    evidence.append(r['id'])
    cut = np.array(r['result']['data'])
    s_xm = np.linalg.norm(np.array(xm) @ np.array([[1.0, 0.0], [0.5, math.sqrt(3) / 2]]))   # |xm| in units of |a1|
    hx['1.54'] = float(cut[np.argmin(abs(cut[:, 0] - s_xm)), 1] * 1000)
    submission = {'dft_chern_K_up': dft, 'continuum_chern_3p89': cont, 'continuum_set_matching_dft': match,
                  'continuum_P_chern_2p14': fixed_214 if match == 'P' else model('P', 2.14).chern(),
                  'hartree_xm_meV': {a: round(hx[a], 4) for a in ('3.15', '1.70', '1.54', '1.25')},
                  'hartree_first_star_xm_meV': {a: round(first[a], 4) for a in ('1.70', '1.25')},
                  'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# Twist-angle evolution of the K-valley spin-up moire Chern bands (reference solution)

**Symmetry indicators.** For each DFT calculation the C3z eigenvalues of the top doublets were bought at gamma and
kappa (kappa' checked at 3.89 deg tMoTe2: identical to kappa, as required by C2y T, not complex conjugate). Only the
spin-up (K-valley) member of each doublet was used; the spin-down member belongs to K' and carries the opposite Chern
number. With spinful C3z (eigenvalues -1, e^(+-i pi/3)) the indicator is exp(2 pi i C/3) = -xi_gamma xi_kappa xi_kappa',
giving C mod 3; the smallest-magnitude representative is reported. Chern pairs (C1, C2): {dft}.
tWSe2: (+1,+1) at 3.15 and 1.70 deg, (0,0) at 1.54 deg (bands 1-2 invert at gamma, 2-3 at kappa), (-1,-1) at 1.25 deg.
tMoTe2: (+1,-1) at 3.89/3.15, (+1,+1) at 2.14 (band-2/3 inversion at gamma), (0,0) at 1.25 deg.

**Continuum model.** The fig1 maps have first-star harmonics (1,0),(0,1),(1,1) (60-deg grid basis, not the 120-deg
vectors of plot-sky.py); with Delta_T = w(1 + e^(-i g2 r) + e^(-i g3 r)) and right-handed axes the basis is b1 = g1,
b2 = g3. Fukui-Hatsugai-Suzuki Chern numbers (18x18 mesh, 91 plane waves per layer; converged): {cont}.
RMS distance of the two top bands to the 3.89 deg SIESTA bands: {rms} meV, so set {match} is the DFT-fitted set; it
reproduces the DFT Chern pair (+1,-1). The other set (local-stacking) gives the opposite top-band Chern number:
its Delta_z poles are reversed (opposite skyrmion number). Keeping set {match} fixed and changing only the twist to
2.14 deg gives {fixed_214}, while DFT gives (+1,+1): a twist-independent continuum model does not capture the
evolution.

**Hartree potential.** XM (top-layer Se above bottom-layer W) is at fractional {tuple(round(x, 3) for x in xm)} of the
archived cells (checked on the 3.15 and 1.25 deg structures). Delta v_H(XM) = {hx} meV: negative at 3.15 deg and
positive at 1.70/1.54/1.25 deg - the sign flip of the layer potential (north/south pole exchange). The first-star
(|G| = g1) part of Delta v_H at XM is {first} meV: at 1.70 deg it nearly vanishes although the full drop is ~+9 meV,
so the first-star continuum expansion fails there (higher harmonics, six-petal MM structure); at 1.25 deg it has the
reversed sign. The archived Fourier file holds magnitudes only and cannot show the sign change.

**Interpretation.** The flip of the layer potential tracks the competition between the out-of-plane moire
ferroelectric polarization (stacking dependent, XM vs MX) and the relaxation-induced in-plane piezoelectric charge
(strain gradients concentrated at domain walls at small angles). It causes band inversions at gamma and kappa that
change the Chern numbers; a fixed local-stacking continuum model cannot describe it.

**Limits.** C3 indicators fix C only mod 3; the DFT bands come in near-degenerate K-up/K'-down doublets and the
spin-up assignment relies on the authors' labels; tMoTe2 has no archived Hartree data, so its mechanism is inferred
by analogy; continuum Chern numbers near 1.25 deg (band 1-2 gap < 0.01 meV in the fixed model) are not reliable.
""")
    lab('submit', submission)
    for _ in range(120):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
