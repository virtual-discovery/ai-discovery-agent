"""Reference calculations for tmd-moire-chern-twist (run inside the runtime image, see run_in_image.sh).

    bash run_in_image.sh reference.py continuum   # -> reference/continuum.json
    bash run_in_image.sh reference.py hartree     # -> reference/hartree.json

continuum: K-valley spin-up Chern numbers of the two archived fig1 continuum parameter maps (Eq. 1 form) at 3.89 deg
and, for the set that reproduces the 3.89 deg DFT bands, at 2.14 deg with the same (twist-independent) parameters;
convergence variants and the wrong-convention variants used by the ablations.
hartree: Delta v_H = v_H(top surface) - v_H(bottom surface) of tWSe2 at XM (fractional (1/3, 1/3) of the archived
moire cell, located from the relaxed structures) and its first-star (six shortest G) contribution, with estimator
variants and wrong paths.
"""
import json
import sys
import time
from pathlib import Path
import numpy as np
from continuum import Continuum, coefficients_from_maps, moire, wu_basis

SRC = Path('/data/x/source-data')
OUT = Path(__file__).resolve().parent/'reference'
# set label -> (archive folder, m*/m_e, a0 in A, reference of the parameter set)
SETS = {'P': ('fig1b', 0.60, 3.52, 'DFT-fitted at 3.89 deg (paper ref. 25)'),
        'Q': ('fig1a', 0.62, 3.472, 'local-stacking approximation (paper ref. 20)')}


def load_maps(folder):
    rd = lambda name: np.loadtxt(SRC/'fig1'/folder/name)
    t, b, T = rd('Delta-top.dat'), rd('Delta-bot.dat'), rd('Delta-T.dat')
    return coefficients_from_maps(t[:, 0] + 1j * t[:, 1], b[:, 0] + 1j * b[:, 1], T[:, 0] + 1j * T[:, 1])


def chern(label, theta, nshell=5, nk=18, variant=None):
    folder, m, a0, _ = SETS[label]
    cb, ct, cT = load_maps(folder)
    basis = None
    if variant == 'conj_tunnelling':          # the plot script reads Delta-T as complex(re, -im)
        cT = {k: np.conj(v) for k, v in cT.items()}
    elif variant == 'layers_swapped':         # Delta-top.dat used as the first (kappa_+) block
        cb, ct = ct, cb
    elif variant == 'plot_script_basis':      # fig1 plot-sky.py lattice vectors (120 deg apart)
        aM = moire(theta, a0)[0]
        A = aM * np.array([[0.5, np.sqrt(3) / 2], [-0.5, np.sqrt(3) / 2]])
        basis = 2 * np.pi * np.linalg.inv(A).T
    elif variant == 'left_handed_basis':      # mirror image: b2 at -120 deg
        basis = wu_basis(theta, a0) * np.array([1, -1])
    t0 = time.time()
    M = Continuum(theta, a0, m, cb, ct, cT, nshell=nshell, basis=basis)
    C, E = M.chern(nk=nk, nb=3)
    gam = M.bands([0.0, 0.0], 4)[0]
    return {'set': label, 'theta': theta, 'nshell': nshell, 'plane_waves': M.n, 'nk': nk, 'variant': variant,
            'chern_raw': [round(float(c), 4) for c in C], 'chern': [int(round(c)) for c in C[:2]],
            'width_band1_meV': round(float(E[..., 0].max() - E[..., 0].min()), 4),
            'min_direct_gap_12_meV': round(float((E[..., 0] - E[..., 1]).min()), 4),
            'min_direct_gap_23_meV': round(float((E[..., 1] - E[..., 2]).min()), 4),
            'gamma_meV': [round(float(x), 3) for x in gam], 'seconds': round(time.time() - t0, 1)}


def path_bands(label, theta, nshell=5, n=20):
    folder, m, a0, _ = SETS[label]
    M = Continuum(theta, a0, m, *load_maps(folder), nshell=nshell)
    pts = [np.zeros(2), (M.kp + M.km) / 2, M.kp, np.zeros(2)]     # gamma - m - kappa - gamma
    ks = [pts[s] + (pts[s + 1] - pts[s]) * t for s in range(3) for t in np.linspace(0, 1, n + 1)[:(n if s < 2 else n + 1)]]
    return np.array([M.bands(k, 3)[0] for k in ks])


def continuum():
    rows = []
    for label in SETS:
        for ns, nk in ((4, 12), (5, 18), (7, 24)):
            rows.append(chern(label, 3.89, ns, nk)); print(rows[-1], flush=True)
        for v in ('conj_tunnelling', 'layers_swapped', 'plot_script_basis', 'left_handed_basis'):
            rows.append(chern(label, 3.89, 5, 18, v)); print(rows[-1], flush=True)
    for ns, nk in ((5, 18), (7, 24)):
        rows.append(chern('P', 2.14, ns, nk)); print(rows[-1], flush=True)
    # band comparison with the archived 3.89 deg SIESTA bands (gamma-m-kappa-gamma, 61 points, top doublets)
    raw = (SRC/'fig2/fig2d/3.89/Es-band.dat').read_text().strip().split('\n\n')
    E = np.array([[float(l.split()[1]) for l in b.strip().split('\n')] for b in raw])[::-1] * 1000
    dft = E[[0, 2]].T                             # one member of the two top doublets (meV), top = 0 at kappa
    cmp = {}
    for label in SETS:
        b = path_bands(label, 3.89)[:, :2]
        b = b - b[:, 0].max()
        cmp[label] = round(float(np.sqrt(np.mean((b - dft) ** 2))), 3)
    print('rms vs DFT 3.89', cmp)
    (OUT/'continuum.json').write_text(json.dumps({'sets': {k: v[3] for k, v in SETS.items()}, 'rows': rows,
                                                  'rms_vs_dft_3p89_meV': cmp}, indent=1) + '\n')


GRIDS = {'3.15': 320, '1.70': 300, '1.54': 320, '1.25': 360}
STAR1 = [(1, 0), (0, 1), (1, 1), (-1, 0), (0, -1), (-1, -1)]


def hartree_map(angle):
    n = GRIDS[angle]
    bot = 'pos_2d_bot_sm.dat' if angle == '3.15' else 'pot_2d_bot_sm.dat'
    top = np.loadtxt(SRC/'fig3/fig3a-d'/angle/'pot_2d_top_sm.dat').reshape(n, n)
    return (top - np.loadtxt(SRC/'fig3/fig3a-d'/angle/bot).reshape(n, n)) * 1000.0


def first_star_at(F, r):
    return float(sum((F[g] * np.exp(2j * np.pi * (g[0] * r[0] + g[1] * r[1]))).real for g in STAR1))


def hartree():
    out = {}
    for angle, n in GRIDS.items():
        d = hartree_map(angle)
        N = n - 1
        Fp = np.fft.fft2(d[:-1, :-1]) / N ** 2              # the grid repeats its first row/column at the end
        Ff = np.fft.fft2(d) / n ** 2                        # wrong-but-close: duplicated endpoint kept
        k = np.fft.fftfreq(N, 1 / N)
        xm, mx = (1 / 3, 1 / 3), (2 / 3, 2 / 3)
        fourier = lambda r: float((Fp * np.exp(2j * np.pi * (k[:, None] * r[0] + k[None, :] * r[1]))).sum().real)
        f = N / 3
        i0, t = int(np.floor(f)), f - np.floor(f)
        bil = float((1 - t) ** 2 * d[i0, i0] + t * (1 - t) * (d[i0 + 1, i0] + d[i0, i0 + 1]) + t * t * d[i0 + 1, i0 + 1])
        cut = np.loadtxt(SRC/'fig3/fig3e'/angle/'pot_diag-top.dat')
        mags = np.loadtxt(SRC/'fig3/fig3f'/angle/'pot_top_bot_G.dat')
        first_mag = float(np.mean(mags[(mags[:, 0] > 0.5) & (mags[:, 0] < 1.5), 1])) * 1000
        out[angle] = {
            'xm_fourier': fourier(xm), 'xm_nearest': float(d[int(round(f)), int(round(f))]), 'xm_bilinear': bil,
            'xm_linecut': float(cut[int(round((len(cut) - 1) / 3)), 1] * 1000),
            'mx_fourier': fourier(mx),
            'first_star_xm': first_star_at(Fp, xm), 'first_star_xm_fullgrid': first_star_at(Ff, xm),
            'first_star_mx': first_star_at(Fp, mx),
            'first_star_coefficients': {str(g): [float(Fp[g].real), float(Fp[g].imag)] for g in STAR1},
            'archived_first_star_magnitude_meV': first_mag,
            'first_star_xm_from_magnitudes_with_3p15_sign': -6 * first_mag * np.sin(2 * np.pi / 3),
            'min': float(d.min()), 'max': float(d.max())}
        print(angle, {k: (round(v, 4) if isinstance(v, float) else '') for k, v in out[angle].items()}, flush=True)
    (OUT/'hartree.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    OUT.mkdir(exist_ok=True)
    {'continuum': continuum, 'hartree': hartree}[sys.argv[1]]()
