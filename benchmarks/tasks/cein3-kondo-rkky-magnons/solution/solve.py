"""Reference solution for cein3-kondo-rkky-magnons; runs as the Harbor oracle and in check_task.

Buys zone-boundary 12-meV spectra, high-resolution spectra and cuts near R, fits linear spin-wave models (with the
archived momentum-integration windows) of increasing exchange range, and predicts the withheld R-M data.
"""
import json
import os
import re
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import fftconvolve

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
EVIDENCE = []
S_EFF = 0.5                      # effective spin of the Gamma7 doublet


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(1)
        r = lab(f"job/{r['id']}")
    if r['status'] != 'complete':
        raise RuntimeError(r.get('error'))
    EVIDENCE.append(r['id'])
    return r['result']


# ---------------- linear spin-wave forward model ----------------
R = np.array([0.5, 0.5, 0.5])
_N = np.array([(i, j, k) for i in range(-6, 7) for j in range(-6, 7) for k in range(-6, 7) if (i, j, k) != (0, 0, 0)])
SHELLS = []
for r2 in sorted(set((_N ** 2).sum(1))):
    vecs = _N[(_N ** 2).sum(1) == r2]
    groups = {}
    for v in vecs:
        groups.setdefault(tuple(sorted(np.abs(v))), []).append(v)
    SHELLS += [np.array(groups[k]) for k in sorted(groups)]
SHELLS = SHELLS[:10]
EF = np.arange(-1.0, 12.0, 0.01)


def design(P):
    return np.stack([np.cos(2 * np.pi * P @ v.T).sum(-1) for v in SHELLS], -1)


CR = design(R[None])[0]


def d2R(P):
    return ((P - np.floor(P) - 0.5) ** 2).sum(1)


def unit(v):
    v = np.array(v, float)
    return v / np.linalg.norm(v)


def perp(u):
    a = np.array([1., 0, 0]) if abs(u[0]) < 0.9 else np.array([0, 1., 0])
    e1 = unit(np.cross(u, a))
    return e1, unit(np.cross(u, e1))


class Bins:
    def __init__(self, centers, dirs, half_along, half_perp, n_perp=14, n_along=3):
        pts = []
        for c, u in zip(centers, dirs):
            e1, e2 = perp(u)
            a = ((np.arange(n_perp) + 0.5) / n_perp * 2 - 1) * half_perp
            t = ((np.arange(n_along) + 0.5) / n_along * 2 - 1) * half_along
            A_, B_, T_ = np.meshgrid(a, a, t, indexing='ij')
            pts.append(np.array(c) + T_.reshape(-1, 1) * u + A_.reshape(-1, 1) * e1 + B_.reshape(-1, 1) * e2)
        self.nb, self.npt = len(pts), len(pts[0])
        P = np.concatenate(pts)
        self.Ca, self.Cb, self.da, self.db = design(P) - CR, design(P + R) - CR, d2R(P), d2R(P + R)

    def spectra(self, p):
        J = np.r_[p['SJ'], np.zeros(len(SHELLS) - len(p['SJ']))]
        a, b = self.Ca @ J, self.Cb @ J
        if p['A']:
            k2 = 2 * p['kappa'] ** 2
            a = a - p['A'] * (np.exp(-self.da / k2) - 1)
            b = b - p['A'] * (np.exp(-self.db / k2) - 1)
        a, b = np.clip(a, 1e-9, None), np.clip(b, 0, None)
        w, A = np.sqrt(a * b), np.sqrt(b / a)
        de = EF[1] - EF[0]
        x = np.clip((w - EF[0]) / de, 0, len(EF) - 1.001)
        k = x.astype(int)
        f = x - k
        base = np.repeat(np.arange(self.nb), self.npt) * len(EF)
        n = self.nb * len(EF)
        H = np.bincount(base + k, A * (1 - f), n) + np.bincount(base + k + 1, A * f, n)
        return H.reshape(self.nb, len(EF)) / (self.npt * de)


def broaden(H, sigma, gamma):
    half = min(4.0, 6 * sigma + 15 * gamma + 0.05)
    x = np.arange(-half, half + 1e-9, EF[1] - EF[0])
    g = np.exp(-0.5 * (x / sigma) ** 2)
    g /= g.sum()
    if gamma > 1e-4:
        g = fftconvolve(gamma / np.pi / (x ** 2 + gamma ** 2) * (EF[1] - EF[0]), g, mode='same')
    return fftconvolve(H, g[None, :], mode='same', axes=1)


def ebins(Hc, Es, ebin):
    c = np.concatenate([np.zeros((Hc.shape[0], 1)), np.cumsum(Hc, 1)], 1)
    lo, hi = np.searchsorted(EF, Es - ebin / 2 - 1e-9), np.searchsorted(EF, Es + ebin / 2 - 1e-9)
    return (c[:, hi] - c[:, lo]) / np.maximum(hi - lo, 1)


HR_SIGMA = 0.106 / 2.3548
HR_HALF_PERP = 0.08
HR_DIR = {'RG': unit([-1, -1, -1]), 'RX': unit([0, -1, -1]), 'RM': unit([0, 0, -1])}


def path_dir(seg):
    return {'RG': unit([-1, -1, -1]), 'GX': unit([1, 0, 0]), 'XM': unit([0, 1, 0]), 'MG': unit([-1, -1, 0])}[seg]


class Data:
    """A purchased data set with model evaluation on its own binning."""
    def __init__(self, kind, y2, mask, e, model):
        self.kind, self.mask, self.model_fn = kind, mask, model
        self.y = y2[mask] if mask is not None else y2
        self.e = e

    def model(self, p):
        m = self.model_fn(p)
        return m[self.mask] if self.mask is not None else m


def project(m, y, e, offset=True):
    w = 1 / e ** 2
    if offset:
        A = np.c_[m, np.ones_like(m)] * np.sqrt(w)[:, None]
        c = np.linalg.lstsq(A, y * np.sqrt(w), rcond=None)[0]
        if c[0] < 0:
            c = np.array([0.0, np.sum(w * y) / np.sum(w)])
        return c[0] * m + c[1]
    return max(np.sum(w * m * y) / np.sum(w * m * m), 0) * m


NAMES = ['A', 'kappa', 'gamma_hr', 'sigma_he', 'gamma_he']
BND = {'A': (0, 30), 'kappa': (0.015, 0.12), 'gamma_hr': (0, 1.5), 'sigma_he': (0.05, 0.6), 'gamma_he': (0, 1.5)}


def fit(dsets, N, x0, fixed):
    names = [f'SJ{i+1}' for i in range(N)] + NAMES
    free = [n for n in names if n not in fixed]
    lo = np.array([BND.get(n, (-3, 3))[0] for n in free])
    hi = np.array([BND.get(n, (-3, 3))[1] for n in free])

    def to_p(x):
        p = dict(fixed)
        p.update(zip(free, x))
        p['SJ'] = [p.pop(f'SJ{i+1}') for i in range(N)]
        return p

    def res(x):
        p = to_p(x)
        return np.concatenate([(project(d.model(p), d.y, d.e) - d.y) / d.e / np.sqrt(len(d.y)) for d in dsets])
    xs = np.clip([x0.get(n, 0.0) for n in free], lo + 1e-6, hi - 1e-6)
    r = least_squares(res, xs, bounds=(lo, hi), x_scale='jac', max_nfev=300, diff_step=2e-2)
    return to_p(r.x), float(2 * r.cost)


def props(p):
    J = np.r_[p['SJ'], np.zeros(len(SHELLS) - len(p['SJ']))]

    def ab(P):
        a, b = (design(P) - CR) @ J, (design(P + R) - CR) @ J
        if p['A']:
            k2 = 2 * p['kappa'] ** 2
            a, b = a - p['A'] * (np.exp(-d2R(P) / k2) - 1), b - p['A'] * (np.exp(-d2R(P + R) / k2) - 1)
        return a, b
    vs = []
    for u in [unit([1, 1, 1]), unit([0, 1, 1]), unit([0, 0, 1])]:
        a, b = ab(np.array([R + 2e-3 * u]))
        vs.append(float(np.sqrt(max(a[0] * b[0], 0)) / 2e-3))
    f = np.linspace(0, 1, 24, endpoint=False)
    a, b = ab(np.array(np.meshgrid(f, f, f, indexing='ij')).reshape(3, -1).T)
    W = float(np.sqrt(np.clip(a * b, 0, None)).max())
    return {'v': float(np.mean(vs)), 'W': W, 'eta': float(np.mean(vs)) / (2 * np.pi * W)}


def grid(spec, n):
    """Evaluate a catalog grid expression such as 'q = -0.33 + 0.02 * index' for index n."""
    expr = spec.split('=', 1)[1].split(',')[0]
    var = [w for w in ('index', 'row', 'column', 'point') if w in expr][0]
    return float(eval(expr, {'__builtins__': {}, 'sqrt': np.sqrt, var: n}))


def count(spec):
    return int(re.search(r'0\.\.(\d+)', spec).group(1)) + 1


def nearest(avail, spec, targets):
    vals = {i: grid(spec, i) for i in avail}
    return sorted(set(min(vals, key=lambda i: abs(vals[i] - t)) for t in targets))


def main():
    req = CATALOG['requests']
    bands = buy({'kind': 'band_structure'})
    # 12-meV spectra: zone boundary and the middle of R-Gamma, where a nearest-neighbour model is tested on the bandwidth
    plan = {'RG': [0.116, 0.216, 0.316, 0.416, 0.516], 'GX': [0.966, 1.066, 1.166, 1.266, 1.366],
            'XM': [1.416, 1.516, 1.616, 1.716, 1.866], 'MG': [1.941, 2.091, 2.241, 2.391, 2.491]}
    av = req['he_spectrum']['available']
    he = {seg: [buy({'kind': 'he_spectrum', 'path': seg, 'index': i}) for i in nearest(av[seg]['indices'], av[seg]['s_of_index'], plan[seg])]
          for seg in plan}
    # high-resolution spectra near R (the long-wavelength magnons) and the constant-energy cuts
    qs = [-0.09, -0.05, -0.03, -0.01, 0.01, 0.03, 0.05, 0.07, 0.09, 0.13, 0.17, 0.25]
    av = req['hr_spectrum']['available']
    hr = {path: [buy({'kind': 'hr_spectrum', 'path': path, 'index': i}) for i in nearest(av[path]['indices'], av[path]['q_of_index'], qs)]
          for path in ['RG', 'RX']}
    cuts = [buy({'kind': 'hr_cut', 'path': path, 'energy_meV': e}) for path in ['RG', 'RX'] for e in [0.6, 1.4]]

    dsets_he, dsets = [], []
    for seg, rows in he.items():
        Es = np.array(rows[0]['energies_meV'])
        y2 = np.array([[np.nan if v is None else v for v in r['chi2']] for r in rows])
        mask = np.isfinite(y2) & (Es[None, :] > 0.3)
        B = Bins([r['q_center_rlu'] for r in rows], [path_dir(seg)] * len(rows), 0.0125, rows[0]['perpendicular_halfwidth_rlu'], 12)
        eb = rows[0]['energy_bin_meV']
        y = y2[mask]
        dsets_he.append(Data('he', y2, mask, np.sqrt(0.1 ** 2 + (0.15 * np.abs(y)) ** 2),
                             lambda p, B=B, Es=Es, eb=eb: ebins(broaden(B.spectra(p), p['sigma_he'], p['gamma_he']), Es, eb)))
    for path, rows in hr.items():
        Es = np.array(rows[0]['energies_meV'])
        q = np.array([r['q_rlu'] for r in rows])
        y2 = np.array([[np.nan if v is None else v for v in r['chi2']] for r in rows])
        mask = np.isfinite(y2) & (Es[None, :] > 0.25)
        pairs = [(i, int(np.argmin(abs(q + q[i])))) for i in range(len(q)) if q[i] > 0]
        diffs = np.array([y2[i] - y2[j] for i, j in pairs if abs(q[j] + q[i]) < 1e-6])[:, Es > 0.25]
        sig = float(np.nanstd(diffs) / np.sqrt(2))
        B = Bins([R + t * HR_DIR[path] for t in q], [HR_DIR[path]] * len(q), 0.01, HR_HALF_PERP)
        dsets.append(Data('hr', y2, mask, np.full(mask.sum(), sig),
                          lambda p, B=B, Es=Es: ebins(broaden(B.spectra(p), HR_SIGMA, p['gamma_hr']), Es, 0.1)))
    for c in cuts:
        q = np.array(c['q_rlu'])
        B = Bins([R + t * HR_DIR[c['path']] for t in q], [HR_DIR[c['path']]] * len(q), 0.005, HR_HALF_PERP)
        m = (EF >= c['energy_meV'] - 0.1 - 1e-9) & (EF < c['energy_meV'] + 0.1 - 1e-9)
        dsets.append(Data('cut', np.array(c['chi2']), None, np.array(c['error']),
                          lambda p, B=B, m=m: broaden(B.spectra(p), HR_SIGMA, p['gamma_hr'])[:, m].mean(1)))
    allsets = dsets_he + dsets

    # 1. nearest-neighbour local-moment model of the 12-meV data alone
    p_he, c_he = fit(dsets_he, 1, {'SJ1': 0.45, 'sigma_he': 0.4}, {'A': 0, 'kappa': 0.05, 'gamma_hr': 0.05, 'gamma_he': 0.0})
    J1 = p_he['SJ'][0] / S_EFF
    cost_j1_all = sum(float(np.sum(((project(d.model(p_he), d.y, d.e) - d.y) / d.e) ** 2) / len(d.y)) for d in dsets)
    # 2. short-range (J1-J3) model of all data, 3. short-range plus a long-range staggered term (sharp minimum of J(q) at R)
    fits = {'J1-J3': fit(allsets, 3, {'SJ1': 0.4, 'gamma_hr': 0.3, 'sigma_he': 0.4}, {'A': 0, 'kappa': 0.05}),
            'J1+long-range': fit(allsets, 1, {'SJ1': 0.3, 'A': 1.0, 'kappa': 0.06, 'gamma_hr': 0.3, 'sigma_he': 0.4}, {}),
            'J1-J3+long-range': fit(allsets, 3, {'SJ1': 0.2, 'A': 1.5, 'kappa': 0.06, 'gamma_hr': 0.3, 'sigma_he': 0.4}, {})}
    best = min(['J1+long-range', 'J1-J3+long-range'], key=lambda k: fits[k][1])
    p = fits[best][0]
    pr, pr_sr = props(p), props(fits['J1-J3'][0])

    # predictions of the withheld data on their own binning
    tp = CATALOG['to_predict']
    sq, se = tp['rm_map']['q_rlu_of_row'], tp['rm_map']['energy_meV_of_column']
    qrm = np.array([grid(sq, i) for i in range(count(sq))])
    erm = np.array([grid(se, j) for j in range(count(se))])
    Brm = Bins([R + t * HR_DIR['RM'] for t in qrm], [HR_DIR['RM']] * len(qrm), 0.01, HR_HALF_PERP)
    rm_map = ebins(broaden(Brm.spectra(p), HR_SIGMA, p['gamma_hr']), erm, erm[1] - erm[0])
    sc = tp['rm_cut_1p4meV']['q_rlu_of_point']
    qc = np.array([grid(sc, i) for i in range(count(sc))])
    Bc = Bins([R + t * HR_DIR['RM'] for t in qc], [HR_DIR['RM']] * len(qc), 0.005, HR_HALF_PERP)
    m = (EF >= 1.3 - 1e-9) & (EF < 1.5 - 1e-9)
    rm_cut = broaden(Bc.spectra(p), HR_SIGMA, p['gamma_hr'])[:, m].mean(1)

    def integrated(pp):
        out = []
        for pt in tp['integrated_chi']['points']:
            B = Bins([pt['q_center_rlu']], [path_dir(pt['segment'])], 0.0125, pt['perpendicular_halfwidth_rlu'], 16)
            out.append(float(B.spectra(pp)[0].sum() * (EF[1] - EF[0])))
        return out
    integ, integ_sr = integrated(p), integrated(fits['J1-J3'][0])
    submission = {'rm_map_prediction': np.round(rm_map, 6).tolist(), 'rm_cut_1p4meV_prediction': np.round(rm_cut, 6).tolist(),
                  'integrated_chi_prediction': [round(v, 6) for v in integ], 'bandwidth_meV': round(pr['W'], 4),
                  'velocity_bandwidth_ratio': round(pr['eta'], 4), 'J1_nearest_neighbour_meV': round(J1, 4),
                  'exchange_range': 'long_range', 'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission) + '\n')
    nf = int((np.abs(np.array(bands['f_bands_eV'])).min(0) < 0.02).sum())
    (APP/'report.md').write_text(f"""# CeIn3 magnons: nearest-neighbour versus long-range exchange (reference solution)

Model: linear spin waves of H = sum_<ij> J_ij S_i.S_j (each bond once) for the Gamma7 doublet (S = 1/2) on the simple-cubic
lattice ordered at R, omega(q) = S sqrt[(J(q)-J(R))(J(q+R)-J(R))] and chi'' ~ sqrt[(J(q+R)-J(R))/(J(q)-J(R))]. Every bin is
averaged over its archived momentum window (12-meV data: +-0.0125 r.l.u. along the path and +-0.09 or +-0.17 r.l.u.
perpendicular; high-resolution data: +-0.08 r.l.u. perpendicular), energy-binned and broadened (0.106 meV FWHM resolution
plus a Lorentzian HWHM fitted separately from the resolution). Each data set has its own scale and constant offset
(the archived background subtraction leaves negative intensities).

1. The 12-meV zone-boundary data alone are fitted well by a nearest-neighbour model: S*J1 = {p_he['SJ'][0]:.3f} meV,
   J1 = {J1:.3f} meV for S = 1/2, bandwidth {6 * p_he['SJ'][0]:.2f} meV. Its velocity at R is fixed by the bandwidth
   (eta = v/(2 pi W) = 1/sqrt3 = 0.577) and it fails the high-resolution data near R (mean reduced chi2 there {cost_j1_all / len(dsets):.2f}).
2. The broad constant-q profiles near R come mostly from integrating a steep magnon cone over the +-0.08 r.l.u. window, not
   from an incoherent continuum: with the window modelled, a dispersive model reproduces both the spectra and the cuts.
3. A short-range J1-J3 model of all data reaches eta = {pr_sr['eta']:.2f} (W = {pr_sr['W']:.2f} meV) and predicts a Q-integrated
   chi'' that varies strongly along the path (ratio max/min {max(integ_sr) / min(integ_sr):.1f}); short-range exchange cannot
   exceed eta ~ 1. Adding a long-range staggered term (a sharp minimum of J(q) at R of width kappa = {p['kappa']:.3f} r.l.u., real-space
   range ~1/(2 pi kappa) ~ {1 / (2 * np.pi * p['kappa']):.1f} a) fits better (fit cost {fits[best][1]:.2f} versus {fits['J1-J3'][1]:.2f}):
   v = {pr['v']:.1f} meV/r.l.u., W = {pr['W']:.2f} meV, eta = {pr['eta']:.2f}, damping HWHM {p['gamma_hr']:.2f} meV beyond resolution.
   Such a sharp minimum of J(q) is what RKKY exchange through the Fermi surface produces; in the archived band structure
   {nf} of the 14 f bands lie within 20 meV of the chemical potential (heavy, hybridised Gamma7 states), so the spectrum
   needs long-range Kondo-lattice exchange rather than nearest-neighbour superexchange between localised moments.
4. Predictions of the withheld R-M map, the R-M 1.4-meV cut and the Q-integrated chi'' use the {best} model with the same binning.

Not supported by these data: a spin gap below ~0.2 meV (the lowest bins are dominated by the subtracted elastic line), a
precise velocity (it depends on the reading of the 0.08 r.l.u. integration window and on the shape of the long-range
term; values of 20-35 meV/r.l.u. are compatible), the absolute high-energy damping (the 12-meV resolution is not archived),
and the sign pattern of individual distant exchange shells.
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
