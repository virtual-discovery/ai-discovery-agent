"""Analysis library for raman-polymorph-orientation (reference solution and construction variants).

Self-contained (numpy/scipy only); make_task.py writes solution/solve.py as this file plus a driver, so the
oracle and the construction-time ablations share one implementation. Works on the catalog's record format.
"""
import itertools
import math
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation

VASP2THZ, THZ2CM, E2, HC_K = 15.633302, 33.35641, 14.399645, 1.4387769
OFFSET = 500          # detector offset per accumulation (instrument description)
RES = 1.2             # instrumental Gaussian FWHM (cm-1)


# ----------------------------------------------------------------------------- database records
def parse(rec):
    e = {'mpid': rec['mpid'], 'cell': np.array(rec['cell_A'], float), 'numbers': np.array(rec['atomic_numbers']),
         'positions': np.array(rec['positions_A'], float), 'masses': np.array(rec['masses_amu'], float),
         'freq': np.array(rec['frequencies_cm1'], float), 'eigvec': np.array(rec['eigenvectors'], float),
         'raman_active': np.array(rec['raman_active']), 'ir_active': np.array(rec['ir_active']),
         'bands': rec['degenerate_sets'], 'labels': rec['irrep_labels']}
    n = len(e['freq'])
    e['E'] = e['eigvec'].reshape(n, -1)
    e['Rsym'] = np.array([np.zeros((3, 3)) if t is None else 0.5 * (np.array(t) + np.array(t).T)
                          for t in rec['raman_tensors']])
    e['born_rows'] = None if rec['born_charges'] is None else np.array(rec['born_charges'], float)
    e['eps'] = None if rec['epsilon_inf'] is None else np.array(rec['epsilon_inf'], float)
    e['volume'] = abs(np.linalg.det(e['cell']))
    e['frac'] = e['positions'] @ np.linalg.inv(e['cell'])
    e['Z'] = expand_born(e) if e['born_rows'] is not None else None
    return e


def space_group_ops(e, tol=1e-3):
    A, frac, num = e['cell'], e['frac'], e['numbers']
    G = A @ A.T
    ops = []
    for flat in itertools.product((-1, 0, 1), repeat=9):
        W = np.array(flat).reshape(3, 3)
        if abs(round(np.linalg.det(W))) != 1 or not np.allclose(W.T @ G @ W, G, atol=1e-3 * np.abs(G).max()):
            continue
        Wx = frac @ W.T
        for j in np.where(num == num[0])[0]:
            y = Wx + frac[j] - Wx[0]
            perm = []
            for k in range(len(frac)):
                dd = y[k] - frac
                dd -= np.round(dd)
                m = np.where((np.linalg.norm(dd @ A, axis=1) < 10 * tol) & (num == num[k]))[0]
                if len(m) != 1:
                    break
                perm.append(int(m[0]))
            if len(perm) == len(frac):
                ops.append((W, np.array(perm)))
    return ops


def expand_born(e, mode='rotate'):
    """The archived rows belong to the symmetry-independent atoms (first appearance in each orbit)."""
    ops = space_group_ops(e)
    A = e['cell']
    n = len(e['numbers'])
    rep = -np.ones(n, int)
    for k in range(n):
        if rep[k] < 0:
            for W, p in ops:
                if rep[p[k]] < 0:
                    rep[p[k]] = k
    reps = sorted(set(rep.tolist()))
    rows = e['born_rows']
    Z = np.zeros((n, 3, 3))
    for k in range(n):
        Zi = rows[reps.index(rep[k])]
        if mode == 'unrotated':
            Z[k] = Zi
            continue
        for W, p in ops:
            if p[rep[k]] == k:
                R = A.T @ W @ np.linalg.inv(A.T)
                Z[k] = R @ Zi @ R.T
                break
    if mode == 'isotropic':
        Z = np.array([np.trace(z) / 3 * np.eye(3) for z in Z])
    return Z


def modes(e, qhat=None, Z=None, eps=None):
    """Optical modes (frequency in cm-1 on the database scale, Raman tensor) for phonon direction qhat
    in the database Cartesian frame; qhat=None returns the archived transverse (analytic) modes."""
    f0, E = e['freq'], e['E']
    Z = e['Z'] if Z is None else Z
    if qhat is None or Z is None:
        idx = np.argsort(f0)[3:]
        return f0[idx], e['Rsym'][idx]
    eps = e['eps'] if eps is None else eps
    w = f0 / THZ2CM / VASP2THZ
    D = E.T @ np.diag(np.sign(w) * w ** 2) @ E
    q = np.asarray(qhat, float) / np.linalg.norm(qhat)
    zq = np.einsum('a,kab->kb', q, Z).reshape(-1) / np.sqrt(np.repeat(e['masses'], 3))
    D = D + 4 * np.pi * E2 / e['volume'] * np.outer(zq, zq) / (q @ eps @ q)
    val, vec = np.linalg.eigh(D)
    f = np.sign(val) * np.sqrt(np.abs(val)) * VASP2THZ * THZ2CM
    return f[3:], np.einsum('jm,jab->mab', E @ vec, e['Rsym'])[3:]


def strongest_band(e):
    vals = [sum(float((e['Rsym'][k] ** 2).sum()) for k in b) for b in e['bands']]
    return e['bands'][int(np.argmax(vals))]


def powder(R):
    a = np.trace(R) / 3
    g2 = 0.5 * ((R[0, 0] - R[1, 1]) ** 2 + (R[0, 0] - R[2, 2]) ** 2 + (R[1, 1] - R[2, 2]) ** 2
                + 6 * (R[0, 1] ** 2 + R[0, 2] ** 2 + R[1, 2] ** 2))
    return (45 * a * a + 4 * g2) / 45, g2 / 15


# ----------------------------------------------------------------------------- spectra
def pol(a):
    a = math.radians(a)
    return np.array([math.cos(a), math.sin(a), 0.0])


def pvoigt(x, wl, wg=RES):
    f = (wg ** 5 + 2.69269 * wg ** 4 * wl + 2.42843 * wg ** 3 * wl ** 2 + 4.47163 * wg ** 2 * wl ** 3
         + 0.07842 * wg * wl ** 4 + wl ** 5) ** 0.2
    r = wl / f
    eta = 1.36603 * r - 0.47719 * r ** 2 + 0.11116 * r ** 3
    g, s = f / 2, f / 2.354820045
    return eta * (g / np.pi) / (x * x + g * g) + (1 - eta) * np.exp(-x * x / (2 * s * s)) / (s * np.sqrt(2 * np.pi))


def thermal(f, laser, T=295.0, prefactors=True):
    """Stokes and anti-Stokes factors omega_s^3 (n+1)/omega and omega_as^3 n/omega (photon counting)."""
    if not prefactors:
        return np.ones_like(f), np.zeros_like(f)
    nb = 1.0 / np.expm1(HC_K * f / T)
    return ((laser - f) / laser) ** 3 * (nb + 1) / f, ((laser + f) / laser) ** 3 * nb / f


def prepare(results):
    """Purchase results -> spectra dicts (counts above the detector offset)."""
    out = []
    for r in results:
        acc = r['accumulations']
        out.append({'ai': r['polarizer_deg'], 'as': r['analyzer_deg'], 'center': r['center_cm1'], 'accum': acc,
                    'nu': np.array(r['raman_shift_cm1']), 'y': np.array(r['counts'], float) - OFFSET * acc})
    return out


def config_strengths(Rlab, sp):
    ei = pol(sp['ai'])
    outs = [pol(0.0), pol(90.0)] if sp['as'] is None else [pol(sp['as'])]
    return sum(np.einsum('i,mij,j->m', es, Rlab, ei) ** 2 for es in outs)


class CandidateFit:
    """Global fit of all purchased spectra with one candidate polymorph.
    Nonlinear parameters: rotation vector (3), ln C, s, delta, g0, g1 (line FWHM g0 + g1*omega).
    Per spectrum, the amplitude a_k (= C times the laser-power drift, prior 1.5 % rms) and a linear
    background are solved linearly inside the residual."""

    def __init__(self, e, spectra, nac=True, prefactors=True, fix_offset=None, laser=1e7 / 532.0, drift_sd=0.015):
        self.e, self.sp, self.nac, self.pref, self.laser = e, spectra, nac, prefactors, laser
        self.fix_offset, self.drift_sd = fix_offset, drift_sd
        self.w = [1.0 / np.sqrt(np.maximum(sp['y'], 1.0) + 25.0 * sp['accum']) for sp in spectra]

    def lines(self, Q):
        if self.nac and self.e['Z'] is not None:
            f, R = modes(self.e, Q.T @ np.array([0.0, 0.0, 1.0]))
        else:
            f, R = modes(self.e)
        keep = (f > 5) & ((R ** 2).sum((1, 2)) > 1e-10)
        return f[keep], np.einsum('ia,mab,jb->mij', Q, R[keep], Q) / math.sqrt(self.e['volume'])

    def unpack(self, p):
        d = self.fix_offset if self.fix_offset is not None else p[5]
        return Rotation.from_rotvec(p[:3]).as_matrix(), p[3], p[4], d, p[6], p[7]

    def signals(self, p):
        Q, lnC, s, d, g0, g1 = self.unpack(p)
        f, Rl = self.lines(Q)
        f = f * s
        w = g0 + g1 * f
        st, an = thermal(f, self.laser - d, prefactors=self.pref)
        out = []
        for sp in self.sp:
            x = (sp['nu'] - d)[:, None]
            prof = st * pvoigt(x - f, w) + an * pvoigt(x + f, w)
            out.append(prof @ config_strengths(Rl, sp))
        return out

    def solve_linear(self, p, sig):
        C = math.exp(p[3])
        res, coefs = [], []
        for sp, w, sg in zip(self.sp, self.w, sig):
            acc = sp['accum']
            A = np.array([acc * sg, acc * np.ones_like(sg), acc * (sp['nu'] - sp['center']) / 100]).T * w[:, None]
            b = sp['y'] * w
            lam = 1.0 / (self.drift_sd * C)
            A2 = np.vstack([A, [lam, 0.0, 0.0]])
            b2 = np.r_[b, lam * C]
            c = np.linalg.lstsq(A2, b2, rcond=None)[0]
            res.append(A2 @ c - b2)
            coefs.append(c)
        return np.concatenate(res), coefs

    def resid(self, p):
        return self.solve_linear(p, self.signals(p))[0]

    def model(self, p):
        sig = self.signals(p)
        _, coefs = self.solve_linear(p, sig)
        return [sp['accum'] * (c[0] * sg + c[1] + c[2] * (sp['nu'] - sp['center']) / 100)
                for sp, sg, c in zip(self.sp, sig, coefs)], coefs

    def start(self, rotvec, s, d):
        p = np.array([*rotvec, 0.0, s, d, 1.0, 0.006])
        sig = self.signals(p)
        num = sum(float((sg * (sp['y'] - np.percentile(sp['y'], 10))).sum()) / sp['accum'] for sg, sp in zip(sig, self.sp))
        den = sum(float((sg * sg).sum()) for sg in sig)
        p[3] = math.log(max(num, 1e-30) / max(den, 1e-30))
        return p

    def fit(self, p0, max_nfev=60):
        lo = [-np.inf] * 4 + [0.93, -10.0, 0.2, 0.0]
        hi = [np.inf] * 4 + [1.07, 10.0, 6.0, 0.03]
        p0 = np.clip(np.array(p0, float), np.array(lo) + 1e-9, np.array(hi) - 1e-9)
        o = least_squares(self.resid, p0, bounds=(lo, hi), x_scale='jac', max_nfev=max_nfev)
        o.chi2 = float(2 * o.cost)
        o.npts = len(o.fun)
        o.Q = Rotation.from_rotvec(o.x[:3]).as_matrix()
        o.offset = self.fix_offset if self.fix_offset is not None else float(o.x[5])
        return o


# ----------------------------------------------------------------------------- calibration and search
def peak_position(sp, guess, half=6.0):
    """Centre of the strongest peak near `guess` (pseudo-Voigt + linear background)."""
    m = np.abs(sp['nu'] - guess) < half
    x, y = sp['nu'][m], sp['y'][m]
    x0 = x[np.argmax(y)]

    def res(p):
        return (p[0] * pvoigt(x - p[1], p[2]) + p[3] + p[4] * (x - x0) - y) / np.sqrt(np.maximum(y, 1) + 25)
    p0 = [max(y.max() - y.min(), 1) * 3, x0, 3.0, y.min(), 0.0]
    o = least_squares(res, p0, bounds=([0, x0 - 3, 0.3, -np.inf, -np.inf], [np.inf, x0 + 3, 15, np.inf, np.inf]))
    return o.x[1]


def calibrate(spectra):
    """Strongest Stokes line above 200 cm-1 and its anti-Stokes partner -> axis offset and true shift."""
    stokes = [sp for sp in spectra if sp['nu'].min() < 250 and sp['nu'].max() > 330]
    anti = [sp for sp in spectra if sp['nu'].min() < -330 and sp['nu'].max() > -250]
    sp = max(stokes, key=lambda s: s['y'][(s['nu'] > 200)].max() / s['accum'])
    sel = sp['nu'] > 200
    guess = sp['nu'][sel][np.argmax(sp['y'][sel])]
    nu_s = peak_position(sp, guess)
    vals = []
    for a in anti:
        sel = (a['nu'] < -200)
        if sel.sum() == 0:
            continue
        vals.append(peak_position(a, -nu_s, half=9.0))
    nu_as = float(np.median(vals))
    return 0.5 * (nu_s + nu_as), 0.5 * (nu_s - nu_as)


def random_rotations(n, seed=12345):
    return Rotation.random(n, random_state=seed)


def coarse_orientation(e, spectra, s, d, n_rot=6000, keep=6, prefactors=True):
    """Grid search over orientations using fitted areas of isolated non-polar lines (unaffected by LO mixing)."""
    f, R = modes(e)
    act = (R ** 2).sum((1, 2)) > 1e-8
    nonpolar = np.array([e['ir_active'][k] == 0 for k in np.argsort(e['freq'])[3:]])
    f, R, nonpolar = f[act], R[act], nonpolar[act]
    pos = s * f + d
    use_sp = [k for k, sp in enumerate(spectra) if sp['nu'].min() < 120 and sp['nu'].max() > 300]
    iso = []
    for m in range(len(f)):
        others = np.abs(pos - pos[m]) < 3.0
        others[m] = False
        if nonpolar[m] and not others.any() and pos[m] > 40:
            iso.append(m)
    st, _ = thermal(s * f, 1e7 / 532.0 - d, prefactors=prefactors)
    areas = []
    for k in use_sp:
        sp = spectra[k]
        inwin = [m for m in range(len(f)) if sp['nu'].min() + 2 < pos[m] < sp['nu'].max() - 2]
        cols = [pvoigt(sp['nu'] - pos[m], 1.0 + 0.0065 * s * f[m]) for m in inwin]
        cols += [np.ones_like(sp['nu']), (sp['nu'] - sp['center']) / 100]
        A = np.array(cols).T
        coef = np.linalg.lstsq(A / np.sqrt(np.maximum(sp['y'], 1) + 25)[:, None],
                               sp['y'] / np.sqrt(np.maximum(sp['y'], 1) + 25), rcond=None)[0]
        areas.append({m: coef[i] / sp['accum'] for i, m in enumerate(inwin)})
    iso = [m for m in iso if all(m in a for a in areas)]
    obs = np.array([[a[m] for a in areas] for m in iso])            # (modes, spectra)
    rots = random_rotations(n_rot).as_matrix()
    Rl = np.einsum('ria,mab,rjb->rmij', rots, R[iso], rots)
    pred = np.zeros((n_rot, len(iso), len(use_sp)))
    for c, k in enumerate(use_sp):
        pred[:, :, c] = config_strengths_batch(Rl, spectra[k])
    pred *= st[iso][None, :, None]
    w = 1.0 / (np.abs(obs) + 0.02 * np.abs(obs).max())
    num = (pred * obs[None] * w[None] ** 2).sum((1, 2))
    den = (pred ** 2 * w[None] ** 2).sum((1, 2))
    cost = ((obs[None] - (num / den)[:, None, None] * pred) ** 2 * w[None] ** 2).sum((1, 2))
    order = np.argsort(cost)
    chosen = []
    for i in order:
        c_ax = rots[i][:, 2]
        if all(abs(float(c_ax @ rots[j][:, 2])) < 0.985 or np.linalg.norm(Rotation.from_matrix(rots[i].T @ rots[j]).as_rotvec()) > 0.35
               for j in chosen):
            chosen.append(i)
        if len(chosen) == keep:
            break
    return [Rotation.from_matrix(rots[i]).as_rotvec() for i in chosen], [float(cost[i]) for i in chosen]


def config_strengths_batch(Rl, sp):
    ei = pol(sp['ai'])
    outs = [pol(0.0), pol(90.0)] if sp['as'] is None else [pol(sp['as'])]
    return sum(np.einsum('i,rmij,j->rm', es, Rl, ei) ** 2 for es in outs)


def fit_candidate(e, spectra, s0, d0, nac=True, prefactors=True, fix_offset=None, n_starts=4, rotvecs=None):
    fit = CandidateFit(e, spectra, nac=nac, prefactors=prefactors, fix_offset=fix_offset)
    if rotvecs is None:
        rotvecs, _ = coarse_orientation(e, spectra, s0, d0 if fix_offset is None else fix_offset,
                                        keep=n_starts, prefactors=prefactors)
    best = None
    for rv in rotvecs:
        o = fit.fit(fit.start(rv, s0, d0), max_nfev=25)
        if best is None or o.chi2 < best.chi2:
            best = o
    best = fit.fit(best.x, max_nfev=150)
    best.fitter = fit
    return best


def minority_test(main_fit, candidates, spectra):
    """Add an orientation-averaged component of each other polymorph to the best main-phase model."""
    fitter, p = main_fit.fitter, main_fit.x
    _, lnC, s, d, g0, g1 = fitter.unpack(p)
    model, coefs = fitter.model(p)
    r0 = np.concatenate([(m - sp['y']) * w for m, sp, w in zip(model, spectra, fitter.w)])
    chi0 = float((r0 ** 2).sum())
    out = {}
    for mpid, e in candidates.items():
        f, R = modes(e)
        keep = (f > 5) & ((R ** 2).sum((1, 2)) > 1e-10)
        f, R = f[keep] * s, R[keep]
        inv = np.array([powder(r) for r in R]) / e['volume']
        st, an = thermal(f, 1e7 / 532.0 - d)
        best = {'weight': 0.0, 'chi2': chi0, 'width_multiplier': None}
        for wmul in (1.0, 1.5, 2.0, 3.0):
            w = wmul * (g0 + g1 * f)
            cols = []
            for sp, c, wt in zip(spectra, coefs, fitter.w):
                x = (sp['nu'] - d)[:, None]
                if sp['as'] is None:
                    I = inv[:, 0] + inv[:, 1]
                else:
                    dd = math.radians(sp['as'] - sp['ai'])
                    I = inv[:, 0] * math.cos(dd) ** 2 + inv[:, 1] * math.sin(dd) ** 2
                prof = st * pvoigt(x - f, w) + an * pvoigt(x + f, w)
                cols.append(sp['accum'] * c[0] * (prof @ I) * wt)
            a = np.concatenate(cols)
            x = max(0.0, -float(a @ r0) / float(a @ a))
            chi = float(((r0 + x * a) ** 2).sum())
            if chi < best['chi2']:
                best = {'weight': x, 'chi2': chi, 'width_multiplier': wmul}
        best['chi2_drop_fraction'] = 1 - best['chi2'] / chi0
        out[mpid] = best
    return chi0, out


def strongest_above(e, fmin=200.0):
    vals = [(sum(float((e['Rsym'][k] ** 2).sum()) for k in b), b) for b in e['bands'] if e['freq'][b[0]] > fmin]
    return float(e['freq'][max(vals, key=lambda v: v[0])[1][0]])


def analyse(catalog, results, use=None, nac=True, prefactors=True, offset='fit', born='rotate', eps='inf',
            candidates=None, minority=True, minority_threshold=0.05, log=print):
    """Full analysis. Variant switches reproduce the wrong paths and the defensible alternatives."""
    spectra = prepare([results[i] for i in (range(len(results)) if use is None else use)])
    recs = catalog['polymorph_records']
    parsed = {m: parse(r) for m, r in recs.items()}
    d0, nu_main = calibrate(spectra)
    fix = 0.0 if offset == 'zero' else None
    log('calibration: axis offset %.3f cm-1, strongest line above 200 cm-1 at %.3f cm-1' % (d0, nu_main))
    fits = {}
    for mpid, e in parsed.items():
        if candidates is not None and mpid not in candidates:
            continue
        s0 = (nu_main + (d0 if fix is not None else 0.0)) / strongest_above(e)
        o = fit_candidate(e, spectra, s0, d0 if fix is None else 0.0, nac=nac, prefactors=prefactors, fix_offset=fix)
        fits[mpid] = o
        log('%s: chi2/point %.3f  s %.5f  offset %.3f  c-axis %s' % (mpid, o.chi2 / o.npts, o.x[4], o.offset,
                                                                    np.round(o.Q @ (e['cell'][2] / np.linalg.norm(e['cell'][2])), 3)))
    main = min(fits, key=lambda m: fits[m].chi2 / fits[m].npts)
    o = fits[main]
    out = {'main_phase': main, 'frequency_scale': float(o.x[4]), 'axis_offset_cm1': o.offset,
           'calibration_offset_cm1': d0, 'chi2_per_point': {m: fits[m].chi2 / fits[m].npts for m in fits},
           'Q': o.Q.tolist(), 'linewidth': [float(o.x[6]), float(o.x[7])], 'fit_params': o.x.tolist()}
    if minority:
        others = {m: e for m, e in parsed.items() if m != main}
        chi0, mt = minority_test(o, others, spectra)
        out['minority_test'] = mt
        best = max(mt, key=lambda m: mt[m]['chi2_drop_fraction'])
        out['minority_phase'] = best if mt[best]['chi2_drop_fraction'] > minority_threshold else 'none'
        log('minority test: %s' % {m: round(v['chi2_drop_fraction'], 4) for m, v in mt.items()})
    pred = face_b_predictions(parsed[main], o.Q, float(o.x[4]), born=born, eps=eps, nac=nac)
    out.update(pred)
    return out


def face_b_predictions(e, Q, s, born='rotate', eps='inf', nac=True, face_deg=30.0):
    b = math.radians(face_deg)
    n = np.array([math.cos(b), math.sin(b), 0.0])
    t, z = np.array([-n[1], n[0], 0.0]), np.array([0.0, 0.0, 1.0])
    band = strongest_band(e)
    Rl = [Q @ e['Rsym'][k] @ Q.T for k in band]
    ratio = sum(float(t @ R @ t) ** 2 for R in Rl) / sum(float(z @ R @ z) ** 2 for R in Rl)
    if not nac or e['Z'] is None:
        lo = s * float(e['freq'].max())
    else:
        Z = e['Z'] if born == 'rotate' else expand_born(e, born)
        f, _ = modes(e, Q.T @ n, Z=Z, eps=np.eye(3) if eps == 'vacuum' else None)
        lo = s * float(f.max())
    c = Q @ (e['cell'][2] / np.linalg.norm(e['cell'][2]))
    return {'ratio_tt_over_zz': ratio, 'highest_optical_cm1': lo, 'c_axis_lab': c.tolist()}


def area_orientation(e, spectra, s, d, rotvec0):
    """Defensible alternative: orientation from fitted areas of isolated non-polar lines only (no LO model)."""
    f, R = modes(e)
    act = (R ** 2).sum((1, 2)) > 1e-8
    nonpolar = np.array([e['ir_active'][k] == 0 for k in np.argsort(e['freq'])[3:]])
    f, R, nonpolar = f[act], R[act], nonpolar[act]
    pos = s * f + d
    use_sp = [k for k, sp in enumerate(spectra) if sp['nu'].min() < 120 and sp['nu'].max() > 300]
    iso = [m for m in range(len(f)) if nonpolar[m] and pos[m] > 40 and (np.abs(pos - pos[m]) < 3.0).sum() == 1]
    st, _ = thermal(s * f, 1e7 / 532.0 - d)
    areas = []
    for k in use_sp:
        sp = spectra[k]
        inwin = [m for m in range(len(f)) if sp['nu'].min() + 2 < pos[m] < sp['nu'].max() - 2]
        cols = [pvoigt(sp['nu'] - pos[m], 1.0 + 0.0065 * s * f[m]) for m in inwin]
        cols += [np.ones_like(sp['nu']), (sp['nu'] - sp['center']) / 100]
        wt = 1.0 / np.sqrt(np.maximum(sp['y'], 1) + 25)
        coef = np.linalg.lstsq(np.array(cols).T * wt[:, None], sp['y'] * wt, rcond=None)[0]
        areas.append({m: coef[i] / sp['accum'] for i, m in enumerate(inwin)})
    iso = [m for m in iso if all(m in a for a in areas)]
    obs = np.array([[a[m] for a in areas] for m in iso])
    w = 1.0 / (np.abs(obs) + 0.02 * np.abs(obs).max())

    def res(p):
        Q = Rotation.from_rotvec(p[:3]).as_matrix()
        Rl = np.einsum('ia,mab,jb->mij', Q, R[iso], Q)
        pred = np.array([config_strengths(Rl, spectra[k]) for k in use_sp]).T * st[iso][:, None]
        return ((obs - math.exp(p[3]) * pred) * w).ravel()
    Q0 = Rotation.from_rotvec(rotvec0).as_matrix()
    Rl = np.einsum('ia,mab,jb->mij', Q0, R[iso], Q0)
    pred = np.array([config_strengths(Rl, spectra[k]) for k in use_sp]).T * st[iso][:, None]
    c0 = float((obs * pred * w * w).sum() / (pred * pred * w * w).sum())
    o = least_squares(res, np.r_[rotvec0, math.log(max(c0, 1e-30))])
    return Rotation.from_rotvec(o.x[:3]).as_matrix()
