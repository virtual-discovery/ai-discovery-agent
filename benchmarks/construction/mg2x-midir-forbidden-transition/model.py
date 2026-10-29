"""Readers and estimators for the mg2x-midir-forbidden-transition construction (numpy only).

The estimators port the archive's Julia programs (MIT, H. Sugawara 2026) to Python:
  analysis/qbox.jl            nested boxes around X, trapezoid weights, Boltzmann <q_perp^2>
  analysis/ir_raman_ratio.jl  Fermi-Dirac <q_perp^2> on the boxes, the 300 K / 77 K ratio
  analysis/effective_masses.jl parabolic masses from the fine scans, thermal average along X-W
and are checked against the archive's JSON outputs in make_task.py.

ARCHIVE (default ~/.cache/ai-theorist/mg2x-midir-forbidden-transition/x) holds the extracted members.
"""
import csv
import glob
import json
import os
from pathlib import Path
import numpy as np

ARCHIVE = Path(os.environ.get('MG2X_ARCHIVE', Path.home()/'.cache/ai-theorist/mg2x-midir-forbidden-transition/x'))
KB_EV = 8.617333262e-5
HARTREE_EV = 27.211386245988
BOHR_ANG = 0.529177210903
RY_EV = 13.605693122994

# compound: (tag, lattice constant in A, number of valence bands in the ecalj numbering, lower, upper band at X,
#            X1 band, X3 band)  -- X1/X3 from the transverse/longitudinal masses (archive: quasiparticle_path.jl)
COMPOUNDS = {
    'Mg2Si': dict(tag='si', a=6.3380, nv=4, lo=5, hi=6, x1=5, x3=6),
    'Mg2Ge': dict(tag='ge', a=6.3933, nv=9, lo=10, hi=11, x1=10, x3=11),
    'Mg2Sn': dict(tag='sn', a=6.7601, nv=9, lo=10, hi=11, x1=11, x3=10),
}
BOXES = {
    ('Mg2Si', 'coarse'): ('qsgw-qbox2-scf8-cont', 0.02, (0.24, 0.16, 0.16)),
    ('Mg2Si', 'fine_0.01'): ('qsgw-qbox-step0.01-scf8-cont', 0.01, (0.20, 0.14, 0.14)),
    ('Mg2Si', 'fine_0.005'): ('qsgw-qbox-step0.005-scf8-cont', 0.005, (0.10, 0.07, 0.07)),
    ('Mg2Ge', 'coarse'): ('mg2ge-qbox-scf8-cont', 0.02, (0.24, 0.16, 0.16)),
    ('Mg2Sn', 'coarse'): ('mg2sn-qbox-scf8-cont', 0.02, (0.24, 0.16, 0.16)),
}


# ----------------------------------------------------------------------------- band files
def read_bnd(path, bands=None):
    """Rows (band, path length, energy eV, k1, k2, k3) of an ecalj bnd*.spin1 file; '#' lines dropped."""
    rows = []
    for ln in open(path):
        s = ln.split()
        if not s or s[0].startswith('#'):
            continue
        b = int(float(s[0]))
        if bands is None or b in bands:
            rows.append((b, float(s[1]), float(s[2]), float(s[4]), float(s[5]), float(s[6])))
    return rows


def path_bands(compound, level):
    """L-Gamma-X-W-Gamma path: list of segments, each a dict band -> list of (x, k, E)."""
    c = COMPOUNDS[compound]
    d = ARCHIVE/f"data/bands/mg2{c['tag']}-{level.lower()}"
    out = []
    for f in sorted(glob.glob(str(d/'bnd00*.spin1'))):
        seg = {}
        for b, x, e, k1, k2, k3 in read_bnd(f):
            seg.setdefault(b, []).append((x, (k1, k2, k3), e))
        for b in seg:
            seg[b].sort()
        out.append(seg)
    return out


def x_scan(compound, level, direction):
    """Fine scan at X: q from X (units of 2 pi / a) and energies per band (eV)."""
    c = COMPOUNDS[compound]
    f = {'X-Gamma': 'bnd001.spin1', 'X-W': 'bnd002.spin1'}[direction]
    rows = read_bnd(ARCHIVE/f"data/mass-scan/mg2{c['tag']}-{level.lower()}/{f}")
    out = {}
    for b, x, e, k1, k2, k3 in rows:
        out.setdefault(b, []).append((x, e))
    for b in out:
        out[b].sort()
        x0 = out[b][0][0]
        out[b] = [(round(x - x0, 6), e) for x, e in out[b]]
    return out


def load_box(compound, box, bands):
    """Points (q_par, q_y, q_z) from X in units of 2 pi / a (one octant) and energies of `bands` (eV)."""
    d = ARCHIVE/'data/qbox'/BOXES[(compound, box)][0]
    seen = {b: {} for b in bands}
    for f in sorted(glob.glob(str(d/'bnd*.spin1'))):
        for b, x, e, k1, k2, k3 in read_bnd(f, set(bands)):
            key = (round(k1 - 1.0, 6), round(k2, 6), round(k3, 6))
            seen[b].setdefault(key, e)
    keys = sorted(seen[bands[0]])
    for b in bands:
        assert sorted(seen[b]) == keys
    return np.array(keys), {b: np.array([seen[b][k] for k in keys]) for b in bands}


def box_trap(K, upper):
    w = np.ones(len(K))
    for j in range(3):
        w *= np.where(np.abs(K[:, j]) < 1e-9, 0.5, 1.0)
        w *= np.where(np.abs(K[:, j] - upper[j]) < 1e-9, 0.5, 1.0)
    return w


def nested(compound, names, band):
    """Nested boxes (coarse to fine): K, energy from the band minimum at X (eV), weights (2 pi / a)^3."""
    data = []
    for name in names:
        K, E = load_box(compound, name, [band])
        _, h, r = BOXES[(compound, name)]
        data.append((K, E[band], h, r))
    i0 = lambda K: int(np.argmin(np.abs(K).sum(1)))
    EX = data[0][1][i0(data[0][0])]
    Ks, Es, Ws = [], [], []
    for i, (K, E, h, r) in enumerate(data):
        assert abs(E[i0(K)] - EX) < 1e-8, 'X energy differs between boxes'
        w = box_trap(K, r) * h ** 3
        if i < len(data) - 1:
            inner = data[i + 1][3]
            inside = np.all(K <= np.array(inner) + 1e-9, axis=1)
            w = w - np.where(inside, box_trap(K, inner), 0.0) * h ** 3
        Ks.append(K); Es.append(E - EX); Ws.append(w)
    return np.vstack(Ks), np.concatenate(Es), np.concatenate(Ws)


def qperp2(K):
    return K[:, 1] ** 2 + K[:, 2] ** 2


def boltzmann_mean(box, T, weight=qperp2):
    K, e, w = box
    f = w * np.exp(-e / (KB_EV * T))
    return float((f * weight(K)).sum() / f.sum())


def fd_mean(box, n_cm3, T, a_ang, weight=qperp2, degeneracy=6.0, octants=8.0):
    """Fermi-Dirac average over the box: n = degeneracy x octants x (2 pi/a)^3/(2 pi)^3 x sum w f."""
    K, e, w = box
    G = 2 * np.pi / (a_ang * 1e-10)
    vol = octants * degeneracy * G ** 3 / (2 * np.pi) ** 3
    kT = KB_EV * T
    occ = lambda mu: w / (1 + np.exp(np.clip((e - mu) / kT, -700, 700)))
    lo, hi = -2.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if vol * occ(mid).sum() - n_cm3 * 1e6 < 0:
            lo = mid
        else:
            hi = mid
    f = occ((lo + hi) / 2)
    return float((f * weight(K)).sum() / f.sum())


# ----------------------------------------------------------------------------- parabolic valley
def parabolic_mass(q, e, a_ang, qfit=0.020):
    k2pi = 2 * np.pi / (a_ang / BOHR_ANG)
    q, e = np.asarray(q), np.asarray(e)
    keep = q <= qfit + 1e-9
    x = (q[keep] * k2pi) ** 2
    y = (e[keep] - e[0]) / HARTREE_EV
    slope = np.polyfit(x, y, 1)[0]
    return 1 / (2 * slope)


def fd_integral(j, eta, n=4000):
    hi = np.sqrt(60.0 + max(eta, 0.0))
    u = np.linspace(0, hi, n + 1)
    f = 2 * u ** (2 * j + 1) / (1 + np.exp(np.clip(u ** 2 - eta, -700, 700)))
    return float(np.sum((f[1:] + f[:-1]) * np.diff(u)) / 2)


def parabolic_ratio(n_cm3, m_dos, valleys=3, T_hot=300.0, T_cold=77.0):
    """<q_perp^2> ratio of a parabolic valley with Fermi-Dirac statistics (ir_raman_ratio.jl)."""
    me, hbar, ev = 9.1093837139e-31, 1.054571817e-34, 1.602176634e-19

    def mean_E(T):
        nc = 2 * (m_dos * me * KB_EV * T * ev / (2 * np.pi * hbar ** 2)) ** 1.5 * valleys * 1e-6
        g = lambda eta: nc * 2 / np.sqrt(np.pi) * fd_integral(0.5, eta) - n_cm3
        lo, hi = -80.0, 80.0
        for _ in range(200):
            mid = (lo + hi) / 2
            if g(mid) < 0:
                lo = mid
            else:
                hi = mid
        eta = (lo + hi) / 2
        return T * fd_integral(1.5, eta) / fd_integral(0.5, eta)
    return mean_E(T_hot) / mean_E(T_cold)


def line_isotropic_ratio(q, e, a_ang, T_hot=300.0, T_cold=77.0):
    """effective_masses.jl thermal_q2: the X-W scan taken as isotropic in the transverse plane (measure q dq)."""
    k2pi = 2 * np.pi / (a_ang / BOHR_ANG)
    k = np.asarray(q) * k2pi
    en = (np.asarray(e) - e[0]) / HARTREE_EV
    outer = k > 0.6 * k.max()
    A = np.c_[np.ones(outer.sum()), k[outer] ** 2, k[outer] ** 4]
    c = np.linalg.lstsq(A, en[outer], rcond=None)[0]
    grid = np.linspace(0, 1.6 * k.max(), 4000)
    ee = np.where(grid <= k.max(), np.interp(grid, k, en), c[0] + c[1] * grid ** 2 + c[2] * grid ** 4)

    def q2(T):
        w = np.exp(-ee / (T * 3.166811563e-6)) * grid
        return np.trapz(grid ** 2 * w, grid) / np.trapz(w, grid)
    return q2(T_hot) / q2(T_cold)


# ----------------------------------------------------------------------------- other archive outputs
def tote_rows(path, kpts=((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)), bands=None):
    """Header line and rows of an ecalj TOTE file at the given q (units of 2 pi / a), as text columns."""
    lines = open(path).read().splitlines()
    rows = []
    for ln in lines[1:]:
        f = ln.replace('D', 'E').split()
        if len(f) < 9:
            continue
        q = tuple(float(x) for x in f[:3])
        if any(np.allclose(q, k, atol=1e-6) for k in kpts) and (bands is None or int(f[3]) in bands):
            rows.append([float(x) for x in f[:3]] + [int(f[3]), int(f[4])] + [float(x) for x in f[5:9]])
    return lines[0].strip(), rows


def tote_split(path, x1, x3, col):
    """X3 - X1 (eV) from a TOTE file; col is the 0-based column of the energy (5: input eigenvalue, 6, 7)."""
    _, rows = tote_rows(path, kpts=((0.0, 1.0, 0.0),))
    E = {r[3]: r[col] for r in rows}
    return E[x3] - E[x1]


def questaal_levels(path):
    """Eigenvalues (Ry) at Gamma and X = (0, 1, 0) from an extracted Questaal .eval.txt file."""
    out, q = {}, None
    for ln in open(path):
        if ln.startswith('#'):
            continue
        if ln.startswith('q '):
            q = tuple(float(x) for x in ln.split()[2:5])
            continue
        if q is not None:
            e = [float(x) for x in ln.split()]
            if np.allclose(q, 0):
                out['Gamma'] = e
            elif np.allclose(q, (0, 1, 0)):
                out['X'] = e
            q = None
    return out


def bgw_rows(path):
    """eqp1.dat of BerkeleyGW: {k-point tuple (crystal): [(band, E_DFT, E_qp)]}."""
    out, k = {}, None
    for ln in open(path):
        f = ln.split()
        if len(f) == 4 and '.' in f[0]:
            k = tuple(round(float(x), 6) for x in f[:3])
            out[k] = []
        elif k is not None and len(f) >= 4:
            out[k].append((int(f[1]), float(f[2]), float(f[3])))
    return out


def pmat_rows():
    rows = list(csv.DictReader(open(ARCHIVE/'data/pmat-pbe/box_pmat.csv')))
    return [{'q': (round(float(r['kx']) - 1, 6), round(float(r['ky']), 6), round(float(r['kz']), 6)),
             'E13': float(r['E13_eV']), 'E14': float(r['E14_eV']), 'P2': float(r['P2_13_14'])} for r in rows]


def eph_json():
    return json.loads((ARCHIVE/'eph/analysis/eph_renormalization.json').read_text())
