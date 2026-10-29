"""Reference solution for mnte-altermagnet-chiral-magnons; runs as the Harbor oracle and in check_task."""
import json
import os
import re
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq, least_squares

APP = Path(os.environ.get('APP_DIR', '/app'))
CAT = json.loads((APP/'catalog.json').read_text())


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    while result.get('status') == 'running':
        time.sleep(1)
        result = lab('job/%d' % result['id'])
    return result


# ---------------------------------------------------------------- structure and bond classes (P6_3/mmc)
A_LAT, C_LAT = CAT['structure']['a_A'], CAT['structure']['c_A']
CELL = np.array([[A_LAT, 0, 0], [-A_LAT / 2, A_LAT * np.sqrt(3) / 2, 0], [0, 0, C_LAT]])
BMAT = 2 * np.pi * np.linalg.inv(CELL).T
MN = np.array([[0, 0, 0], [0, 0, 0.5]])
GEN = ['x,y,z', '-y,x-y,z', '-x+y,-x,z', '-x,-y,z+1/2', 'y,-x+y,z+1/2', 'x-y,x,z+1/2',
       'y,x,-z', 'x-y,-y,-z', '-x,-x+y,-z', '-y,-x,-z+1/2', '-x+y,y,-z+1/2', 'x,x-y,-z+1/2']


def parse_op(s):
    R, t = np.zeros((3, 3)), np.zeros(3)
    for r, expr in enumerate(s.split(',')):
        for sign, term in re.findall(r'([+-]?)([xyz]|\d+/\d+)', expr):
            v = -1.0 if sign == '-' else 1.0
            if term in 'xyz':
                R[r, 'xyz'.index(term)] += v
            else:
                a, b = term.split('/')
                t[r] += v * float(a) / float(b)
    return R, t


OPS = [parse_op(s) for s in GEN]
OPS += [(-R, -t) for R, t in OPS]


def site(frac):
    for k, p in enumerate(MN):
        d = frac - p
        if np.allclose(d, np.round(d), atol=1e-6):
            return k, np.round(d).astype(int)


def same_sublattice_classes(dist, tol=1e-3):
    """Symmetry classes of the undirected same-sublattice bonds of length dist, as {class: [(i, delta)]}."""
    bonds = []
    for i in range(2):
        for R in np.array(np.meshgrid(*[range(-4, 5)] * 3)).reshape(3, -1).T:
            if abs(np.linalg.norm(R @ CELL) - dist) < tol:
                bonds.append((i, tuple(R)))
    cls = {}
    for b in bonds:
        if b in cls:
            continue
        n = len(set(cls.values()))
        for Rot, t in OPS:
            ri = Rot @ MN[b[0]] + t
            rj = Rot @ (MN[b[0]] + np.array(b[1])) + t
            ki, Ti = site(ri)
            kj, Tj = site(rj)
            for img in ((ki, tuple(Tj - Ti)), (ki, tuple(Ti - Tj))):
                cls.setdefault(img, n)
    out = {}
    for (i, R), c in cls.items():
        out.setdefault(c, []).append((i, np.array(R, float)))
    return out


# ---------------------------------------------------------------- linear spin waves
S = 2.5
PRIOR = CAT['prior_hamiltonian']
D_ANI = PRIOR['D_meV']


def bonds_model(dJ):
    """Bond list (i, j, delta_frac, J), every undirected pair once: prior shells plus the two 9.83 A classes."""
    out = []
    for R in np.array(np.meshgrid(*[range(-3, 4)] * 3)).reshape(3, -1).T:
        dl = R + MN[1] - MN[0]                       # A -> B pairs, each once
        r = np.linalg.norm(dl @ CELL)
        for dist, J in ((3.3555, 3.99), (5.3353, 0.472)):
            if abs(r - dist) < 1e-3:
                out.append((0, 1, dl, J))
        if tuple(R) > (0, 0, 0) and abs(np.linalg.norm(R @ CELL) - 4.148) < 1e-3:
            for i in range(2):
                out.append((i, i, R.astype(float), -0.120))
    for k, c in enumerate(sorted(CLASSES)):
        Jc = JBAR + (dJ / 2 if k == 0 else -dJ / 2)
        for i, R in CLASSES[c]:
            if tuple(R) > (0, 0, 0):
                out.append((i, i, R, Jc))
    return out


JBAR = 0.0230
CLASSES = same_sublattice_classes(9.8313)


def frames():
    t = np.radians(30.0)
    n = np.array([np.cos(t), np.sin(t), 0.0])
    zc = np.array([0, 0, 1.0])
    out = []
    for s in (1, -1):
        e = s * n
        u = np.cross(zc, e)
        w = np.cross(e, u)
        out.append((e, u - 1j * w))
    return out


FR = frames()


def blocks(k, bonds):
    A = np.zeros((2, 2), complex)
    B = np.zeros((2, 2), complex)
    for i, j, dl, J in bonds:
        (ei, ci), (ej, cj) = FR[i], FR[j]
        p = np.exp(2j * np.pi * np.dot(dl, k))
        K = J * S / 2
        A[i, j] += K * np.dot(np.conj(ci), cj) * p
        A[j, i] += K * np.dot(ci, np.conj(cj)) * np.conj(p)
        b = K * np.dot(np.conj(ci), np.conj(cj))
        B[i, j] += b * p
        B[j, i] += b * np.conj(p)
        A[i, i] -= J * S * np.dot(ei, ej)
        A[j, j] -= J * S * np.dot(ei, ej)
    for s in range(2):
        cz = FR[s][1][2]
        A[s, s] += S * D_ANI * abs(cz) ** 2
        B[s, s] += S * D_ANI * np.conj(cz) ** 2
    return A, B


def energies(k, bonds):
    A, B = blocks(np.asarray(k, float), bonds)
    Am, _ = blocks(-np.asarray(k, float), bonds)
    M = np.block([[A, B], [B.conj().T, Am.T]])
    M = (M + M.conj().T) / 2
    g = np.diag([1, 1, -1, -1.0])
    ev = np.linalg.eigvals(g @ M).real
    return np.sort(ev[ev > 0])[:2]


# ---------------------------------------------------------------- data
def fit_two_gauss(x, y, e, guess):
    f = lambda p: p[0] + p[6] * (x - x.mean()) + p[1] * np.exp(-0.5 * ((x - p[2]) / (p[5] / 2.3548)) ** 2) \
        + p[3] * np.exp(-0.5 * ((x - p[4]) / (p[5] / 2.3548)) ** 2)
    return least_squares(lambda p: (f(p) - y) / e, list(guess) + [0.0]).x


def main():
    evidence = []
    got = {}
    for req in [{'kind': 'hyspec_cut', 'cut': c} for c in 'ghij'] + [{'kind': 'pnd', 'instrument': 'HYSPEC'},
                                                                     {'kind': 'pnd', 'instrument': 'IN20', 'state': 'A'},
                                                                     {'kind': 'pnd', 'instrument': 'IN20', 'state': 'B'}]:
        r = lab('evaluate', req)
        evidence.append(r['id'])
        got[r['result']['item']] = r['result']
    eg = CAT['hyspec']['energy_grid_meV']
    Emap = RegularGridInterpolator((np.array(eg['u']), np.array(eg['v'])), np.array(eg['E']), bounds_error=False, fill_value=None)
    Ehys = lambda u, v: float(Emap([[u, v]])[0])
    uv = lambda u, v: np.array([u + v, 0.0, -2.23 * u + 1.6 * v])
    # 1. peak positions in the symmetric cuts (two Gaussians on a sloping background)
    cuts = {'g': ('u', 0.45, (0.80, 0.95), [0.5, 0.7, 0.845, 0.45, 0.88, 0.03]),
            'h': ('v', 1.00, (0.52, 0.72), [0.7, 0.5, 0.60, 0.6, 0.635, 0.03])}
    pos = {}
    for c, (ax, fixed, rng, guess) in cuts.items():
        d = got['hyspec_cut:' + c]
        x, y, e = (np.array(d[k], float) for k in ('x', 'intensity', 'error'))
        p = fit_two_gauss(x, y * 1e4, e * 1e4, guess)
        pos[c] = sorted([p[2], p[4]])

    def crossings(bonds, c):
        ax, fixed, rng, _ = cuts[c]
        xs = np.linspace(rng[0], rng[1], 31)
        out = []
        for m in range(2):
            f = lambda x: energies(uv(*((x, fixed) if ax == 'u' else (fixed, x))), bonds)[m] - Ehys(*((x, fixed) if ax == 'u' else (fixed, x)))
            v = [f(x) for x in xs]
            out += [(brentq(f, xs[i], xs[i + 1]), m) for i in range(30) if v[i] * v[i + 1] < 0]
        return sorted(out)
    # 2. exchange difference of the two 9.83 A classes from the splitting seen in both cuts
    grid = np.arange(0.03, 0.161, 0.01)
    est = {}
    for c in cuts:
        seps = []
        for dJ in grid:
            cr = crossings(bonds_model(dJ), c)
            seps.append(cr[-1][0] - cr[0][0])
        est[c] = float(np.interp(pos[c][1] - pos[c][0], seps, grid))
    dJ = (est['g'] + est['h']) / 2
    model = bonds_model(dJ)
    # 3. chirality of the branches in HYSPEC: sign of the chiral cut at the two peak positions of cut g
    ci = got['hyspec_cut:i']
    x, y, e = (np.array(ci[k], float) for k in ('x', 'intensity', 'error'))
    w = 0.03 / 2.3548
    G = np.c_[np.exp(-0.5 * ((x - pos['g'][0]) / w) ** 2), np.exp(-0.5 * ((x - pos['g'][1]) / w) ** 2)]
    amp = np.linalg.lstsq(G / e[:, None], y / e, rcond=None)[0]
    branch_g = [m for _, m in crossings(model, 'g')]          # branch index (0 lower, 1 upper) of each peak along u
    upper_sign_hyspec = np.sign(amp[branch_g.index(1)])
    pnd = got['pnd:HYSPEC']
    nmi_hyspec = np.sign(sum(pnd['I_y_plus0']) - sum(pnd['I_y_minus0']))

    def nmi_in20(raw):
        lines = raw.splitlines()
        k = next(i for i, l in enumerate(lines) if l.startswith('DATA_:'))
        cols = lines[k + 1].split()
        rows = np.array([l.split() for l in lines[k + 2:] if len(l.split()) == len(cols)], float)
        py, i9, i10, cnt, mon = (rows[:, cols.index(c)] for c in ('PY1', 'I9', 'I10', 'CNTS', 'M1'))
        y = py > 0.5
        pp = cnt[y & (i9 > 0) & (i10 > 0)].sum() / mon[y & (i9 > 0) & (i10 > 0)].sum()
        mm = cnt[y & (i9 < 0) & (i10 < 0)].sum() / mon[y & (i9 < 0) & (i10 < 0)].sum()
        return np.sign(pp - mm)
    # 4. IN20: instrument lattice in the scan headers -> true l of the nominal (h, 0, 1.25) scans
    head = CAT['in20']['withheld_inelastic_scans']['scans'][0]['header']
    cs = float(re.search(r'CS=\s*([-0-9.]+)', head).group(1))
    l_true = 1.25 * C_LAT / cs
    hs = np.linspace(1.2, 1.5, 31)
    cr = []
    for m in range(2):
        f = lambda h: energies([h, 0, l_true], model)[m] - 32.0
        v = [f(h) for h in hs]
        cr += [(brentq(f, hs[i], hs[i + 1]), m) for i in range(30) if v[i] * v[i + 1] < 0]
    cr = sorted(cr)
    peaks = {'low': cr[0][0], 'high': cr[-1][0]}
    low_is_upper = cr[0][1] == 1
    chir = {}
    for st in ('A', 'B'):
        s_upper = upper_sign_hyspec * nmi_in20(got['pnd:IN20:' + st]['raw_file']) * nmi_hyspec
        s_low = s_upper if low_is_upper else -s_upper
        chir[st] = {'low': '+' if s_low > 0 else '-', 'high': '-' if s_low > 0 else '+'}
    # 5. wave vector of maximal altermagnetic splitting: exchange-only form factor of the split classes
    best, qbest = -1, None
    cl = [CLASSES[c] for c in sorted(CLASSES)]
    for h in np.linspace(0, 1, 61):
        for k in np.linspace(0, 1, 61):
            for l in np.linspace(0, 0.5, 21):
                q = np.array([h, k, l])
                val = abs(sum(np.cos(2 * np.pi * q @ R) for i, R in cl[0] if i == 0) - sum(np.cos(2 * np.pi * q @ R) for i, R in cl[1] if i == 0))
                if val > best + 1e-9:
                    best, qbest = val, q
    sub = {'split_exchange_difference_meV': dJ, 'in20_peak_qh': peaks, 'in20_chirality': chir,
           'max_splitting_q_rlu': [float(x) for x in qbest], 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(sub, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# MnTe chiral magnons (reference solution)

Model: H = sum_pairs J S_i.S_j + D sum (S^z)^2 with the prior J1-J3, D and the 9.83 A same-sublattice bonds. Spin-group
symmetry splits the 12 bonds at 9.83 A into two classes of six that are exchanged between the sublattices (the Mn
sites are related by the 6_3 screw, not by translation or inversion); all shorter same-sublattice shells stay single
classes. Only J_a - J_b lifts the degeneracy of the two magnon branches, which then carry opposite chirality.
The splitting form factor is proportional to sin(2 pi l) [sin 2pi(2h+k) + sin 2pi(k-h) - sin 2pi(h+2k)]: nodal on l = 0, 1/2
and on the (h h l) planes (g-wave), maximal at {sub['max_splitting_q_rlu']}. The 4 meV difference of the branches at the zone centre
is the single-ion (easy-plane) anisotropy gap, not altermagnetic splitting.
J_a - J_b: from the peak separations in the HYSPEC cuts g and h (branch crossings with the kinematic energy of each
pixel): {est} meV, mean {dJ:.4f} meV.
IN20: the scan headers use a = 4.17, c = {cs} A instead of 4.148 / 6.711 A, so nominal l = 1.25 is l = {l_true:.4f}.
Predicted 32 meV peaks at h = {peaks['low']:.4f} and {peaks['high']:.4f} (the low-h peak is the upper branch).
Chirality: in HYSPEC the upper branch has chiral sign {upper_sign_hyspec:+.0f} with NMI sign {nmi_hyspec:+.0f} at (1 0 -1). The NMI sign
and the splitting form factor both change sign with l, so the rule (upper-branch chirality = HYSPEC sign x NMI ratio) carries
over to IN20, whose PND (1 0 1) and inelastic scans share the instrument frame. Only the y-channel NMI measures the
a*-projection of the net Neel vector that the chiral term at (h 0 l) needs; the z channel measures the b projection.
Field cooling reverses the antiphase-domain balance and therefore every chiral sign. Relativistic (DM) terms are
forbidden on the Mn inversion centres; magnon-phonon mixing cannot produce opposite-sign chiral peaks that follow the
domain state. The chiral amplitudes are fractions of the symmetric intensity set by the domain imbalance, which these
data constrain only roughly.
""")
    lab('submit', sub)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
