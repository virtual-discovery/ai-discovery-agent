"""Reference solution for co3sn2s2-dmi-diffuse; runs as the Harbor oracle and in check_task."""
import itertools
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, minimize

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
STEP = 0.015
EDGES_PIX = CATALOG['profile_annuli_pixels']
QC = 0.5*STEP*(np.array(EDGES_PIX[1:]) + np.array(EDGES_PIX[:-1]))
A_LAT, C_LAT = CATALOG['structure']['a_A'], CATALOG['structure']['c_A']
ASTAR, CSTAR = 4*np.pi/(np.sqrt(3)*A_LAT), 2*np.pi/C_LAT


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(1)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(r)
    return r['id'], r['result']


# ---------------- diffuse scattering ----------------
def profile(res, h, l):
    """Annular profile: pooled counts over pooled exposure monitor*normalization/normcount, Poisson errors.
    Peak pixel = brightest valid pixel within 0.05 1/A of the nominal (h0l); its 7x7 block is excluded."""
    c, m, nc = (np.array(res[k], float) for k in ('counts', 'monitor', 'normcount'))
    n = np.array([[np.nan if v is None else v for v in row] for row in res['normalization']], float)
    cx = (np.array(res['qx_edges'])[:-1] + STEP/2)[:, None]*np.ones(c.shape)
    cy = (np.array(res['qy_edges'])[:-1] + STEP/2)[None, :]*np.ones(c.shape)
    valid = (nc > 0) & (m > 0) & np.isfinite(n) & (n > 0)
    with np.errstate(all='ignore'):
        inten = np.where(valid, c*nc/(m*n), -np.inf)
        expo = np.where(valid, m*n/nc, 0.0)
    near = valid & (np.hypot(cx - h*ASTAR, cy - l*CSTAR) < 0.05)
    i, j = np.unravel_index(np.argmax(np.where(near, inten, -np.inf)), c.shape)
    ii, jj = np.indices(c.shape)
    r2 = (ii - i)**2 + (jj - j)**2
    keep = valid & ~((np.abs(ii - i) <= 3) & (np.abs(jj - j) <= 3))
    vals, errs = [], []
    for lo, hi in zip(EDGES_PIX[:-1], EDGES_PIX[1:]):
        s = keep & (r2 >= lo*lo) & (r2 < hi*hi)
        vals.append(1e6*c[s].sum()/expo[s].sum())
        errs.append(1e6*np.sqrt(c[s].sum())/expo[s].sum())
    return np.array(vals), np.array(errs)


def lorentz(q, A, k):
    return A/(1 + (q/k)**2)


def global_fit(profiles):
    """y_T(q) = b0 + b1 (T - 200) + s q + A_T/(1 + (q/kappa_T)^2), fitted jointly to all (003) profiles."""
    Ts = sorted(profiles)
    n = len(Ts)

    def resid(p):
        b0, b1, s, A, k = p[0], p[1], p[2], p[3:3+n], p[3+n:]
        return np.concatenate([(profiles[T][0] - (b0 + b1*(T - 200) + s*QC + lorentz(QC, A[i], k[i])))/profiles[T][1]
                               for i, T in enumerate(Ts)])
    p0 = np.r_[8.2, 0.002, 0.0, np.full(n, 0.5), np.full(n, 0.1)]
    lb = np.r_[5, -1, -20, np.zeros(n), np.full(n, 0.02)]
    ub = np.r_[12, 1, 20, np.full(n, 30), np.full(n, 0.5)]
    p = least_squares(resid, p0, bounds=(lb, ub), x_scale='jac').x
    return {'Ts': Ts, 'b0': p[0], 'b1': p[1], 's': p[2], 'A': dict(zip(Ts, p[3:3+n])), 'k': dict(zip(Ts, p[3+n:]))}


def critical(g, T):
    if T in g['A']:
        return lorentz(QC, g['A'][T], g['k'][T])
    lo = max(t for t in g['Ts'] if t < T)
    hi = min(t for t in g['Ts'] if t > T)
    w = (T - lo)/(hi - lo)
    return lorentz(QC, (1 - w)*g['A'][lo] + w*g['A'][hi], (1 - w)*g['k'][lo] + w*g['k'][hi])


def linear_q(values, errors):
    return np.polyval(np.polyfit(QC, values, 1, w=1/errors), QC)


# ---------------- spin Hamiltonian ----------------
LAT = np.array([[A_LAT, 0, 0], [-A_LAT/2, A_LAT*np.sqrt(3)/2, 0], [0, 0, C_LAT]])
OPS_TXT = ['x,y,z', '-y,x-y,z', '-x+y,-x,z', 'y,x,-z', 'x-y,-y,-z', '-x,-x+y,-z',
           '-x,-y,-z', 'y,-x+y,-z', 'x-y,x,-z', '-y,-x,z', '-x+y,y,z', 'x,x-y,z']


def parse(op):
    W = np.zeros((3, 3))
    for row, expr in enumerate(op.split(',')):
        for term in expr.replace('-', '+-').split('+'):
            if term:
                W[row, 'xyz'.index(term.lstrip('-'))] = -1 if term.startswith('-') else 1
    return W


OPS = [(parse(s), t) for s in OPS_TXT for t in (np.zeros(3), np.array([2, 1, 1])/3, np.array([1, 2, 2])/3)]
SITES = []
for W_, t_ in OPS:
    p_ = W_ @ np.array([0.5, 0, 0.5]) + t_
    p_ = p_ - np.floor(p_ + 1e-9)
    if not any(np.allclose(p_, q_, atol=1e-6) for q_ in SITES):
        SITES.append(p_)
SITES = np.array(SITES)


def site_of(f):
    d = np.abs(((f - SITES) + 0.5) % 1 - 0.5).sum(1)
    return int(np.argmin(d))


def bond_image(W, t, b):
    i, j, n = b
    ri, rj = W @ SITES[i] + t, W @ (SITES[j] + np.array(n)) + t
    i2, j2 = site_of(ri), site_of(rj)
    shift = np.round(ri - SITES[i2])
    return (i2, j2, tuple(int(x) for x in np.round(rj - shift - SITES[j2])))


def rev(b):
    return (b[1], b[0], tuple(-x for x in b[2]))


def build_hamiltonian():
    """SpinW_Co3Sn2S2_DMI.m: J1 = -56 and J1dmi on bond 1 (2.6785 A), Jc1 = -3 on bond 3 (4.6409 A; bond 2 is the
    in-plane 4.6393 A shell), S = 1, E = sum over bonds of S_i^T J S_j. The J1dmi matrix belongs to the bond
    Co(0,1/2,1/2) -> Co(1/2,1/2,1/2); images carry R J R^T, reversed bonds the transpose."""
    shells = {}
    for i, j in itertools.product(range(9), repeat=2):
        for n in itertools.product(range(-2, 3), repeat=3):
            L = np.linalg.norm((SITES[j] + np.array(n) - SITES[i]) @ LAT)
            if 1e-6 < L < 7:
                shells.setdefault(round(L, 4), []).append((i, j, n))
    lengths = sorted(shells)
    groups = [[lengths[0]]]
    for L in lengths[1:]:
        (groups[-1].append(L) if L - groups[-1][-1] < 1e-3 else groups.append([L]))

    def unique(group):
        out = []
        for L in group:
            for b in shells[L]:
                if b not in out and rev(b) not in out:
                    out.append(b)
        return out
    nn, inter = unique(groups[0]), unique(groups[2])
    D = 15.0
    Mdm = np.array([[0, 0, D], [0, 0, 0], [-D, 0, 0]]) - 0.4*np.array([[0, D, 0], [-D, 0, 0], [0, 0, 0]])
    gen = (site_of(np.array([0, 0.5, 0.5])), site_of(np.array([0.5, 0.5, 0.5])), (0, 0, 0))
    mats = {}
    for W, t in OPS:
        b = bond_image(W, t, gen)
        R = LAT.T @ W @ np.linalg.inv(LAT.T)
        if b not in mats and rev(b) not in mats:
            mats[b] = R @ Mdm @ R.T
    terms = []
    for b in nn:
        M = -56.0*np.eye(3) + (mats[b] if b in mats else mats[rev(b)].T)
        terms.append((b[0], b[1], (SITES[b[1]] + np.array(b[2]) - SITES[b[0]]) @ LAT, M))
    for b in inter:
        terms.append((b[0], b[1], (SITES[b[1]] + np.array(b[2]) - SITES[b[0]]) @ LAT, -3.0*np.eye(3)))
    return terms


def ground_state(terms):
    I = np.array([t[0] for t in terms])
    J = np.array([t[1] for t in terms])
    Ms = np.array([t[3] for t in terms])

    def spins(x):
        th, ph = x[:9], x[9:]
        return np.c_[np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)]

    def energy(x):
        S = spins(x)
        return np.einsum('ba,bac,bc->', S[I], Ms, S[J])/9
    rng = np.random.default_rng(7)
    best = None
    for _ in range(12):
        r = minimize(energy, np.r_[rng.uniform(0, np.pi, 9), rng.uniform(0, 2*np.pi, 9)], method='BFGS',
                     options={'gtol': 1e-11})
        if best is None or r.fun < best.fun - 1e-10:
            best = r
    S = spins(best.x)
    return best.fun, (S if S[:, 2].sum() > 0 else -S)


def lswt_gamma(S, terms, spin=1.0):
    """Linear spin waves at k = 0 (Holstein-Primakoff, paraunitary eigenproblem)."""
    N = len(S)
    u, v = [], []
    for s in S:
        t = np.array([0, 0, 1.0]) if abs(s[2]) < 0.9 else np.array([1.0, 0, 0])
        e1 = np.cross(t, s)
        e1 /= np.linalg.norm(e1)
        u.append(e1 + 1j*np.cross(s, e1))
        v.append(s)
    u, v = np.array(u), np.array(v)
    A = np.zeros((N, N), complex)
    B = np.zeros((N, N), complex)
    C = np.zeros(N)
    for i, j, d, M in terms:
        A[i, j] += spin/2*(u[i] @ M @ np.conj(u[j]))
        A[j, i] += spin/2*(np.conj(u[i]) @ M @ u[j])
        B[i, j] += spin/2*(u[i] @ M @ u[j])
        B[j, i] += spin/2*(u[i] @ M @ u[j])
        C[i] += spin*(v[i] @ M @ v[j])
        C[j] += spin*(v[i] @ M @ v[j])
    A -= np.diag(C)
    Hk = np.block([[A, B], [B.conj().T, A.T]])
    w = np.linalg.eigvals(np.diag(np.r_[np.ones(N), -np.ones(N)]) @ Hk)
    return np.sort(w.real)[N:]


def main():
    evidence, prof003, prof101 = [], {}, {}
    for T in (140, 178, 185, 192, 220):
        i, res = buy({'kind': 'narrow_scan', 'temperature_K': T})
        evidence.append(i)
        prof003[T] = profile(res, 0, 3)
    for T, (h, l) in ((165, (0, 3)), (180, (0, 3)), (165, (1, 1)), (300, (1, 1))):
        i, res = buy({'kind': 'wide_region', 'temperature_K': T, 'hkl': [h, 0, l]})
        evidence.append(i)
        (prof003 if (h, l) == (0, 3) else prof101)[T] = profile(res, h, l)
    # Narrow (003) profiles: the critical term rises steeply between 165 K (wide scan) and 178 K, so the
    # withheld 170 and 175 K scans are bracketed with the 165 K wide-scan point; the background is linear in T.
    g = global_fit(prof003)
    narrow = {str(T): [float(x) for x in g['b0'] + g['b1']*(T - 200) + g['s']*QC + critical(g, T)] for T in (170, 175, 200)}
    # (101) at 180 K: same |Q| and form factor as (003); ferromagnetic correlations have |F|^2 = 1 at (101) and 9 at
    # (003), and for isotropic spin components the polarisation sums are equal (2 at both reflections).
    w = (180 - 165)/(300 - 165)
    bg101 = (1 - w)*linear_q(*prof101[165]) + w*linear_q(*prof101[300])
    wide101 = [float(x) for x in bg101 + critical(g, 180)/9]
    # Spin Hamiltonian of the archived SpinW script.
    terms = build_hamiltonian()
    E, S = ground_state(terms)
    tilt = float(np.degrees(np.arccos(np.clip(np.abs(S[:, 2]), 0, 1))).mean())
    order = [site_of(np.array([0.5, 0, 0.5])), site_of(np.array([0, 0.5, 0.5])), site_of(np.array([0.5, 0.5, 0.5]))]
    rep = {}
    for k in range(9):
        rep[k] = next(o for o in range(3) if any(np.allclose((SITES[k] - SITES[order[o]] + 1e-9) % 1, c, atol=1e-6)
                                                 for c in (np.zeros(3), np.array([2, 1, 1])/3, np.array([1, 2, 2])/3)))
    assert all(np.allclose(S[k], S[order[rep[k]]], atol=1e-4) for k in range(9))   # k = 0 in the primitive cell
    gamma = lswt_gamma(S[order], [(rep[i], rep[j], d, M/3) for i, j, d, M in terms])
    gap = float(gamma.min())
    submission = {'narrow_003_profiles': narrow, 'wide_101_profile_180K': wide101,
                  'canting_angle_deg': tilt, 'spin_wave_gap_meV': gap, 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    amp = {str(T): round(float(g['A'][T]), 2) for T in g['Ts']}
    (APP/'report.md').write_text(f"""# Co3Sn2S2: critical diffuse scattering and the DMI-canted state (reference solution)

Data. Annular profiles around (0 0 3) were reduced with the archived normalisation (counts*normcount/(monitor*
normalization)), pooling counts and exposure per annulus and excluding the 7x7 Bragg block. Purchases: narrow scans at
140, 178, 185, 192 and 220 K; wide-scan (003) regions at 165 and 180 K; wide-scan (101) regions at 165 and 300 K.

Critical scattering. A Lorentzian centred on (003) appears near T_C and dominates the first annuli; its amplitude
(fit: {amp}) is largest at 178-180 K, falls slowly above T_C and collapses within ~10 K below it: the wide 165 K scan
shows no excess at (003). The narrow scans alone jump from 140 to 178 K, so the 165 K wide-scan point is needed to
predict 170 and 175 K; the scattering is not symmetric in |T - T_C|. The background rises smoothly with T
(b1 = {g['b1']:.4f} per K, in 1e-6 units) and was modelled as linear in T.

Which spins fluctuate. At (003) Q is parallel to c, so only in-plane spin components scatter: the critical signal
shows that in-plane ferromagnetic fluctuations become critical at T_C, which a strongly Ising (single-ion) easy axis
would suppress; the effective anisotropy is weak. (101) has the same |Q| (1.436 1/A) and form factor, but the
ferromagnetic structure factor per kagome triangle is 1 there versus 9 at (003), while an in-plane q = 0 120-degree
(antiferromagnetic) component has |F|^2 = 4 at (101) and 0 at (003). The 180 K (101) prediction is therefore the
smooth (101) background (interpolated between 165 and 300 K) plus 1/9 of the (003) critical term; an incipient
in-plane AFM order would instead put a strong peak at (101). The 165 and 1.5 K (101) data show no such signal, and
there is a single critical anomaly at T_C, so these data give no evidence for a separate intermediate-temperature AFM
phase. Bragg intensities in these scans change non-monotonically even between 180 and 300 K (above T_C), so they
were not used for magnetic structure.

Hamiltonian. From the SpinW script: J1 = -56 meV (FM), DM matrix with in-plane component 15 meV perpendicular to
the bond and out-of-plane 6 meV (|D|/|J1| = 0.29), Jc1 = -3 meV on the 4.641 A interlayer bond (bond 3; bond 2 is
the 4.639 A in-plane shell), S = 1, each bond counted once. The bond-1 matrix was propagated with the R-3m operations
(R J R^T, transpose for reversed bonds). The classical ground state is a k = 0 umbrella: moments tilted
{tilt:.2f} deg from c with in-plane components forming a 120-degree pattern (energy {E:.3f} meV per Co). The
in-plane DM component produces the canting (linear in D); the out-of-plane part alone would leave a collinear state.
LSWT about this state gives a k = 0 gap of {gap:.2f} meV, produced entirely by the DMI (no single-ion term).
The measured canting (about 17 deg) and gap (2-3 meV) are not reproduced exactly by this effective model.

Not supported by these data: a value of D/J1 from the diffuse scattering itself (the critical scattering near
(003) does not separate D and J); absolute ordered moments or canting angles from the Bragg peaks; a critical
exponent with useful precision (it depends on T_C and on the background model).
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
