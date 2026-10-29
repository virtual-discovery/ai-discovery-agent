"""Classical ground state and linear spin waves of the authors' SpinW Hamiltonian (construction side).

SpinW_Co3Sn2S2_DMI.m: R-3m, a = 5.357 A, c = 13.127 A, Co (S = 1) at (1/2, 0, 1/2); gencoupling('maxDistance', 7)
numbers bonds by length (tolerance 1e-3 A): bond 1 = 2.6785 A (in-plane NN), bond 2 = 4.6393 A (in-plane second
neighbour), bond 3 = 4.6409 A (interlayer). Couplings: J1 = -56 (bond 1), J1dmi (bond 1), Jc1 = -3 (bond 3), with
E = sum over bonds (each once) of S_i^T J S_j, Cartesian x || a, z || c.

The task fixes the generator of bond 1: the J1dmi matrix belongs to the bond from Co(0, 1/2, 1/2) to
Co(1/2, 1/2, 1/2) (vector +a/2). Symmetry images carry R J R^T (R the Cartesian rotation of the space-group
operation) and a reversed bond carries the transpose.
"""
import itertools
import numpy as np
from scipy.optimize import minimize

A_LAT, C_LAT = 5.357, 13.127
LAT = np.array([[A_LAT, 0, 0], [-A_LAT/2, A_LAT*np.sqrt(3)/2, 0], [0, 0, C_LAT]])   # rows a, b, c (Cartesian)
RECIP = 2*np.pi*np.linalg.inv(LAT).T
OPS_TXT = ['x,y,z', '-y,x-y,z', '-x+y,-x,z', 'y,x,-z', 'x-y,-y,-z', '-x,-x+y,-z',
           '-x,-y,-z', 'y,-x+y,-z', 'x-y,x,-z', '-y,-x,z', '-x+y,y,z', 'x,x-y,z']
CENTRING = [np.zeros(3), np.array([2, 1, 1])/3, np.array([1, 2, 2])/3]
J1, D, DZ_RATIO, JC = -56.0, 15.0, 0.4, -3.0


def parse(op):
    W = np.zeros((3, 3))
    for row, expr in enumerate(op.split(',')):
        for term in expr.replace('-', '+-').split('+'):
            if term:
                W[row, 'xyz'.index(term.lstrip('-'))] = -1 if term.startswith('-') else 1
    return W


OPS = [(parse(s), t) for s in OPS_TXT for t in CENTRING]


def wrap(f):
    return f - np.floor(f + 1e-9)


def co_sites():
    sites = []
    for W, t in OPS:
        p = wrap(W @ np.array([0.5, 0, 0.5]) + t)
        if not any(np.allclose(p, q, atol=1e-6) for q in sites):
            sites.append(p)
    return np.array(sites)


CO = co_sites()
GEN = (3, 6, (0, 0, 0))     # Co(0,1/2,1/2) -> Co(1/2,1/2,1/2)
assert np.allclose(CO[3], [0, 0.5, 0.5]) and np.allclose(CO[6], [0.5, 0.5, 0.5])


def cart(f):
    return np.asarray(f) @ LAT


def cart_rotation(W):
    return LAT.T @ W @ np.linalg.inv(LAT.T)


def site_index(f):
    d = np.abs(((wrap(f) - CO) + 0.5) % 1 - 0.5).sum(1)
    k = int(np.argmin(d))
    assert d[k] < 1e-5
    return k


def image(W, t, bond):
    i, j, n = bond
    ri, rj = W @ CO[i] + t, W @ (CO[j] + np.array(n)) + t
    i2, j2 = site_index(ri), site_index(rj)
    shift = np.round(ri - CO[i2])
    return (i2, j2, tuple(int(x) for x in np.round(rj - shift - CO[j2])))


def reverse(b):
    i, j, n = b
    return (j, i, tuple(-x for x in n))


def vec(b):
    i, j, n = b
    return cart(CO[j] + np.array(n) - CO[i])


def shells(maxd=7.0):
    out = {}
    for i, j in itertools.product(range(len(CO)), repeat=2):
        for n in itertools.product(range(-2, 3), repeat=3):
            L = np.linalg.norm(cart(CO[j] + np.array(n) - CO[i]))
            if 1e-6 < L < maxd:
                out.setdefault(round(L, 4), []).append((i, j, n))
    keys = sorted(out)
    groups, cur = [], [keys[0]]
    for L in keys[1:]:                       # SpinW-style grouping with 1e-3 A tolerance
        if L - cur[-1] < 1e-3:
            cur.append(L)
        else:
            groups.append(cur)
            cur = [L]
    groups.append(cur)
    result = []
    for g in groups:
        bonds = []
        for L in g:
            for b in out[L]:
                if b not in bonds and reverse(b) not in bonds:
                    bonds.append(b)
        result.append((g, bonds))
    return result


SHELLS = shells()
NN = SHELLS[0][1]
INTERLAYER = SHELLS[2][1]
assert abs(SHELLS[0][0][0] - 2.6785) < 1e-3 and abs(SHELLS[2][0][0] - 4.6409) < 1e-3


def dm_authors(d=D, ratio=DZ_RATIO):
    return np.array([[0, 0, d], [0, 0, 0], [-d, 0, 0]]) - ratio*np.array([[0, d, 0], [-d, 0, 0], [0, 0, 0]])


def dm_vector(M):
    A = (M - M.T)/2
    return np.array([A[1, 2], -A[0, 2], A[0, 1]])


def propagate(M, gen=GEN):
    mats = {}
    for W, t in OPS:
        b2 = image(W, t, gen)
        R = cart_rotation(W)
        M2 = R @ M @ R.T
        if b2 in mats:
            assert np.allclose(mats[b2], M2, atol=1e-9), 'matrix not compatible with bond symmetry'
        elif reverse(b2) in mats:
            assert np.allclose(mats[reverse(b2)], M2.T, atol=1e-9), 'matrix not compatible with bond symmetry'
        else:
            mats[b2] = M2
    return mats


def terms(dm=None, gen=GEN, j1=J1, jc=JC, copy=False):
    """List of (i, j, d_cart, M) with E = sum S_i^T M S_j over bonds, each once (9-site hexagonal cell)."""
    mats = propagate(dm, gen) if (dm is not None and not copy) else {}
    out = []
    for b in NN:
        M = j1*np.eye(3)
        if dm is not None:
            if copy:                                 # same matrix on every bond, oriented along +x / +y
                v = vec(b)
                bb = b if (v[0] > 1e-6 or (abs(v[0]) < 1e-6 and v[1] > 0)) else reverse(b)
                out.append((bb[0], bb[1], vec(bb), M + dm))
                continue
            M = M + (mats[b] if b in mats else mats[reverse(b)].T)
        out.append((b[0], b[1], vec(b), M))
    for b in INTERLAYER:
        out.append((b[0], b[1], vec(b), jc*np.eye(3)))
    return out


def energy(S, T):
    return sum(S[i] @ M @ S[j] for i, j, d, M in T)/len(S)


def angles_to_spins(x):
    n = len(x)//2
    th, ph = x[:n], x[n:]
    return np.c_[np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)]


def ground_state(T, nrun=16, seed=1):
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(nrun):
        x0 = np.r_[rng.uniform(0, np.pi, 9), rng.uniform(0, 2*np.pi, 9)]
        r = minimize(lambda x: energy(angles_to_spins(x), T), x0, method='BFGS', options={'gtol': 1e-11})
        if best is None or r.fun < best.fun - 1e-10:
            best = r
    S = angles_to_spins(best.x)
    if S[:, 2].sum() < 0:          # bilinear energy: time reversal is a symmetry
        S = -S
    return best.fun, S


def frame(s):
    s = s/np.linalg.norm(s)
    t = np.array([0, 0, 1.0]) if abs(s[2]) < 0.9 else np.array([1.0, 0, 0])
    e1 = np.cross(t, s)
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(s, e1), s


def lswt(S_dirs, T, kvecs, spin=1.0):
    """Linear spin-wave energies (Holstein-Primakoff, paraunitary diagonalisation)."""
    N = len(S_dirs)
    fr = [frame(s) for s in S_dirs]
    u = np.array([f[0] + 1j*f[1] for f in fr])
    v = np.array([f[2] for f in fr])

    def amat(k):
        A = np.zeros((N, N), complex)
        for i, j, d, M in T:
            ph = np.exp(1j*np.dot(k, d))
            A[i, j] += spin/2*(u[i] @ M @ np.conj(u[j]))*ph
            A[j, i] += spin/2*(np.conj(u[i]) @ M @ u[j])*np.conj(ph)
        return A
    eta = np.diag(np.r_[np.ones(N), -np.ones(N)])
    out = []
    for k in kvecs:
        k = np.asarray(k, float)
        Dm = np.zeros((N, N), complex)
        C = np.zeros(N)
        for i, j, d, M in T:
            ph = np.exp(1j*np.dot(k, d))
            cs = spin/2*(u[i] @ M @ u[j])
            Dm[i, j] += cs*ph
            Dm[j, i] += cs*np.conj(ph)
            C[i] += spin*(v[i] @ M @ v[j])
            C[j] += spin*(v[i] @ M @ v[j])
        H = np.block([[amat(k) - np.diag(C), Dm], [Dm.conj().T, (amat(-k) - np.diag(C)).T]])
        w = np.linalg.eigvals(eta @ H)
        out.append(np.sort(w.real)[N:])
    return np.array(out)


PRIM = {0: 0, 1: 0, 2: 0, 3: 1, 4: 1, 5: 1, 6: 2, 7: 2, 8: 2}


def gap_and_tilt(T, spin=1.0, nrun=16):
    E, S = ground_state(T, nrun=nrun)
    tilt = np.degrees(np.arccos(np.clip(np.abs(S[:, 2]), 0, 1)))
    k0 = all(np.allclose(S[a], S[b], atol=1e-4) for a, b in [(0, 1), (0, 2), (3, 4), (3, 5), (6, 7), (6, 8)])
    if k0:
        w = lswt(S[[0, 3, 6]], [(PRIM[i], PRIM[j], d, M/3) for i, j, d, M in T], [np.zeros(3)], spin)[0]
    else:
        w = lswt(S, T, [np.zeros(3)], spin)[0]
    return {'energy_meV_per_spin': float(E), 'tilt_deg': float(np.mean(tilt)), 'tilt_spread_deg': float(np.ptp(tilt)),
            'gap_meV': float(max(np.min(w), 0.0)), 'k0_primitive': bool(k0), 'mz': float(abs(S[:, 2].mean())),
            'spins': S.round(6).tolist()}


def landau_lifshitz_k0(S9, T):
    """Independent check: k = 0 frequencies of the linearised classical spin dynamics (S = 1)."""
    Sp = S9[[0, 3, 6]]
    Tp = [(PRIM[i], PRIM[j], d, M/3) for i, j, d, M in T]
    fr = [frame(s) for s in Sp]

    def spins(x):
        out = []
        for n, (e1, e2, s) in enumerate(fr):
            w = s + x[2*n]*e1 + x[2*n+1]*e2
            out.append(w/np.linalg.norm(w))
        return np.array(out)

    def rhs(x):
        S = spins(x)
        h = np.zeros_like(S)
        for i, j, d, M in Tp:
            h[i] -= M @ S[j]
            h[j] -= M.T @ S[i]
        dS = np.cross(S, h)
        return np.array([[dS[n] @ fr[n][0], dS[n] @ fr[n][1]] for n in range(3)]).ravel()
    eps = 1e-6
    J = np.array([(rhs(eps*e) - rhs(-eps*e))/(2*eps) for e in np.eye(6)]).T
    return np.sort(np.abs(np.linalg.eigvals(J).imag))
