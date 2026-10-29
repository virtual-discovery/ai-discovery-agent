"""Linear spin-wave theory for NiAs-type MnTe (two Mn sublattices, collinear in-plane Neel order).

H = sum_{pairs} J_ij S_i.S_j + D sum_i (S_i^z)^2   (each pair once, S = 5/2, z || c, Cartesian x || a).
Exchange shells are numbered by distance with sym.bond_classes(); shells 10 and 11 are the two
symmetry-inequivalent 9.831 A same-sublattice classes (altermagnetic pair).
"""
import numpy as np
import sym

S = 2.5
BONDS, CLASSES = sym.bond_classes(10.5)
_order = sorted(set(CLASSES), key=lambda c: min(b[3] for b, cc in zip(BONDS, CLASSES) if cc == c))
SHELL_OF_CLASS = {c: n + 1 for n, c in enumerate(_order)}
SHELL_BONDS = {}
for b, c in zip(BONDS, CLASSES):
    i, j, R, dist = b
    if i > j or (i == j and R <= (0, 0, 0)):
        continue                      # each undirected pair once
    SHELL_BONDS.setdefault(SHELL_OF_CLASS[c], []).append((i, j, np.array(R) + sym.MN[j] - sym.MN[i]))
SHELL_ARR = {s: (np.array([b[0] for b in v]), np.array([b[1] for b in v]), np.array([b[2] for b in v]))
             for s, v in SHELL_BONDS.items()}
SHELL_DIST = {SHELL_OF_CLASS[c]: min(b[3] for b, cc in zip(BONDS, CLASSES) if cc == c) for c in set(CLASSES)}


def frames(neel_angle_deg=30.0, sign=+1):
    """Spin directions e_s and complex transverse vectors c_s = u_s - i w_s for sublattices A (z=0), B (z=1/2).

    The Neel vector n lies in the ab plane at neel_angle_deg from Cartesian x (x || a; a* is at 30 deg).
    Sublattice A spin = sign * n."""
    t = np.radians(neel_angle_deg)
    n = np.array([np.cos(t), np.sin(t), 0.0])
    zc = np.array([0, 0, 1.0])
    out = []
    for s in (sign, -sign):
        e = s * n
        u = np.cross(zc, e)
        w = np.cross(e, u)
        out.append((e, u - 1j * w))
    return out


def _blocks(k, J, D, fr):
    A = np.zeros((2, 2), complex)
    B = np.zeros((2, 2), complex)
    k = np.asarray(k, float)
    for shell, Jv in J.items():
        if shell == 'extra' or Jv == 0:
            continue
        I, Jx, dl = SHELL_ARR[shell]
        ph = np.exp(2j * np.pi * dl @ k)
        for (i, j) in ((0, 0), (0, 1), (1, 0), (1, 1)):
            m = (I == i) & (Jx == j)
            if not m.any():
                continue
            (ei, ci), (ej, cj) = fr[i], fr[j]
            K = Jv * S / 2
            p = ph[m].sum()
            A[i, j] += K * np.dot(np.conj(ci), cj) * p
            A[j, i] += K * np.dot(ci, np.conj(cj)) * np.conj(p)
            b = K * np.dot(np.conj(ci), np.conj(cj))
            B[i, j] += b * p
            B[j, i] += b * np.conj(p)
            A[i, i] += -Jv * S * np.dot(ei, ej) * m.sum()
            A[j, j] += -Jv * S * np.dot(ei, ej) * m.sum()
    for (i, j, dl, Jv) in J.get('extra', []):        # explicit extra bonds (i, j, delta_frac, J), each pair once
        (ei, ci), (ej, cj) = fr[i], fr[j]
        K = Jv * S / 2
        p = np.exp(2j * np.pi * np.dot(dl, k))
        A[i, j] += K * np.dot(np.conj(ci), cj) * p
        A[j, i] += K * np.dot(ci, np.conj(cj)) * np.conj(p)
        b = K * np.dot(np.conj(ci), np.conj(cj))
        B[i, j] += b * p
        B[j, i] += b * np.conj(p)
        A[i, i] += -Jv * S * np.dot(ei, ej)
        A[j, j] += -Jv * S * np.dot(ei, ej)
    for s in range(2):
        cz = fr[s][1][2]
        A[s, s] += S * D * abs(cz) ** 2
        B[s, s] += S * D * np.conj(cz) ** 2
    return A, B


def bdg(k, J, D, fr):
    A, B = _blocks(k, J, D, fr)
    Am, _ = _blocks(-np.asarray(k, float), J, D, fr)
    M = np.block([[A, B], [B.conj().T, Am.T]])
    return (M + M.conj().T) / 2


def modes(k, J, D, fr, eps=1e-9):
    """Colpa diagonalisation: energies (2,), ascending, and T with X = T Gamma."""
    M = bdg(k, J, D, fr)
    g = np.diag([1, 1, -1, -1.0])
    K = np.linalg.cholesky(M + eps * np.eye(4)).conj().T
    L, U = np.linalg.eigh(K @ g @ K.conj().T)
    pos, neg = np.where(L > 0)[0], np.where(L <= 0)[0]
    pos = pos[np.argsort(L[pos])]
    neg = neg[np.argsort(-L[neg])]
    idx = np.r_[pos, neg]
    L, U = L[idx], U[:, idx]
    E = np.abs(L)
    T = np.linalg.solve(K, U @ np.diag(np.sqrt(E)))
    return E[:2], T


def correlations(k, J, D, fr):
    """Energies (2,) and S^{ab}_n (2,3,3) for magnon creation at momentum k (r.l.u.)."""
    E, T = modes(k, J, D, fr)
    W = np.zeros((3, 4), complex)
    for s in range(2):
        ph = np.exp(-2j * np.pi * np.dot(k, sym.MN[s]))
        c = fr[s][1]
        W[:, s] = np.sqrt(S / 2) * c * ph
        W[:, 2 + s] = np.sqrt(S / 2) * np.conj(c) * ph
    V = W @ T[:, :2]
    Sab = np.einsum('an,bn->nab', np.conj(V), V)
    return E, Sab


BMAT = 2 * np.pi * np.linalg.inv(sym.CELL).T


def qcart(hkl):
    return BMAT @ np.asarray(hkl, float)
