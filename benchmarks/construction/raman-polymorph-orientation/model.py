"""Physics utilities for the raman-polymorph-orientation construction (numpy/scipy only).

Database conventions (Bagheri & Komsa, Sci. Data 10, 80 (2023), Eqs. 4-6):
  * eigenvectors e_j(k) are the normalized eigenvectors of the mass-weighted dynamical matrix (phonopy);
  * the Raman tensor of mode j is R_j = (V_c / 4 pi) d chi / d xi_j, with xi_j the normal coordinate along
    e_j(k)/sqrt(M_k); the Stokes cross-section is proportional to omega_S^4 (n+1)/omega_j |e_s . R_j . e_i|^2;
  * Born charges are archived in phonopy BORN form: one 3x3 tensor per symmetry-independent atom
    (first index = electric field, second = atomic displacement), in order of first appearance.
Longitudinal-optical (LO) mixing at q -> 0 along a direction q is the standard non-analytic term
  D_NA(k a, k' b) = (4 pi e^2 / V_c) (q.Z_k)_a (q.Z_k')_b / (q . eps_inf . q) / sqrt(M_k M_k')
added to the analytic (TO) dynamical matrix D = sum_j omega_j^2 e_j e_j^T. The Raman tensors of the mixed
modes are projections of the archived ones on the new eigenvectors (deformation-potential term; the
electro-optic contribution is not available from the database and is neglected).
"""
import itertools
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
VASP2THZ = 15.633302          # sqrt(eV/A^2/amu) -> THz (phonopy)
THZ2CM = 33.35641
E2 = 14.399645                # e^2/(4 pi eps0) in eV A
HC_K = 1.4387769              # cm K
MASSES = {9: 18.998403163, 50: 118.71}   # standard atomic weights used by phonopy
POLYMORPHS = {'mp-2843': '3216', 'mp-7456': '2564', 'mp-7457': '2390'}


def _arr(x):
    if isinstance(x, dict) and '__ndarray__' in x:
        shape, _, data = x['__ndarray__']
        return np.array(data).reshape(shape)
    return np.array(x)


def load_pool():
    return json.loads((HERE/'pool.json').read_text())


def record(rec):
    """Parse one archived ASE-database record into numpy arrays (no change of content)."""
    d, kv = rec['data'], rec['key_value_pairs']
    e = {'mpid': kv['mpid'], 'spacegroup': kv['spacegroup'], 'pointgroup': str(d['pointgroup']),
         'cell': _arr(rec['cell']['array']), 'numbers': _arr(rec['numbers']).astype(int),
         'positions': _arr(rec['positions']),
         'freq': np.array([float(v) for v in kv['frequencies_cm'].strip('[]').split()]),
         'eigvec': np.array(d['eigenvec'], float)}
    nmode = len(e['freq'])
    tensors = list(d['raman_tensors']) + [None] * (nmode - len(d['raman_tensors']))
    e['raman_tensors'] = [None if t is None else np.array(t, float) for t in tensors]
    e['raman_active'] = _arr(d['Ramanactive']).astype(int)
    e['ir_active'] = _arr(d['IRactive']).astype(int)
    e['born_rows'] = (np.array([[float(v) for v in row.split()] for row in d['born']]).reshape(-1, 3, 3)
                      if isinstance(d['born'], list) else None)
    e['eps_inf'] = np.array([float(v) for v in d['diel']]).reshape(3, 3) if isinstance(d['diel'], list) else None
    e['labels'] = eval(kv['IRlabels']) if kv.get('IRlabels') else None   # archived python-literal strings
    e['bands'] = eval(kv['IRbands']) if kv.get('IRbands') else None
    e['flags'] = {k: kv.get(k) for k in ('mp_e_above_hull', 'negative_freq_Gamma', 'negative_freq_path', 'bandgap')}
    e['masses'] = np.array([MASSES[z] for z in e['numbers']])
    e['volume'] = abs(np.linalg.det(e['cell']))
    e['frac'] = e['positions'] @ np.linalg.inv(e['cell'])
    e['Rsym'] = np.array([np.zeros((3, 3)) if t is None else 0.5 * (t + t.T) for t in e['raman_tensors']])
    e['E'] = e['eigvec'].reshape(nmode, -1)
    e['Z'] = expand_born(e) if e['born_rows'] is not None else None
    return e


def load(mpid):
    return record(load_pool()[POLYMORPHS[mpid]])


def space_group_ops(e, tol=1e-3):
    """All (W, t, perm) with W integer in fractional coordinates mapping the structure onto itself."""
    A, frac, num = e['cell'], e['frac'], e['numbers']
    G = A @ A.T
    ops = []
    for flat in itertools.product((-1, 0, 1), repeat=9):
        W = np.array(flat).reshape(3, 3)
        if abs(round(np.linalg.det(W))) != 1 or not np.allclose(W.T @ G @ W, G, atol=1e-3 * np.abs(G).max()):
            continue
        Wx = frac @ W.T
        for j in np.where(num == num[0])[0]:
            t = frac[j] - Wx[0]
            y = Wx + t
            perm = []
            for k in range(len(frac)):
                dd = y[k] - frac
                dd -= np.round(dd)
                m = np.where((np.linalg.norm(dd @ A, axis=1) < 10 * tol) & (num == num[k]))[0]
                if len(m) != 1:
                    break
                perm.append(int(m[0]))
            if len(perm) == len(frac):
                ops.append((W, t - np.round(t), np.array(perm)))
    return ops


def cart(e, W):
    A = e['cell']
    return A.T @ W @ np.linalg.inv(A.T)


def expand_born(e, rotate=True):
    """Born charges of all atoms from the archived rows (independent atoms in order of first appearance)."""
    ops = space_group_ops(e)
    n = len(e['numbers'])
    rep = -np.ones(n, int)
    for k in range(n):
        if rep[k] < 0:
            for W, t, p in ops:
                if rep[p[k]] < 0:
                    rep[p[k]] = k
    reps = sorted(set(rep.tolist()))
    rows = e['born_rows']
    assert len(reps) == len(rows), (len(reps), len(rows))
    Z = np.zeros((n, 3, 3))
    for k in range(n):
        i = rep[k]
        Zi = rows[reps.index(i)]
        if not rotate:
            Z[k] = Zi
            continue
        for W, t, p in ops:
            if p[i] == k:
                R = cart(e, W)
                Z[k] = R @ Zi @ R.T
                break
    if rotate:
        assert np.abs(Z.sum(0)).max() < 1e-3, 'acoustic sum rule violated after expansion'
        for idx, i in enumerate(reps):        # site-symmetry invariance of each archived row
            for W, t, p in ops:
                if p[i] == i:
                    R = cart(e, W)
                    assert np.abs(R @ rows[idx] @ R.T - rows[idx]).max() < 1e-3
    return Z


def modes(e, qhat=None, Z=None, eps=None):
    """Optical-mode frequencies (cm-1, database scale) and Raman tensors for phonon direction qhat
    (Cartesian frame of the database cell). qhat=None gives the archived (TO) modes."""
    f0, E = e['freq'], e['E']
    Z = e['Z'] if Z is None else Z
    if qhat is None or Z is None:
        idx = np.argsort(f0)[3:]
        return f0[idx], e['Rsym'][idx]
    eps = e['eps_inf'] if eps is None else eps
    w = f0 / THZ2CM / VASP2THZ
    D = E.T @ np.diag(np.sign(w) * w ** 2) @ E
    q = np.asarray(qhat, float) / np.linalg.norm(qhat)
    zq = np.einsum('a,kab->kb', q, Z).reshape(-1) / np.sqrt(np.repeat(e['masses'], 3))
    D = D + 4 * np.pi * E2 / e['volume'] * np.outer(zq, zq) / (q @ eps @ q)
    val, vec = np.linalg.eigh(D)
    f = np.sign(val) * np.sqrt(np.abs(val)) * VASP2THZ * THZ2CM
    R = np.einsum('jm,jab->mab', E @ vec, e['Rsym'])
    return f[3:], R[3:]


def strongest_mode(e):
    """Index (in the archived mode list) of the Raman mode with the largest sum R_ab^2, degenerate sets summed."""
    best, val = None, -1
    for band in e['bands']:
        s = sum(float((e['raman_tensors'][k] ** 2).sum()) for k in band if e['raman_tensors'][k] is not None)
        if s > val:
            best, val = band, s
    return best, val


def powder_invariants(R):
    """Orientation-averaged backscattering intensities (parallel, crossed) of a symmetric tensor."""
    a = np.trace(R) / 3
    g2 = 0.5 * ((R[0, 0] - R[1, 1]) ** 2 + (R[0, 0] - R[2, 2]) ** 2 + (R[1, 1] - R[2, 2]) ** 2
                + 6 * (R[0, 1] ** 2 + R[0, 2] ** 2 + R[1, 2] ** 2))
    return (45 * a * a + 4 * g2) / 45, g2 / 15


def euler(phi, theta, psi):
    """Q = Rz(phi) Ry(theta) Rz(psi) (degrees): crystal (database Cartesian frame) -> lab."""
    def rz(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    def ry(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    return rz(np.radians(phi)) @ ry(np.radians(theta)) @ rz(np.radians(psi))


def axis_error_deg(c, c_true):
    """Angle between axes, modulo sign and the z-reflection that face-A backscattering cannot resolve."""
    c = np.asarray(c, float) / np.linalg.norm(c)
    ct = np.asarray(c_true, float)
    return min(np.degrees(np.arccos(min(1.0, abs(float(c @ x))))) for x in (ct, ct * np.array([1, 1, -1])))


FACE_B_DEG = 30.0


def face_b():
    b = np.radians(FACE_B_DEG)
    n = np.array([np.cos(b), np.sin(b), 0.0])
    return n, np.array([-n[1], n[0], 0.0]), np.array([0.0, 0.0, 1.0])


def targets(e, Q, scale, Z=None, eps=None):
    """Face-B predictions for polymorph e with orientation Q (crystal -> lab) and frequency scale."""
    n, t, z = face_b()
    band, _ = strongest_mode(e)
    Rl = [Q @ e['Rsym'][k] @ Q.T for k in band]
    ratio = sum((t @ R @ t) ** 2 for R in Rl) / sum((z @ R @ z) ** 2 for R in Rl)
    f, _ = modes(e, Q.T @ n, Z=Z, eps=eps)
    return {'ratio_tt_over_zz': float(ratio), 'highest_optical_cm1': float(scale * f.max()),
            'c_axis_lab': (Q @ (e['cell'][2] / np.linalg.norm(e['cell'][2]))).tolist()}
