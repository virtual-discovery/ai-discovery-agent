"""Readers and reference estimators for the inmnti2-weyl-census construction (numpy only).

The archive is read from DATA/x/<member> as extracted by fetch.py (DATA is /data inside run_in_image.sh).
"""
import os
import re
from pathlib import Path
import numpy as np

DATA = Path(os.environ.get('INMNTI2_DATA', '/data'))
X = DATA/'x'
A_PRIM = 3.1335920710                      # Angstrom; CELL_PARAMETERS entries of every electronic-structure run
ALAT = A_PRIM * np.sqrt(2)                 # pw.x alat for ibrav=0 with CELL_PARAMETERS angstrom (|a1|)
TPIBA = 2 * np.pi / ALAT                   # 1/Angstrom per 2pi/alat
EF = {'hubbard': 16.3261, 'standard': 15.6841}   # scf Fermi energies (eV), asserted below
OFFSET = {'hubbard': 36, 'standard': 34}         # DFT band index = Wannier band index + OFFSET (asserted)
NOCC_DFT = 52                                    # electrons per cell -> crossing between DFT bands 52 and 53


def read_hr(path):
    with open(path) as f:
        f.readline()
        nw = int(f.readline()); nr = int(f.readline())
        deg = []
        while len(deg) < nr:
            deg += [int(x) for x in f.readline().split()]
        data = np.loadtxt(f)
    R = data[::nw * nw, :3].astype(int)
    m = data[:, 3].astype(int) - 1
    n = data[:, 4].astype(int) - 1
    H = np.zeros((nr, nw, nw), complex)
    H[np.repeat(np.arange(nr), nw * nw), m, n] = data[:, 5] + 1j * data[:, 6]
    return R, np.array(deg), H


def qe_bands(path):
    """k-points (2pi/alat Cartesian, as printed) and eigenvalues (eV) of the last band listing of a pw.x output."""
    txt = Path(path).read_text()
    for marker in ('End of band structure calculation', 'End of self-consistent calculation'):
        if marker in txt:
            txt = txt[txt.rfind(marker):]
            break
    blocks = re.split(r'\n\s+k =', txt)[1:]
    ks, es = [], []
    for b in blocks:
        head, rest = b.split('\n', 1)
        ks.append([float(x) for x in re.findall(r'-?\d+\.\d+', head.replace('-', ' -'))[:3]])
        body = re.split(r'occupation numbers|the Fermi energy|highest occupied|Writing', rest)[0]
        es.append([float(x) for x in re.findall(r'-?\d+\.\d+', body.replace('-', ' -'))])
    n = min(len(e) for e in es)
    return np.array(ks), np.array([e[:n] for e in es])


def qe_header(path):
    """Selected lines of a pw.x output: program, cell, electrons, k-points, smearing, Fermi energy, timing."""
    keys = ('Program PWSCF', 'lattice parameter', 'unit-cell volume', 'number of electrons', 'Kohn-Sham states',
            'number of k points', 'the Fermi energy', 'highest occupied', 'kinetic-energy cutoff', 'charge density cutoff',
            'Exchange-correlation', 'No symmetry found', 'Sym. Ops.', 'Noncollinear calculation', 'JOB DONE',
            'Hubbard', 'a(1) =', 'a(2) =', 'a(3) =')
    return [line.rstrip() for line in Path(path).read_text().splitlines() if any(k in line for k in keys)][:60]


def tpiba_b_path(nscf_in_text):
    """Exact Cartesian k (1/Angstrom) of a K_POINTS tpiba_b list as pw.x generates it, and the corner list."""
    lines = nscf_in_text.split('K_POINTS')[1].split('\n')
    npts = int(lines[1])
    pts = [(np.array([float(x) for x in lines[2 + i].split()[:3]]), int(lines[2 + i].split()[3])) for i in range(npts)]
    ks = []
    for i, (p, n) in enumerate(pts):
        if i + 1 < len(pts):
            q = pts[i + 1][0]
            ks += [p + (q - p) * j / n for j in range(n)]
        else:
            ks.append(p)
    return np.array(ks) * TPIBA, pts


def xyz_segments(pts):
    """Slices of the three lines (x, y, z) of an nscfXYZ path with corners [P1 n1, P2 1, P3 n3, P4 1, P5 n5, P6 1]."""
    n1, n3, n5 = pts[0][1], pts[2][1], pts[4][1]
    return [('x', slice(0, n1 + 1)), ('y', slice(n1 + 1, n1 + n3 + 2)), ('z', slice(n1 + n3 + 2, n1 + n3 + n5 + 3))]


def crossing_on_line(s, lo, hi, i0, window=0.03, skip=0.004):
    """Crossing of two straight bands near index i0 of a line (s: arc length in 1/A; lo, hi: the two bands in eV).

    The four half-branches within `window` of the gap minimum (excluding |s - s0| < skip) are fitted linearly,
    joined into the two straight bands (upper-left with lower-right and vice versa) and intersected.
    Returns (s_cross, E_cross, (slope_a, slope_b))."""
    left = [i for i in range(len(s)) if skip < s[i0] - s[i] < window]
    right = [i for i in range(len(s)) if skip < s[i] - s[i0] < window]
    fit = lambda idx, y: np.polyfit(s[idx], y[idx], 1)
    a = (fit(left, hi) + fit(right, lo)) / 2
    b = (fit(left, lo) + fit(right, hi)) / 2
    sc = (b[1] - a[1]) / (a[0] - b[0])
    return float(sc), float(np.polyval(a, sc)), (float(a[0]), float(b[0]))


def symmetry_ops(scf_out_text):
    """Point operations of a pw.x run: list of (cartesian 3x3 matrix, time_reversal flag)."""
    ops = []
    blocks = re.split(r'\n\s+isym =\s*\d+', scf_out_text.split('Sym. Ops.')[1].split('Cartesian axes')[0])[1:]
    for b in blocks:
        tr = re.search(r'Time Reversal\s+(\d)', b)
        cart = re.split(r'cart\.\s+s\(', b)[1]
        rows = re.findall(r'\(\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*\)', cart)[:3]
        ops.append((np.array(rows, float), int(tr.group(1)) if tr else 0))
    return ops


def kspace_images(k, ops, nonmagnetic=False):
    """Images of Cartesian k (and chirality factor) under the operations; antiunitary ones map k -> -R k.
    For a non-magnetic run time reversal itself is added (k -> -k, chirality unchanged)."""
    out = []
    for R, t in ops:
        sign = -1 if t else 1
        out.append((sign * R @ k, round(np.linalg.det(R))))
        if nonmagnetic:
            out.append((-R @ k, round(np.linalg.det(R))))
    return out


def orbits(kcart, chi, ops, B, tol=0.03, nonmagnetic=False):
    """Partition nodes into symmetry orbits (union-find on images matched modulo reciprocal lattice vectors)."""
    kcart = np.asarray(kcart); n = len(kcart)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    Binv = np.linalg.inv(B)
    for i in range(n):
        for img, det in kspace_images(kcart[i], ops, nonmagnetic):
            d = (kcart - img) @ Binv
            d = (d - np.round(d)) @ B
            j = int(np.argmin(np.linalg.norm(d, axis=1)))
            if np.linalg.norm(d[j]) < tol:
                assert chi[j] == chi[i] * det, 'chirality inconsistent with symmetry'
                parent[find(i)] = find(j)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return sorted(groups.values(), key=len)
