"""Numpy utilities for the construction of polar-phonons-1d-loto-collapse.

- parse_dyn: Quantum ESPRESSO ph.x dynamical-matrix files (fildyn format).
- Fourier interpolation along the wire axis (QE q2r/matdyn conventions, Wigner-Seitz weights along z),
  with an optional long-range term: the 1D dipole-dipole kernel of the authors' modified matdyn
  (filled dielectric cylinder of radius t, Rivano et al. 2023) or nothing.
- The 3D q->0 nonanalytic term of a periodic supercell (standard QE).

Units: Rydberg atomic units (e^2 = 2, masses in units of m_e/2 as in the dyn files), Ry/bohr^2 for the
force-constant matrices, cm^-1 for frequencies. The 1D kernel reproduces the authors' archived
interpolated dispersions of BN(5,5), BN(6,6), the BN chain and the GaAs nanowire (make_task.py asserts it).
"""
import glob
import os
import re
import numpy as np
from scipy.special import i0, i1, k0, k1

RY_TO_CMM1 = 109737.31570111268
E2 = 2.0


def parse_dyn(path):
    lines = open(path).read().split('\n')
    assert 'Dynamical matrix file' in lines[0], path
    h = lines[2].split()
    ntyp, nat, ibrav = int(h[0]), int(h[1]), int(h[2])
    celldm = [float(x) for x in h[3:9]]
    i = 3
    lattice = None
    if ibrav == 0:
        i += 1
        lattice = np.array([[float(x) for x in lines[i + k].split()] for k in range(3)])
        i += 3
    species = {}
    for _ in range(ntyp):
        m = re.match(r"\s*(\d+)\s+'([^']*)'\s+([\d.Ee+-]+)", lines[i])
        species[int(m.group(1))] = (m.group(2).strip(), float(m.group(3)))
        i += 1
    ityp, tau = [], []
    for _ in range(nat):
        p = lines[i].split()
        ityp.append(int(p[1]))
        tau.append([float(x) for x in p[2:5]])
        i += 1
    out = {'nat': nat, 'ibrav': ibrav, 'celldm': celldm, 'lattice': lattice, 'tau_alat': np.array(tau),
           'mass': np.array([species[t][1] for t in ityp]), 'symbols': [species[t][0] for t in ityp], 'blocks': []}
    n = len(lines)
    while i < n:
        line = lines[i]
        if 'Dynamical  Matrix in cartesian axes' in line:
            q = [float(x) for x in lines[i + 2].split('(')[1].split(')')[0].split()]
            i += 4
            D = np.zeros((nat, 3, nat, 3), complex)
            for _ in range(nat * nat):
                a, b = [int(x) for x in lines[i].split()]
                for r in range(3):
                    v = [float(x) for x in lines[i + 1 + r].split()]
                    D[a - 1, r, b - 1, :] = np.array(v[0::2]) + 1j * np.array(v[1::2])
                i += 4
            out['blocks'].append({'q': np.array(q), 'D': D})
            continue
        if 'Dielectric Tensor' in line:
            out['epsilon'] = np.array([[float(x) for x in lines[i + 2 + r].split()] for r in range(3)])
            i += 5
            continue
        if 'Effective Charges E-U' in line:
            Z = np.zeros((nat, 3, 3))
            i += 2
            for a in range(nat):
                i += 1
                for r in range(3):
                    Z[a, r] = [float(x) for x in lines[i].split()]
                    i += 1
            out['zeu'] = Z          # [atom][field alpha][displacement beta]
            continue
        if 'Diagonalizing the dynamical matrix' in line:
            break
        i += 1
    return out


def cell_of(d):
    alat = d['celldm'][0]
    if d['ibrav'] == 0:
        return d['lattice'] * alat
    return np.diag([alat, alat, alat * d['celldm'][2]])


def load_system(folder, grid_files=None):
    """Grid dynamical matrices of one archived folder (dyn0 lists the irreducible grid)."""
    g = parse_dyn(f'{folder}/dyn1')
    head = open(f'{folder}/dyn0').read().split()
    nq, nirr = int(head[2]), int(head[3])
    lat = cell_of(g)
    alat, c = g['celldm'][0], lat[2, 2]
    Dq = {}
    for k in range(1, nirr + 1):
        d = parse_dyn(f'{folder}/dyn{k}')
        for b in d['blocks']:
            m = int(round(b['q'][2] * 2 * np.pi / alat * c / (2 * np.pi) * nq)) % nq
            Dq[m] = b['D']
    assert len(Dq) == nq, (folder, sorted(Dq))
    return {'nq': nq, 'alat': alat, 'lat': lat, 'c': c, 'Dq': Dq, 'tau': g['tau_alat'] * alat, 'Z': g['zeu'],
            'eps': g['epsilon'], 'mass': g['mass'], 'nat': g['nat'], 'area': lat[0, 0] * lat[1, 1]}


def delta_full(x, eps):
    return 2 * i1(x) * k1(x) / (x * (i0(x) * k1(x) + eps * i1(x) * k0(x)))


def delta_expansion(x, eps):
    """Small-qt limit of delta_full (the paper's asymptotic form)."""
    return 1 - x * x / 2 * (eps * (np.log(2) - np.euler_gamma) + 0.25 - eps * np.log(x))


def na_1d(s, qz, t=None, eps1d=None, ewald=True, alph=1.0, gmax=14.0, nG=3, kernel=delta_full, scale=1.0):
    """1D long-range dynamical matrix (Ry/bohr^2) at q = qz (bohr^-1) along the axis.

    D_ai,bj = (4 pi e^2 / L) sum_G W(q+G) (q+G) Z_a,zi (q+G) Z_b,zj e^{i(q+G)(tau_a - tau_b)},
    W(g) = [1 - Delta(|g| t)] / (pi t^2 eps1d g^2) [x exp(-g^2/(4 alph)), g in 2pi/alat units, if ewald];
    defaults: t = alat/4 and eps1d = A (eps_zz - 1)/(pi t^2) (the authors' choices; A = cell cross-section).
    """
    nat, Z, tau, c, alat = s['nat'], s['Z'], s['tau'], s['c'], s['alat']
    t = alat / 4 if t is None else t
    eps1d = s['area'] / (np.pi * t ** 2) * (s['eps'][2, 2] - 1) if eps1d is None else eps1d
    fac = 4 * np.pi * E2 / c * scale
    D = np.zeros((nat, 3, nat, 3), complex)
    for m in range(-nG, nG + 1):
        g = qz + m * 2 * np.pi / c
        gt = g * alat / (2 * np.pi)
        if abs(g) < 1e-12 or gt ** 2 / alph / 4 > gmax:
            continue
        W = (1 - kernel(abs(g) * t, eps1d)) / (np.pi * t ** 2 * eps1d * g ** 2)
        if ewald:
            W *= np.exp(-gt ** 2 / alph / 4)
        zg = g * Z[:, 2, :]
        ph = np.exp(1j * g * (tau[:, None, 2] - tau[None, :, 2]))
        D += fac * W * np.einsum('ai,bj,ab->aibj', zg, zg, ph)
    return D


def na_3d(s, eps=None, omega=None):
    """3D nonanalytic term for q -> 0 along z of the periodic supercell (Ry/bohr^2)."""
    eps = s['eps'] if eps is None else eps
    omega = abs(np.linalg.det(s['lat'])) if omega is None else omega
    zq = s['Z'][:, 2, :]
    return 4 * np.pi * E2 / omega * np.einsum('ai,bj->aibj', zq, zq) / eps[2, 2]


def ws_weights(s):
    nq, c, tau, nat = s['nq'], s['c'], s['tau'], s['nat']
    L = nq * c
    ns = np.arange(-2 * nq, 2 * nq + 1)
    w = np.zeros((nat, nat, len(ns)))
    for a in range(nat):
        for b in range(nat):
            r = ns * c + tau[a, 2] - tau[b, 2]
            d0 = np.abs(r)
            dm, dp = np.abs(r - L), np.abs(r + L)
            inside = (dm >= d0 - 1e-6) & (dp >= d0 - 1e-6)
            cnt = 1 + (np.abs(dm - d0) < 1e-6) + (np.abs(dp - d0) < 1e-6)
            w[a, b] = np.where(inside, 1.0 / cnt, 0.0)
    return ns, w


class Interpolator:
    """q2r + matdyn along z with an optional 1D long-range term (lr=None or dict of na_1d options)."""

    def __init__(self, s, lr=None):
        self.s, self.lr = s, lr
        nq, c = s['nq'], s['c']
        Dsr = {}
        for m, D in s['Dq'].items():
            qz = m * 2 * np.pi / (nq * c)
            if qz > np.pi / c + 1e-9:
                qz -= 2 * np.pi / c
            Dsr[m] = D - (na_1d(s, qz, **lr) if lr is not None else 0)
        self.C = [sum(Dsr[m] * np.exp(2j * np.pi * m * n / nq) for m in range(nq)) / nq for n in range(nq)]
        self.ns, self.w = ws_weights(s)

    def D(self, qz):
        s = self.s
        D = np.zeros((s['nat'], 3, s['nat'], 3), complex)
        for i, n in enumerate(self.ns):
            if not self.w[:, :, i].any():
                continue
            D += np.einsum('ab,aibj->aibj', self.w[:, :, i], self.C[n % s['nq']]) * np.exp(-1j * qz * n * s['c'])
        if self.lr is not None:
            D += na_1d(s, qz, **self.lr)
        return D


def modes(s, D):
    nat, M = s['nat'], s['mass']
    Dm = (D / np.sqrt(M[:, None, None, None] * M[None, None, :, None])).reshape(3 * nat, 3 * nat)
    w2, v = np.linalg.eigh(0.5 * (Dm + Dm.conj().T))
    f = np.sign(w2) * np.sqrt(np.abs(w2)) * RY_TO_CMM1
    ev = v.reshape(nat, 3, -1) / np.sqrt(M)[:, None, None]
    Pz = np.abs(np.einsum('aj,ajn->n', s['Z'][:, 2, :], ev))
    return f, Pz, v


def lo(s, D):
    """Frequency of the mode with the largest longitudinal mode effective charge."""
    f, Pz, _ = modes(s, D)
    return float(f[np.argmax(Pz)])


def atomic_radius(s):
    xy = s['tau'][:, :2]
    return float(np.mean(np.linalg.norm(xy - xy.mean(0), axis=1)))


def parse_modes(path, nat):
    """matdyn.modes -> list of (q (2pi/alat), freqs cm^-1, eigenvectors [mode, atom, 3])."""
    txt = open(path).read().split('\n')
    out, i = [], 0
    while i < len(txt):
        if txt[i].strip().startswith('q ='):
            q = np.array([float(x) for x in txt[i].split('=')[1].split()])
            fr, ev, j = [], [], i + 2
            for _ in range(3 * nat):
                mm = re.search(r'=\s*([-\d.]+) \[THz\] =\s*([-\d.]+) \[cm-1\]', txt[j])
                fr.append(float(mm.group(2)))
                vec = []
                for a in range(nat):
                    v = [float(x) for x in txt[j + 1 + a].replace('(', ' ').replace(')', ' ').split()]
                    vec.append(np.array(v[0::2]) + 1j * np.array(v[1::2]))
                ev.append(np.array(vec))
                j += 1 + nat
            out.append((q, np.array(fr), np.array(ev)))
            i = j
        else:
            i += 1
    return out
