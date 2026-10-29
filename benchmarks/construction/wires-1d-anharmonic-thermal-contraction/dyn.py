"""Reader for Quantum ESPRESSO / CellConstructor dynamical-matrix files and simple 1D lattice dynamics.

Units: force constants in Ry/bohr^2, masses in Ry atomic units (amu * 911.444), lengths in bohr.
"""
import re
import numpy as np

RY_TO_CM = 109737.31570111268


def read_dyn(path):
    """Return dict(alat, cell (bohr), masses (per atom), species, pos (bohr), blocks=[(q_2pi/alat, Phi[nat,nat,3,3] complex)])."""
    lines = open(path).read().splitlines()
    head = lines[2].split()
    ntyp, nat, ibrav = int(head[0]), int(head[1]), int(head[2])
    alat = float(head[3])
    i = 3
    cell = None
    if ibrav == 0:
        assert lines[i].strip().startswith('Basis')
        cell = np.array([[float(x) for x in lines[i + k].split()] for k in (1, 2, 3)]) * alat
        i += 4
    tmass, tname = {}, {}
    for k in range(ntyp):
        m = re.match(r"\s*(\d+)\s+'([^']*)'\s+(\S+)", lines[i + k])
        tname[int(m.group(1))] = m.group(2).strip()
        tmass[int(m.group(1))] = float(m.group(3))
    i += ntyp
    species, masses, pos = [], [], []
    for k in range(nat):
        p = lines[i + k].split()
        t = int(p[1])
        species.append(tname[t]); masses.append(tmass[t]); pos.append([float(x) for x in p[2:5]])
    i += nat
    pos = np.array(pos) * alat
    blocks = []
    while i < len(lines):
        if 'Dynamical  Matrix in cartesian axes' in lines[i] or 'Dynamical Matrix in cartesian axes' in lines[i]:
            j = i + 1
            while 'q = (' not in lines[j]:
                j += 1
            q = np.array([float(x) for x in lines[j].split('(')[1].split(')')[0].split()])
            j += 1
            phi = np.zeros((nat, nat, 3, 3), complex)
            for _ in range(nat * nat):
                while not lines[j].strip():
                    j += 1
                a, b = (int(x) for x in lines[j].split())
                for r in range(3):
                    v = [float(x) for x in lines[j + 1 + r].split()]
                    phi[a - 1, b - 1, r] = np.array(v[0::2]) + 1j * np.array(v[1::2])
                j += 4
            blocks.append((q, phi))
            i = j
        else:
            i += 1
    return {'alat': alat, 'cell': cell, 'masses': np.array(masses), 'species': species, 'pos': pos, 'blocks': blocks}


def read_set(prefix, nq):
    """Read files prefix1..prefixN (CellConstructor naming); return header of first and list of (q, Phi)."""
    first = read_dyn(f'{prefix}1')
    blocks = list(first['blocks'])
    for k in range(2, nq + 1):
        blocks += read_dyn(f'{prefix}{k}')['blocks']
    first['blocks'] = blocks
    return first


def frequencies(d, phi):
    """Eigen-frequencies (cm^-1, negative = imaginary) of one q block."""
    nat = len(d['masses'])
    m = np.sqrt(np.repeat(d['masses'], 3))
    D = phi.transpose(0, 2, 1, 3).reshape(3 * nat, 3 * nat) / np.outer(m, m)
    D = 0.5 * (D + D.conj().T)
    w2 = np.linalg.eigvalsh(D)
    return np.sign(w2) * np.sqrt(np.abs(w2)) * RY_TO_CM


def real_space_ifc(d, ncell):
    """Inverse Fourier transform of the commensurate q blocks of a 1 x 1 x ncell grid along z.

    Returns R (bohr, along z) and Phi_R [nR, nat, nat, 3, 3] (real), with Wigner-Seitz weights for |n| = ncell/2.
    Convention (QE): Phi(q) = sum_R Phi(0, R) exp(-i q R).
    """
    c = d['cell'][2, 2]
    ns = list(range(-(ncell // 2), ncell // 2 + 1))
    out = []
    for n in ns:
        acc = 0
        for q, phi in d['blocks']:
            qz = q[2] * 2 * np.pi / d['alat']
            acc = acc + phi * np.exp(1j * qz * n * c)
        acc = acc / len(d['blocks'])
        w = 0.5 if (ncell % 2 == 0 and abs(n) == ncell // 2) else 1.0
        out.append(acc.real * w)
    return np.array(ns) * c, np.array(out)


def dynmat_at(d, R, PhiR, qz):
    """Interpolated mass-reduced dynamical matrix at qz (1/bohr)."""
    nat = len(d['masses'])
    phi = sum(P * np.exp(-1j * qz * r) for r, P in zip(R, PhiR))
    m = np.sqrt(np.repeat(d['masses'], 3))
    D = phi.transpose(0, 2, 1, 3).reshape(3 * nat, 3 * nat) / np.outer(m, m)
    return 0.5 * (D + D.conj().T)
