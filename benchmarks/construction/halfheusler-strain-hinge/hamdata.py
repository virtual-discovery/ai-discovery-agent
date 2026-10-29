"""Reader for the FPLO '+hamdata' Wannier Hamiltonians of Zenodo record 16927281 (numpy only).

File layout (FPLO 18.57, full-relativistic, 8 Wannier functions): 'lattice_vectors' (rows, Bohr), 'wannames',
'wancenters' (Bohr), then one block per orbital pair 'Tij, Hij:' / 'i j' / rows 'Tx Ty Tz ReH ImH' / 'end Tij, Hij:'.
T is the vector (Bohr) from the centre of orbital i in the home cell to the centre of orbital j in cell R,
T = R + w_j - w_i, as in the authors' process_hamdata() (code/strained_hH.py).  Energies are in eV.
"""
import numpy as np

BOHR = 0.529177210903  # Angstrom


def read(path):
    lines = open(path).read().splitlines()
    i = lines.index('lattice_vectors:')
    lat = np.array([[float(x) for x in lines[i + 1 + j].split()] for j in range(3)])
    nw = int(lines[lines.index('nwan:') + 1])
    j = lines.index('wannames:')
    names = [lines[j + 1 + n].strip() for n in range(nw)]
    j = lines.index('wancenters:')
    cen = np.array([[float(x) for x in lines[j + 1 + n].split()] for n in range(nw)])
    inv = np.linalg.inv(lat.T)
    blocks = {}
    n = 0
    while n < len(lines):
        if lines[n].startswith('Tij, Hij:'):
            a, b = [int(x) - 1 for x in lines[n + 1].split()]
            m = n + 2
            rows = []
            while not lines[m].startswith('end Tij, Hij:'):
                rows.append([float(x) for x in lines[m].split()])
                m += 1
            rows = np.array(rows)
            R = inv @ (rows[:, :3] - cen[b] + cen[a]).T
            Ri = np.round(R).astype(int)
            assert np.abs(R - Ri).max() < 1e-4, path
            for r, h in zip(map(tuple, Ri.T), rows[:, 3] + 1j * rows[:, 4]):
                blocks.setdefault(r, np.zeros((nw, nw), complex))[a, b] += h
            n = m
        n += 1
    Rs = sorted(blocks)
    return {'lattice_bohr': lat, 'centers_bohr': cen, 'names': names, 'R': np.array(Rs, int),
            'H': np.array([blocks[r] for r in Rs])}
