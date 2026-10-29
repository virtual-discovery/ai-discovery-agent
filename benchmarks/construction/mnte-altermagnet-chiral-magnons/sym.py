"""P6_3/mmc (No. 194) operations and Mn-Mn bond classes of NiAs-type MnTe."""
import numpy as np

A_LAT, C_LAT = 4.148, 6.711
CELL = np.array([[A_LAT, 0, 0], [-A_LAT / 2, A_LAT * np.sqrt(3) / 2, 0], [0, 0, C_LAT]])
MN = np.array([[0, 0, 0], [0, 0, 0.5]])
TE = np.array([[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]])

_G = ['x,y,z', '-y,x-y,z', '-x+y,-x,z', '-x,-y,z+1/2', 'y,-x+y,z+1/2', 'x-y,x,z+1/2',
      'y,x,-z', 'x-y,-y,-z', '-x,-x+y,-z', '-y,-x,-z+1/2', '-x+y,y,-z+1/2', 'x,x-y,-z+1/2']


def _op(s):
    rows, t = [], []
    for e in s.split(','):
        r = np.zeros(3)
        e = e.replace('-', '+-').lstrip('+')
        tr = 0.0
        for term in e.split('+'):
            if not term:
                continue
            sign = -1 if term.startswith('-') else 1
            term = term.lstrip('-')
            if term in 'xyz':
                r['xyz'.index(term)] += sign
            else:
                tr += sign * eval(term)
        rows.append(r)
        t.append(tr)
    return np.array(rows), np.array(t)


OPS = [_op(s) for s in _G]
OPS += [(-R, -t) for R, t in OPS]


def site_index(frac):
    for k, p in enumerate(MN):
        d = frac - p
        if np.allclose(d, np.round(d), atol=1e-6):
            return k, np.round(d).astype(int)
    raise ValueError(frac)


def bond_classes(rmax=10.5):
    """Classes of directed bonds (i, j, R): Mn_i at r_i -> Mn_j at r_j + R."""
    rng = range(-4, 5)
    bonds = []
    for i in range(2):
        for j in range(2):
            for a in rng:
                for b in rng:
                    for cc in range(-3, 4):
                        R = np.array([a, b, cc])
                        d = (MN[j] + R - MN[i]) @ CELL
                        dist = np.linalg.norm(d)
                        if 1e-6 < dist < rmax:
                            bonds.append((i, j, tuple(R), round(dist, 4)))
    key = {b[:3]: n for n, b in enumerate(bonds)}
    cls = [-1] * len(bonds)
    nc = 0
    for n, (i, j, R, dist) in enumerate(bonds):
        if cls[n] >= 0:
            continue
        for Rot, t in OPS:
            for (p, q, RR) in [(i, j, R), (j, i, tuple(-np.array(R)))]:   # bonds are undirected exchange
                ri = Rot @ MN[p] + t
                rj = Rot @ (MN[q] + np.array(RR)) + t
                ki, Ti = site_index(ri)
                kj, Tj = site_index(rj)
                k2 = (ki, kj, tuple(Tj - Ti))
                if k2 in key:
                    if cls[key[k2]] < 0:
                        cls[key[k2]] = nc
        nc += 1
    return bonds, cls


if __name__ == '__main__':
    bonds, cls = bond_classes()
    import collections
    info = collections.OrderedDict()
    for b, c in zip(bonds, cls):
        info.setdefault(c, []).append(b)
    rows = sorted(info.items(), key=lambda kv: kv[1][0][3])
    for n, (c, bl) in enumerate(rows, 1):
        i, j, R, dist = bl[0]
        kind = 'AB' if i != j else 'same'
        fromA = sorted(set(b[2] for b in bl if b[0] == 0 and b[1] == 0))
        fromB = sorted(set(b[2] for b in bl if b[0] == 1 and b[1] == 1))
        alt = kind == 'same' and set(fromA) != set(fromB)
        print(n, 'class', c, 'dist %.4f' % dist, kind, 'mult/site', sum(1 for b in bl if b[0] == 0), 'ALTERMAGNETIC' if alt else '',
              fromA[:6] if alt else '')
