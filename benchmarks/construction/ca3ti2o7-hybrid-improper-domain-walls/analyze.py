"""Local order parameters of Ca3Ti2O7 structures: octahedral rotation (psi_z), tilts (phi_x, phi_y)
and polar displacements per Ti-O octahedron, binned along a chosen axis."""
import numpy as np


def mic(d, cell):
    f = d @ np.linalg.inv(cell)
    f -= np.round(f)
    return f @ cell


def octahedra(symbols, pos, cell, theta0=0.0):
    sym = np.array(symbols)
    ti = np.nonzero(sym == 'Ti')[0]
    ox = np.nonzero(sym == 'O')[0]
    res = []
    for t in ti:
        d = mic(pos[ox] - pos[t], cell)
        r = np.linalg.norm(d, axis=1)
        nb = np.argsort(r)[:6]
        v = d[nb]
        ap = nb[np.argsort(-np.abs(v[:, 2]))[:2]]
        eq = [k for k in nb if k not in ap]
        dap = d[ap]; deq = d[eq]
        # rotation about z: angle of equatorial vectors relative to nearest of +-x, +-y
        ang = []
        for w in deq:
            a = np.arctan2(w[1], w[0]) - theta0
            a0 = np.round(a / (np.pi / 2)) * (np.pi / 2)
            ang.append(a - a0)
        psi = np.degrees(np.mean(ang))
        # tilt: apical O offsets (top minus bottom) in x and y
        top, bot = (dap[0], dap[1]) if dap[0][2] > 0 else (dap[1], dap[0])
        axis = top - bot
        tilt_x = np.degrees(np.arctan2(axis[1], axis[2]))   # rotation about x tips the axis toward -y... sign kept raw
        tilt_y = np.degrees(np.arctan2(axis[0], axis[2]))
        cen = pos[t] + np.mean(d[nb], axis=0)
        off = -np.mean(d[nb], axis=0)                          # Ti off-centring relative to O6 centroid
        res.append({'ti': int(t), 'pos': pos[t], 'psi': psi, 'tx': tilt_x, 'ty': tilt_y, 'off': off})
    return res


def reference_sites(symbols, pos, cell, T, ref_struct):
    """Match each atom to an ideal I4/mmm site; returns displacement u (N,3).
    T: (2,3) fractional coordinates (in this cell) of the two in-plane tetragonal translations.
    ref_struct: (symbols, fractional positions) of the 24-atom I4/mmm conventional cell."""
    rs, rf = ref_struct
    sym = np.array(symbols)
    T = np.asarray(T, float)
    if len(T) == 2:
        T = np.vstack([T, [0, 0, 1.0]])
    base = np.array(rf)
    bfrac = base @ T
    n = int(round(1 / np.abs(T).max(0)[0:2].min())) + 1 if False else 0
    span = int(round(max(1 / max(abs(T[0, 0]), 1e-9), 1 / max(abs(T[1, 1]), 1e-9)))) + 2
    sites, ssym = [], []
    nz = int(round(1 / T[2, 2]))
    for i in range(-span, span + 1):
        for j in range(-span, span + 1):
            for k in range(nz):
                sites.append(bfrac + i * T[0] + j * T[1] + k * T[2]); ssym += list(rs)
    sites = np.concatenate(sites) @ cell; ssym = np.array(ssym)
    ti = np.nonzero(sym == 'Ti')[0]
    tsites = sites[ssym == 'Ti']
    best = None
    for k in range(len(tsites)):
        shift = pos[ti[0]] - tsites[k]
        d = mic(pos[ti][:, None, :] - (tsites[None, :, :] + shift), cell)
        err = np.linalg.norm(d, axis=2).min(1).mean()
        if best is None or err < best[0] - 1e-9:
            best = (err, shift)
        if err < 0.05:
            break
    shift = best[1]
    u = np.zeros_like(pos)
    for s in set(symbols):
        m = sym == s
        st = sites[ssym == s] + shift
        d = mic(pos[m][:, None, :] - st[None, :, :], cell)
        r = np.linalg.norm(d, axis=2)
        u[m] = d[np.arange(m.sum()), r.argmin(1)]
    # remove the arbitrary uniform shift of the Ti sublattice origin choice
    return u, best[0]
