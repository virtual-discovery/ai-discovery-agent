"""Defect thermodynamics utilities for zno-hydrogen-shallow-donor (numpy/scipy only; runs in the lab image)."""
import json
import math
import numpy as np
from pathlib import Path
from scipy.special import erfc

COUL = 14.399645  # e^2/(4 pi eps0) in eV A
KB = 8.617333262e-5


def madelung_point_charge(lattice, q=1.0, tol=1e-10):
    """Ewald energy of a lattice of point charges q with neutralizing background (eV), like pymatgen EwaldSummation."""
    A = np.array(lattice, float)
    V = abs(np.linalg.det(A))
    alpha = math.sqrt(math.pi) * (1.0 / V) ** (1 / 3)  # splitting (1/A)
    rmax = math.sqrt(-math.log(tol)) / alpha
    gmax = 2 * alpha * math.sqrt(-math.log(tol))
    B = 2 * math.pi * np.linalg.inv(A).T
    nr = [int(math.ceil(rmax / np.linalg.norm(A[i]) * 1.0)) + 1 for i in range(3)]
    ng = [int(math.ceil(gmax / np.linalg.norm(B[i]))) + 1 for i in range(3)]
    real = 0.0
    for i in range(-nr[0], nr[0] + 1):
        for j in range(-nr[1], nr[1] + 1):
            for k in range(-nr[2], nr[2] + 1):
                if i == j == k == 0:
                    continue
                r = np.linalg.norm(i * A[0] + j * A[1] + k * A[2])
                real += erfc(alpha * r) / r
    recip = 0.0
    for i in range(-ng[0], ng[0] + 1):
        for j in range(-ng[1], ng[1] + 1):
            for k in range(-ng[2], ng[2] + 1):
                if i == j == k == 0:
                    continue
                g2 = np.sum((i * B[0] + j * B[1] + k * B[2]) ** 2)
                recip += math.exp(-g2 / (4 * alpha ** 2)) / g2
    e = 0.5 * q * q * (real + 4 * math.pi / V * recip - 2 * alpha / math.sqrt(math.pi) - math.pi / (V * alpha ** 2))
    return COUL * e


def lz_image_charge(lattice, q, eps):
    return -2.0 / 3.0 * madelung_point_charge(lattice, q) / eps


def frac_dist(lattice, f1, f2):
    d = np.array(f2) - np.array(f1)
    d -= np.round(d)
    return float(np.linalg.norm(d @ np.array(lattice)))


def atom_mapping(s1, s2, tol=0.5):
    """For each atom of s1, the index of the nearest same-species atom of s2 within tol (else -1)."""
    L = np.array(s1['lattice_A'])
    f1, f2 = np.array(s1['frac']), np.array(s2['frac'])
    out = []
    for i, sp in enumerate(s1['species']):
        best, bi = 1e9, -1
        for j, sp2 in enumerate(s2['species']):
            if sp2 != sp:
                continue
            d = f2[j] - f1[i]
            d -= np.round(d)
            dist = np.linalg.norm(d @ L)
            if dist < best:
                best, bi = dist, j
        out.append(bi if best <= tol else -1)
    return out


def defect_center(host, defect):
    nh, nd = len(host['species']), len(defect['species'])
    if nh == nd + 1:
        m = atom_mapping(host, defect)
        return [host['frac'][i] for i, s in enumerate(m) if s == -1][-1]
    if nh == nd - 1:
        m = atom_mapping(defect, host)
        return [defect['frac'][i] for i, s in enumerate(m) if s == -1][-1]
    if nh == nd:
        m = atom_mapping(defect, host)
        return [defect['frac'][i] for i, s in enumerate(m) if s == -1][-1]
    # complex: reproduces the authors' averaging (their code indexes host coordinates for both passes)
    c, n = np.zeros(3), 0
    for i, s in enumerate(atom_mapping(host, defect)):
        if s == -1:
            c += host['frac'][i]; n += 1
    for i, s in enumerate(atom_mapping(defect, host)):
        if s == -1:
            c += host['frac'][i]; n += 1
    return (c / n).tolist()


def potential_alignment(host, hpot, defect, dpot, rexcl=5.0, center=None):
    center = defect_center(host, defect) if center is None else center
    L = host['lattice_A']
    m = atom_mapping(host, defect)
    diffs, rows = [], []
    for i, s in enumerate(m):
        if s == -1:
            continue
        d = frac_dist(L, host['frac'][i], center)
        rows.append((d, dpot[s] - hpot[i]))
        if d >= rexcl:
            diffs.append(dpot[s] - hpot[i])
    return float(np.mean(diffs)), rows


def band_edges(run):
    vb, cb = -1e9, 1e9
    for sp, ks in run['eigenvalues'].items():
        for b in ks:
            for e, o in b:
                if o > 0.5:
                    vb = max(vb, e)
                else:
                    cb = min(cb, e)
    return vb, cb


# ---------------------------------------------------------------- thermodynamics
H_PLANCK = 4.135667696e-15  # eV s
M0C2 = 0.51099895e6         # eV
HBARC = 1.973269804e-5      # eV cm
EPS = 10.4
M_E, M_H = 0.24, 0.59
DEFECTS = {'V_O': [0, 1, 2], 'Zn_i': [0, 1, 2], 'OH_O': [0, 1], 'H_O': [0, 1, 2], 'H_Zn': [-1, -2],
           'Zn_i-H_O': [0, 1, 2], 'H_Zn-OH_O': [0, -1]}


def key(defect, q, level):
    return 'supercell_3_3_2/%s/q_%s%d/%s' % (defect, '+' if q > 0 else '', q, level)


def count(run, el):
    return run['final']['species'].count(el)


def chem_pots(P, level, h_limit='all', fu_errors=False):
    E = lambda k: P['chem_pot/%s/%s' % (k, level)]['energy_sigma0_eV']
    run = lambda k: P['chem_pot/%s/%s' % (k, level)]
    eZn = E('Zn') / count(run('Zn'), 'Zn')
    eZnO = E('ZnO') / count(run('ZnO'), 'Zn')
    eO2, eH2, eH2O = E('O2') / 2, E('H2') / 2, E('H2O')
    nfu = 1 if fu_errors else count(run('ZnOH2'), 'Zn')
    eZOH = E('ZnOH2') / nfu
    muZn, muO = eZn, eZnO - eZn
    lims = {'H2': eH2, 'H2O': (eH2O - muO) / 2, 'ZnOH2': (eZOH - muZn - 2 * muO) / 2}
    use = lims if h_limit == 'all' else {k: lims[k] for k in h_limit}
    muH = min(use.values())
    return {'Zn': muZn, 'O': muO, 'H': muH, 'dmu_O': muO - eO2, 'dmu_H': muH - eH2, 'limits_dmuH': {k: v - eH2 for k, v in lims.items()}}


def band_filling(run, e_cbm, e_vbm, weights=None):
    """Lany-Zunger/Persson band-filling correction: electrons above the aligned CBM, holes below the aligned VBM."""
    w = np.array(run['kweights'] if weights is None else weights, float)
    w = w / w.sum()
    occmax = 1.0 if run['ispin'] == 2 else 2.0
    corr = 0.0
    for sp, ks in run['eigenvalues'].items():
        for ik, b in enumerate(ks):
            for e, o in b:
                if e > e_cbm and o > 1e-6:
                    corr -= w[ik] * o * (e - e_cbm)
                if e < e_vbm and o < occmax - 1e-6:
                    corr -= w[ik] * (occmax - o) * (e_vbm - e)
    return corr


def formation_energies(P, level='HSE06', mu=None, lz=True, align=True, bf=True, rexcl=5.0, bf_align=True,
                       bf_equal_weights=False, eps=EPS, lz_factor=2.0 / 3.0):
    """Formation energies at E_F = VBM (eV) for every archived defect/charge."""
    host = P['supercell_3_3_2/' + level]
    vbm, cbm = band_edges(host)
    mu = mu or chem_pots(P, level)
    out = {}
    for d, qs in DEFECTS.items():
        for q in qs:
            r = P[key(d, q, level)]
            dn = {el: count(r, el) - count(host, el) for el in ('Zn', 'O', 'H')}
            ef = r['energy_sigma0_eV'] - host['energy_sigma0_eV'] - sum(dn[el] * mu[el] for el in dn)
            have_pot = 'core_potentials_eV' in r
            pa = potential_alignment(host['final'], host['core_potentials_eV'], r['final'], r['core_potentials_eV'], rexcl)[0] \
                if (align and have_pot) else 0.0
            ic = lz_factor * (-madelung_point_charge(r['final']['lattice_A'], q) / eps) if (lz and q) else 0.0
            b = 0.0
            if bf and 'eigenvalues' in r:
                shift = pa if bf_align else 0.0
                b = band_filling(r, cbm + shift, vbm + shift, [1] * len(r['kweights']) if bf_equal_weights else None)
            out[(d, q)] = {'E0': ef + ic + b + q * (vbm + pa), 'raw': ef + q * vbm, 'pa': pa, 'ic': ic, 'bf': b}
    return out, vbm, cbm


def carriers(ef, gap, T, me=M_E, mh=M_H):
    """Free electron and hole densities (cm^-3), parabolic bands, Fermi-Dirac statistics."""
    kT = KB * T
    def nc(m):
        return 2.0 * (m * M0C2 * kT / (2 * math.pi * HBARC ** 2)) ** 1.5
    x = np.linspace(0, 40, 8001)
    def F(eta):
        return 2 / math.sqrt(math.pi) * np.trapz(np.sqrt(x) / (1 + np.exp(np.clip(x - eta, -700, 700))), x)
    return nc(me) * F((ef - gap) / kT), nc(mh) * F(-ef / kT)


def sc_fermi(E0, gap, T, nsite, include=None):
    kT = KB * T
    items = [(dq, v) for dq, v in E0.items() if include is None or dq[0] in include]
    def net(ef):
        n, p = carriers(ef, gap, T)
        return p - n + sum(dq[1] * nsite * math.exp(-(v + dq[1] * ef) / kT) for dq, v in items)
    lo, hi = -0.5, gap + 0.5
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if net(mid) > 0:
            lo = mid
        else:
            hi = mid
    ef = 0.5 * (lo + hi)
    conc = {dq: nsite * math.exp(-(v + dq[1] * ef) / kT) for dq, v in items}
    n, p = carriers(ef, gap, T)
    return ef, conc, n, p


def stable_levels(E0q, gap):
    """Thermodynamic transition levels between consecutive stable charge states in [0, gap]."""
    qs = sorted(E0q)
    grid = np.linspace(-1.0, gap + 1.0, 20001)
    stable = [max(qs, key=lambda q: -(E0q[q] + q * x)) for x in grid]
    out = []
    for i in range(1, len(grid)):
        if stable[i] != stable[i - 1]:
            q1, q2 = stable[i - 1], stable[i]
            out.append({'from': q1, 'to': q2, 'level_eV': (E0q[q2] - E0q[q1]) / (q1 - q2)})
    return out, sorted(set(s for s, x in zip(stable, grid) if 0 <= x <= gap))
