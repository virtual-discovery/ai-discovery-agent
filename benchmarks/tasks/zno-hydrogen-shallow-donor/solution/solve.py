"""Reference solution for zno-hydrogen-shallow-donor; runs as the Harbor oracle and in check_task."""
import json
import math
import os
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.special import erfc

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
COUL, KB = 14.399645, 8.617333262e-5
M0C2, HBARC = 0.51099895e6, 1.973269804e-5
EPS, ME, MH, T, REXCL = 10.4, 0.24, 0.59, 1000.0, 5.0
DEFECTS = {'V_O': [0, 1, 2], 'Zn_i': [0, 1, 2], 'OH_O': [0, 1], 'H_O': [0, 1, 2], 'H_Zn': [-1, -2],
           'Zn_i-H_O': [0, 1, 2], 'H_Zn-OH_O': [0, -1]}
EVIDENCE = []


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(kind, rid):
    r = lab('evaluate', {'kind': kind, 'id': rid})
    while r['status'] == 'running':
        time.sleep(1)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(r)
    EVIDENCE.append(r['id'])
    return r['result']


def did(defect, q):
    return 'supercell_3_3_2/%s/q_%s%d/HSE06' % (defect, '+' if q > 0 else '', q)


# ------------------------------------------------------------------ electrostatics and alignment
def madelung(lattice, q):
    A = np.array(lattice, float)
    V = abs(np.linalg.det(A))
    a = math.sqrt(math.pi) / V ** (1 / 3)
    B = 2 * math.pi * np.linalg.inv(A).T
    rmax, gmax = math.sqrt(23.0) / a, 2 * a * math.sqrt(23.0)
    nr = [int(rmax / np.linalg.norm(A[i])) + 2 for i in range(3)]
    ng = [int(gmax / np.linalg.norm(B[i])) + 2 for i in range(3)]
    idx = lambda n: np.array([(i, j, k) for i in range(-n[0], n[0] + 1) for j in range(-n[1], n[1] + 1)
                              for k in range(-n[2], n[2] + 1) if (i, j, k) != (0, 0, 0)], float)
    r = np.linalg.norm(idx(nr) @ A, axis=1)
    g2 = np.sum((idx(ng) @ B) ** 2, axis=1)
    e = 0.5 * q * q * (np.sum(erfc(a * r) / r) + 4 * math.pi / V * np.sum(np.exp(-g2 / (4 * a * a)) / g2)
                       - 2 * a / math.sqrt(math.pi) - math.pi / (V * a * a))
    return COUL * e


def mic(L, d):
    d = np.array(d) - np.round(d)
    return np.linalg.norm(d @ np.array(L), axis=-1)


def mapping(s1, s2, tol=0.5):
    f2 = np.array(s2['frac'])
    sp2 = np.array(s2['species'])
    out = []
    for f, sp in zip(s1['frac'], s1['species']):
        cand = np.where(sp2 == sp)[0]
        if len(cand) == 0:
            out.append(-1)
            continue
        d = mic(s1['lattice_A'], f2[cand] - np.array(f))
        out.append(int(cand[np.argmin(d)]) if d.min() <= tol else -1)
    return out


def alignment(host, hpot, defect, dpot):
    m = mapping(host['final_structure'], defect['final_structure'])
    back = mapping(defect['final_structure'], host['final_structure'])
    centres = [host['final_structure']['frac'][i] for i, s in enumerate(m) if s == -1]
    centres += [defect['final_structure']['frac'][i] for i, s in enumerate(back) if s == -1]
    ref = np.array(centres[0])
    c = ref + np.mean([np.array(x) - np.round(np.array(x) - ref) - ref for x in centres], axis=0)
    L = host['final_structure']['lattice_A']
    diffs = [dpot[s] - hpot[i] for i, s in enumerate(m)
             if s != -1 and mic(L, np.array(host['final_structure']['frac'][i]) - c) >= REXCL]
    return float(np.mean(diffs))


def edges(eig):
    vb, cb = -1e9, 1e9
    for ks in eig['eigenvalues'].values():
        for b in ks:
            for n, e, o in b:
                if o > 0.5:
                    vb = max(vb, e)
                else:
                    cb = min(cb, e)
    return vb, cb


def band_filling(eig, ecbm, evbm):
    w = np.array(eig['kpoint_weights'], float)
    w = w / w.sum()   # the k-points of the 2x2x2 mesh carry unequal weights
    corr = 0.0
    for ks in eig['eigenvalues'].values():
        for ik, b in enumerate(ks):
            for n, e, o in b:
                if e > ecbm and o > 1e-6:
                    corr -= w[ik] * o * (e - ecbm)
                if e < evbm and o < 1 - 1e-6:
                    corr -= w[ik] * (1 - o) * (evbm - e)
    return corr


# ------------------------------------------------------------------ statistics
def carriers(ef, gap):
    kT = KB * T
    x = np.linspace(0, 40, 8001)
    F = lambda eta: 2 / math.sqrt(math.pi) * np.trapz(np.sqrt(x) / (1 + np.exp(np.clip(x - eta, -700, 700))), x)
    nc = lambda m: 2.0 * (m * M0C2 * kT / (2 * math.pi * HBARC ** 2)) ** 1.5
    return nc(ME) * F((ef - gap) / kT), nc(MH) * F(-ef / kT)


def sc_fermi(E0, gap, nsite):
    kT = KB * T
    def net(ef):
        n, p = carriers(ef, gap)
        return p - n + sum(q * nsite * math.exp(-(v + q * ef) / kT) for (d, q), v in E0.items())
    lo, hi = -0.5, gap + 0.5
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if net(mid) > 0 else (lo, mid)
    ef = 0.5 * (lo + hi)
    return ef, {k: nsite * math.exp(-(v + k[1] * ef) / kT) for k, v in E0.items()}


def ladder(E0q, gap):
    grid = np.linspace(-1.0, gap + 1.0, 20001)
    st = [min(E0q, key=lambda q: E0q[q] + q * x) for x in grid]
    out = []
    for i in range(1, len(grid)):
        if st[i] != st[i - 1]:
            q1, q2 = st[i - 1], st[i]
            out.append({'from': q1, 'to': q2, 'level_eV': (E0q[q2] - E0q[q1]) / (q1 - q2)})
    return out, sorted(set(s for s, x in zip(st, grid) if 0 <= x <= gap))


def classify_site(run):
    """Bond-centre vs antibonding, parallel vs perpendicular to c, from the relaxed structure."""
    s = run['final_structure']
    L = np.array(s['lattice_A'])
    f = np.array(s['frac'])
    sp = np.array(s['species'])
    h = int(np.where(sp == 'H')[0][0])
    vec = lambda i, j: (lambda d: (d - np.round(d)) @ L)(f[j] - f[i])
    o = min(np.where(sp == 'O')[0], key=lambda j: np.linalg.norm(vec(h, j)))
    oh = vec(o, h)
    zns = [j for j in np.where(sp == 'Zn')[0] if np.linalg.norm(vec(o, j)) < 3.4]
    cosines = [float(np.dot(oh, vec(o, j)) / np.linalg.norm(oh) / np.linalg.norm(vec(o, j))) for j in zns]
    # bond centre: H lies on the segment from O towards one of its Zn neighbours; antibonding: on the far side of O
    kind = 'BC' if any(c > 0.8 and np.linalg.norm(vec(o, j)) > np.linalg.norm(oh) for c, j in zip(cosines, zns)) else 'AB'
    caxis = L[2] / np.linalg.norm(L[2])
    orient = 'parallel' if abs(np.dot(oh, caxis)) / np.linalg.norm(oh) > 0.9 else 'perpendicular'
    return f'{kind}_{orient}', float(np.linalg.norm(oh)), [round(c, 2) for c in cosines]


def main():
    ref_phases = ['Zn', 'ZnO', 'O2', 'H2', 'H2O', 'ZnOH2']
    phases = {p: buy('run', 'chem_pot/%s/HSE06' % p) for p in ref_phases}
    host = buy('run', 'supercell_3_3_2/HSE06')
    host_eig = buy('eigenvalues', 'supercell_3_3_2/HSE06')
    vbm, cbm = edges(host_eig)
    gap = cbm - vbm
    runs = {(d, q): buy('run', did(d, q)) for d, qs in DEFECTS.items() for q in qs}
    # Cells whose magnetization does not match a localized state (odd electron count but |m| != 1, or fractional m)
    # hold carriers in host bands: buy their eigenvalues to quantify the band filling.
    nel_host = host['NELECT']
    suspicious = []
    for (d, q), r in runs.items():
        m = abs(r['total_magnetization_muB'])
        odd = int(round(r['NELECT'])) % 2 == 1
        if (odd and abs(m - 1) > 0.1) or (abs(m - round(m)) > 0.1):
            suspicious.append((d, q))
    eig = {}
    for key in [('OH_O', 0), ('H_O', 0), ('Zn_i', 1)]:
        eig[key] = buy('eigenvalues', did(*key))
    for key in [('OH_O', 1), ('V_O', 1)]:
        eig[key] = buy('eigenvalues', did(*key))
    # chemical potentials: Zn-rich, maximal mu_H allowed by every H-containing reference phase
    per = lambda p, el: phases[p]['final_structure']['species'].count(el)
    E = lambda p: phases[p]['final_energy_sigma0_eV']
    muZn = E('Zn') / per('Zn', 'Zn')
    muO = E('ZnO') / per('ZnO', 'Zn') - muZn
    nfu = per('ZnOH2', 'Zn')
    hlim = {'H2': E('H2') / 2, 'H2O': (E('H2O') - muO) / 2, 'ZnOH2': (E('ZnOH2') / nfu - muZn - 2 * muO) / 2}
    muH = min(hlim.values())
    mu = {'Zn': muZn, 'O': muO, 'H': muH}
    dmu = {'O': muO - E('O2') / 2, 'H': muH - E('H2') / 2}
    hs = host['final_structure']['species']
    E0, parts = {}, {}
    for (d, q), r in runs.items():
        ds = r['final_structure']['species']
        ef = r['final_energy_sigma0_eV'] - host['final_energy_sigma0_eV'] - sum(
            (ds.count(el) - hs.count(el)) * mu[el] for el in mu)
        pa = alignment(host, host['core_potentials_eV'], r, r['core_potentials_eV'])
        ic = -2.0 / 3.0 * madelung(r['final_structure']['lattice_A'], q) / EPS if q else 0.0
        bf = band_filling(eig[(d, q)], cbm + pa, vbm + pa) if (d, q) in eig else 0.0
        E0[(d, q)] = ef + ic + bf + q * (vbm + pa)
        parts['%s %+d' % (d, q)] = {'raw': round(ef + q * vbm, 3), 'image': round(ic, 3), 'align': round(pa, 3),
                                    'band_filling': round(bf, 3), 'E_f(VBM)': round(E0[(d, q)], 3),
                                    'magnetization': r['total_magnetization_muB']}
    nsite = 36 / (abs(np.linalg.det(host['final_structure']['lattice_A'])) * 1e-24)
    ef, conc = sc_fermi(E0, gap, nsite)
    levels, stable = {}, {}
    for d in DEFECTS:
        levels[d], stable[d] = ladder({q: E0[(d, q)] for q in DEFECTS[d]}, gap)
    in_gap = lambda d: [x for x in levels[d] if 0.0 <= x['level_eV'] <= gap]
    # EPR candidates: odd-electron (paramagnetic) charge states that are thermodynamically stable somewhere in the gap
    para = []
    for (d, q), r in runs.items():
        if int(round(r['NELECT'])) % 2 == 1 and q in stable[d]:
            para.append((conc[(d, q)], d, q))
    para.sort(reverse=True)
    site, doh, cosines = classify_site(runs[('OH_O', 0)])
    lg = math.log10
    submission = {
        'band_gap_eV': gap,
        'delta_mu_eV': dmu,
        'transition_levels': {d: in_gap(d) for d in ('OH_O', 'H_O', 'V_O')},
        'V_O_plus_ever_stable': 1 in stable['V_O'],
        'epr_centre': {'defect': para[0][1], 'charge': para[0][2]},
        'oh_o_site': site,
        'fermi_level_eV': ef,
        'log10_concentration_cm3': {'V_O': lg(sum(conc[('V_O', q)] for q in (0, 1, 2))),
                                    'OH_O^+': lg(conc[('OH_O', 1)]), 'OH_O^0': lg(conc[('OH_O', 0)])},
        'evidence_ids': EVIDENCE,
    }
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    table = '\n'.join(f'| {k} | {v["raw"]} | {v["image"]} | {v["align"]} | {v["band_filling"]} | {v["E_f(VBM)"]} | {v["magnetization"]} |'
                      for k, v in parts.items())
    (APP/'report.md').write_text(f"""# Hydrogen and the g = 1.96 shallow donor in ZnO (reference solution)

## Hamiltonian and band edges
The defect energetics use the hybrid single points (HSE06 folders; INCAR: LHFCALC, AEXX = 0.37, ENCUT 400 eV) on the
PBE-relaxed structures. PBE gives a 0.72 eV gap and puts every level in the wrong place; the hybrid supercell gives
E_VBM = {vbm:.3f} eV, E_CBM = {cbm:.3f} eV, gap {gap:.3f} eV. All reference phases are taken at the same level.

## Chemical potentials
Zn-rich: mu_Zn = E(Zn metal), mu_O = E(ZnO) - mu_Zn, i.e. Delta mu_O = {dmu['O']:.3f} eV. Hydrogen limits
(Delta mu_H relative to H2): {json.dumps({k: round(v - E('H2') / 2, 3) for k, v in hlim.items()})}. Zn(OH)2 (4 formula
units per cell) is the binding phase, so Delta mu_H = {dmu['H']:.3f} eV rather than the H2 limit (0).

## Finite-size treatment
Localized charges: Lany-Zunger image charge (2/3 of the point-charge Madelung energy, eps = 10.4) and alignment of the
core potentials beyond 5 A. Cells flagged by their magnetization ({suspicious}) were inspected. The neutral OH_O and
H_O cells (odd electron count, total moment -0.25 and 0) do not hold a localized spin: the extra electron fills
conduction-band states at several k-points of the 2x2x2 mesh (weights 1/8 and 3/8), up to ~1 eV above the CBM. That
is a perturbed-host (effective-mass) state; I removed the band filling with the OUTCAR k-weights (Lany-Zunger/Persson),
instead of an image-charge correction that would treat it as localized. Without it OH_O(+/0) lies 0.5 eV above the
CBM and OH_O^0 is never stable, which contradicts a shallow donor. Zn_i^+ and Zn_i^0 are likewise Zn_i^2+ with
electrons in the conduction band.

| state | raw | image | align | band filling | E_f(E_F = VBM) | magnetization |
|---|---|---|---|---|---|---|
{table}

## Levels and the EPR centre
OH_O(+/0) = {in_gap('OH_O')}, H_O(+/0) = {in_gap('H_O')}: both are shallow donors just below the CBM (an effective-
mass binding of tens of meV is below what a 72-atom cell resolves). V_O has a negative-U (2+/0) level at
{in_gap('V_O')}: V_O^+ is never the stable state, so V_O cannot give an equilibrium EPR signal; its levels are deep.
H in OH_O sits at a {site} site (O-H {doh:.2f} A, cosines to the O-Zn bonds {cosines}).
Paramagnetic stable states at the 1000 K Fermi level, by concentration: {[(d, q, f'{c:.2e}') for c, d, q in para[:4]]}.
The g = 1.96 centre is therefore neutral OH_O (H at the bond centre, an effective-mass donor), introduced during synthesis.

## Equilibrium at 1000 K
Charge neutrality with parabolic bands (m_e 0.24, m_h 0.59, Fermi-Dirac) and Boltzmann defect statistics gives
E_F = {ef:.3f} eV above the VBM, set by ionized OH_O^+ against free electrons: [OH_O^+] = {conc[('OH_O', 1)]:.2e},
[OH_O^0] = {conc[('OH_O', 0)]:.2e}, [V_O] = {sum(conc[('V_O', q)] for q in (0, 1, 2)):.2e} cm^-3.

## What the data do not support
The OH_O^0 concentration depends exponentially on the 0.1 eV-scale position of a shallow level that the supercell
cannot resolve, so it is only good to within about half a decade. The equilibrium values describe the synthesis
temperature; nothing in the archive describes hydrogen out-diffusion during the 873 K vacuum anneal or the quench, so
the room-temperature EPR spin density is not predicted. Only one H configuration per defect was archived (no AB site,
no vibrational zero-point energies), the chemical potentials are 0 K DFT values, and eps = 10.4 is the authors' choice.
""")
    lab('submit', submission)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
