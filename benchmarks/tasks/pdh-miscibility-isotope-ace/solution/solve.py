"""Reference solution for pdh-miscibility-isotope-ace; runs as the Harbor oracle and in check_task."""
import json
import math
import os
import subprocess
import time
from pathlib import Path
import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import eigsh

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
KB = 8.617333262e-5
HBAR = 0.0646541          # hbar * sqrt(eV / (amu A^2)) in eV
M_H, M_D = 1.00782503, 2.01410178
T300 = 300.0
GAS = CATALOG['gas_phase_300K_1bar']
BASE = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


EVIDENCE = []


def evaluate(req):
    r = lab('evaluate', req)
    while r['status'] == 'running':
        time.sleep(3)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(r.get('error'))
    EVIDENCE.append(r['id'])
    return r['result']


def fcc(a, n):
    pts = np.array([b + [i, j, k] for i in range(n) for j in range(n) for k in range(n) for b in BASE]) * a
    return (np.eye(3) * a * n).tolist(), pts


def structure(cell, symbols, positions):
    return {'cell': np.asarray(cell).tolist(), 'symbols': list(symbols), 'positions': np.asarray(positions).tolist()}


def schroedinger(V, h, mass, k=30):
    """Eigenvalues of -hbar^2/2m Lap + V on a cubic grid (Dirichlet walls)."""
    n = V.shape[0]
    D = sps.diags([np.ones(n - 1), -2 * np.ones(n), np.ones(n - 1)], [-1, 0, 1]) / h ** 2
    I = sps.identity(n)
    L = sps.kron(sps.kron(D, I), I) + sps.kron(sps.kron(I, D), I) + sps.kron(sps.kron(I, I), D)
    H = -HBAR ** 2 / (2 * mass) * L + sps.diags(V.ravel())
    return np.sort(eigsh(H.tocsc(), k=k, sigma=V.min() - 0.01, which='LM', return_eigenvectors=False))


def free_energy(levels, T):
    w = np.asarray(levels)
    return w[0] - KB * T * math.log(np.exp(-(w - w[0]) / (KB * T)).sum())


def main():
    # ---- Pd host: lattice parameter and elastic constants
    r = evaluate({'kind': 'relax', 'structure': structure(np.eye(3) * 3.9, ['Pd'] * 4, BASE * 3.9), 'mode': 'cell', 'fmax': 1e-4})
    a0 = float(np.cbrt(r['volume_A3']))
    cell1, pos1 = fcc(a0, 1)
    d = 0.005
    strains = [np.diag([d, 0, 0]), np.diag([-d, 0, 0]), np.diag([d, d, 0]), np.diag([-d, -d, 0])]
    sh = np.zeros((3, 3))
    sh[0, 1] = sh[1, 0] = d / 2
    strains += [sh, -sh]
    structs = [structure(np.array(cell1) @ (np.eye(3) + e).T, ['Pd'] * 4, pos1 @ (np.eye(3) + e).T) for e in strains]
    cell3, pos3 = fcc(a0, 3)
    structs.append(structure(cell3, ['Pd'] * 108, pos3))
    structs.append(structure(cell1, ['Pd'] * 4, pos1))
    res = evaluate({'kind': 'energy', 'structures': structs})['results']
    E = [x['energy_eV'] for x in res]
    e0, V0 = E[7], a0 ** 3
    c11 = (E[0] + E[1] - 2 * e0) / d ** 2 / V0
    c12 = ((E[2] + E[3] - 2 * e0) / d ** 2 / V0 - 2 * c11) / 2
    c44 = (E[4] + E[5] - 2 * e0) / d ** 2 / V0
    E_Pd108 = E[6]
    E_Pd_atom = e0 / 4
    # ---- H2 reference (static) from the same potential
    h2 = evaluate({'kind': 'relax', 'structure': structure(np.eye(3) * 12, ['H', 'H'], [[0, 0, 0], [0.75, 0, 0]]), 'fmax': 1e-4})
    E_H2 = h2['energy_eV']
    # ---- dilute interstitials: 108-Pd cell at zero pressure (the 4-Pd cell carries a large elastic image interaction)
    dil = {}
    for site, frac in (('O', [0.5, 0, 0]), ('T', [0.25, 0.25, 0.25])):
        s = structure(cell3, ['Pd'] * 108 + ['H'], np.vstack([pos3, np.array(frac) * a0]))
        rr = evaluate({'kind': 'relax', 'structure': s, 'mode': 'volume'})
        dil[site] = {'E_sol': rr['energy_eV'] - E_Pd108 - E_H2 / 2, 'dV': rr['volume_A3'] - V0 * 27, 'relaxed': rr}
    # size check in a 256-Pd cell
    cell4, pos4 = fcc(a0, 4)
    rr = evaluate({'kind': 'relax', 'structure': structure(cell4, ['Pd'] * 256 + ['H'], np.vstack([pos4, [a0 / 2, 0, 0]])), 'mode': 'volume'})
    E_sol_256 = rr['energy_eV'] - 256 * E_Pd_atom - E_H2 / 2
    dV_256 = rr['volume_A3'] - V0 * 64
    E_sol_O = 0.5 * (dil['O']['E_sol'] + E_sol_256)
    dV = 0.5 * (dil['O']['dV'] + dV_256)
    dE_OT = dil['T']['E_sol'] - dil['O']['E_sol']
    # ---- quantum H/D motion in the frozen relaxed host (grid scans; harmonic ZPE badly underestimates the soft O site)
    levels = {}
    for site, half, h in (('O', 0.70, 0.07), ('T', 0.45, 0.045)):
        rel = dil[site]['relaxed']
        c = np.array(rel['positions'][-1])
        g = np.arange(-int(round(half / h)), int(round(half / h)) + 1) * h
        sc = evaluate({'kind': 'scan', 'structure': structure(rel['cell'], rel['symbols'], rel['positions']),
                       'index': len(rel['symbols']) - 1, 'grid': {'center': c.tolist(), 'half_width': half, 'spacing': h}})
        V = np.array(sc['energies_eV']).reshape(len(g), len(g), len(g))
        V -= V[len(g) // 2, len(g) // 2, len(g) // 2]
        for iso, m in (('H', M_H), ('D', M_D)):
            levels[(site, iso)] = schroedinger(V, h, m)
    zpe = {k: v[0] for k, v in levels.items()}
    F = {k: free_energy(v, T300) for k, v in levels.items()}
    sol0K = E_sol_O + zpe[('O', 'H')] - 0.5 * GAS['H2']['zpe_eV']
    dE_OT_zpe = dE_OT + zpe[('T', 'H')] - zpe[('O', 'H')]
    # dilute equilibrium with the gas: x = Z_int exp(mu_gas/2kT) sqrt(p); Z_int sums one O and two T sites per Pd
    lnx, fT = {}, {}
    for iso, gas in (('H', 'H2'), ('D', 'D2')):
        GO = E_sol_O + F[('O', iso)]
        GT = E_sol_O + dE_OT + F[('T', iso)]
        Z = math.exp(-GO / (KB * T300)) + 2 * math.exp(-GT / (KB * T300))
        fT[iso] = 2 * math.exp(-GT / (KB * T300)) / Z
        lnx[iso] = math.log(Z) + 0.5 * GAS[gas]['mu_minus_Estatic_eV'] / (KB * T300)
    sieverts_ratio = math.exp(lnx['H'] - lnx['D'])
    # ---- configurational energetics: random O-site occupations, fully relaxed (zero stress), 108-Pd cells
    rng = np.random.default_rng(7)
    octa = pos3 + np.array([0.5, 0, 0]) * a0
    xs, es = [0.0], [0.0]
    for nh in (4, 8, 13, 13, 19, 19, 27, 27, 40, 54, 70, 85):
        sel = rng.choice(108, nh, replace=False)
        s = structure(cell3, ['Pd'] * 108 + ['H'] * nh, np.vstack([pos3, octa[sel]]))
        rr = evaluate({'kind': 'relax', 'structure': s, 'mode': 'cell'})
        xs.append(nh / 108)
        es.append((rr['energy_eV'] - E_Pd108 - nh * E_H2 / 2) / 108)
    rr = evaluate({'kind': 'relax', 'structure': structure(cell1, ['Pd'] * 4 + ['H'] * 4, np.vstack([pos1, pos1 + [a0 / 2, 0, 0]])), 'mode': 'cell'})
    e1 = (rr['energy_eV'] - 4 * E_Pd_atom - 2 * E_H2) / 4
    xs.append(1.0)
    es.append(e1)
    xs, es = np.array(xs), np.array(es)
    dmix = es - xs * e1
    m = (xs > 0) & (xs < 1)
    deg = 4
    A = np.array([xs[m] * (1 - xs[m]) * (1 - 2 * xs[m]) ** j for j in range(deg)]).T
    L = np.linalg.lstsq(A, dmix[m], rcond=None)[0]

    def mix(z):
        return z * (1 - z) * sum(L[j] * (1 - 2 * z) ** j for j in range(deg))
    z = np.linspace(0.002, 0.998, 4000)
    hh = 1e-4
    d2 = (mix(z + hh) + mix(z - hh) - 2 * mix(z)) / hh ** 2
    Tc_inc = float(np.max(-d2 * z * (1 - z) / KB))
    # coherent: Cahn, cubic; directional modulus minimised over directions; Vegard coefficient from the dilute limit
    eta = dV / (3 * V0 / 4)
    J = np.linspace(0, 1 / 3, 200)       # l^2 m^2 + m^2 n^2 + n^2 l^2 spans [0, 1/3]
    Yn = 0.5 * (c11 + 2 * c12) * (3 - (c11 + 2 * c12) / (c11 + 2 * (2 * c44 - c11 + c12) * J))
    Ymin = float(Yn.min())
    coh = 2 * eta ** 2 * Ymin * (V0 / 4)
    Tc_coh = float(np.max(-(d2 + coh) * z * (1 - z) / KB))
    # incoherent plateau at 300 K: common tangent of f(x) = e(x) + x F_vib(H, O) - T s(x)
    G = F[('O', 'H')]
    zz = np.linspace(1e-6, 1 - 1e-6, 200001)
    f = mix(zz) + zz * (e1 + G) + KB * T300 * (zz * np.log(zz) + (1 - zz) * np.log(1 - zz))
    hull = []
    for p in zip(zz, f):
        p = np.array(p)
        while len(hull) >= 2 and np.cross(hull[-1] - hull[-2], p - hull[-2]) <= 0:
            hull.pop()
        hull.append(p)
    hull = np.array(hull)
    j = int(np.argmax(np.diff(hull[:, 0])))
    mu = (hull[j + 1, 1] - hull[j, 1]) / (hull[j + 1, 0] - hull[j, 0])
    log10p = 2 * (mu - 0.5 * GAS['H2']['mu_minus_Estatic_eV']) / (KB * T300) / math.log(10)
    xa, xb = hull[j, 0], hull[j + 1, 0]

    submission = {'solution_energy_0K_eV': sol0K, 'OT_difference_static_eV': dE_OT, 'OT_difference_with_zpe_eV': dE_OT_zpe,
                  'lattice_expansion_A3_per_H': dV, 'tetrahedral_fraction_H_300K': fT['H'],
                  'sieverts_ratio_H_over_D_300K': sieverts_ratio, 'log10_plateau_pressure_H_300K_bar': log10p,
                  'Tc_incoherent_K': Tc_inc, 'Tc_coherent_K': Tc_coh, 'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    harm = {}
    (APP/'report.md').write_text(f"""# Pd-H miscibility gap and isotope effect from the released ACE (reference solution)

Host: a0 = {a0:.4f} A, C11/C12/C44 = {c11*160.2177:.0f}/{c12*160.2177:.0f}/{c44*160.2177:.0f} GPa. H2 reference: relaxed molecule of the same
potential (static energy {E_H2:.4f} eV); gas-phase free energies and ZPE from the supplied thermochemistry.

Dilute hydrogen: interstitials were relaxed at zero pressure in 108- and 256-Pd cells; small cells are not usable because the
elastic image interaction of the H dilatation lowers the energy per H strongly at zero pressure (and raises it at fixed volume).
Static solution energy (O) {E_sol_O:.4f} eV/H, O->T static {dE_OT:.4f} eV, lattice expansion {dV:.2f} A^3/H (256 cell {dV_256:.2f}).

Quantum H motion: the octahedral site is very soft and anharmonic. Harmonic small-displacement curvature would give a far too small
zero-point energy; instead the nuclear Schroedinger equation was solved on 3D frozen-host energy grids around each site. ZPE(O,H) =
{zpe[('O','H')]*1000:.1f} meV, ZPE(T,H) = {zpe[('T','H')]*1000:.1f} meV, ZPE(O,D) = {zpe[('O','D')]*1000:.1f} meV; first O-site excitation
{(levels[('O','H')][1]-levels[('O','H')][0])*1000:.1f} meV (H). With ZPE the octahedral preference grows to {dE_OT_zpe*1000:.0f} meV; T-site fraction at 300 K:
H {fT['H']:.4f}, D {fT['D']:.4f}. Thermal (excited-state) free energies at 300 K matter for the isotope effect. Sieverts ratio K_H/K_D =
{sieverts_ratio:.2f}: Pd dissolves more H than D at equal pressure (inverse isotope effect), because the in-metal ZPE difference H-D
is smaller than half the gas-phase one.

Phase separation: fixed-volume pair energies are only weakly attractive and the fixed-volume mixing energy is convex (no gap);
the attraction that drives separation is the volume-mediated (elastic) interaction, visible only when the cell relaxes. At high x the
H-H interaction turns repulsive (filling of the Pd d-band holes), which closes the gap near x ~ 0.3 at 300 K. Random-occupation
zero-stress energies (108-Pd cells) fitted with a Redlich-Kister series give the incoherent mean-field Tc = {Tc_inc:.0f} K; the coherent
spinodal (Cahn, Ymin = {Ymin*160.2177:.0f} GPa along <100>, eta = {eta:.4f}) lowers it to {Tc_coh:.0f} K. Incoherent plateau at 300 K:
10^{log10p:.2f} bar between x = {xa:.3f} and {xb:.3f}.

Not supported: the absolute Tc depends on the mean-field (ideal-entropy) model and on PBE; short-range order, H vibrational changes
with x and anharmonic Pd phonons are neglected; the isotope dependence of Tc is not resolved by this model.
""")
    lab('submit', submission)
    for _ in range(120):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
