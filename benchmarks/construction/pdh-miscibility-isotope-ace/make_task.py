"""Generate catalog, lab data, hidden targets, ablation and robustness values for pdh-miscibility-isotope-ace.

Inputs: ref/*.json and ref/grid_*.npy (reference.py), $CACHE/work/PdH_ace_ctilde.json (convert_potential.py) and
$CACHE/work/dft_subset.json (extract_dft.py). Run inside the runtime image (needs scipy):
  docker run --rm --user $(id -u):$(id -g) -v $PWD:/c -v $CACHE/work:/w -v $TASK:/t ai-theorist-benchmark-runtime:20261007 \
      python3 -B /c/make_task.py
"""
import glob
import json
import math
import shutil
import sys
from pathlib import Path
import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import eigsh

HERE = Path(__file__).resolve().parent
TASK = Path('/t') if Path('/t').exists() else HERE.parents[1]/'tasks/pdh-miscibility-isotope-ace'
WORK = Path('/w')
REF = HERE/'ref'
KB = 8.617333262e-5
HBAR = 0.0646541
M_H, M_D = 1.00782503, 2.01410178
T = 300.0
EVA3_GPA = 160.21766208
# Gas phase at 300 K, 1 bar from experimental spectroscopic constants (rigid rotor with nuclear-spin statistics,
# anharmonic ZPE, ideal-gas translation): mu - E_static per molecule, and ZPE.
GAS = {'H2': {'mu_minus_Estatic_eV': -0.0483, 'zpe_eV': 0.2701},
       'D2': {'mu_minus_Estatic_eV': -0.1700, 'zpe_eV': 0.1918}}
BUDGET = 30

base = json.load(open(REF/'base.json'))
dil = json.load(open(REF/'dilute.json'))
grids = json.load(open(REF/'grids.json'))
rand = sum([json.load(open(f)) for f in sorted(glob.glob(str(REF/'random_*.json')))], [])
ordered = json.load(open(REF/'ordered.json'))


def schro(V, h, mass, k=30):
    n = V.shape[0]
    D = sps.diags([np.ones(n - 1), -2 * np.ones(n), np.ones(n - 1)], [-1, 0, 1]) / h ** 2
    I = sps.identity(n)
    L = sps.kron(sps.kron(D, I), I) + sps.kron(sps.kron(I, D), I) + sps.kron(sps.kron(I, I), D)
    H = -HBAR ** 2 / (2 * mass) * L + sps.diags(V.ravel())
    return np.sort(eigsh(H.tocsc(), k=k, sigma=V.min() - 0.01, which='LM', return_eigenvectors=False))


def Fvib(levels):
    w = np.asarray(levels)
    return w[0] - KB * T * math.log(np.exp(-(w - w[0]) / (KB * T)).sum())


def harmonic(kdiag, mass):
    w = [HBAR * math.sqrt(k / mass) for k in kdiag]
    zpe = sum(0.5 * x for x in w)
    return zpe, zpe + sum(KB * T * math.log(1 - math.exp(-x / (KB * T))) for x in w)


def vib_sets():
    """Vibrational energies per (site, isotope): anharmonic levels from the grids, plus variants."""
    out = {}
    for variant, sub, crops in (('ref', 1, (0, 0)), ('coarse_h', 2, (0, 0)), ('small_box', 1, (3, 1))):
        lv = {}
        for site, crop in zip(('O', 'T'), crops):
            V = np.load(REF/f'grid_{site}.npy')
            h = grids[site]['h']
            if crop:
                V = V[crop:-crop, crop:-crop, crop:-crop]
            V = V[::sub, ::sub, ::sub]
            for iso, m in (('H', M_H), ('D', M_D)):
                lv[(site, iso)] = schro(V, h * sub, m)
        out[variant] = lv
    return out


def isotope_quantities(E_sol, dE_OT, E0, F):
    """E0[(site, iso)] ground-state vibrational energy; F[(site, iso)] free energy at 300 K."""
    res = {'solution_energy_0K_eV': E_sol + E0[('O', 'H')] - 0.5 * GAS['H2']['zpe_eV'],
           'OT_difference_with_zpe_eV': dE_OT + E0[('T', 'H')] - E0[('O', 'H')]}
    lnx, fT = {}, {}
    for iso, gas in (('H', 'H2'), ('D', 'D2')):
        GO, GT = E_sol + F[('O', iso)], E_sol + dE_OT + F[('T', iso)]
        Z = math.exp(-GO / (KB * T)) + 2 * math.exp(-GT / (KB * T))
        fT[iso] = 2 * math.exp(-GT / (KB * T)) / Z
        lnx[iso] = math.log(Z) + 0.5 * GAS[gas]['mu_minus_Estatic_eV'] / (KB * T)
    res['tetrahedral_fraction_H_300K'] = fT['H']
    res['sieverts_ratio_H_over_D_300K'] = math.exp(lnx['H'] - lnx['D'])
    return res, fT


def rk_fit(xs, es, deg):
    xs, es = np.asarray(xs), np.asarray(es)
    e1 = es[xs == 1].mean()
    d = es - xs * e1
    m = (xs > 0) & (xs < 1)
    A = np.array([xs[m] * (1 - xs[m]) * (1 - 2 * xs[m]) ** j for j in range(deg)]).T
    L = np.linalg.lstsq(A, d[m], rcond=None)[0]
    return L, e1


def mixf(L):
    return lambda z: z * (1 - z) * sum(L[j] * (1 - 2 * z) ** j for j in range(len(L)))


def tc(L, coh=0.0):
    mix = mixf(L)
    z = np.linspace(0.002, 0.998, 4000)
    h = 1e-4
    d2 = (mix(z + h) + mix(z - h) - 2 * mix(z)) / h ** 2
    return float(max(0.0, np.max(-(d2 + coh) * z * (1 - z) / KB)))


def plateau(L, e1, G):
    mix = mixf(L)
    zz = np.linspace(1e-6, 1 - 1e-6, 200001)
    f = mix(zz) + zz * (e1 + G) + KB * T * (zz * np.log(zz) + (1 - zz) * np.log(1 - zz))
    hull = []
    for p in zip(zz, f):
        p = np.array(p)
        while len(hull) >= 2 and np.cross(hull[-1] - hull[-2], p - hull[-2]) <= 0:
            hull.pop()
        hull.append(p)
    hull = np.array(hull)
    j = int(np.argmax(np.diff(hull[:, 0])))
    if hull[j + 1, 0] - hull[j, 0] < 1e-3:       # no two-phase region: chemical potential at x = 0.3
        mu = float((f[60001] - f[59999]) / (zz[60001] - zz[59999]))
    else:
        mu = (hull[j + 1, 1] - hull[j, 1]) / (hull[j + 1, 0] - hull[j, 0])
    return 2 * (mu - 0.5 * GAS['H2']['mu_minus_Estatic_eV']) / (KB * T) / math.log(10)


def Ymin(c11, c12, c44, direction='min'):
    J = np.linspace(0, 1 / 3, 400)
    Y = 0.5 * (c11 + 2 * c12) * (3 - (c11 + 2 * c12) / (c11 + 2 * (2 * c44 - c11 + c12) * J))
    return float({'min': Y.min(), '111': Y[-1], '100': Y[0]}[direction])


def main():
    a0 = base['a0']
    v0 = a0 ** 3 / 4
    c11, c12, c44 = base['Cij_eV_A3_d0.005']
    E_sol = 0.5 * (dil['3_O_volume']['E_sol'] + dil['4_O_volume']['E_sol'])
    dV = 0.5 * (dil['3_O_volume']['dV'] + dil['4_O_volume']['dV'])
    dE_OT = 0.5 * (dil['3_T_volume']['E_sol'] - dil['3_O_volume']['E_sol'] + dil['4_T_volume']['E_sol'] - dil['4_O_volume']['E_sol'])
    vs = vib_sets()
    lv = vs['ref']
    E0 = {k: v[0] for k, v in lv.items()}
    F = {k: Fvib(v) for k, v in lv.items()}
    iso, fT = isotope_quantities(E_sol, dE_OT, E0, F)
    # configurational energetics
    xs = np.array([0.0] + [r['x'] for r in rand])
    ec = np.array([0.0] + [r['e_cell'] for r in rand])
    ef = np.array([0.0] + [r['e_fixed'] for r in rand])
    ev = np.array([0.0] + [r['e_volume'] for r in rand])
    Ns = np.array([0] + [r['N'] for r in rand])
    eta = dV / (3 * v0)
    coh = 2 * eta ** 2 * Ymin(c11, c12, c44) * v0
    L, e1 = rk_fit(xs, ec, 4)
    G = F[('O', 'H')]
    ref = {'solution_energy_0K_eV': iso['solution_energy_0K_eV'], 'OT_difference_static_eV': dE_OT,
           'OT_difference_with_zpe_eV': iso['OT_difference_with_zpe_eV'], 'lattice_expansion_A3_per_H': dV,
           'tetrahedral_fraction_H_300K': iso['tetrahedral_fraction_H_300K'],
           'sieverts_ratio_H_over_D_300K': iso['sieverts_ratio_H_over_D_300K'],
           'log10_plateau_pressure_H_300K_bar': plateau(L, e1, G), 'Tc_incoherent_K': tc(L), 'Tc_coherent_K': tc(L, coh)}
    tol = {'solution_energy_0K_eV': ['abs', 0.012], 'OT_difference_static_eV': ['abs', 0.010],
           'OT_difference_with_zpe_eV': ['abs', 0.015], 'lattice_expansion_A3_per_H': ['rel', 0.08],
           'tetrahedral_fraction_H_300K': ['log', 2.0], 'sieverts_ratio_H_over_D_300K': ['rel', 0.20],
           'log10_plateau_pressure_H_300K_bar': ['abs', 0.4], 'Tc_incoherent_K': ['rel', 0.15], 'Tc_coherent_K': ['rel', 0.20]}

    # ---- robustness (defensible alternatives)
    R = {}
    for name in ('coarse_h', 'small_box'):
        lvv = vs[name]
        r_iso, _ = isotope_quantities(E_sol, dE_OT, {k: v[0] for k, v in lvv.items()}, {k: Fvib(v) for k, v in lvv.items()})
        r_iso['log10_plateau_pressure_H_300K_bar'] = plateau(L, e1, Fvib(lvv[('O', 'H')]))
        R['grid_' + name] = r_iso
    e108 = (Ns == 108) | (Ns == 0) | (xs == 1)
    e32 = (Ns == 32) | (Ns == 0)
    for name, sel, deg, energies in (('rk3_all', None, 3, ec), ('rk5_all', None, 5, ec), ('rk4_108_only', e108, 4, ec),
                                     ('rk4_32_only', e32, 4, ec), ('rk4_isotropic_volume', None, 4, ev)):
        xx, yy = (xs, energies) if sel is None else (xs[sel], energies[sel])
        Lr, e1r = rk_fit(xx, yy, deg)
        R[name] = {'Tc_incoherent_K': tc(Lr), 'Tc_coherent_K': tc(Lr, coh), 'log10_plateau_pressure_H_300K_bar': plateau(Lr, e1r, G)}
    R['dilute_108_only'] = {'lattice_expansion_A3_per_H': dil['3_O_volume']['dV'],
                            'OT_difference_static_eV': dil['3_T_volume']['E_sol'] - dil['3_O_volume']['E_sol'],
                            'solution_energy_0K_eV': dil['3_O_volume']['E_sol'] + E0[('O', 'H')] - 0.5 * GAS['H2']['zpe_eV']}
    R['dilute_256_fixed_volume'] = {'solution_energy_0K_eV': dil['4_O_positions']['E_sol'] + E0[('O', 'H')] - 0.5 * GAS['H2']['zpe_eV'],
                                    'lattice_expansion_A3_per_H': dil['4_O_positions']['p_GPa'] / EVA3_GPA * 256 * v0 / ((c11 + 2 * c12) / 3)}
    R['elastic_constants_d0.01'] = {'Tc_coherent_K': tc(L, 2 * eta ** 2 * Ymin(*base['Cij_eV_A3_d0.01']) * v0)}

    # ---- ablations (realistic wrong paths computed from the reference data)
    A = {}
    Lf, e1f = rk_fit(xs, ef, 4)
    A['fixed_volume_configurational'] = {'Tc_incoherent_K': tc(Lf), 'Tc_coherent_K': tc(Lf, coh), 'log10_plateau_pressure_H_300K_bar': plateau(Lf, e1f, G)}
    A['no_coherency'] = {'Tc_coherent_K': ref['Tc_incoherent_K']}
    A['coherency_111'] = {'Tc_coherent_K': tc(L, 2 * eta ** 2 * Ymin(c11, c12, c44, '111') * v0)}
    sm = {k: dil[f'1_{k}_volume']['E_sol'] for k in ('O', 'T')}
    A['small_cell_dilute'] = {'solution_energy_0K_eV': sm['O'] + E0[('O', 'H')] - 0.5 * GAS['H2']['zpe_eV'],
                              'OT_difference_static_eV': sm['T'] - sm['O'],
                              'OT_difference_with_zpe_eV': sm['T'] - sm['O'] + E0[('T', 'H')] - E0[('O', 'H')],
                              'lattice_expansion_A3_per_H': dil['1_O_volume']['dV']}
    smf = {k: dil[f'1_{k}_positions']['E_sol'] for k in ('O', 'T')}
    A['small_cell_fixed_volume_dilute'] = {'solution_energy_0K_eV': smf['O'] + E0[('O', 'H')] - 0.5 * GAS['H2']['zpe_eV'],
                                           'OT_difference_static_eV': smf['T'] - smf['O'],
                                           'OT_difference_with_zpe_eV': smf['T'] - smf['O'] + E0[('T', 'H')] - E0[('O', 'H')]}
    hz, hf = {}, {}
    for site in ('O', 'T'):
        for isot, m in (('H', M_H), ('D', M_D)):
            hz[(site, isot)], hf[(site, isot)] = harmonic(grids[site]['k_diag'], m)
    h_iso, _ = isotope_quantities(E_sol, dE_OT, hz, hf)
    h_iso['log10_plateau_pressure_H_300K_bar'] = plateau(L, e1, hf[('O', 'H')])
    A['harmonic_zpe'] = h_iso
    z_iso, _ = isotope_quantities(E_sol, dE_OT, E0, E0)
    z_iso['log10_plateau_pressure_H_300K_bar'] = plateau(L, e1, E0[('O', 'H')])
    A['zpe_only_no_thermal'] = z_iso
    zero = {k: 0.0 for k in E0}
    s_iso, _ = isotope_quantities(E_sol, dE_OT, zero, zero)
    A['static_site_preference'] = {'tetrahedral_fraction_H_300K': s_iso['tetrahedral_fraction_H_300K'],
                                   'OT_difference_with_zpe_eV': dE_OT}
    u_el = (c11 + 2 * c12) / 3 * dV ** 2 / v0
    A['elastic_only_wagner_horner'] = {'Tc_incoherent_K': u_el / (4 * KB), 'Tc_coherent_K': (u_el - coh) / (4 * KB)}
    xo = np.array([0.0] + [r['x'] for r in ordered])
    eo = np.array([0.0] + [r['e_cell'] for r in ordered])
    Lo, e1o = rk_fit(xo, eo, 3)
    A['ordered_small_cells'] = {'Tc_incoherent_K': tc(Lo), 'Tc_coherent_K': tc(Lo, coh), 'log10_plateau_pressure_H_300K_bar': plateau(Lo, e1o, G)}
    A['tetrahedral_site_as_dilute'] = {'solution_energy_0K_eV': dil['3_T_volume']['E_sol'] + E0[('T', 'H')] - 0.5 * GAS['H2']['zpe_eV'],
                                       'lattice_expansion_A3_per_H': dil['3_T_volume']['dV']}

    def check(vals):
        out = {}
        for f, v in vals.items():
            kind, t = tol[f]
            r = ref[f]
            ok = abs(v - r) <= t if kind == 'abs' else (abs(v / r - 1) <= t if kind == 'rel' else (v > 0 and abs(math.log(v / r)) <= math.log(t)))
            out[f] = (round(v, 5), bool(ok))
        return out
    print('reference', json.dumps({k: round(v, 5) for k, v in ref.items()}, indent=1))
    for n, v in R.items():
        print('ROBUST', n, check(v))
    for n, v in A.items():
        print('ABLATE', n, check(v))

    # ---- write the task files
    lab = TASK/'environment/lab'
    for f in ('ace_numpy.py', 'relax.py'):
        shutil.copy(HERE/f, lab/f)
    shutil.copy(WORK/'PdH_ace_ctilde.json', lab/'PdH_ace_ctilde.json')
    dft = json.load(open(WORK/'dft_subset.json'))
    # drop the few records whose coordinate text happens to contain a graded number (release leak scan)
    needles = set()
    for v in ref.values():
        needles |= {t for t in {f'{v:.4g}', f'{v:.5g}', f'{v:.6g}', repr(v)} if len(t.replace('-', '').replace('.', '')) >= 4}

    def rec_text(r):
        return json.dumps(r['cell']) + json.dumps(np.round(np.array(r['positions']) @ np.linalg.inv(np.array(r['cell'])), 7).tolist())
    dropped = [r['id'] for r in dft if any(n in rec_text(r) for n in needles)]
    dft = [r for r in dft if r['id'] not in dropped]
    print('DFT records dropped for leak-scan collisions:', dropped)
    json.dump(dft, open(lab/'dft_subset.json', 'w'))
    hidden = {'reference': ref, 'tolerance': tol}
    (lab/'private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    (lab/'private/hidden_values.json').write_text(json.dumps(list(ref.values())) + '\n')
    catalog = {
        'potential': 'ACE PdH_PBE_v03a (atomic cluster expansion fitted to VASP PBE data; Qamar, Ngoipala, Mrovec, Drautz et al. 2026, '
                     'Zenodo 20591788). Evaluated exactly by the lab.',
        'units': {'energy': 'eV', 'length': 'Angstrom', 'force': 'eV/Angstrom', 'stress': 'eV/Angstrom^3 (ASE sign: pressure = -trace/3)'},
        'energy_zero': 'Isolated atoms have zero energy (the potential is fitted to DFT energies minus free-atom reference energies). '
                       'DFT records below use the same zero.',
        'structure_format': {'symbols': ['Pd', 'H'], 'positions': [[0.0, 0.0, 0.0], [1.97, 0.0, 0.0]],
                             'cell': [[3.94, 0.0, 0.0], [0.0, 3.94, 0.0], [0.0, 0.0, 3.94]],
                             'notes': 'Cartesian Angstrom, periodic in all directions, at most 300 atoms '
                                      '(DFT records below list fractional positions: Cartesian = fractional @ cell)'},
        'requests': {
            'energy': {'example': {'kind': 'energy', 'structures': ['<structure>', '...']}, 'cost': 1,
                       'returns': 'energy, forces, stress, pressure and volume of up to 8 structures'},
            'relax': {'example': {'kind': 'relax', 'structure': '<structure>', 'mode': 'positions | volume | cell', 'fixed': [], 'fmax': 0.001},
                      'cost': 1, 'returns': 'relaxed structure and its properties; volume = isotropic cell scaling to zero pressure, '
                                            'cell = full cell relaxation to zero stress, fixed = atom indices held in place'},
            'scan': {'example': {'kind': 'scan', 'structure': '<structure>', 'index': 0, 'points': [[0.0, 0.0, 0.0]]},
                     'example_grid': {'kind': 'scan', 'structure': '<structure>', 'index': 0,
                                      'grid': {'center': [1.97, 0.0, 0.0], 'half_width': 0.5, 'spacing': 0.05}}, 'cost': 2,
                     'returns': 'total energy with atom `index` moved to each point, every other atom frozen. Either an explicit '
                                'list of points within a 4 A box, or a cubic grid center + (i, j, k) * spacing with |i|, |j|, |k| <= '
                                'floor(half_width / spacing), returned in C order (x slowest). At most 10000 points; requests are '
                                'limited to 64 kB, so use the grid form for large scans'},
            'dft_reference': {'example': {'kind': 'dft_reference', 'ids': ['DFT-001']}, 'cost': 1,
                              'returns': 'archived DFT (VASP PBE) energy, forces and stress of up to 4 records listed below'}},
        'budget': BUDGET,
        'gas_phase_300K_1bar': {g: dict(v, note='mu(300 K, 1 bar) - E_static(molecule) per molecule, and zero-point energy, from experimental '
                                               'spectroscopic data (ideal gas, rotor with nuclear-spin statistics)') for g, v in GAS.items()},
        'masses_amu': {'H': M_H, 'D': M_D},
        'dft_reference_records': [dict({k: r[k] for k in ('id', 'label', 'split', 'formula', 'cell', 'symbols')},
                                       fractional_positions=np.round(np.array(r['positions']) @ np.linalg.inv(np.array(r['cell'])), 7).tolist())
                                  for r in dft],
    }
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog) + '\n')
    for f in ('environment/lab/public.json', 'tests/public.json'):
        p = json.loads((TASK/f).read_text())
        p['solver_budget'] = BUDGET
        if f.startswith('environment'):
            p['requests'] = {k: {'kind': k, 'cost': v['cost']} for k, v in catalog['requests'].items()}
            p['catalog'] = '/app/catalog.json'
        (TASK/f).write_text(json.dumps(p, indent=1) + '\n')
    json.dump({'reference': ref, 'robustness': R, 'ablations': A, 'tolerance': tol,
               'extra': {'a0': a0, 'Cij_GPa': [c * EVA3_GPA for c in (c11, c12, c44)], 'E_sol_static': E_sol, 'eta': eta,
                         'coh_eV': coh, 'Ymin_GPa': Ymin(c11, c12, c44) * EVA3_GPA, 'u_el_eV': u_el,
                         'E0': {f'{a}_{b}': v for (a, b), v in E0.items()}, 'F300': {f'{a}_{b}': v for (a, b), v in F.items()},
                         'fT': fT, 'harmonic_zpe': {f'{a}_{b}': v for (a, b), v in hz.items()}}},
              open(TASK/'solution/ablation_values.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
