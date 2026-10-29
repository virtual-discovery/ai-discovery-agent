"""Generate the staged, lab, private and ablation data for tb2j-groundstate-cuo-multiferroic.

    python3 -B benchmarks/construction/tb2j-groundstate-cuo-multiferroic/make_task.py     (inside the runtime image)

Reads the extracted archive members (fetch.py) and writes:
  environment/agent/catalog.json, environment/lab/data.json, environment/lab/private/hidden.json,
  environment/lab/private/hidden_values.json, solution/ablation_values.json and calibration.json (this folder).
It asserts that the authors' archived h(k) spectra of the FM/ligand and AFM/ligand calculations are reproduced
from the archived pair lists, and that the FM/ligand k_min equals the paper's SI value.
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cuomodel as M  # noqa: E402

TASK = HERE.parents[1]/'tasks/tb2j-groundstate-cuo-multiferroic'
WINDOWS = [(0, 4), (4, 6), (6, 8), (8, 10), (10, 12), (12, 14), (14, 17), (17, 20), (20, 25), (25, 33)]
COSTS = [1, 1, 1, 1, 1, 1, 2, 2, 3, 3]
BUDGET = 24
SPIN = 0.5
THRESHOLDS = {'kmin_ligand_rlu': 0.006, 'kmin_uncorrected_rlu': 0.003, 'converged_q_rlu': 0.006,
              'chain_J_rel': 0.04, 'magnon_mean_rel': 0.05, 'magnon_max_rel': 0.04}
LABEL = {'FM/ligand': 'FM reference, spin-ligand correction on', 'AFM/ligand': 'AFM reference, spin-ligand correction on',
         'FM/no_ligand': 'FM reference, spin-ligand correction off', 'AFM/no_ligand': 'AFM reference, spin-ligand correction off'}


def wname(w):
    return f'{w[0]}-{w[1]}'


def authors_path(A, Y, X):
    ks = []
    for r in range(300):
        if r < 144:
            ks.append(A + (Y - A)*r/144)
        elif r < 219:
            ks.append(Y*(1 - (r - 144)/75))
        else:
            ks.append(X*(r - 219)/80)
    return np.array(ks)


def check_author_dispersions(data):
    """Fit the unarchived path nodes A and X on the FM/ligand file, then reproduce FM/ligand and AFM/ligand."""
    fm = M.Native('FM/ligand', data=data['FM/ligand'])
    _, f = M.dispersion('FM/ligand')
    Y = np.array([0.5, 0, -0.5])            # paper SI 'Y = (1/2, 0, 1/2)' in the authors' axes

    def res(p):
        ks = authors_path(p[:3], Y, p[3:])
        return (np.sort(np.linalg.eigvalsh(fm.h(ks[::3]))/1000, axis=1) - np.sort(f[::3], axis=1)).ravel()
    p = least_squares(res, [-0.39, 0.5, -0.39, 0.39, 0, 0.39]).x
    afm = M.Native('AFM/ligand', data=data['AFM/ligand'])
    _, g = M.dispersion('AFM/ligand')
    # full-precision cells: FM from ligand/FM/TB2J.yaml (exchange.out rounds to 1e-3 A), AFM from CuO_AFM.xsf
    fm_cell = np.array([[4.9791366, -7.6e-06, -0.1677558], [-5.5e-06, 3.2000538, -2.1e-06], [-1.0605033, -1.8e-06, 5.0657527]])
    xsf = (M.RAW/'ligand/AFM/CuO_AFM.xsf').read_text().split('\n')
    afm_cell = np.array([[float(x) for x in xsf[xsf.index('PRIMVEC') + k].split()] for k in (1, 2, 3)])
    assert np.abs(afm_cell - afm.cell).max() < 1e-3
    to_sc = np.linalg.inv(fm_cell).T @ afm_cell.T
    # the archived AFM spectrum uses the FM path's Cartesian k-vectors in the AFM calculation's own frame

    def res2(p):
        ks = authors_path(p[:3], Y, p[3:])
        a = np.sort(np.linalg.eigvalsh(fm.h(ks[::3]))/1000, axis=1) - np.sort(f[::3], axis=1)
        b = np.sort(np.linalg.eigvalsh(afm.h(ks[::10] @ to_sc))/1000, axis=1) - np.sort(g[::10], axis=1)
        return np.r_[a.ravel(), b.ravel()]
    p = least_squares(res2, p).x
    path = authors_path(p[:3], Y, p[3:])
    dev_fm = np.abs(np.sort(np.linalg.eigvalsh(fm.h(path))/1000, axis=1) - np.sort(f, axis=1)).max()
    k_sc = path @ to_sc
    dev_afm = np.abs(np.sort(np.linalg.eigvalsh(afm.h(k_sc))/1000, axis=1) - np.sort(g, axis=1)).max()
    dev_afm_yx = np.abs(np.sort(np.linalg.eigvalsh(afm.h(k_sc[144:]))/1000, axis=1) - np.sort(g[144:], axis=1)).max()
    # floor ~1e-5 eV from the 4-decimal J in exchange.out; the A-Y segment of the AFM file deviates by <= 2e-4 eV
    # (0.1 % of the bandwidth) because the unarchived node A is only fitted
    assert dev_fm < 3e-5 and dev_afm_yx < 3e-5 and dev_afm < 2e-4, (dev_fm, dev_afm_yx, dev_afm)
    return {'A_rlu': p[:3].round(4).tolist(), 'X_rlu': p[3:].round(4).tolist(), 'max_dev_eV_FM_ligand': float(dev_fm),
            'max_dev_eV_AFM_ligand': float(dev_afm), 'max_dev_eV_AFM_ligand_Y_Gamma_X': float(dev_afm_yx)}


def moment_wavevector(info):
    """Fourier components of the reference moment pattern on the supercell's reciprocal points (conventional r.l.u.)."""
    G = np.array(np.meshgrid(range(-2, 3), range(-2, 3), range(-2, 3), indexing='ij')).reshape(3, -1).T
    qs = np.unique(np.round((G @ np.linalg.inv(info['Mi']).T) % 1, 6), axis=0)
    out = []
    for q in qs:
        amp = max(abs(np.sum(info['spins'][info['types'] == t]*np.exp(2j*np.pi*info['fconv'][info['types'] == t] @ q)))
                  for t in (0, 1))
        if amp > 1e-3:
            out.append((q.tolist(), float(amp)))
    return out


def main():
    data = {k: M.load(k) for k in M.SETS}
    path_check = check_author_dispersions(data)

    # ---- references ----
    kmin_lig, lam_lig = M.Folded(M.folded('FM/ligand', data=data['FM/ligand'])).lt()
    kmin_nol, lam_nol = M.Folded(M.folded('FM/no_ligand', data=data['FM/no_ligand'])).lt()
    assert M.q_distance(kmin_lig, [0.491, 0.0, -0.474]) < 0.002          # paper SI k_min, authors' c-axis sign
    assert M.q_distance(kmin_nol, [0.0, 0.5, 0.0]) < 0.006               # paper SI k_min without the correction
    info_afm = M.cu_info(data['AFM/ligand'])
    comps = moment_wavevector(info_afm)
    assert len(comps) == 1 and M.q_distance(comps[0][0], [0.5, 0, 0.5]) < 1e-9, comps
    native = M.Native('AFM/ligand', data=data['AFM/ligand'])
    K0 = native.K([0, 0, 0])[0].real
    w, v = np.linalg.eigh(K0)
    assert abs(abs(v[:, -1] @ native.spins)/4 - 1) < 1e-4    # AF1 is the LT optimum of its own pair list in the supercell
    stability = {}
    min_eigs = {}
    for key in M.SETS:
        e = M.min_h_eig(M.Native(key, data=data[key]))
        min_eigs[key] = e
        stability[key] = bool(e > -1e-3)
    assert stability == {'FM/ligand': False, 'AFM/ligand': True, 'FM/no_ligand': False, 'AFM/no_ligand': False}, min_eigs
    cj = M.chain_J('AFM/ligand', data=data['AFM/ligand'])
    chain = -8*cj[0]                                    # S = 1/2, each bond once, J > 0 antiferromagnetic
    mag = M.magnon_stats(native, n=8, spin=SPIN)
    assert mag['max_imag'] < 1e-8 and mag['min_h_eig'] > -1e-6
    hidden = {
        'targets': {'kmin_FM_reference_rlu': {'ligand_corrected': M.canonical(kmin_lig).tolist(), 'uncorrected': M.canonical(kmin_nol).tolist()},
                    'stable_at_own_reference': stability, 'converged_propagation_vector_rlu': [0.5, 0.0, 0.5],
                    'chain_exchange_meV': chain, 'mean_magnon_energy_meV': mag['mean'],
                    'magnon_band_maximum_meV': mag['max']},
        'thresholds': THRESHOLDS,
        'reference_details': {'chain_bond_TB2J_meV': cj, 'lt_top_eigenvalue_meV': {'FM/ligand': lam_lig, 'FM/no_ligand': lam_nol},
                              'min_h_eigenvalue_meV': min_eigs, 'AF1_energy_meV_per_Cu': native.energy(),
                              'author_path_check': path_check, 'magnon_details': mag}}

    # ---- ablations (realistic wrong paths, from the same data) ----
    sites = M.af1_sites()
    abl, rob = {}, {}
    fm_af1 = M.magnon_stats(M.FoldedOnSupercell(M.folded('FM/ligand', data=data['FM/ligand']), sites), spin=SPIN)
    abl['fm_reference_as_final'] = {'converged_q': kmin_lig.tolist(), 'chain': -8*M.chain_J('FM/ligand', data=data['FM/ligand'])[0],
                                    'mean': fm_af1['mean'], 'max': fm_af1['max']}
    abl['uncorrected_exchange_as_final'] = {'chain': -8*M.chain_J('FM/no_ligand', data=data['FM/no_ligand'])[0],
                                            'stability': {**stability, 'AFM/no_ligand': True}}
    k6l = M.Folded(M.folded('FM/ligand', rmax=6.0, data=data['FM/ligand'])).lt()[0]
    k6n = M.Folded(M.folded('FM/no_ligand', rmax=6.0, data=data['FM/no_ligand'])).lt()[0]
    abl['exchange_truncated_at_6A'] = {'ligand_corrected': k6l.tolist(), 'uncorrected': k6n.tolist()}
    abl['commensurate_approximant_as_kmin'] = {'ligand_corrected': [0.5, 0.0, 0.5], 'uncorrected': [0.0, 0.5, 0.0]}
    sym = M.folded('AFM/ligand', data=data['AFM/ligand'])
    q_sym, _ = M.Folded(sym).lt()
    sym_af1 = M.magnon_stats(M.FoldedOnSupercell(sym, sites), spin=SPIN)
    abl['symmetrized_afm_exchange'] = {'converged_q': q_sym.tolist(), 'mean': sym_af1['mean'], 'max': sym_af1['max'],
                                       'af1_minus_spiral_meV_per_Cu': M.Folded(sym).top(q_sym) - M.Folded(sym).top([0.5, 0, 0.5])}
    abl['hamiltonian_eigenvalues_as_magnons'] = {'mean': mag['h_mean'], 'max': mag['h_max']}
    abl['spin_length_one'] = {'mean': mag['mean']*SPIN, 'max': mag['max']*SPIN}
    abl['bond_counted_twice_convention'] = {'chain': chain/2}
    abl['unit_vector_bond_convention'] = {'chain': chain/4}
    abl['memorised_experimental_AF2_as_kmin'] = {'ligand_corrected': [0.506, 0.0, -0.483]}
    abl['paper_SI_kmin_verbatim'] = {'ligand_corrected': [0.491, 0.0, 0.474], 'uncorrected': [0.0, 0.5, 0.0]}
    # ---- robustness (defensible alternatives) ----
    for r in (8.0, 10.0, 14.0, 17.0):
        rob[f'kmin_r{r:g}'] = {'ligand_corrected': M.Folded(M.folded('FM/ligand', rmax=r, data=data['FM/ligand'])).lt()[0].tolist(),
                               'uncorrected': M.Folded(M.folded('FM/no_ligand', rmax=r, data=data['FM/no_ligand'])).lt()[0].tolist()}
    for r in (6.5, 10.0, 14.0):
        s = M.magnon_stats(M.Native('AFM/ligand', rmax=r, data=data['AFM/ligand']), spin=SPIN)
        rob[f'magnons_native_r{r:g}'] = {'mean': s['mean'], 'max': s['max']}
    coarse = M.magnon_stats(native, n=4, spin=SPIN)
    rob['magnons_grid4_unrefined_max'] = {'mean': coarse['mean'], 'max': coarse['max_grid_only']}
    rob['chain_single_bond_extremes'] = {'chain_low': -8*cj[3], 'chain_high': -8*cj[2]}
    s6 = M.magnon_stats(M.Native('AFM/ligand', rmax=6.0, data=data['AFM/ligand']), spin=SPIN)
    rob['magnons_native_r6_info'] = {'mean': s6['mean'], 'max': s6['max']}

    # ---- staged data ----
    def calc_meta(key):
        d = data[key]
        return {'calculation': key, 'description': LABEL[key], 'tb2j_header': d['header'], 'cell_A': d['cell_A'],
                'atom_columns': d['atom_columns'], 'atoms': d['atoms']}
    lab_sets = {}
    for key in M.SETS:
        wins = {wname(w): [] for w in WINDOWS}
        for p in data[key]['pairs']:
            for w in WINDOWS:
                if w[0] <= p[9] < w[1]:
                    wins[wname(w)].append(p)
                    break
            else:
                raise AssertionError(p)
        lab_sets[key] = {'meta': calc_meta(key), 'windows': wins}
    lab_disp = {}
    for key in M.SETS:
        head, arr = M.dispersion(key)
        lab_disp[key] = {'archived_header': head, 'values': arr.tolist()}
    fm = data['FM/ligand']
    catalog = {
        'material': 'CuO (tenorite), monoclinic C2/c, 4 Cu per conventional cell. Released calculations of Tellez-Mora, He, '
                    'Bousquet, Wirtz and Romero (2024): SIESTA DFT+U (U = 7.14 eV on Cu 3d), exchange from the magnetic '
                    'force theorem (Green-function method) with TB2J 0.7.4.5.',
        'exchange_convention': 'Archived TB2J convention: E = -sum_{i != j} J_ij e_i.e_j over ordered pairs (each bond '
                               'appears twice, as (i,j,R) and (j,i,-R)), e_i unit vectors along the local moments, '
                               'J > 0 ferromagnetic, meV. R is in lattice vectors of the calculation\'s own cell; '
                               'vector = r_j + R.cell - r_i in Angstrom.',
        'isotropic_only': 'Only isotropic exchange is served: the archive README states that only the isotropic exchange '
                          'constants are valid.',
        'conventional_cell': {'cell_A': fm['cell_A'], 'atoms': [{'name': a['name'], 'position_A': a['position_A']}
                                                                for a in fm['atoms']],
                              'note': 'Cell and Cartesian frame of the FM-reference calculations. Wave vectors in this task '
                                      'are in reciprocal-lattice units of this cell.'},
        'calculations': {key: LABEL[key] for key in M.SETS},
        'calculation_notes': 'Each calculation was done in its own cell and magnetic configuration; the cell, atoms, '
                             'charges and moments of a calculation are returned with every exchange purchase from it. '
                             '"Spin-ligand correction" is TB2J\'s correction for the ligand (O) spin polarisation.',
        'requests': {
            'exchange': {'example': {'kind': 'exchange', 'calculation': 'FM/ligand', 'window': '0-4'},
                         'windows_A': {wname(w): {'distance_range_A': list(w), 'cost': c} for w, c in zip(WINDOWS, COSTS)},
                         'returns': 'all isotropic pairs of that calculation with lo <= distance < hi, as rows '
                                    '[i, j, Rx, Ry, Rz, J_iso_meV, dx_A, dy_A, dz_A, distance_A], plus the calculation '
                                    'metadata (TB2J header lines, cell, atoms with charges and moments).'},
            'author_dispersion': {'example': {'kind': 'author_dispersion', 'calculation': 'AFM/ligand'}, 'cost': 1,
                                  'returns': 'the archived file of eigenvalues of the paper\'s bosonic matrix h(k) '
                                             '(Eq. 9) for that calculation\'s reference configuration, 300 k-points along '
                                             'the path named in its header line. The k-point coordinates and units are '
                                             'not archived.'}},
        'budget': BUDGET,
        'catalog_value': 4*sum(COSTS) + 4,
    }
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps({'windows': [wname(w) for w in WINDOWS],
                                                              'costs': dict(zip([wname(w) for w in WINDOWS], COSTS)),
                                                              'sets': lab_sets, 'dispersions': lab_disp}) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    t = hidden['targets']
    cl, cn = M.canonical(kmin_lig), M.canonical(kmin_nol)
    values = [float(cl[0]), float(cl[2]), float(1 + cl[2]), float(1 - cl[0]), float(cn[1]), float(1 + cn[1]),
              chain, mag['mean'], mag['max']]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1) + '\n')
    for f in ('environment/lab/public.json', 'tests/public.json'):
        p = json.loads((TASK/f).read_text())
        p['solver_budget'] = BUDGET
        (TASK/f).write_text(json.dumps(p, indent=1) + '\n')

    # ---- calibration ----
    def qd(field, q):
        return round(M.q_distance(q, t['kmin_FM_reference_rlu'][field]), 4)
    cal = {
        'kmin_ligand_corrected': {'reference': [round(x, 4) for x in kmin_lig], 'tolerance_rlu': THRESHOLDS['kmin_ligand_rlu'],
                                  'robustness_distance': {k: qd('ligand_corrected', v['ligand_corrected']) for k, v in rob.items() if k.startswith('kmin')},
                                  'ablation_distance': {'truncated_6A': qd('ligand_corrected', k6l), 'commensurate': qd('ligand_corrected', [0.5, 0, 0.5]),
                                                        'AF2_memorised': qd('ligand_corrected', [0.506, 0, -0.483]),
                                                        'paper_verbatim': qd('ligand_corrected', [0.491, 0, 0.474]),
                                                        'uncorrected_set': qd('ligand_corrected', kmin_nol)}},
        'kmin_uncorrected': {'reference': [round(x, 4) for x in kmin_nol], 'tolerance_rlu': THRESHOLDS['kmin_uncorrected_rlu'],
                             'robustness_distance': {k: qd('uncorrected', v['uncorrected']) for k, v in rob.items() if k.startswith('kmin')},
                             'ablation_distance': {'truncated_6A': qd('uncorrected', k6n), 'commensurate': qd('uncorrected', [0, 0.5, 0])}},
        'converged_q': {'reference': [0.5, 0, 0.5], 'tolerance_rlu': THRESHOLDS['converged_q_rlu'],
                        'ablation_distance': {'symmetrized_afm': round(M.q_distance(q_sym, [0.5, 0, 0.5]), 4),
                                              'fm_reference_final': round(M.q_distance(kmin_lig, [0.5, 0, 0.5]), 4)}},
        'chain_exchange_meV': {'reference': round(chain, 3), 'tolerance': '+-4%',
                               'robustness': {'single_bond_low': round(-8*cj[3], 3), 'single_bond_high': round(-8*cj[2], 3)},
                               'ablations': {'fm_reference': round(abl['fm_reference_as_final']['chain'], 2),
                                             'uncorrected': round(abl['uncorrected_exchange_as_final']['chain'], 2),
                                             'bond_twice': round(chain/2, 2), 'unit_vector': round(chain/4, 2)}},
        'mean_magnon_energy_meV': {'reference': round(mag['mean'], 3), 'tolerance': '+-5%',
                                   'robustness': {k: round(v['mean'], 3) for k, v in rob.items() if k.startswith('magnons')},
                                   'ablations': {k: round(v['mean'], 2) for k, v in abl.items() if 'mean' in v}},
        'magnon_band_maximum_meV': {'reference': round(mag['max'], 3), 'tolerance': '+-4%',
                                    'robustness': {k: round(v['max'], 3) for k, v in rob.items() if k.startswith('magnons')},
                                    'ablations': {k: round(v['max'], 2) for k, v in abl.items() if 'max' in v}},
        'stability': {'reference': stability, 'min_h_eigenvalue_meV': {k: round(v, 3) for k, v in min_eigs.items()}},
        'author_path_check': path_check,
        'symmetrized_afm_af1_minus_spiral_meV_per_Cu': abl['symmetrized_afm_exchange']['af1_minus_spiral_meV_per_Cu'],
        'AF1_energy_native_vs_symmetrized_meV_per_Cu': [native.energy(), -M.Folded(sym).top([0.5, 0, 0.5])],
    }
    (HERE/'calibration.json').write_text(json.dumps(cal, indent=1) + '\n')
    print(json.dumps(cal, indent=1))


if __name__ == '__main__':
    main()
