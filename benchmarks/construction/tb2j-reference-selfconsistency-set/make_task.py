"""Generate the staged, lab, private and ablation data for tb2j-reference-selfconsistency-set.

    python3 -B benchmarks/construction/tb2j-reference-selfconsistency-set/make_task.py     (inside the runtime image)

Reads the extracted archive members (fetch.py) and writes:
  environment/agent/catalog.json, environment/lab/data.json, environment/lab/private/hidden.json,
  environment/lab/private/hidden_values.json, solution/ablation_values.json, both public.json budgets,
  and calibration.json (this folder).
It asserts that the paper's FM-reference results and the authors' archived h(k) files are reproduced from the archived
pair lists, and that the FePS3 DFT+U FM pair list duplicates the U = 0 one.
"""
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import tb2jmodel as M  # noqa: E402

TASK = HERE.parents[1]/'tasks/tb2j-reference-selfconsistency-set'
WINDOWS = [(0, 5), (5, 8), (8, 11), (11, 14), (14, 18), (18, 23), (23, 30), (30, 40), (40, 50)]
COSTS = [1, 1, 1, 1, 2, 2, 2, 3, 3]
BUDGET = 68
THRESHOLDS = {'kmin_fm_FeP_rlu': 0.005, 'kmin_converged_FeP_rlu': 0.012, 'energy_rel': 0.04}
MATERIAL_TEXT = {
    'NiO': 'rock-salt NiO, one Ni per primitive fcc cell in the FM reference',
    'FeCl2': '1T-FeCl2 (bulk, P-3m1 layers in rhombohedral stacking), three Fe layers per hexagonal cell',
    'FePS3': 'monoclinic FePS3 (C2/m van der Waals layers, Fe honeycomb), four Fe per cell; DFT (U = 0) and DFT+U (U = 3 eV on Fe 3d per the paper) workflows',
    'FeP': 'orthorhombic FeP (MnP type, Pnma; metallic), four Fe per cell'}


def wname(w):
    return f'{w[0]}-{w[1]}'


def window_of(dist):
    for w in WINDOWS:
        if w[0] <= dist < w[1]:
            return wname(w)
    raise ValueError(dist)


def upto(d, hi):
    return M.subset(d, lambda p: p[9] < hi)


def meta(d, label):
    return {'calculation': label, 'archive_member': M.CALCS[label], 'tb2j_header': d['header'], 'cell_A': d['cell_A'],
            'atom_columns': d['atom_columns'], 'atoms': d['atoms'],
            'n_pairs_total': len(d['pairs'])}


def path(nodes, idx, n):
    ks = []
    for r in range(n):
        s = max(i for i in range(len(idx) - 1) if idx[i] <= r) if r < idx[-1] else len(idx) - 2
        t = (r - idx[s])/(idx[s + 1] - idx[s])
        ks.append(np.array(nodes[s]) + t*(np.array(nodes[s + 1]) - np.array(nodes[s])))
    return np.array(ks)


def validate(data):
    """Reproduce the paper's FM-reference results and the authors' archived h(k) files from the archived pair lists."""
    out = {}
    # FePS3 (U = 0) FM reference: Z = (0, 0, 1/2) below Gamma by 2.0 meV (paper text); archived file at Gamma and Z
    ps = M.Model(data['FePS3/DFT/FM'])
    hz, hg = (np.sort(np.linalg.eigvalsh(ps.h([k]))[0]) for k in ([0, 0, 0.5], [0, 0, 0]))
    f = np.sort(np.loadtxt(M.RAW/'FePS3/DFT/FM/dispersion_FM'), axis=1)*1000
    out['FePS3_FM_Z_minus_Gamma_meV'] = float(hz[0] - hg[0])
    out['FePS3_FM_archived_dev_at_Gamma_Z_meV'] = float(max(np.abs(f[0] - hg).max(), np.abs(f[206] - hz).max()))
    assert abs(out['FePS3_FM_Z_minus_Gamma_meV'] + 2.0) < 0.05 and out['FePS3_FM_archived_dev_at_Gamma_Z_meV'] < 0.01
    # FeCl2 FM/AFM archived files at the path nodes M, K, Gamma, A, L (hexagonal r.l.u. of each calculation's cell)
    nodes = [[.5, 0, 0], [1/3, 1/3, 0], [0, 0, 0], [0, 0, .5], [.5, 0, .5]]
    for key, rel in (('FeCl2/FM', 'FeCl2/FM/dispersion_FM'), ('FeCl2/AFM', 'FeCl2/AFM/dispersion_AFM')):
        a = np.sort(np.loadtxt(M.RAW/rel), axis=1)*1000
        mine = np.sort(np.linalg.eigvalsh(M.Model(data[key]).h(nodes)), axis=1)
        out[key + '_archived_dev_at_nodes_meV'] = float(np.abs(mine - a[[0, 138, 416, 456, 699]]).max())
        assert out[key + '_archived_dev_at_nodes_meV'] < 0.02
    # FeP AFM (spiral) archived file: reproduced by the pair list as listed and the archived moments; flat along Gamma-Z
    fa = M.Model(data['FeP/AFM'])
    a = np.sort(np.loadtxt(M.RAW/'FeP/AFM/dispersion_AFM'), axis=1)*1000
    mine = np.sort(np.linalg.eigvalsh(fa.h(path([[0, .5, 0], [0, 0, 0], [0, 0, .5]], [0, 157, 300], 300))), axis=1)
    out['FeP_AFM_archived_dev_meV'] = float(np.abs(mine - a).max())
    out['FeP_AFM_archived_GammaZ_flatness_meV'] = float(np.abs(a[157:] - a[157]).max())
    out['FeP_AFM_archived_min_meV'] = float(a.min())
    assert out['FeP_AFM_archived_dev_meV'] < 0.15 and out['FeP_AFM_archived_GammaZ_flatness_meV'] < 1e-6
    # paper Tables: FeP NN/NNN (AFM reference), NiO J1/J2 (PBE+U rows), FeCl2 J1/J2 (SI Table 1), FePS3 Table 2
    def shell(key, lo, hi):
        return float(np.mean([p[5] for p in data[key]['pairs'] if lo < p[9] < hi]))
    out['paper_checks'] = {'FeP_AFM_NN': shell('FeP/AFM', 2.6, 2.7), 'NiO_FM_J2': shell('NiO/FM', 4.2, 4.3),
                           'NiO_AFM_J2': shell('NiO/AFM', 4.1, 4.15), 'FeCl2_AFM_J1': shell('FeCl2/AFM', 3.5, 3.7),
                           'FePS3_U_AFM_J3': shell('FePS3/DFT+U/AFM', 6.95, 7.0)}
    assert abs(out['paper_checks']['FeP_AFM_NN'] + 4.642) < 0.002 and abs(out['paper_checks']['NiO_FM_J2'] + 11.58) < 0.01
    assert abs(out['paper_checks']['FeCl2_AFM_J1'] - 7.330) < 0.002 and abs(out['paper_checks']['FePS3_U_AFM_J3'] + 3.35) < 0.01
    # the FePS3 DFT+U FM list duplicates the U = 0 list (same multiset of pairs, cell and moments)
    key = lambda p: (p[0], p[1], p[2], p[3], p[4], p[5])
    out['FePS3_DFT+U_FM_duplicates_DFT_FM'] = (sorted(map(key, data['FePS3/DFT/FM']['pairs'])) ==
                                               sorted(map(key, data['FePS3/DFT+U/FM']['pairs'])) and
                                               data['FePS3/DFT/FM']['atoms'] == data['FePS3/DFT+U/FM']['atoms'] and
                                               data['FePS3/DFT/FM']['cell_A'] == data['FePS3/DFT+U/FM']['cell_A'])
    assert out['FePS3_DFT+U_FM_duplicates_DFT_FM']
    # supercell lists computed with R = 0 along the long axis only (FeCl2 and FeP AFM references)
    out['R_long_axis_values'] = {k: sorted({p[4] for p in data[k]['pairs']}) for k in ('FeCl2/AFM', 'FeP/AFM', 'FePS3/DFT/AFM', 'NiO/AFM')}
    assert out['R_long_axis_values']['FeCl2/AFM'] == [0] and out['R_long_axis_values']['FeP/AFM'] == [0]
    return out


def main():
    data = {k: M.load(k) for k in M.CALCS}
    checks = validate(data)
    fepF, fepA = data['FeP/FM'], data['FeP/AFM']
    fm = M.Model(fepF)
    starts = ([0, 0.03, 0.22], [0, 0, 0.2], [0, 0.5, 0])

    # ---- graded targets ----
    k_fm, low_fm = M.kmin(fm, starts=starts)
    assert M.qdist(k_fm, [0, 0.03, 0.21]) < 0.01                      # paper: k_min = (0.00, 0.03, 0.21)
    acc_mi = M.fold(fm, M.minimage(fepA))
    k_conv, low_conv = M.kmin(M.FoldedModel(acc_mi, 4), starts=starts)
    energies, conf = {}, {}
    for mat, (f, a) in M.MATERIALS.items():
        ef, ea = M.stabilization(data[f], data[a])
        energies[mat] = {'fm_reference': ef, 'converged_reference': ea}
    energies['FePS3_DFT+U']['fm_reference'] = None                    # the archived FM file is the U = 0 calculation
    for mat, (f, a) in M.MATERIALS.items():
        assert M.Model(data[a]).n > M.Model(data[f]).n or mat == 'NiO'
    # stability of the converged references at their own (archived) moments, and of FeP after relaxation
    stab = {}
    for mat, (f, a) in M.MATERIALS.items():
        am = M.Model(data[a])
        stab[a] = float(np.linalg.eigvalsh(am.h(M.grid(6))).min())
        stab[f] = float(np.linalg.eigvalsh(M.Model(data[f]).h(M.grid(6))).min())
    e_rel, E_rel = M.relax(M.Model(fepA).K0(), M.Model(fepA).e)
    fepA_mi = M.Model(M.minimage(fepA))
    stab['FeP/AFM relaxed, as listed'] = float(np.linalg.eigvalsh(M.Model(fepA).h(M.grid(8), e=e_rel)).min())
    stab['FeP/AFM relaxed, minimum image'] = float(np.linalg.eigvalsh(fepA_mi.h(M.grid(8), e=e_rel)).min())
    torque = np.linalg.norm(np.cross(M.Model(fepA).e, M.Model(fepA).K0() @ M.Model(fepA).e), axis=1)

    # ---- ablations: realistic wrong paths computed from the same data ----
    abl, rob = {}, {}
    abl['kmin_FeP_from_pairs_below_14A'] = M.kmin(M.Model(upto(fepF, 14)), starts=starts)[0].tolist()
    abl['kmin_FeP_from_pairs_below_11A'] = M.kmin(M.Model(upto(fepF, 11)), starts=starts)[0].tolist()
    abl['kmin_FeP_paper_verbatim'] = [0.0, 0.03, 0.21]
    abl['kmin_FeP_commensurate_approximant'] = [0.0, 0.0, 0.2]
    abl['FeP_fm_energy_pairs_below_18A'] = M.stabilization(upto(fepF, 18), fepA)[0]
    abl['kmin_FeP_pairs_below_18A'] = M.kmin(M.Model(upto(fepF, 18)), starts=starts)[0].tolist()
    abl['insulators_truncated'] = {
        'FeCl2_fm_below_8A': M.stabilization(upto(data['FeCl2/FM'], 8), data['FeCl2/AFM'])[0],
        'FePS3_DFT_fm_below_11A': M.stabilization(upto(data['FePS3/DFT/FM'], 11), data['FePS3/DFT/AFM'])[0],
        'FePS3_DFT_conv_below_11A': M.stabilization(data['FePS3/DFT/FM'], upto(data['FePS3/DFT/AFM'], 11))[1]}
    abl['FeCl2_conv_long_windows_skipped'] = M.stabilization(data['FeCl2/FM'], upto(data['FeCl2/AFM'], 30))[1]
    abl['FeP_conv_long_windows_skipped'] = M.stabilization(fepF, upto(fepA, 23))[1]
    abl['kmin_conv_FeP_folded_as_listed'] = M.kmin(M.FoldedModel(M.fold(fm, fepA), 4), starts=starts)[0].tolist()
    abl['kmin_conv_FeP_folded_as_listed_sum5'] = M.kmin(M.FoldedModel(M.fold(fm, fepA), 4, average=5), starts=starts)[0].tolist()
    abl['FePS3_DFT+U_fm_from_archived_FM_file'] = M.stabilization(data['FePS3/DFT+U/FM'], data['FePS3/DFT+U/AFM'])[0]
    w = np.linalg.eigvalsh(fm.K([k_fm])[0])
    abl['FeP_fm_energy_of_its_own_LT_helix'] = float(w[-1] - fm.K0().sum(1).mean())
    abl['kmin_conv_FeP_assumed_reference_independent'] = k_fm.tolist()
    # robustness: defensible alternatives
    rob['kmin_FeP_pairs_below_18A'] = abl['kmin_FeP_pairs_below_18A']
    rob['kmin_FeP_pairs_below_23A'] = M.kmin(M.Model(upto(fepF, 23)), starts=starts)[0].tolist()
    ef_r, ea_r = M.stabilization(fepF, fepA, e=e_rel)
    rob['FeP_relaxed_spiral'] = {'fm_reference': ef_r, 'converged_reference': ea_r}
    rob['insulators_pairs_below_14A'] = {
        'FePS3_DFT': dict(zip(('fm_reference', 'converged_reference'), M.stabilization(upto(data['FePS3/DFT/FM'], 14), upto(data['FePS3/DFT/AFM'], 14)))),
        'FeCl2_fm_below_11A': M.stabilization(upto(data['FeCl2/FM'], 11), data['FeCl2/AFM'])[0],
        'NiO': dict(zip(('fm_reference', 'converged_reference'), M.stabilization(upto(data['NiO/FM'], 5), upto(data['NiO/AFM'], 5))))}
    rob['kmin_conv_FeP_sum5'] = M.kmin(M.FoldedModel(acc_mi, 4, average=5), starts=starts)[0].tolist()
    rob['kmin_conv_FeP_drop_half_box'] = M.kmin(M.FoldedModel(M.fold(fm, M.subset(M.minimage(fepA), lambda p: abs(p[8]) < 13.0)), 4), starts=starts)[0].tolist()
    rob['kmin_conv_FeP_bought_below_30A'] = M.kmin(M.FoldedModel(M.fold(fm, M.minimage(upto(fepA, 30))), 4), starts=starts)[0].tolist()

    # ---- staged and lab files ----
    sets = {}
    for label, d in data.items():
        win = {wname(w): [] for w in WINDOWS}
        for p in d['pairs']:
            win[window_of(p[9])].append(p)
        sets[label] = {'meta': meta(d, label), 'windows': win}
    (TASK/'environment/lab/data.json').write_text(json.dumps({'sets': sets, 'costs': {wname(w): c for w, c in zip(WINDOWS, COSTS)}}) + '\n')
    catalog = {
        'source': 'Released TB2J output of Tellez-Mora, He, Bousquet, Wirtz and Romero (2024): SIESTA DFT (PBE; DFT+U where stated), '
                  'isotropic exchange from the magnetic force theorem (Green-function method) with TB2J 0.7.4.x.',
        'materials': MATERIAL_TEXT,
        'calculations': {label: ('FM reference' if label.endswith('/FM') else 'AFM reference (the authors\' updated reference)')
                         for label in M.CALCS},
        'calculation_notes': 'Calculation labels are the archive folder names. Each calculation was done in its own cell and '
                             'magnetic configuration; its TB2J header lines, cell, atoms (positions, Wannier charges and moments) '
                             'are returned with every exchange purchase from it.',
        'exchange_convention': 'Archived TB2J convention: E = -sum_{i != j} J_ij e_i.e_j over ordered pairs (each bond appears '
                               'twice, as (i,j,R) and (j,i,-R)), e_i unit vectors along the local moments, J > 0 ferromagnetic, meV. '
                               'R is in lattice vectors of the calculation\'s own cell; vector = r_j + R.cell - r_i in Angstrom; '
                               'distance is |vector|, as archived.',
        'isotropic_only': 'Only isotropic exchange is served (the archive README states that only the isotropic constants are valid '
                          'for NiO, FeCl2 and FePS3).',
        'requests': {'exchange': {
            'example': {'kind': 'exchange', 'calculation': 'FeP/FM', 'window': '0-5'},
            'windows_A': {wname(w): {'distance_range_A': list(w), 'cost': c} for w, c in zip(WINDOWS, COSTS)},
            'returns': 'all isotropic pairs of that calculation with lo <= archived distance < hi, as rows '
                       '[i, j, Rx, Ry, Rz, J_iso_meV, dx_A, dy_A, dz_A, distance_A], plus the calculation metadata.'}},
        'budget': BUDGET, 'catalog_value': len(M.CALCS)*sum(COSTS)}
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    for rel in ('environment/lab/public.json', 'tests/public.json'):
        pub = json.loads((TASK/rel).read_text())
        pub['solver_budget'] = BUDGET
        (TASK/rel).write_text(json.dumps(pub, indent=1) + '\n')

    hidden = {'targets': {'fm_reference_kmin_FeP_rlu': k_fm.tolist(), 'converged_exchange_kmin_FeP_rlu': k_conv.tolist(),
                          'stabilization_energy_meV_per_ion': energies},
              'thresholds': THRESHOLDS,
              'reference_details': {'validation': checks, 'lowest_h_FM_kmin_meV': {'FeP_FM': low_fm, 'FeP_converged_folded': low_conv},
                                    'min_h_eigenvalue_6grid_meV': stab, 'FeP_AFM_torque_meV': [float(torque.min()), float(torque.max())],
                                    'FeP_AFM_relaxed_energy_meV_per_ion': E_rel, 'FeP_AFM_archived_energy_meV_per_ion': M.Model(fepA).energy()}}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = [abs(k_fm[1]), abs(k_fm[2]), abs(k_conv[1]), abs(k_conv[2])]
    values += [v for e in energies.values() for v in e.values() if v is not None]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps([float(v) for v in values]) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1) + '\n')

    # ---- calibration ----
    def dist(q, ref):
        return round(M.qdist(q, ref), 4)
    def rel(v, ref):
        return round(v/ref - 1, 4)
    cal = {
        'fm_reference_kmin_FeP': {'reference': np.round(M.qabs(k_fm), 4).tolist(), 'tolerance_rlu': THRESHOLDS['kmin_fm_FeP_rlu'],
                                  'robustness_distance': {k: dist(rob[k], k_fm) for k in ('kmin_FeP_pairs_below_18A', 'kmin_FeP_pairs_below_23A')},
                                  'ablation_distance': {k: dist(abl[k], k_fm) for k in ('kmin_FeP_from_pairs_below_14A', 'kmin_FeP_from_pairs_below_11A',
                                                                                         'kmin_FeP_paper_verbatim', 'kmin_FeP_commensurate_approximant')}},
        'converged_exchange_kmin_FeP': {'reference': np.round(M.qabs(k_conv), 4).tolist(), 'tolerance_rlu': THRESHOLDS['kmin_converged_FeP_rlu'],
                                        'robustness_distance': {k: dist(rob[k], k_conv) for k in ('kmin_conv_FeP_sum5', 'kmin_conv_FeP_drop_half_box', 'kmin_conv_FeP_bought_below_30A')},
                                        'ablation_distance': {k: dist(abl[k], k_conv) for k in ('kmin_conv_FeP_folded_as_listed', 'kmin_conv_FeP_folded_as_listed_sum5',
                                                                                                'kmin_conv_FeP_assumed_reference_independent', 'kmin_FeP_commensurate_approximant')}},
        'stabilization_energy_meV_per_ion': {
            'reference': {m: {k: (None if v is None else round(v, 4)) for k, v in e.items()} for m, e in energies.items()},
            'tolerance': 'relative 4% per entry; FePS3_DFT+U fm_reference must be null',
            'robustness_relative': {
                'FeP_relaxed_spiral': {k: rel(rob['FeP_relaxed_spiral'][k], energies['FeP'][k]) for k in ('fm_reference', 'converged_reference')},
                'FePS3_DFT_below_14A': {k: rel(rob['insulators_pairs_below_14A']['FePS3_DFT'][k], energies['FePS3_DFT'][k]) for k in ('fm_reference', 'converged_reference')},
                'FeCl2_fm_below_11A': rel(rob['insulators_pairs_below_14A']['FeCl2_fm_below_11A'], energies['FeCl2']['fm_reference']),
                'NiO_below_5A': {k: rel(rob['insulators_pairs_below_14A']['NiO'][k], energies['NiO'][k]) for k in ('fm_reference', 'converged_reference')}},
            'ablation_relative': {
                'FeP_fm_pairs_below_18A': rel(abl['FeP_fm_energy_pairs_below_18A'], energies['FeP']['fm_reference']),
                'FeP_fm_own_LT_helix': rel(abl['FeP_fm_energy_of_its_own_LT_helix'], energies['FeP']['fm_reference']),
                'FeCl2_fm_below_8A': rel(abl['insulators_truncated']['FeCl2_fm_below_8A'], energies['FeCl2']['fm_reference']),
                'FePS3_DFT_fm_below_11A': rel(abl['insulators_truncated']['FePS3_DFT_fm_below_11A'], energies['FePS3_DFT']['fm_reference']),
                'FePS3_DFT_conv_below_11A': rel(abl['insulators_truncated']['FePS3_DFT_conv_below_11A'], energies['FePS3_DFT']['converged_reference']),
                'FeCl2_conv_long_windows_skipped': rel(abl['FeCl2_conv_long_windows_skipped'], energies['FeCl2']['converged_reference']),
                'FeP_conv_long_windows_skipped': rel(abl['FeP_conv_long_windows_skipped'], energies['FeP']['converged_reference']),
                'FePS3_DFT+U_fm_from_duplicate': round(abl['FePS3_DFT+U_fm_from_archived_FM_file'], 4),
                'bond_counted_once': -0.5}},
        'reference_ratio_converged_over_fm': {m: (None if e['fm_reference'] is None else round(e['converged_reference']/e['fm_reference'], 4))
                                              for m, e in energies.items()},
        'stability': stab, 'validation': checks}
    (HERE/'calibration.json').write_text(json.dumps(cal, indent=1, default=float) + '\n')
    print(json.dumps(cal, indent=1, default=float))


if __name__ == '__main__':
    main()
