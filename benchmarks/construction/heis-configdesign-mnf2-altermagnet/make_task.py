"""Generate the staged, lab, private and ablation data for heis-configdesign-mnf2-altermagnet.

    cd benchmarks/construction/heis-configdesign-mnf2-altermagnet && python3 -B make_task.py

Inputs: pool.json (extract.py), j_reference.json (reference.py), tn_reference.json and tn_coarse.json (tn_mc.py).
Writes the catalog, lab data, private targets and ablation values, and prints the calibration numbers.
"""
import json
import numpy as np
from pathlib import Path
import model as M

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/heis-configdesign-mnf2-altermagnet'
P = M.POOL['spirals']
Q = np.array([s['q_reciprocal_lattice'] for s in P])
Y = np.array([s['energy_meV_per_cell'] for s in P])            # = archive CSV energies, meV per 6-atom cell
YCSV = np.array([s['energy_csv_meV'] for s in P])
D = np.array([s['direction'] for s in P])
T = np.round(np.array([s['t'] for s in P]), 2)
CLEAN = np.array([s['enmax_needed_eV'] <= s['encut_eV'] for s in P])
AVAILABLE_DIRS = ['[001]', '[121]', '[214]']
HIDDEN_DIRS = ['[111]', '[113]']
BUDGET = 24
REF = json.loads((HERE/'j_reference.json').read_text())
J_REF = np.array(REF['J_meV'])
THRESHOLDS = {'spiral_rmse_meV': 0.40, 'ground_state_tol_meV': 0.35, 'J1_window_meV': [-1.25, -1.00], 'J2_window_meV': [-3.10, -2.86],
              'J3_window_meV': [-0.09, 0.03], 'TN_rel_tol': 0.12, 'split_shell_distance_tol_A': 0.05,
              'split_shell_mean_window_meV': [-0.11, -0.025]}
SPLIT_SHELL = next(i for i, s in enumerate(M.SHELLS) if M.altermagnetic(s))     # 0-based index of the 7.007 A shell
ORACLE_PLAN = {'G': [0.0], '[001]': [0.15, 0.35, 0.55, 0.75, 0.9, 1.0], '[121]': [0.1, 0.2, 0.3, 0.4, 0.5],
               '[214]': [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85]}


def key_of(i):
    return 'G-00' if D[i] == 'gamma' else f'{D[i]}-{int(round(T[i] / 0.05)):02d}'


def rows(plan):
    out = []
    for d, ts in plan.items():
        for t in ts:
            dd = 'gamma' if d == 'G' else d
            out.append(int(np.where((D == dd) & (np.abs(T - t) < 1e-6))[0][0]))
    return np.array(out)


HID = np.where(np.isin(D, HIDDEN_DIRS))[0]
HID_IDS = [f'P{D[i][1:4]}-{int(round(T[i] / 0.05)):02d}' for i in HID]
# hidden spirals that realise the Neel state: integer q with h + k + l odd (uniform spiral, tau = (1/2,1/2,1/2))
GROUND_IDS = [pid for pid, i in zip(HID_IDS, HID) if np.allclose(Q[i], np.round(Q[i])) and int(round(Q[i].sum())) % 2 == 1]


def fit_rows(r, n, convention='uniform'):
    return M.fit(2 * M.X(Q[r], n, convention), Y[r], intercept=True)


def predict(a, J, n, convention='uniform'):
    return {pid: float(v) for pid, v in zip(HID_IDS, a + 2 * M.X(Q[HID], n, convention) @ J)}


def rmse(pred):
    return float(np.sqrt(np.mean([(pred[pid] - Y[i]) ** 2 for pid, i in zip(HID_IDS, HID)])))


def ground_error(pred):
    return float(np.mean([pred[pid] - YCSV[HID[HID_IDS.index(pid)]] for pid in GROUND_IDS]))


def loo_choice(r, nmax=12):
    """Shell count by leave-one-out cross-validation (closed form for linear least squares)."""
    cv = {}
    for n in range(3, nmax + 1):
        A = np.c_[np.ones(len(r)), 2 * M.X(Q[r], n)]
        H = A @ np.linalg.pinv(A)
        res = Y[r] - H @ Y[r]
        lev = np.clip(1 - np.diag(H), 1e-9, None)
        cv[n] = float(np.sqrt(np.mean((res / lev) ** 2)))
    return min(cv, key=cv.get), cv


def orthorhombic_difference(r, n):
    """Wrong-symmetry fit: the split shell gets separate [110] and [1-10] constants, equal on both sublattices."""
    X = 2 * M.X(Q[r], n)
    v = M.SHELLS[SPLIT_SHELL]['vectors']
    pos = v[(v[:, 0] * v[:, 1]) > 0]
    neg = v[(v[:, 0] * v[:, 1]) < 0]
    cpos = 2 * (1 - np.cos(2 * np.pi * Q[r] @ pos.T)).sum(1)
    cneg = 2 * (1 - np.cos(2 * np.pi * Q[r] @ neg.T)).sum(1)
    Xo = np.c_[np.delete(X, SPLIT_SHELL, axis=1), cpos, cneg]
    a, c = M.fit(Xo, Y[r], intercept=True)
    return float(c[-2] - c[-1])


def tn_ref():
    tn = json.loads((HERE/'tn_reference.json').read_text())
    cross = [c['T_cross_K'] for c in tn['crossings'] if c['T_cross_K']]
    coarse = json.loads((HERE/'tn_coarse.json').read_text())
    return float(np.mean(cross)), tn, coarse


def unit(m):
    m = np.array(m, float)
    return (m / np.linalg.norm(m)).round(6).tolist()


def lab_record(i):
    s = P[i]
    return {'q_reciprocal_lattice': s['q_reciprocal_lattice'],
            'vasp_total_energy_eV': s['total_energy_eV'],
            'cell_content': '6-atom primitive cell: Mn1 at (1/2,1/2,1/2), Mn2 at (0,0,0), 4 F',
            'archive_csv_energy': s['energy_csv_meV'],
            'archive_csv_energy_units_as_documented': 'meV / atom (archive README: "DFT energies: meV / atom")',
            'constrained_moments_muB': {'Mn1': s['mn_moments_muB'][0], 'Mn2': s['mn_moments_muB'][1],
                                        'F': s['moments_all_atoms_muB'][2:]},
            'incar': s['incar'],
            'kpoints': 'Monkhorst-Pack ' + ' x '.join(s['kmesh']),
            'outcar_warnings': s['outcar_warnings'],
            'electronic_steps': s['scf_steps']}


def main():
    t_ref, tn_data, tn_coarse = tn_ref()
    # ---------------- catalog, lab data, hidden data
    catalog_spirals, lab = [], {}
    for d in ['gamma'] + AVAILABLE_DIRS:
        for i in np.where(D == d)[0][np.argsort(T[D == d])]:
            k = key_of(i)
            catalog_spirals.append({'key': k, 'direction': 'Gamma' if d == 'gamma' else d, 't': float(T[i]),
                                    'q_reciprocal_lattice': Q[i].tolist()})
            lab[k] = lab_record(i)
    predict_list = [{'id': pid, 'direction': D[i], 't': float(T[i]), 'q_reciprocal_lattice': Q[i].tolist()}
                    for pid, i in zip(HID_IDS, HID)]
    s0 = P[0]
    catalog = {
        'material': 'rutile MnF2 (P4_2/mnm), noncollinear PBE spin spirals (VASP 6.4.2, constrained moments, '
                    'generalized Bloch theorem), released by Li, von Malottki and Rignanese (2026)',
        'structure': {'cell_A': M.POOL['cell_A'], 'species': M.POOL['species'],
                      'fractional_positions': M.POOL['fractional_positions'],
                      'labels': ['Mn1', 'Mn2', 'F', 'F', 'F', 'F']},
        'spiral_definition': ('LSPIRAL (generalized Bloch theorem): the magnetization in the cell displaced by a '
                              'lattice vector R equals that of the reference cell rotated about the spin z axis by '
                              '+2*pi*q.R; the in-cell moments are those of MAGMOM / M_CONSTR in each calculation\'s INCAR. '
                              'q (QSPIRAL) is in reciprocal-lattice coordinates of the tetragonal cell. All released '
                              'spirals are flat (moments in the spin xy plane).'),
        'archive_documentation': {'energies': 'meV / atom', 'exchange parameters': 'meV',
                                  'q': 'reciprocal lattice coordinates', 'source': 'archive README.md'},
        'common_settings': {k: s0['incar'][k] for k in ('ENCUT', 'ENINI', 'PREC', 'LREAL', 'I_CONSTRAINED_M', 'LAMBDA', 'RWIGS', 'ISMEAR', 'SIGMA')},
        'purchase_costs': {'spiral': 1}, 'budget': BUDGET,
        'spirals_available': catalog_spirals, 'spirals_to_predict': predict_list}
    hidden_targets = {pid: float(YCSV[i]) for pid, i in zip(HID_IDS, HID)}
    split = M.SHELLS[SPLIT_SHELL]
    hidden = {'spiral_targets': hidden_targets, 'ground_state_ids': GROUND_IDS, 'thresholds': THRESHOLDS,
              'reference': {'J_meV': J_REF.tolist(), 'J1_meV': float(J_REF[0]), 'J2_meV': float(J_REF[1]),
                            'J3_meV': float(J_REF[2]), 'TN_K': t_ref, 'split_shell_distance_A': split['distance'],
                            'split_shell_mean_meV': float(J_REF[SPLIT_SHELL])}}
    # ---------------- oracle replica, ablations, robustness
    r_or = rows(ORACLE_PLAN)
    n_or, cv = loo_choice(r_or)
    a_or, J_or = fit_rows(r_or, n_or)
    pred_or = predict(a_or, J_or, n_or)
    base_alt = {'distance_A': split['distance'], 'mean_meV': float(J_or[SPLIT_SHELL]), 'difference_meV': None}
    abl, rob, cal = {}, {}, {}
    # dp1 energy normalisation: trust the documented 'meV / atom' label (6 atoms per cell) or read it as per Mn
    for name, f in (('per_atom_label', 6.0), ('per_mn_label', 2.0)):
        abl[name] = {'predictions': {k: f * v for k, v in pred_or.items()}, 'J': (f * J_or[:3]).tolist(),
                     'TN': f * t_ref, 'split_mean': f * float(J_or[SPLIT_SHELL])}
    # dp2 basis-truncated spirals ([214] t >= 0.90) kept in the fit
    r_bad = np.r_[r_or, rows({'[214]': [0.9, 0.95, 1.0]})]
    n_b = loo_choice(r_bad)[0]
    a_b, J_b = fit_rows(r_bad, n_b)
    abl['include_truncated_214'] = {'predictions': predict(a_b, J_b, n_b), 'J': J_b[:3].tolist(), 'n': n_b}
    r_bad2 = np.r_[rows({'G': [0.0], '[001]': [1.0]}), rows({'[214]': [round(0.05 * k, 2) for k in range(1, 21)]})]
    a_b2, J_b2 = fit_rows(r_bad2, 10)
    abl['dense_214_with_truncated'] = {'predictions': predict(a_b2, J_b2, 10), 'J': J_b2[:3].tolist()}
    # dp3 direction aliasing: high-symmetry / [214]-free designs
    r_001_121 = rows({'G': [0.0], '[001]': [round(0.05 * k, 2) for k in range(1, 21)], '[121]': [round(0.05 * k, 2) for k in range(1, 11)]})
    a3, J3 = fit_rows(r_001_121, 10)
    abl['no_214_ten_shells'] = {'predictions': predict(a3, J3, 10), 'J': J3[:3].tolist()}
    r_001 = rows({'G': [0.0], '[001]': [round(0.05 * k, 2) for k in range(1, 21)]})
    a1, J1 = fit_rows(r_001, 10)
    abl['001_only'] = {'predictions': predict(a1, J1, 10), 'J': J1[:3].tolist()}
    # dp4 commensurate coarse q grid (t multiples of 0.2) with a long-range model
    grid = [0.2, 0.4, 0.6, 0.8]
    r_grid = rows({'G': [0.0], '[001]': grid + [1.0], '[121]': grid, '[214]': grid})
    a4, J4 = fit_rows(r_grid, 10)
    abl['coarse_commensurate_grid'] = {'predictions': predict(a4, J4, 10), 'J': J4[:3].tolist()}
    # dp5 model order
    for n in (3, 6):
        a5, J5 = fit_rows(r_or, n)
        abl[f'{n}_shells'] = {'predictions': predict(a5, J5, n), 'J': J5[:3].tolist(),
                              'split_mean': float(J5[SPLIT_SHELL]) if n > SPLIT_SHELL else 0.0}
    # dp6 spiral reference: q = 0 taken as the Neel state (gauge-equivalent fit, J2 sign flipped)
    a6, J6 = fit_rows(r_or, n_or, 'neel_reference')
    abl['neel_reference'] = {'predictions': predict(a6, J6, n_or, 'neel_reference'), 'J': J6[:3].tolist()}
    a6b, J6b = fit_rows(r_or, n_or, 'lattice_phase')
    abl['lattice_phase_only'] = {'predictions': predict(a6b, J6b, n_or, 'lattice_phase'), 'J': J6b[:3].tolist()}
    # dp7 T_N estimator
    abl['TN_mean_field'] = M.neel_mean_field_K(J_or)
    # dp8 altermagnetic identifiability
    abl['alt_zero_difference'] = dict(base_alt, difference_meV=0.0)
    abl['alt_orthorhombic_fit'] = dict(base_alt, difference_meV=orthorhombic_difference(r_bad, n_or))
    abl['alt_wrong_shell'] = dict(base_alt, distance_A=M.SHELLS[2]['distance'])
    # robustness
    for n in (8, 12):
        ar, Jr = fit_rows(r_or, n)
        rob[f'{n}_shells'] = {'predictions': predict(ar, Jr, n), 'J': Jr[:3].tolist(), 'split_mean': float(Jr[SPLIT_SHELL])}
    r_all = np.where(np.isin(D, ['gamma'] + AVAILABLE_DIRS) & CLEAN)[0]
    ar, Jr = fit_rows(r_all, 10)
    rob['all_clean_available'] = {'predictions': predict(ar, Jr, 10), 'J': Jr[:3].tolist(), 'split_mean': float(Jr[SPLIT_SHELL])}
    coarse_cross = [c['T_cross_K'] for c in tn_coarse['crossings']]
    rob['TN_coarse_L8_L12'] = float(coarse_cross[0])
    # ---------------- write files
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps({'spirals': lab}) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = list(hidden_targets.values()) + [float(J_REF[0]), float(J_REF[1]), float(J_REF[2]), t_ref,
                                              float(J_REF[SPLIT_SHELL])]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob,
                                                                  'oracle_replica': {'n': n_or, 'J': J_or.tolist(), 'alt': base_alt}},
                                                                 indent=1) + '\n')
    cal = {'J_ref': J_REF[:3].tolist(), 'J7_ref': float(J_REF[SPLIT_SHELL]), 'TN_ref': t_ref,
           'ground_ids': GROUND_IDS,
           'oracle': {'n': n_or, 'cv': cv, 'rmse': rmse(pred_or), 'ground': ground_error(pred_or), 'J': J_or[:3].tolist(), 'J7': float(J_or[SPLIT_SHELL]),
                      'TN_MF': M.neel_mean_field_K(J_or)},
           'robust': {k: {'rmse': rmse(v['predictions']), 'ground': ground_error(v['predictions']), 'J': v['J'], 'J7': v['split_mean']} for k, v in rob.items() if isinstance(v, dict)},
           'TN_coarse': rob['TN_coarse_L8_L12'],
           'ablations': {k: ({'rmse': rmse(v['predictions']), 'ground': ground_error(v['predictions']), 'J': v.get('J')} if isinstance(v, dict) and 'predictions' in v else v)
                         for k, v in abl.items()}}
    print(json.dumps(cal, indent=1, default=float))


if __name__ == '__main__':
    main()
