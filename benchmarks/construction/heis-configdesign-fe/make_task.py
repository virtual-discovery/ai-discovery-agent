"""Generate the staged, lab, private and ablation data for heis-configdesign-fe from pool.json.

    python3 -B benchmarks/construction/heis-configdesign-fe/make_task.py
"""
import json
import numpy as np
from pathlib import Path
import model as M

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/heis-configdesign-fe'
P = M.POOL
SP, RS = P['spirals'], P['random_states']
DIRS = np.array([M.direction_of(s['q_reciprocal_lattice']) for s in SP])
Q = np.array([s['q_reciprocal_lattice'] for s in SP])
TH = np.array([s['cone_angle_deg'] for s in SP])
Y = np.array([s['energy_meV_per_atom'] for s in SP])
B = 2 * np.pi * np.linalg.inv(M.A).T
QA = np.linalg.norm(Q @ B, axis=1)
AVAILABLE_DIRS = ['[001]', '[111]', '[113]', '[121]']
NOMINAL_CONE = 30.0
THRESHOLDS = {'spiral_rmse_meV': 2.5, 'J1_window_meV': [18.0, 24.5], 'stiffness_rel_tol': 0.15,
              'random_mean_tol_meV': 10.0, 'random_rmse_meV': 12.0}


def by_q(direction):
    idx = np.where(DIRS == direction)[0]
    return idx[np.argsort(QA[idx])]


def stiffness_from(idx_by_dir, kmax=2):
    vals = []
    for idx in idx_by_dir:
        idx = idx[:kmax]
        s2 = np.sin(np.radians(TH[idx])) ** 2
        vals.append(np.sum(Y[idx] / s2 * QA[idx] ** 2) / np.sum(QA[idx] ** 4))
    return float(np.mean(vals))


def unit(m):
    m = np.array(m)
    return (m / np.linalg.norm(m, axis=1, keepdims=True)).round(6).tolist()


def main():
    spiral_keys = {}
    catalog_spirals = []
    lab_spirals = {}
    for d in AVAILABLE_DIRS:
        for k, i in enumerate(by_q(d)):
            key = f'{d}-{k:02d}'
            spiral_keys[key] = i
            catalog_spirals.append({'key': key, 'direction': d, 'q_index': k, 'q_reciprocal_lattice': Q[i].tolist(),
                                    'cone_angle_target_deg': NOMINAL_CONE})
            lab_spirals[key] = {'q_reciprocal_lattice': Q[i].tolist(), 'energy_meV_per_atom': float(Y[i]),
                                'cone_angle_deg': float(TH[i]), 'moment_muB': SP[i]['moment_muB']}
    hidden_spiral = {}
    predict_spirals = []
    for k, i in enumerate(by_q('[214]')):
        sid = f'S214-{k + 1:02d}'
        hidden_spiral[sid] = float(Y[i])
        predict_spirals.append({'id': sid, 'direction': '[214]', 'q_reciprocal_lattice': Q[i].tolist(),
                                'cone_angle_target_deg': NOMINAL_CONE})
    catalog_random, lab_random, predict_random, hidden_random = [], {}, [], {}
    for s in RS:
        geometry = {'supercell': s['supercell'], 'cell_A': s['cell_A'], 'fractional_positions': s['fractional_positions'],
                    'spin_directions': unit(s['moments_muB'])}
        if s['supercell'] == '334':
            rid = f'R334-{len(predict_random) + 1:02d}'
            hidden_random[rid] = s['energy_meV_per_atom']
            predict_random.append({'id': rid, **geometry})
        else:
            catalog_random.append({'label': s['label'], **geometry})
            lab_random[s['label']] = {'energy_meV_per_atom': s['energy_meV_per_atom'], 'moments_muB': s['moments_muB']}
    catalog = {
        'material': 'bcc Fe, noncollinear PBE (VASP, constrained moments), released by Li, von Malottki and Rignanese (2026)',
        'energy_reference': 'All energies are meV/atom relative to the collinear ferromagnet computed with the same protocol.',
        'spiral_cell': {'lattice_constant_A': 2.87, 'primitive_cell_A': P['primitive_cell_A'],
                        'q_coordinates': 'reciprocal-lattice coordinates of this primitive cell'},
        'purchase_costs': {'spiral': 1, 'random_state': 1}, 'budget': 60,
        'spirals_available': catalog_spirals, 'spirals_to_predict': predict_spirals,
        'random_states_available': catalog_random, 'random_states_to_predict': predict_random}
    # Reference answers and ablation values
    tr = np.r_[np.where(DIRS == 'gamma')[0], *[by_q(d) for d in AVAILABLE_DIRS]]
    c0, J = M.fit(M.spiral_design(Q[tr], TH[tr], 8), Y[tr], intercept=True)
    t214 = by_q('[214]')
    d_ref = stiffness_from([by_q(d) for d in AVAILABLE_DIRS])
    Xr = np.array([M.random_design(s, 8) for s in RS])
    yr = np.array([s['energy_meV_per_atom'] for s in RS])
    sc = np.array([s['supercell'] for s in RS])
    hidden = {'spiral_targets': hidden_spiral, 'random_targets': hidden_random, 'thresholds': THRESHOLDS,
              'reference': {'J_meV': J.tolist(), 'J1_meV': float(J[0]), 'stiffness_meV_A2': d_ref,
                            'random_mean_meV': float(np.mean(list(hidden_random.values())))}}

    def spiral_pred(train_idx, n, cone=True):
        th_train = TH[train_idx] if cone else np.full(len(train_idx), 90.0)
        a, Jf = M.fit(M.spiral_design(Q[train_idx], th_train, n), Y[train_idx], intercept=True)
        pred = a + M.spiral_design(Q[t214], np.full(len(t214), NOMINAL_CONE if cone else 90.0), n) @ Jf
        return {f'S214-{k + 1:02d}': float(v) for k, v in enumerate(pred)}, Jf

    def random_pred(mask, n, X=None):
        X = Xr if X is None else X
        a, Jf = M.fit(X[mask, :n], yr[mask], intercept=True)
        return {f'R334-{k + 1:02d}': float(v) for k, v in enumerate(a + X[sc == '334', :n] @ Jf)}, Jf

    gamma = np.where(DIRS == 'gamma')[0]
    abl = {}
    abl['spiral_001_only'], J001 = spiral_pred(np.r_[gamma, by_q('[001]')], 8)
    abl['J1_001_only'] = float(J001[0])
    _, Jflat = spiral_pred(tr, 8, cone=False)
    abl['J1_flat_spiral'] = float(Jflat[0])
    abl['stiffness_flat_spiral'] = d_ref * 0.25
    abl['spiral_four_shells'], _ = spiral_pred(tr, 4)
    abl['random_from_spiral_model'] = {f'R334-{k + 1:02d}': float(v) for k, v in enumerate(Xr[sc == '334', :8] @ J)}
    abl['random_341_only'], _ = random_pred(sc == '341', 4)
    abl['random_221_three_shells'], _ = random_pred(sc == '221', 3)
    abl['random_wrong_lattice'] = {k: float(yr[np.isin(sc, ['236', '341'])].mean()) for k in hidden_random}
    abl['stiffness_no_2pi'] = d_ref * (2 * np.pi) ** 2
    _, Jrand = random_pred(np.isin(sc, ['236', '341']), 6)
    abl['J1_from_random_states'] = float(Jrand[0])
    abl['stiffness_from_truncated_model'] = float(M.stiffness(spiral_pred(tr, 10)[1]))
    robust = {'spiral_ten_shells': spiral_pred(tr, 10)[0], 'J1_ten_shells': float(spiral_pred(tr, 10)[1][0]),
              'random_eight_shells': random_pred(np.isin(sc, ['221', '341', '236']), 8)[0],
              'stiffness_kmax3': stiffness_from([by_q(d) for d in AVAILABLE_DIRS], kmax=3)}
    # write files
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps({'spirals': lab_spirals, 'random_states': lab_random}) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = list(hidden_spiral.values()) + list(hidden_random.values()) + [float(J[0]), d_ref, hidden['reference']['random_mean_meV']]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': robust}, indent=1) + '\n')
    rmse = lambda pred, ref: float(np.sqrt(np.mean([(pred[k] - ref[k]) ** 2 for k in ref])))
    print(json.dumps({'J1_ref': J[0], 'D_ref': d_ref, 'random_mean_ref': hidden['reference']['random_mean_meV'],
                      'abl_spiral_001_rmse': rmse(abl['spiral_001_only'], hidden_spiral), 'abl_J1_001': abl['J1_001_only'],
                      'abl_J1_flat': abl['J1_flat_spiral'], 'abl_four_shell_rmse': rmse(abl['spiral_four_shells'], hidden_spiral),
                      'abl_random_spiral_model_mean_err': np.mean(list(abl['random_from_spiral_model'].values())) - hidden['reference']['random_mean_meV'],
                      'abl_random_341_mean_err': np.mean(list(abl['random_341_only'].values())) - hidden['reference']['random_mean_meV'],
                      'abl_random_221_3_mean_err': np.mean(list(abl['random_221_three_shells'].values())) - hidden['reference']['random_mean_meV'],
                      'abl_random_wrong_lattice_mean_err': np.mean(list(abl['random_wrong_lattice'].values())) - hidden['reference']['random_mean_meV'],
                      'abl_J1_random': abl['J1_from_random_states'], 'abl_D_truncated_model': abl['stiffness_from_truncated_model'],
                      'rob_spiral10_rmse': rmse(robust['spiral_ten_shells'], hidden_spiral),
                      'rob_random8_mean_err': np.mean(list(robust['random_eight_shells'].values())) - hidden['reference']['random_mean_meV'],
                      'rob_random8_rmse': rmse(robust['random_eight_shells'], hidden_random), 'rob_D_kmax3': robust['stiffness_kmax3']},
                     indent=1, default=float))


if __name__ == '__main__':
    main()
