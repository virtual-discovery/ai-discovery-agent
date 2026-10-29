"""Write the staged, lab, private and ablation files of mnwo4-halfmoon-frustration.

    python3 -B make_task.py CACHE      (after extract.py and calibrate.py; run in the runtime image)

Inputs: CACHE/extracted.npz, CACHE/model.json, calibration_fits.json. Prints the calibration table.
"""
import json
import sys
from pathlib import Path
import numpy as np
from model import truth, map14_blocks, ko_blocks, CELL, TAU, NAX, Model

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/mnwo4-halfmoon-frustration'
CACHE = Path(sys.argv[1])
X = np.load(CACHE/'extracted.npz')
CAL = json.loads((HERE/'calibration_fits.json').read_text())
COSTS = {'spectrum': 1, 'map': 4, 'reduced_model_map': 4, 'zone_slice': 1}
BUDGET = 33

# ---- thresholds (see README calibration table) ----
TH = {'k_tolerance': 0.020,
      'map14_min_r': 0.85,
      'reduced_min_r': {'5d': 0.68, '5e': 0.88}}


def r4(a):
    a = np.asarray(a, float)
    return np.vectorize(lambda x: float('%.4g' % x))(a)


def pearson(a, b):
    a = np.asarray(a, float).ravel() - np.mean(a); b = np.asarray(b, float).ravel() - np.mean(b)
    return float(a @ b / np.sqrt((a @ a) * (b @ b)))


def kdist(k1, k2):
    return min(float(np.abs((np.asarray(k1) - s * np.asarray(k2)) - np.round(np.asarray(k1) - s * np.asarray(k2))).max())
               for s in (1, -1))


def zone_cube(q, S, n1, n2, n3, o1, o2, o3):
    cube = np.full((n1, n2, n3, 3), np.nan)
    i = np.round(q[:, 0] * n1).astype(int) + o1
    j = np.round(q[:, 1] * n2).astype(int) + o2
    m = np.round(q[:, 2] * n3).astype(int) + o3
    cube[i, j, m] = S[:, :3]
    assert not np.isnan(cube).any()
    return cube


def main():
    shells, J0, D0, dist = truth(CACHE/'model.json')
    M = Model(shells)
    kref, lam = M.lt(J0)
    dEref = (-lam + D0 / 2) - M.af1_energy(J0, D0)
    lab = TASK/'environment/lab'
    # ---- lab data (archived values, 4 significant digits, float32 storage) ----
    data = {'spectra_q': X['spectra_q'], 'spectra_S': r4(X['spectra_S']).astype(np.float32)}
    for T in (20, 30, 60):
        data['map_T%d' % T] = r4(X['map_T%d' % T][186:, :, :3]).astype(np.float32)   # k >= 0 half
    for m in ('5b', '5c'):
        data['ko_' + m] = r4(X['ko_' + m][:, :, :3]).astype(np.float32)
    data['zone_36_24_36'] = r4(zone_cube(X['bz_36_24_36_q'], X['bz_36_24_36_S'], 36, 24, 36, 17, 11, 17)).astype(np.float32)
    data['zone_48_24_24'] = r4(zone_cube(X['bz_48_24_24_q'], X['bz_48_24_24_S'], 48, 24, 24, 23, 11, 11)).astype(np.float32)
    np.savez_compressed(lab/'data.npz', **data)
    dE = float(X['spectra_E'][1] - X['spectra_E'][0])
    meta = {'costs': COSTS, 'spectra_dE': dE,
            'map_grid': 'q = (-l/2, k, l) r.l.u.; k_i = i * %.10f (i = 0 .. 186, the k >= 0 half; S(-q) = S(q)); '
                        'l_j = -1.776 + j * 0.0111 (j = 0 .. 320)' % (4.66 / 372),
            'ko_grid': 'q = (-l/2, k, l) r.l.u.; k_i = -2 + 0.025 i (i = 0 .. 160); l_j = j / 90 (j = 0 .. 160)',
            'zone_grids': {'36x24x36': 'h_i = (i - 17)/36 (i = 0 .. 35), l_m = (m - 17)/36 (m = 0 .. 35), k = (k_index - 11)/24',
                           '48x24x24': 'h_i = (i - 23)/48 (i = 0 .. 47), l_m = (m - 11)/24 (m = 0 .. 23), k = (k_index - 11)/24'}}
    (lab/'data_meta.json').write_text(json.dumps(meta, indent=1) + '\n')
    # ---- private targets ----
    ref14 = map14_blocks(X['map_T14'][:, :, 3])
    refred = {m: ko_blocks(X['ko_' + m][:, :, 3]) for m in ('5d', '5e')}
    hidden = {'map14_shape': [27, 29], 'reduced_shape': [32, 32],
              'map14_ref': ref14.tolist(), 'reduced_refs': {m: v.tolist() for m, v in refred.items()},
              'propagation_vector_rlu': np.asarray(kref).tolist(), 'half_moon_origin': 'competing_3d_exchange',
              'thresholds': TH,
              'reference': {'J_meV': J0.tolist(), 'D_meV': D0, 'spiral_minus_af1_meV_per_Mn': dEref,
                            'note': 'J, D: archived jfile/kfile (mRy) x 13.605693; k and energy difference computed from them'}}
    (lab/'private/hidden.json').write_text(json.dumps(hidden) + '\n')
    values = [*J0.tolist(), D0, *np.asarray(kref).tolist()[::2], dEref] + ref14.ravel().tolist() + \
        [v for m in refred for v in refred[m].ravel().tolist()]
    (lab/'private/hidden_values.json').write_text(json.dumps(values) + '\n')

    # ---- ablations and robustness from calibration fits ----
    def sub(name):
        c = CAL[name]
        out = {}
        if 'J' in c:
            out['exchange_meV'] = {'J%d' % (i + 1): float(v) for i, v in enumerate(c['J'])}
            out['anisotropy_D_meV'] = float(c['D'])
            out['propagation_vector_rlu'] = [float(x) for x in c['k']]
            out['spiral_minus_af1_meV_per_Mn'] = float(c['dE'])
        out['map_14K_prediction'] = r4(c['map14']).tolist()
        out['reduced_model_predictions'] = {m: r4(v).tolist() for m, v in c['reduced'].items()}
        return out
    abl = {n: sub(n) for n in CAL if n.startswith('abl_')}
    abl['k_from_48_box_peak'] = {'propagation_vector_rlu': [-11 / 48, 0.5, 11 / 24]}
    abl['k_experimental_kic'] = {'propagation_vector_rlu': [-0.214, 0.5, 0.457]}
    rob = {n: sub(n) for n in CAL if n.startswith('robust_')}
    nonrecall = {n: sub(n) for n in CAL if n.startswith('nonrecall_')}
    rob['reference_model'] = sub('reference_model')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob, 'nonrecall_seeds': nonrecall}) + '\n')
    # ---- public budget ----
    for p, extra in ((lab/'public.json', {'reserved_validation_charge': 0, 'catalog': '/app/catalog.json'}),
                     (TASK/'tests/public.json', {'min_report_chars': 400})):
        p.write_text(json.dumps({'task_id': 'mnwo4-halfmoon-frustration', 'solver_budget': BUDGET, **extra}, indent=1) + '\n')
    write_catalog(shells, dist)
    # ---- calibration table ----
    table = {}
    for name, c in CAL.items():
        if 'map14' not in c:
            continue
        row = {}
        if 'J' in c:
            row.update({'J1': c['J'][0], 'J6': c['J'][5], 'J9': c['J'][8], 'J3': c['J'][2], 'J11': c['J'][10],
                        'D': c['D'], 'k': [round(x, 4) for x in c['k']], 'k_dist': kdist(c['k'], kref), 'dE': c['dE']})
        row['r14'] = pearson(c['map14'], ref14)
        for m in refred:
            row['r_' + m] = pearson(c['reduced'][m], refred[m])
        table[name] = row
    table['reference_values'] = {'J1': J0[0], 'J6': J0[5], 'J9': J0[8], 'J3': J0[2], 'J11': J0[10], 'D': D0,
                                 'k': np.asarray(kref).tolist(), 'dE': dEref}
    (HERE/'calibration_table.json').write_text(json.dumps(table, indent=1, default=float) + '\n')
    for n, r in table.items():
        print(n, json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, default=float))


def write_catalog(shells, dist):
    paths = []
    for k, bonds in enumerate(shells):
        paths.append({'label': 'J%d' % (k + 1), 'distance_A': round(float(dist[k]), 3),
                      'bonds': [{'from_site': s + 1, 'to_site': s2 + 1, 'cell_offset': [n1, n2, n3]} for s, s2, n1, n2, n3 in bonds]})
    cat = {
        'material': 'MnWO4 (wolframite, P2/c); two Mn2+ per cell forming zigzag chains along c',
        'model_source': 'Atomistic spin-dynamics model of Natsch, Tosic, Soh and Spaldin (Materials Cloud, record p8297-g9n39; '
                        'arXiv:2503.23704), simulated with UppASD 6.0.1. All purchasable data are the authors\' archived outputs '
                        'of that model; its exchange constants and anisotropy are not given.',
        'structure': {'cell_vectors_A': CELL.round(10).tolist(),
                      'sites_fractional': {'1': TAU[0].tolist(), '2': TAU[1].tolist()},
                      'moment_muB': 5.0, 'g_factor': 2.0,
                      'easy_axis_cartesian': [0.81915, 0.0, 0.57358]},
        'hamiltonian': 'H = - sum_{i != j} J_ij e_i . e_j + D sum_i (e_i . n)^2, unit vectors e, every ordered pair counted '
                       '(each bond twice), J and D in meV, n = easy_axis_cartesian.',
        'exchange_paths': paths,
        'exchange_path_note': 'bond from site i in cell R to site j in cell R + cell_offset; vector = (cell_offset + r_j - r_i) . cell; '
                              'every ordered pair is listed; labels follow increasing Mn-Mn distance.',
        'af1_state': {'description': 'collinear up-up-down-down state along n with k = (-1/4, 1/2, 1/2); spin of site s in cell '
                                     '(n1, n2, n3) = sigma * n',
                      'sigma': 'site 1: [+1, +1, -1, -1][n1 mod 4] * (-1)^(n2 + n3); site 2: [+1, -1, -1, +1][n1 mod 4] * (-1)^(n2 + n3)'},
        'simulations': {
            'dynamics': 'stochastic Landau-Lifshitz-Gilbert spin dynamics, damping 0.05, time step 1 fs, classical spins of 5 muB',
            'spectra': '1.5 K, 48x48x48 periodic supercell initialised in AF1 (30 000 initial steps); dynamical structure factor '
                       'S_aa(q, E) along the path below, components in the Cartesian frame of cell_vectors_A, arbitrary units',
            'maps': '36x24x36 periodic supercell cooled in spin-dynamics steps from 60 K (2 K steps to 40 K, then 1 K steps), '
                    '30 000 initial + 40 000 measurement steps per temperature; static structure factor S_aa(q) in the plane '
                    'q = (-l/2, k, l), arbitrary units',
            'reduced_model_maps': '20 K, same supercell and protocol, for reduced Hamiltonians (see reduced_models)',
            'zone_slices': '11 K, full commensurate q grid of a 36x24x36 or a 48x24x24 periodic supercell, static S_aa(q)'},
        'spectrum_path': {'points': 101, 'q_of_index': 'index 0-50: q = (0.5 - 0.005 i, 0.5, 1 - 0.01 i); index 51-100: '
                                                       'q = (0.25, 0.5 + 0.01 (i - 50), 0.5)'},
        'reduced_models': {'5b': 'all 12 exchange paths, D = 0', '5c': 'J2, J5, J8, J12 removed (D kept)',
                           '5d': 'J2, J5, J8, J12 and J9 removed (D kept)  [withheld: predict]',
                           '5e': 'J2, J5, J8, J12 and J3, J7, J10, J11 removed (D kept)  [withheld: predict]'},
        'requests': {
            'spectrum': {'request': {'kind': 'spectrum', 'index': '<0-100>'}, 'cost': COSTS['spectrum'],
                         'returns': 'q, energy grid and archived Sxx, Syy, Szz (E < 5 meV)'},
            'map': {'request': {'kind': 'map', 'temperature_K': '<20 | 30 | 60>'}, 'cost': COSTS['map'],
                    'returns': 'archived Sxx, Syy, Szz on the k >= 0 half of the plane grid (187 x 321)'},
            'reduced_model_map': {'request': {'kind': 'reduced_model_map', 'model': '<5b | 5c>'}, 'cost': COSTS['reduced_model_map'],
                                  'returns': 'archived 20 K Sxx, Syy, Szz on the reduced-model grid (161 x 161)'},
            'zone_slice': {'request': {'kind': 'zone_slice', 'supercell': '<36x24x36 | 48x24x24>', 'k_index': '<0-23>'},
                           'cost': COSTS['zone_slice'], 'returns': 'archived 11 K Sxx, Syy, Szz at k = (k_index - 11)/24 on the full (h, l) grid'}},
        'budget': BUDGET,
        'to_predict': {
            'map_14K': '14 K map of the full model in the plane q = (-l/2, k, l): 27 x 29 values; row r, column c is the average of '
                       'the archived |S| over the 7 x 11 grid points centred at k = 7 r dk, l = 11 (c - 14) dl of the map grid '
                       '(dk = 4.66/372, dl = 0.0111), i.e. k from 0 to 2.2795 and l from -1.7094 to 1.7094',
            'reduced_models': '20 K maps of reduced models 5d and 5e: 32 x 32 values each; row r, column c is the average of the '
                              'archived |S| over the 5 x 5 grid points k_i, l_j with i = 5r .. 5r+4, j = 5c .. 5c+4 of the reduced-model grid'},
        'structure_factor_definition': '|S| = sqrt(Sxx^2 + Syy^2 + Szz^2) of the archived UppASD components (as in the purchasable maps)'}
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(cat, indent=1) + '\n')


if __name__ == '__main__':
    main()
