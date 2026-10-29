"""Generate the catalog, lab data, private targets and ablation values from pool.json.

    python3 -B make_task.py
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/wurtzite-ferroelectric-switching'
POOL = json.loads((HERE/'pool.json').read_text())
SYSTEMS = ['AlN', 'AlBN']
RUNS = ['NEB', 'refined_NEB', 'refined_NEB2']
N_FU = 16            # 32-atom supercells: 16 cations + 16 N
BUDGET = 28
COSTS = {'neb_run': 2, 'neb_image': 1, 'berry_phase': 1, 'born_charges': 1, 'vc_relax': 1}
THRESHOLDS = {'barrier_rel_tol': 0.08, 'polarization_abs_tol_C_m2': 0.015, 'saddle_dP_abs_tol_C_m2': 0.025,
              'intermediate_energy_rel_tol': 0.10}
E_CHARGE = 1.602176634e-19


def energies(system, run):
    return np.array([row['energy_eV'] for row in POOL['neb_runs'][system][run]['final_table']])


def raw_p(system):
    b = POOL['berry'][system]
    return np.array([b[str(k)]['P_C_m2'] for k in range(len(b))]), b['0']['mod_C_m2']


def track(raw, q, idx):
    out = [raw[idx[0]]]
    for a, b in zip(idx[:-1], idx[1:]):
        d = raw[b] - out[-1]
        out.append(out[-1] + d - q * np.round(d / q))
    return np.array(out)


def cation_offsets(cell, species, pos):
    """z of each cation minus the mean z of its three basal (in-plane) N neighbours, in A."""
    cell, pos = np.array(cell), np.array(pos)
    shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]) @ cell
    an = [i for i, s in enumerate(species) if s == 'N']
    out = []
    for c in [i for i, s in enumerate(species) if s != 'N']:
        d = (pos[an][None] + shifts[:, None] - pos[c]).reshape(-1, 3)
        near = d[np.argsort(np.linalg.norm(d, axis=1))[:5]]
        basal = near[np.argsort(-np.linalg.norm(near[:, :2], axis=1))[:3]]
        out.append(-basal[:, 2].mean())
    return np.array(out)


def local_minima(e):
    return [k for k in range(1, len(e) - 1) if e[k] < e[k - 1] and e[k] < e[k + 1]]


def main():
    # ---------------- reference values ----------------
    e_aln, e_albn, e_col = energies('AlN', 'refined_NEB2'), energies('AlBN', 'refined_NEB2'), energies('AlBN', 'NEB')
    barrier = {'AlN': (e_aln.max() - e_aln[0]) / N_FU * 1000, 'AlBN': (e_albn.max() - e_albn[0]) / N_FU * 1000}
    collective = (e_col.max() - e_col[0]) / N_FU * 1000
    saddle = int(np.argmax(e_albn))                       # 0-based = Berry folder index
    mins = local_minima(e_albn)
    inter = min(mins, key=lambda k: e_albn[k])
    inter_e = (e_albn[inter] - e_albn[0]) / N_FU * 1000
    img = POOL['neb_runs']['AlBN']['refined_NEB2']['images']
    cell = POOL['neb_runs']['AlBN']['refined_NEB2']['cell_A']
    off0 = cation_offsets(cell, img[0]['species'], img[0]['positions_A'])
    off_i = cation_offsets(cell, img[inter]['species'], img[inter]['positions_A'])
    inverted = int(np.sum(np.sign(off_i) != np.sign(off0)))
    assert np.abs(off_i).min() > 0.3, 'intermediate must have no near-planar cation'
    ps, dP_full = {}, {}
    for s in SYSTEMS:
        raw, q = raw_p(s)
        t = track(raw, q, list(range(len(raw))))
        dP_full[s] = t - t[0]
        ps[s] = -(t[-1] - t[0]) / 2
    # AlN: the layered-hexagonal midpoint image is the reference; must agree with half the change
    raw, q = raw_p('AlN')
    t = track(raw, q, list(range(len(raw))))
    assert abs((t[0] - t[12]) - ps['AlN']) < 1e-3 and abs(raw[12]) < 1e-3
    saddle_dP = float(dP_full['AlBN'][saddle])
    ref = {'barrier_meV_per_fu': barrier, 'AlBN_collective_barrier_meV_per_fu': collective,
           'spontaneous_polarization_C_m2': ps, 'AlBN_saddle_polarization_change_C_m2': saddle_dP,
           'AlBN_intermediate': {'energy_meV_per_fu': inter_e, 'inverted_cations': inverted},
           'AlBN_switching_mode': 'non-collective'}
    hidden = {'reference': ref, 'thresholds': THRESHOLDS,
              'notes': {'saddle_image_0based': saddle, 'intermediate_image_0based': inter,
                        'AlBN_local_minima_0based': mins}}

    # ---------------- ablation values (real wrong paths) ----------------
    A = {}
    A['barrier_AlBN_first_stage'] = collective
    A['collective_from_AlN'] = barrier['AlN']
    A['per_atom'] = {'AlN': barrier['AlN'] / 2, 'AlBN': barrier['AlBN'] / 2, 'collective': collective / 2,
                     'intermediate': inter_e / 2}
    A['per_supercell'] = {'AlN': barrier['AlN'] * N_FU, 'AlBN': barrier['AlBN'] * N_FU, 'collective': collective * N_FU,
                          'intermediate': inter_e * N_FU}
    raw_n, q_n = raw_p('AlN')
    raw_b, q_b = raw_p('AlBN')
    A['raw_berry'] = {'AlN': float(raw_n[0] - raw_n[12]), 'AlBN': float(-(raw_b[-1] - raw_b[0]) / 2),
                      'saddle': float(raw_b[saddle] - raw_b[0])}
    t2 = track(raw_n, q_n, list(range(0, 13, 2)))
    t6 = track(raw_b, q_b, list(range(0, 75, 6)) + [74])
    A['sparse_tracking'] = {'AlN_stride2_half_path': float(t2[0] - t2[-1]), 'AlBN_stride6': float(-(t6[-1] - t6[0]) / 2)}
    zs = {}
    for r in POOL['born_charges']['AlBN']:
        zs.setdefault(r['atom'], []).append(r['Z33_e'])
    zav = {k: float(np.mean(v)) for k, v in zs.items()}

    def bec_full(s):
        b = POOL['berry'][s]
        n = len(b)
        vol = abs(np.linalg.det(np.array(b['0']['cell_A'])))
        dz = np.array(b[str(n - 1)]['positions_A'])[:, 2] - np.array(b['0']['positions_A'])[:, 2]
        guess = sum(zav[x] * d for x, d in zip(b['0']['species'], dz)) * E_CHARGE * 1e20 / vol
        raw, q = raw_p(s)
        d = raw[-1] - raw[0]
        d += q * np.round((guess - d) / q)
        return float(-d / 2), float(guess)
    A['bec_full_path'] = {s: bec_full(s)[0] for s in SYSTEMS}
    # the metastable half-switched intermediate taken as the nonpolar reference
    A['half_switched_reference_AlBN'] = float(-dP_full['AlBN'][inter])
    A['off_by_one_saddle'] = {'plus': float(dP_full['AlBN'][saddle + 1]), 'minus': float(dP_full['AlBN'][saddle - 1])}
    first_max = [k for k in range(1, 74) if e_albn[k] > e_albn[k - 1] and e_albn[k] > e_albn[k + 1]][0]
    A['first_subbarrier'] = {'barrier': (e_albn[first_max] - e_albn[0]) / N_FU * 1000, 'saddle_dP': float(dP_full['AlBN'][first_max]),
                             'first_minimum_energy': (e_albn[mins[0]] - e_albn[0]) / N_FU * 1000}
    A['collective_saddle_dP'] = float(-ps['AlBN'])      # midpoint of a collective path: half the reversal
    # Robustness: defensible variants
    from_stage1_aln = (energies('AlN', 'NEB').max() - energies('AlN', 'NEB')[0]) / N_FU * 1000
    # spline-located saddle: cubic through the three images around the maximum
    k = saddle
    c = np.polyfit([k - 1, k, k + 1], e_albn[k - 1:k + 2], 2)
    xs = -c[1] / (2 * c[0])
    dp_interp = float(np.interp(xs, np.arange(75), dP_full['AlBN']))
    R = {'barrier_AlN_first_stage_run': from_stage1_aln, 'saddle_dP_parabolic_position': dp_interp,
         'Ps_AlN_midpoint_reference': float(t[0] - t[12]),
         'barrier_AlBN_vs_relaxed_endpoint': (e_albn.max() - e_col.min()) / N_FU * 1000}

    # ---------------- catalog (agent-visible) ----------------
    neb_runs = {s: {r: {'images': len(POOL['neb_runs'][s][r]['final_table'])} for r in RUNS} for s in SYSTEMS}
    catalog = {
        'source': 'Released DFT data of Calderon et al., Science 380, 1034 (2023), Zenodo record 7916565 (CC BY 4.0). '
                  'Quantum ESPRESSO 6.7, PBE, PseudoDojo pseudopotentials, 32-atom wurtzite supercells.',
        'systems': {'AlN': 'pure AlN supercell', 'AlBN': 'B-substituted AlN supercell, (Al,B)N'},
        'budget': BUDGET, 'purchase_costs': COSTS,
        'requests': {
            'neb_run': {'example': {'kind': 'neb_run', 'system': 'AlBN', 'run': 'refined_NEB'},
                        'returns': 'Settings and the final iteration of one archived neb.x run: path and engine settings, '
                                   'the image energy/error table (images numbered 1..N as printed by neb.x), recent '
                                   'activation-energy lines, the convergence message if any, and the end of the log and '
                                   'scheduler output.'},
            'neb_image': {'example': {'kind': 'neb_image', 'system': 'AlBN', 'run': 'NEB', 'image': 13},
                          'returns': 'Final structure of one image of one run (cell, species, Cartesian positions in A).'},
            'berry_phase': {'example': {'kind': 'berry_phase', 'system': 'AlN', 'folder': '00_image'},
                            'returns': 'One archived Berry-phase calculation (scf + nscf with lberry, gdir = 3): run settings, '
                                       'the scf total energy, the structure, and the POLARIZATION CALCULATION section '
                                       'of the nscf output as printed by pw.x.'},
            'born_charges': {'example': {'kind': 'born_charges', 'system': 'AlBN'},
                             'returns': "The authors' file 'Born effective charge (AlBN).csv' (Z33 per atom)."},
            'vc_relax': {'example': {'kind': 'vc_relax', 'key': 'AlN_primitive'},
                         'returns': 'Final cell and coordinates, energy and pressure of one archived vc-relax run.'}},
        'neb_runs': neb_runs,
        'neb_run_note': "Run names are the authors' folder names under <system>/NEB/.",
        'berry_phase_folders': {s: sorted(f'{int(k):02d}_image' for k in POOL['berry'][s]) for s in SYSTEMS},
        'berry_phase_note': "Folder names are the authors' names under <system>/polarization/nscf/; the archive does not "
                            'say which NEB run they were taken from.',
        'vc_relax_keys': {'AlN_primitive': 'AlN/vc-relax (4-atom cell)',
                          'AlN_32atom_initial_state': 'AlN/NEB/initial_state',
                          'AlBN_32atom_initial_state': 'AlBN/NEB/initial_state',
                          'AlBN_192atom': 'AlBN/443_vc-relax (larger B-substituted supercell)'},
        'born_charges_systems': ['AlBN']}

    # ---------------- lab data ----------------
    lab = {'neb_runs': {}, 'neb_images': {}, 'berry': {}, 'vc_relax': POOL['vc_relax'], 'born_charges': POOL['born_charges']}
    for s in SYSTEMS:
        lab['neb_runs'][s], lab['neb_images'][s] = {}, {}
        for r in RUNS:
            info = POOL['neb_runs'][s][r]
            view = {k: info[k] for k in ('settings', 'engine', 'final_table', 'n_iteration_blocks', 'last_iteration_label',
                                         'activation_energy_forward_history_eV', 'convergence_message', 'tail',
                                         'scheduler_out', 'cell_A')}
            if 'earlier_log_activation_history_eV' in info:
                view['earlier_log_activation_history_eV'] = info['earlier_log_activation_history_eV']
            lab['neb_runs'][s][r] = view
            lab['neb_images'][s][r] = [{'species': im['species'], 'positions_A': im['positions_A']} for im in info['images']]
        lab['berry'][s] = {}
        for k, b in POOL['berry'][s].items():
            lab['berry'][s][f'{int(k):02d}_image'] = {
                'nscf_settings': b['nscf_settings'], 'scf_settings': b['scf_settings'],
                'scf_total_energy_Ry': b['scf_total_energy_Ry'], 'cell_A': b['cell_A'], 'species': b['species'],
                'positions_A': b['positions_A'], 'polarization_output': b['text']}

    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps(lab) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = [barrier['AlN'], barrier['AlBN'], collective, ps['AlN'], ps['AlBN'], saddle_dP, inter_e, inverted]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': A, 'robustness': R}, indent=1) + '\n')
    print(json.dumps({'reference': ref, 'notes': hidden['notes'], 'ablations': A, 'robustness': R,
                      'bec_full_guess': {s: bec_full(s)[1] for s in SYSTEMS}, 'Z_avg': zav}, indent=1, default=float))


if __name__ == '__main__':
    main()
