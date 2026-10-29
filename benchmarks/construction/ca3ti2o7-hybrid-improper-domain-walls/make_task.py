"""Write the catalog, private targets and ablation values for ca3ti2o7-hybrid-improper-domain-walls.

Inputs: the Zenodo bulk structures (fetch.sh) and the construction relaxations in runs/*.json
(produced by run_wall.py through queue.sh with jobs1-3.txt, i.e. the lab's own NEP code path).
"""
import json
import os
from pathlib import Path
import numpy as np

HERE = Path(__file__).parent
TASK = HERE.parents[1]/'tasks'/HERE.name
RAW = Path(os.environ.get('RAW', Path.home()/'.cache/ai-theorist'/HERE.name/'raw'))
import sys
sys.path.insert(0, str(HERE))
from xyzio import read_xyz   # noqa: E402

R = {p.stem: json.loads(p.read_text()) for p in (HERE/'runs').glob('*.json')}
E = {k: v['E_DW_mJm2'] for k, v in R.items()}

A_ORTH = 5.456412031232137          # polar axis of the bulk cell (charged-wall normal)
D_PER_CELL_NM = A_ORTH / 2 / 10     # domain width added per orthorhombic cell (two walls per cell)


def slope(prefix, ns):
    x = np.array(ns, float) * D_PER_CELL_NM
    y = np.array([E['%s_n%d' % (prefix, n)] if '%s_n%d' % (prefix, n) in E else E[prefix % n] for n in ns])
    return float(np.polyfit(x, y, 1)[0])


def slope_tags(tags, ns):
    x = np.array(ns, float) * D_PER_CELL_NM
    return float(np.polyfit(x, [E[t] for t in tags], 1)[0])


ref = {
    'E_neutral_psi_wall_mJ_m2': E['psi110neu_n16'],
    'E_neutral_phi_wall_mJ_m2': E['phi110neu_n16'],
    'E_layer_psi_wall_mJ_m2': 0.5 * (E['psi001_rs'] + E['psi001_rs_nz8']),
    'charged_psi_wall_slope_mJ_m2_per_nm': slope('psi110chg', [8, 12, 16, 20, 24]),
    'psi_wall_profile': 'ising',
    'phi_wall_profile': 'rotational',
}
rel = {'E_neutral_psi_wall_mJ_m2': 0.05, 'E_neutral_phi_wall_mJ_m2': 0.05, 'charged_psi_wall_slope_mJ_m2_per_nm': 0.15}
tol = {k: [round(ref[k] * (1 - r), 2), round(ref[k] * (1 + r), 2)] for k, r in rel.items()}
tol['E_layer_psi_wall_mJ_m2'] = [round(ref['E_layer_psi_wall_mJ_m2'] - 4.0, 2), round(ref['E_layer_psi_wall_mJ_m2'] + 4.0, 2)]

ablations = {
    'swap_psi_phi_construction': {'E_neutral_psi_wall_mJ_m2': E['phi110neu_n8'], 'E_neutral_phi_wall_mJ_m2': E['psi110neu_n8'],
                                  'E_layer_psi_wall_mJ_m2': E['phi001_rs'],
                                  'charged_psi_wall_slope_mJ_m2_per_nm': slope('phi110chg', [8, 12, 16, 20])},
    'antiphase_boundary_as_psi_wall': {'E_neutral_psi_wall_mJ_m2': E['apb110neu_n8']},
    'polar_axis_swapped': {'E_neutral_psi_wall_mJ_m2': E['psi110chg_n8'], 'E_neutral_phi_wall_mJ_m2': E['phi110chg_n8'],
                           'charged_psi_wall_slope_mJ_m2_per_nm': slope('psi110neu', [4, 8, 12, 16])},
    'layer_wall_in_bilayer': {'E_layer_psi_wall_mJ_m2': E['psi001_bl']},
    'charged_slope_from_short_cells': {'charged_psi_wall_slope_mJ_m2_per_nm': slope('psi110chg', [4, 8])},
    'charged_wall_treated_as_converged': {'charged_psi_wall_slope_mJ_m2_per_nm': 0.0},
    'slope_per_supercell_length': {'charged_psi_wall_slope_mJ_m2_per_nm': slope('psi110chg', [8, 12, 16, 20, 24]) / 2},
    'short_range_NEP': {'E_neutral_psi_wall_mJ_m2': E['psi110neu_n8_NEP'], 'E_neutral_phi_wall_mJ_m2': E['phi110neu_n8_NEP'],
                        'E_layer_psi_wall_mJ_m2': E['psi001_rs_NEP'],
                        'charged_psi_wall_slope_mJ_m2_per_nm': slope_tags(['psi110chg_n8_NEP', 'psi110chg_n16_NEP', 'psi110chg_n24_NEP'], [8, 16, 24])},
    'one_wall_per_cell_normalization': {'E_neutral_psi_wall_mJ_m2': 2 * E['psi110neu_n8'], 'E_neutral_phi_wall_mJ_m2': 2 * E['phi110neu_n8'],
                                        'E_layer_psi_wall_mJ_m2': 2 * E['psi001_rs'],
                                        'charged_psi_wall_slope_mJ_m2_per_nm': 2 * slope('psi110chg', [8, 12, 16, 20, 24])},
    'profile_from_polarization': {'phi_wall_profile': 'ising'},
}
robustness = {
    'normal_relaxed_cells': {'E_neutral_psi_wall_mJ_m2': E['psi110neu_n8_norm'], 'E_neutral_phi_wall_mJ_m2': E['phi110neu_n8_norm'],
                             'E_layer_psi_wall_mJ_m2': E['psi001_rs_norm'],
                             'charged_psi_wall_slope_mJ_m2_per_nm': slope_tags(['psi110chg_n12_norm', 'psi110chg_n24_norm'], [12, 24])},
    'small_converged_cells_and_long_domain_fit': {'E_neutral_psi_wall_mJ_m2': E['psi110neu_n8'], 'E_neutral_phi_wall_mJ_m2': E['phi110neu_n8'],
                                                   'E_layer_psi_wall_mJ_m2': E['psi001_rs_nz6'],
                                                   'charged_psi_wall_slope_mJ_m2_per_nm': slope('psi110chg', [16, 20, 24])},
    'moderate_cells_fit_8_to_16': {'E_neutral_psi_wall_mJ_m2': E['psi110neu_n12'], 'E_neutral_phi_wall_mJ_m2': E['phi110neu_n12'],
                                   'E_layer_psi_wall_mJ_m2': E['psi001_bl2'],
                                   'charged_psi_wall_slope_mJ_m2_per_nm': slope('psi110chg', [8, 12, 16])},
}

hidden = {'reference': ref, 'tolerance': tol,
          'notes': 'qNEP-w-BECs relaxations with the lab code path (runs/*.json); in-plane lattice of the bulk Cmc2_1 cell'}
priv = TASK/'environment/lab/private'
(priv/'hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
hv = [ref[k] for k in ('E_neutral_psi_wall_mJ_m2', 'E_neutral_phi_wall_mJ_m2', 'E_layer_psi_wall_mJ_m2',
                       'charged_psi_wall_slope_mJ_m2_per_nm')]
(priv/'hidden_values.json').write_text(json.dumps(hv) + '\n')
(TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': ablations, 'robustness': robustness}, indent=1) + '\n')

# agent catalog: bulk structures (qNEP-w-BECs relaxed) and the potentials on offer
def struct(fn):
    s = read_xyz(RAW/'Bulk-structures'/fn)[0]
    cell = s['lattice'].copy()
    cell[np.abs(cell) < 1e-6] = 0.0
    return {'symbols': s['symbols'], 'cell': np.round(cell, 6).tolist(), 'positions': np.round(s['pos'], 6).tolist()}


catalog = {
    'material': 'Ca3Ti2O7 (n = 2 Ruddlesden-Popper titanate), formula unit Ca3Ti2O7 = 12 atoms',
    'bulk_structures': {
        'Cmc2_1_ground_state': {**struct('Cmc2_1-qNEP-w-BECs.xyz'),
                                'note': 'polar ground state relaxed with the qNEP-w-BECs potential (authors\' file); 4 formula units; '
                                        'orthorhombic cell with the stacking axis along z'},
        'I4mmm_reference': {**struct('I4mmm-qNEP-w-BECs.xyz'),
                            'note': 'high-symmetry parent (no rotations, tilts or polarization), conventional body-centred cell, '
                                    'relaxed with qNEP-w-BECs; 2 formula units; tetragonal axes along x and y'},
    },
    'potentials': {
        'qNEP-w-BECs': 'reference model: NEP with machine-learned charges and Ewald electrostatics, trained on energies, forces, '
                       'virials and Born effective charges',
        'qNEP': 'NEP with charges and electrostatics, not trained on Born charges',
        'NEP': 'short-range NEP without electrostatics',
    },
    'source': 'Skogvoll, Fransson, Williamson, Bristowe, Selbach, Erhart (2026), Zenodo 10.5281/zenodo.22801877, CC BY 4.0',
}
(TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')

# calibration print
print('reference', json.dumps(ref, indent=1))
print('tolerance', json.dumps(tol))
for name, vals in list(ablations.items()) + list(robustness.items()):
    flags = []
    for k, v in vals.items():
        if k in tol:
            flags.append('%s=%.2f %s' % (k.split('_')[1] + k.split('_')[2], v, 'in' if tol[k][0] <= v <= tol[k][1] else 'OUT'))
        else:
            flags.append('%s=%s %s' % (k, v, 'in' if v == ref[k] else 'OUT'))
    print('%-45s %s' % (name, '; '.join(flags)))
