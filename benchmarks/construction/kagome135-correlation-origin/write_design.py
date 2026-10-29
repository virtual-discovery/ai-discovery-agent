"""Write the task's design.json from calibration.json (and oracle_run.json when present)."""
import json
from pathlib import Path
HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/kagome135-correlation-origin'
C = json.loads((HERE/'calibration.json').read_text())
ref, ab, rob, tol = C['reference'], C['ablations'], C['robustness'], C['tolerances']
oracle = json.loads((HERE/'oracle_run.json').read_text()) if (HERE/'oracle_run.json').exists() else {}
dps = [
    {'id': 'dp1', 'fork': 'Pristine electron count of each model: from the archive (charge-neutral DFT, Cs 6s donated) or nominal valence',
     'wrong_choice': 'Nominal / paper counts 27, 30, 33 (Cs electron ignored)',
     'consequence': 'Counts wrong; CsCr3Sb5 Zbar %.3f and Zmin %.3f instead of %.3f, %.3f; CsCr3Bi5 Zbar %.3f; Hund ratio %.2f instead of %.2f'
     % (ab['CrSb_nominal']['Z_mean'], ab['CrSb_nominal']['Z_min'], ref['CrSb_Z_mean'], ref['CrSb_Z_min'],
        ab['CrBi_nominal']['Sb_lattice_Z_mean'], ab['hund_nominal'], ref['hund_ratio']),
     'ablation': 'nominal_electron_counts'},
    {'id': 'dp2', 'fork': 'Which archive folder is CsCr3Bi5 at the CsCr3Sb5 lattice (from the .wout lattice vectors)',
     'wrong_choice': 'Take CrBi_9.2 as the CsCr3Sb5-lattice model (name read as c = 9.2 A)',
     'consequence': 'Bi-lattice Zbar values exchanged (%.3f vs %.3f); sigma %.3f instead of %.3f eV'
     % (ref['CrBi_U4p5_Z_mean'], ref['CrBi92_U4p5_Z_mean'], ab['sigma_swapped_CrBi']['CsCr3Bi5_Sb_lattice'], ref['sigma_CrBi']),
     'ablation': 'swapped_CsCr3Bi5_models'},
    {'id': 'dp3', 'fork': 'Order of the d Wannier functions: archive projection list (mr = 5,3,1,2,4) or Wannier90 default',
     'wrong_choice': 'Assume the default order dz2, dxz, dyz, dx2-y2, dxy', 'consequence': 'Most renormalized orbital reported as dx2-y2',
     'ablation': 'w90_default_orbital_order'},
    {'id': 'dp4', 'fork': 'Local frame of the label: the archive projection frame (x radial) or the paper frame (B3g = dyz)',
     'wrong_choice': 'Recall the paper: B3g = dyz', 'consequence': 'Most renormalized orbital reported as dyz instead of dxz',
     'ablation': 'recalled_paper_frame_label'},
    {'id': 'dp5', 'fork': 'Separating filling from Hamiltonian: both counterfactual orders vs a single path or a cross-compound regression',
     'wrong_choice': 'One-at-a-time path or regression of ln Z on filling and bandwidth across compounds',
     'consequence': 'Filling share %.2f / %.2f (one-sided) or %.2f (regression) instead of %.2f'
     % (ab['share_one_sided_Ti'], ab['share_one_sided_Cr'], ab['share_regression'], ref['filling_share']),
     'ablation': 'cross_compound_regression'},
    {'id': 'dp6', 'fork': 'Kinetic-energy scale: full PDOS second moment (includes d-p hybridization) vs windowed PDOS, d-d block or orbital average',
     'wrong_choice': 'Second moment of the PDOS within the plotted +-4 eV window',
     'consequence': 'sigma %.3f instead of %.3f eV (CsCr3Bi5)' % (ab['sigma_windowed']['CsCr3Bi5_Sb_lattice'], ref['sigma_CrBi']),
     'ablation': 'sigma_windowed_pdos'},
    {'id': 'dp7', 'fork': 'Hund versus Mott: the J dependence at fixed U, at the right filling',
     'wrong_choice': 'Read the J dependence at the nominal count (cascade from dp1)',
     'consequence': 'Hund ratio %.2f instead of %.2f' % (ab['hund_nominal'], ref['hund_ratio']),
     'ablation': 'hund_ratio_nominal_count'},
]
calib = {
    'electron_counts': {'reference': ref['electron_counts'], 'nearest_ablation': ab['nominal_counts'], 'tolerance': 'exact integers (+-0.5)'},
    'CsCr3Sb5_Z_mean': {'reference': ref['CrSb_Z_mean'], 'robustness': {'mesh_12x12x8': rob['CrSb_mesh12']['Z_mean']},
                        'nearest_ablation': {'nominal_count': ab['CrSb_nominal']['Z_mean']}, 'tolerance': '+-%g' % tol['Z_abs']},
    'CsCr3Sb5_Z_min': {'reference': ref['CrSb_Z_min'], 'robustness': {'mesh_12x12x8': rob['CrSb_mesh12']['Z_min']},
                       'nearest_ablation': {'nominal_count': ab['CrSb_nominal']['Z_min']}, 'tolerance': '+-%g' % tol['Z_abs']},
    'most_renormalized_orbital': {'reference': ref['most_renormalized_orbital'],
                                  'nearest_ablation': {'w90_default_order': ab['orbital_w90_default_order'], 'paper_frame': ab['orbital_paper_frame']},
                                  'tolerance': 'exact'},
    'filling_share': {'reference': ref['filling_share'], 'robustness': {'mesh_12x12x8': rob['share_mesh12'], 'nominal_counts (not an ablation: insensitive)': ab['share_nominal']},
                      'nearest_ablation': {'regression': ab['share_regression'], 'one_sided_Cr': ab['share_one_sided_Cr'], 'one_sided_Ti': ab['share_one_sided_Ti']},
                      'tolerance': '+-%g' % tol['share_abs']},
    'hund_ratio': {'reference': ref['hund_ratio'], 'robustness': {'mesh_12x12x8': rob['hund_mesh12']},
                   'nearest_ablation': {'nominal_count': ab['hund_nominal']}, 'tolerance': '+-%d%%' % (100 * tol['hund_rel'])},
    'CsCr3Bi5_Sb_lattice_Z_mean_U4.5': {'reference': ref['CrBi_U4p5_Z_mean'], 'robustness': {'mesh_12x12x8': rob['CrBi_mesh12']['Sb_lattice_Z_mean']},
                                        'nearest_ablation': {'nominal_count': ab['CrBi_nominal']['Sb_lattice_Z_mean'], 'swapped': ab['CrBi_swapped']['Sb_lattice_Z_mean']},
                                        'tolerance': '+-%g' % tol['Z_abs']},
    'CsCr3Bi5_TiBi_lattice_Z_mean_U4.5': {'reference': ref['CrBi92_U4p5_Z_mean'], 'robustness': {'mesh_12x12x8': rob['CrBi_mesh12']['TiBi_lattice_Z_mean']},
                                          'nearest_ablation': {'nominal_count': ab['CrBi_nominal']['TiBi_lattice_Z_mean'], 'swapped': ab['CrBi_swapped']['TiBi_lattice_Z_mean']},
                                          'tolerance': '+-%g' % tol['Z_abs']},
    'sigma_CsCr3Bi5_eV': {'reference': ref['sigma_CrBi'], 'robustness': {'coarse_mesh_pdos': rob['sigma_mesh_histogram']['CsCr3Bi5_Sb_lattice']},
                          'nearest_ablation': {'swapped_model': ab['sigma_swapped_CrBi']['CsCr3Bi5_Sb_lattice'], 'windowed': ab['sigma_windowed']['CsCr3Bi5_Sb_lattice'],
                                               'd_only': ab['sigma_d_only']['CsCr3Bi5_Sb_lattice'], 'orbital_mean': ab['sigma_mean_orbital']['CsCr3Bi5_Sb_lattice']},
                          'tolerance': '+-%d%%' % (100 * tol['sigma_rel'])},
    'sigma_CsTi3Sb5_eV': {'reference': ref['sigma_TiSb'], 'robustness': {'coarse_mesh_pdos': rob['sigma_mesh_histogram']['CsTi3Sb5']},
                          'nearest_ablation': {'windowed': ab['sigma_windowed']['CsTi3Sb5'], 'd_only': ab['sigma_d_only']['CsTi3Sb5']},
                          'tolerance': '+-%d%%' % (100 * tol['sigma_rel'])},
}
design = {
    'task_id': 'kagome135-correlation-origin', 'catalog_id': 'kagome135-correlation-origin', 'start': 'released',
    'source_group': 'crispino_2025_kagome135_correlations',
    'title': '135 kagome metals: origin of correlations from Wannier models and slave spins',
    'decision_points': dps,
    'graded_fields': ['electron_counts (exact)', 'CsCr3Sb5_U6.Z_mean, Z_min (+-0.025)', 'most_renormalized_orbital (exact)',
                      'filling_share (+-0.08)', 'hund_ratio (+-10%)', 'CsCr3Bi5_U4p5 Sb- and TiBi-lattice Z_mean (+-0.025)',
                      'sigma_eV for CsCr3Bi5 (Sb lattice) and CsTi3Sb5 (+-3%)'],
    'length_estimate': {'oracle_lab_evaluations': oracle.get('evaluations', 13), 'budget': 20,
                        'oracle_lab_minutes': oracle.get('minutes'), 'tool_calls': '40-70', 'gpt6astra_cost_usd': '2-4',
                        'why': 'Archive forensics (electron counts, lattices, projection codes), local tight-binding analysis, '
                               'counterfactual SSMF design under a 20-unit budget; evaluations take 0.5-9 min each.'},
    'runtime_dependencies': ['numpy, scipy (runtime image); kag_ssmf.py evaluator vendored in the lab'],
    'data_provenance': {'archive_url': 'https://zenodo.org/records/17951035', 'doi': '10.5281/zenodo.17951035',
                        'files_used': ['CsTiSb.zip', 'CsVSb.zip', 'CsCrSb.zip', 'CsTiBi.zip', 'CsCrBi.zip', 'CsCrBi_9.2.zip', 'DFT.zip'],
                        'md5': {'CsTiSb.zip': '3022b11da82d2a3dd91b57d2fe10d75b', 'CsCrSb.zip': 'ba28567cb06c56faaf9188c06ce5a1ca',
                                'CsCrBi.zip': '95d1e17aaeb1c4b9392c096e27afd76f', 'DFT.zip': 'f8b8de95fa310d206ef01bdd27577501',
                                'CsVSb.zip': '35a8b537aaf029def98e427f23785bed', 'CsCrBi_9.2.zip': '9c4a37b89801300f30eb2bd5d86167d3',
                                'CsTiBi.zip': '1d92f3a6309e431f6debbbf1ed61df8f'},
                        'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': ('Released start with a computing lab: targets come from our numpy SSMF evaluator (kag_ssmf.py) run on the authors\' '
              'archived Wannier Hamiltonians, validated against the paper (see construction README). U_c of CsCr3Bi5 is not graded: '
              'the end of the metallic branch is path dependent within ~0.3 eV and costs up to 9 min per probe.'),
    'calibration': calib,
}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print('ok')
