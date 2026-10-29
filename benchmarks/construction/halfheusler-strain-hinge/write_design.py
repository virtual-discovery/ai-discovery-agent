"""Write benchmarks/tasks/halfheusler-strain-hinge/design.json from calibration.json (make_task.py) and
oracle_run.json (record_oracle.py, optional).

    python3 -B benchmarks/construction/halfheusler-strain-hinge/write_design.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/halfheusler-strain-hinge'
CAL = json.loads((HERE/'calibration.json').read_text())
ORACLE = json.loads((HERE/'oracle_run.json').read_text()) if (HERE/'oracle_run.json').exists() else {}
T, TH, A, R = CAL['targets'], CAL['thresholds'], CAL['ablations'], CAL['robustness']
r = lambda x, n=1: round(x, n)
G, EK, W = T['fundamental_gap_meV'], T['hinge_energy_k0_meV'], T['hinge_band_window_meV']
O, S = 'LiSbZn_orig', 'LiSbZn_strained_setup'

DPS = [
    {'id': 'dp1',
     'fork': 'What the archived "strained_setup" models contain. Their bct lattice vectors give c = sqrt(2) a_t = c0 and the '
             'volume of the orig cell, the Gamma8 quartet is intact, and the LiBiZn bands match the orig model to a few meV. '
             'They are the equilibrium cubic structures in a tetragonal setting, and the kwant_ham "*_strain.ham" pickles hold '
             'the same data. The archive contains no strained Hamiltonian.',
     'wrong_choice': 'Trust the file names and the paper (delta = +-2 %)',
     'consequence': 'strain +-2 % instead of 0 for both strained-setup models; LiSbZn calculation differences read as strain effects',
     'ablation': 'strain_from_file_names'},
    {'id': 'dp2',
     'fork': 'Fundamental gap of LiSbZn from a full-zone search. The conduction-band minimum lies next to L (about 1 meV below L), '
             'not at Gamma6, so the gap is indirect.',
     'wrong_choice': 'Gamma6 - Gamma8 at Gamma',
     'consequence': f"{r(A['gap_at_gamma'][O])} and {r(A['gap_at_gamma'][S])} meV instead of {r(G[O])} and {r(G[S])} meV",
     'ablation': 'gap_at_gamma'},
    {'id': 'dp3',
     'fork': 'Whether LiBiZn is a zero-gap semiconductor with E_F at Gamma8, or a compensated semimetal. Band 7 near L lies '
             f"{r(-G['LiBiZn_orig'], 0)} meV below the Gamma8 maximum, and E_F follows from the six-electron count "
             '(Li 1 + Zn 2 + X-p 3) on a converged mesh.',
     'wrong_choice': 'E_F at the Gamma8 touching point (half-filled quartet, as the paper describes it), or a coarse mesh',
     'consequence': f"E_F 0, {r(A['fermi_coarse_mesh_8'], 0)} or {r(A['fermi_mesh_12'], 0)} meV instead of "
                    f"{r(T['libizn_fermi_level_meV'], 0)} meV relative to Gamma8; gap 0 instead of {r(G['LiBiZn_orig'], 0)} meV",
     'ablation': 'fermi_at_gamma8'},
    {'id': 'dp4',
     'fork': 'Energy reference for LiSbZn. The FPLO zero sits at an arbitrary point inside the gap, and at a different point in the '
             'two LiSbZn calculations, so energies must be referred to each model\'s own VBM.',
     'wrong_choice': 'Use the model zero as the reference (taken as E_F)',
     'consequence': f"hinge energies {r(A['hinge_energy_model_zero'][O])} and {r(A['hinge_energy_model_zero'][S])} meV instead of "
                    f"{r(EK[O])} and {r(EK[S])} meV",
     'ablation': 'hinge_energy_model_zero'},
    {'id': 'dp5',
     'fork': 'Which hinge carries the in-gap state. In the service geometry the Zn site sits at +1/4 along both sides, so Zn is '
             'outermost on the t = N-1 facets; the state lies at the Zn/Zn hinge with about two-thirds Zn-4s weight. The paper '
             'puts it at the Bi/Sb hinge with about 95 % Bi/Sb weight. Its real-space picture is mirrored by the authors\' kwant '
             'HoppingKind assignment, which reverses each hopping vector relative to the site positions.',
     'wrong_choice': 'Adopt the paper\'s Bi/Sb-hinge assignment and orbital character',
     'consequence': f"termination X/X and Zn-4s weight 0.05 instead of Zn/Zn and {r(T['hinge_zn_s_weight_k0'], 2)}",
     'ablation': 'paper_termination'},
    {'id': 'dp6',
     'fork': 'Finite-size convergence of the k = 0 hinge level. It sits close to the conduction band and converges slowly: '
             f"+{r(A['hinge_energy_N12'][O] - EK[O])} meV at N = 12 and "
             f"+{r(R['hinge_energy_N20_direct'][O] - EK[O])} meV at N = 20.",
     'wrong_choice': 'Read the level from a small prism (N <= 12)',
     'consequence': f"{r(A['hinge_energy_small_prism_N10'][O])} (N = 10) or {r(A['hinge_energy_N12'][O])} meV (N = 12) instead of "
                    f"{r(EK[O])} meV",
     'ablation': 'hinge_energy_small_prism_N10'},
    {'id': 'dp7',
     'fork': 'Hinge-band window. The minimum is at k = pi, not k = 0. The upper of the two spin-split branches peaks near k = 1.0, '
             'inside the projected gap but far above the global conduction-band minimum.',
     'wrong_choice': 'A coarse k grid, clipping to the global gap, assuming the minimum at k = 0, or keeping one branch only',
     'consequence': f"windows {[r(x) for x in A['window_grid5']]}, {[r(x) for x in A['window_clipped_to_bulk_gap']]}, "
                    f"{[r(x) for x in A['window_bottom_at_k0']]} or {[r(x) for x in A['window_lower_branch_only']]} instead of "
                    f"{[r(x) for x in W]} meV",
     'ablation': 'window_grid5'},
    {'id': 'dp8',
     'fork': 'The same physical prism in the strained-setup file. Its Cartesian frame is rotated by -45 deg about z and its '
             'lattice basis differs, so P_N needs the basis [[1,0,0],[1,1,1],[1,0,1]], which puts the Zn site at reduced '
             '(1/4, 1/4, 1/4).',
     'wrong_choice': 'Reuse the orig integer basis, as the authors\' shape function does in the file\'s own basis',
     'consequence': f"a prism with (100) and (111) facets whose gap is full of facet states; the most corner-weighted level is "
                    f"{r(A['str_hinge_identity_basis'])} meV instead of {r(EK[S])} meV",
     'ablation': 'str_hinge_identity_basis'},
    {'id': 'dp9',
     'fork': 'Trivial or higher-order origin. The decisive test is a local potential at the hinge that carries the state. An '
             f"attractive potential of {r(T['hinge_removal_potential_eV'], 3)} eV on its 2x2 cells pushes the level into the "
             'valence band and pulls new bound states out of the conduction band.',
     'wrong_choice': 'Apply the test at the paper\'s Bi/Sb (X/X) hinge, see no effect and conclude that the state is protected',
     'consequence': f"removal potential null instead of {r(T['hinge_removal_potential_eV'], 3)} eV (and "
                    f"{r(A['removal_to_model_zero'], 3)} eV if the level is driven to the FPLO zero instead of the VBM, dp4)",
     'ablation': 'removal_potential_wrong_corner'},
]

cal = {
    'strain_percent': {'reference': T['strain_percent'], 'nearest_ablation': {'file_names': 2.0, 'compressive': -2.0},
                       'tolerance': f"+-{TH['strain_percent_abs']}"},
    'fundamental_gap_meV': {'reference': {k: r(v, 2) for k, v in G.items()},
                            'robustness': {'grid_L_point': {k: r(v, 2) for k, v in R['gap_grid_L_point'].items()}},
                            'nearest_ablation': {'gap_at_gamma': {k: r(v, 2) for k, v in A['gap_at_gamma'].items()}},
                            'tolerance': {k: f'+-{v}' for k, v in TH['gap_abs_meV'].items()}},
    'libizn_fermi_level_meV': {'reference_96cube': r(T['libizn_fermi_level_meV'], 2),
                               'robustness': {'mesh_48': r(R['fermi_mesh_48'], 2), 'mesh_64_offset': r(R['fermi_mesh_64_lab_max'], 2)},
                               'nearest_ablation': {'mesh_12': r(A['fermi_mesh_12'], 2), 'mesh_8': r(A['fermi_coarse_mesh_8'], 2),
                                                    'at_gamma8': 0.0},
                               'tolerance': f"+-{TH['fermi_level_abs_meV']}"},
    'hinge_energy_k0_meV': {'reference_extrapolated': {k: r(v, 2) for k, v in EK.items()},
                            'robustness': {'N28_direct': {k: r(v, 2) for k, v in R['hinge_energy_N28_direct'].items()},
                                           'N20_direct': {k: r(v, 2) for k, v in R['hinge_energy_N20_direct'].items()}},
                            'nearest_ablation': {'N12': {k: r(v, 2) for k, v in A['hinge_energy_N12'].items()},
                                                 'N10': {k: r(v, 2) for k, v in A['hinge_energy_small_prism_N10'].items()},
                                                 'model_zero': {k: r(v, 2) for k, v in A['hinge_energy_model_zero'].items()},
                                                 'strained_setup_identity_basis': r(A['str_hinge_identity_basis'], 2)},
                            'tolerance': f"+-{TH['hinge_energy_abs_meV']}"},
    'hinge_band_window_meV': {'reference': [r(x, 2) for x in W],
                              'robustness': {k: [r(x, 2) for x in R[k]] for k in ('window_N12_refined', 'window_N12_grid25')},
                              'nearest_ablation': {k: [r(x, 2) for x in A[k]] for k in ('window_grid5', 'window_clipped_to_bulk_gap',
                                                                                         'window_bottom_at_k0', 'window_lower_branch_only')},
                              'tolerance': f"bottom +-{TH['band_bottom_abs_meV']}, top +-{TH['band_top_abs_meV']}"},
    'hinge_termination': {'reference': T['hinge_termination'], 'nearest_ablation': {'paper': 'X/X (Bi/Sb hinge)'}},
    'hinge_zn_s_weight_k0': {'reference_N28': r(T['hinge_zn_s_weight_k0'], 4), 'robustness': {'N16': R['zn_weight_N16']},
                             'nearest_ablation': {'paper_95_percent_on_X': 0.05}, 'tolerance': f"+-{TH['zn_weight_abs']}"},
    'hinge_removal_potential_eV': {'reference_N20_brentq': r(T['hinge_removal_potential_eV'], 4),
                                   'robustness': {'N16_linear': r(R['removal_N16_linear'], 4), 'N16_brentq': r(R['removal_N16_brentq'], 4)},
                                   'nearest_ablation': {'to_model_zero': r(A['removal_to_model_zero'], 4), 'wrong_corner': None},
                                   'tolerance': f"+-{TH['removal_potential_abs_eV']}"},
    'band_inverted': {'reference': T['band_inverted']},
}

design = {
    'task_id': 'halfheusler-strain-hinge', 'catalog_id': 'halfheusler-strain-hinge', 'start': 'released',
    'source_group': 'das_2025_halfheusler_hinge',
    'title': 'Half-Heusler LiSbZn and LiBiZn: hinge states and the archived strain models',
    'decision_points': DPS,
    'graded_fields': [
        f"strain_percent of all four models within +-{TH['strain_percent_abs']} of 0",
        'band_inverted: LiSbZn false, LiBiZn true',
        f"fundamental_gap_meV: LiSbZn_orig {r(G[O])} +-{TH['gap_abs_meV'][O]}, LiSbZn_strained_setup {r(G[S])} "
        f"+-{TH['gap_abs_meV'][S]}, LiBiZn_orig {r(G['LiBiZn_orig'])} +-{TH['gap_abs_meV']['LiBiZn_orig']}",
        f"libizn_fermi_level_meV {r(T['libizn_fermi_level_meV'])} +-{TH['fermi_level_abs_meV']} (relative to Gamma8)",
        'hinge.termination == Zn/Zn',
        f"hinge.energy_k0_meV: LiSbZn_orig {r(EK[O])} +-{TH['hinge_energy_abs_meV']}, LiSbZn_strained_setup {r(EK[S])} "
        f"+-{TH['hinge_energy_abs_meV']} (relative to each VBM)",
        f"hinge.band_window_meV [{r(W[0])} +-{TH['band_bottom_abs_meV']}, {r(W[1])} +-{TH['band_top_abs_meV']}]",
        f"hinge.zn_s_weight_k0 {r(T['hinge_zn_s_weight_k0'], 3)} +-{TH['zn_weight_abs']}",
        f"hinge.removal_potential_eV {r(T['hinge_removal_potential_eV'], 3)} +-{TH['removal_potential_abs_eV']} (null fails)",
    ],
    'length_estimate': {
        'oracle_lab_evaluations': ORACLE.get('evaluations'), 'oracle_charge': ORACLE.get('charge_used'),
        'budget': ORACLE.get('budget'),
        'catalog_units': 'unbounded: every request is a new evaluation (bulk bands, band-edge meshes, Fermi levels, slabs, '
                         'prisms of any size up to 28 with optional potentials)',
        'tool_calls': '50-80', 'gpt6astra_cost_usd': '3-5',
        'why': ('Bulk classification of four models (inversion, full-zone gaps, a semimetal Fermi level), discovery that the '
                '"strained" files are cubic and use another frame, prism size convergence, assignment of the hinge from the cell '
                'geometry against the paper, a k-resolved hinge-band scan, the same prism rebuilt in the rotated file, and a '
                'potential scan; prisms of size 21-28 cost 3 units, so the 30-unit budget forces choices.'),
        'oracle_errors': ORACLE.get('errors'), 'oracle_seconds_per_request': ORACLE.get('seconds')},
    'runtime_dependencies': ['numpy and scipy (runtime image); no extra layer'],
    'data_provenance': {'archive_url': 'https://zenodo.org/records/16927281', 'doi': '10.5281/zenodo.16927281',
                        'paper_doi': '10.1103/PhysRevB.112.125132', 'files_used': ['strained_hH_zenodo.zip'],
                        'md5': {'strained_hH_zenodo.zip': 'f66b568ce17e4055122facfe02e01e24'},
                        'members': 'benchmarks/construction/halfheusler-strain-hinge/fetch.py (MEMBERS)', 'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': ('Adapted from the candidate design (see the construction README). The archive holds no strained Hamiltonian: both '
              '"eq_strained_setup" files are the cubic equilibrium structures in a tetragonal setting. Strain-dependent targets '
              'were therefore replaced by what the data contain: the identification of the archived models, bulk band edges and '
              'the LiBiZn Fermi level, and the LiSbZn hinge state (assignment, converged energy, dispersion window, orbital '
              'character, removal by a local potential) in two independent LiSbZn calculations. Strong Z2 is not graded: it is '
              'trivially 0 for normal-ordered LiSbZn and undefined for cubic LiBiZn, where the Gamma8 touching is unremoved.'),
    'calibration': cal,
}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print(json.dumps({k: design[k] for k in ('graded_fields',)}, indent=1))
