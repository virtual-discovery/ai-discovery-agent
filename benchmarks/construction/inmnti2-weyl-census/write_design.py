"""Write benchmarks/tasks/inmnti2-weyl-census/design.json from calibration.json (make_task.py) and
oracle_run.json (record_oracle.py, optional).

    python3 -B benchmarks/construction/inmnti2-weyl-census/write_design.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/inmnti2-weyl-census'
CAL = json.loads((HERE/'calibration.json').read_text())
ORACLE = json.loads((HERE/'oracle_run.json').read_text()) if (HERE/'oracle_run.json').exists() else {}
ABL = json.loads((TASK/'solution/ablation_values.json').read_text())
A = ABL['ablations']
r = lambda x, n=2: round(x, n)

DPS = [
    {'id': 'dp1', 'fork': 'Which archived calculation describes the magnetic state: `standard` has no starting magnetization and no '
                          'magnetization in its output (non-magnetic PBE+SOC, time-reversal symmetric); `hubbard` is the DFT+U magnetic state',
     'wrong_choice': 'Use the standard (PBE) model as the magnetic Weyl semimetal',
     'consequence': f"node energy {r(A['standard_model_as_magnetic']['node_energy_meV'])} meV (above E_F), separation "
                    f"{r(A['standard_model_as_magnetic']['pair_separation_invA'], 3)} 1/A, kz offset 0, a single orbit of 24",
     'ablation': 'standard_model_as_magnetic'},
    {'id': 'dp2', 'fork': 'Which Wannier bands carry DFT bands 52/53: hubbard 16/17 (exclude_bands 1-34 plus two bands below '
                          'dis_win_min), standard 18/19 (exclude_bands 1-34, no disentanglement)',
     'wrong_choice': 'Offset 34 from exclude_bands for the hubbard model (bands 18/19), or the same Wannier index in both models',
     'consequence': f"{A['band_index_from_exclude_bands']['node_count']['dftu']} DFT+U 'nodes' at +0.66 to +1.44 eV, or "
                    f"{A['same_band_index_both_models']['node_count']['pbe']} PBE nodes at "
                    f"{r(A['same_band_index_both_models']['pbe_node_energy_meV'], 0)} meV",
     'ablation': 'band_index_from_exclude_bands'},
    {'id': 'dp3', 'fork': 'Whether the DFT+U Wannier model can be trusted for node energies and momenta: it reproduces DFT to 0.4 meV '
                          'on its 9x9x9 fitting grid (the scf k-points) but breaks the magnetic symmetry between grid points; '
                          'the archived DFT lines through one node pair are the reliable evidence',
     'wrong_choice': 'Take energy and separation from the Wannier nodes (validated on the scf/nscf grid)',
     'consequence': f"energy {r(A['wannier_nodes_trusted']['node_energy_meV'])} meV (family mean) or "
                    f"{r(A['grid_validated_nearest_pair']['node_energy_meV'])} meV (pair at the DFT line), separation "
                    f"{r(A['wannier_nodes_trusted']['pair_separation_invA'], 3)} or {r(A['grid_validated_nearest_pair']['pair_separation_invA'], 3)} 1/A",
     'ablation': 'wannier_nodes_trusted'},
    {'id': 'dp4', 'fork': 'How the magnetic symmetry (8 operations, 4 combined with time reversal; k -> -Rk for those) organises the '
                          'nodes, including the C2T stabiliser of nodes on the kx = 0 / ky = 0 planes and the Wannier scatter of positions',
     'wrong_choice': 'Cubic non-magnetic group (one orbit of 24, as the paper states), unitary operations only, planar families, '
                     'or in-plane versus out-of-plane split',
     'consequence': 'orbit sizes [24], [4,4,4,4,4,4], [8,8,8] or [8,16] instead of [4,4,4,4,8]',
     'ablation': 'cubic_single_orbit'},
    {'id': 'dp5', 'fork': 'The kz displacement of the near-kz=0 family and its tie to chirality (chirality of bands 1-52 from Berry flux, '
                          'position from the DFT z-line)',
     'wrong_choice': 'Chirality of the upper band (sign flipped), or a symmetric node in the kz = 0 plane',
     'consequence': f"kz offset {r(A['chirality_of_upper_band']['kz_offset_plus_invA'], 4)} or 0 instead of "
                    f"{r(CAL['kz_offset_plus_invA']['reference'], 4)} 1/A",
     'ablation': 'chirality_of_upper_band'},
    {'id': 'dp6', 'fork': 'What the Kubo folder names 10-300 are: k-mesh sizes (berry_kmesh in kubo/InMnTi2.win, fixed 2 meV smearing), '
                          'not broadenings; convergence across meshes is non-monotonic',
     'wrong_choice': "Read '10' as the narrowest broadening, or accept a ~50^3 mesh",
     'consequence': f"sigma_xy {r(A['kubo_folder_as_broadening']['ahc_S_per_cm']['xy'])} or {r(A['kubo_mesh_58']['ahc_S_per_cm']['xy'])} "
                    f"instead of {r(CAL['ahc_S_per_cm']['reference_xy_300'])} S/cm",
     'ablation': 'kubo_folder_as_broadening'},
    {'id': 'dp7', 'fork': 'Whether the archived nonzero sigma_yz, sigma_zx are physical: they are forbidden for moments along z (C2z, C2x T) '
                          'and come from the symmetry-broken Wannier model',
     'wrong_choice': 'Report the archived converged values of all three components',
     'consequence': f"sigma_zx {r(A['raw_hall_components']['ahc_S_per_cm']['zx'])} S/cm instead of 0",
     'ablation': 'raw_hall_components'},
    {'id': 'dp8', 'fork': 'Fermi-energy reference: the standard/nscfXYZ output prints a Fermi energy (15.7968 eV) computed from its '
                          'line k-points; only the scf value (15.6841 eV) is meaningful',
     'wrong_choice': 'Reference the PBE node energies to the latest printed Fermi energy',
     'consequence': f"PBE node energy {r(A['nscf_path_fermi_level']['pbe_node_energy_meV'])} instead of "
                    f"{r(CAL['pbe_node_energy_meV']['reference'])} meV",
     'ablation': 'nscf_path_fermi_level'},
    {'id': 'dp9', 'fork': 'True nodes versus gap minima: the PBE pairs are 0.155 1/A apart (merged on grids up to 24^3) and six '
                          'minima keep a ~6.6 meV gap (avoided crossings)',
     'wrong_choice': 'Count gap minima of a coarse grid, or count every minimised point with a gap of a few meV',
     'consequence': f"PBE count {A['pbe_count_coarse_grid']['node_count']['pbe']} or "
                    f"{A['pbe_count_with_avoided_crossings']['node_count']['pbe']} instead of 24",
     'ablation': 'pbe_count_with_avoided_crossings'},
]

design = {
    'task_id': 'inmnti2-weyl-census', 'catalog_id': 'inmnti2-weyl-census', 'start': 'released',
    'source_group': 'grassano_2024_inmnti2_weyl', 'title': 'InMnTi2: Weyl-node census of a magnetic Heusler',
    'decision_points': DPS,
    'graded_fields': ['node_count.dftu == 24 and node_count.pbe == 24',
                      'orbit_sizes == [4, 4, 4, 4, 8] (any order)',
                      f"node_energy_meV within +-3 of {r(CAL['node_energy_meV']['reference'])}",
                      f"pair_separation_invA within 2% of {r(CAL['pair_separation_invA']['reference'], 4)}",
                      f"kz_offset_plus_invA within +-0.002 of {r(CAL['kz_offset_plus_invA']['reference'], 5)}",
                      f"pbe_node_energy_meV within +-4 of {r(CAL['pbe_node_energy_meV']['reference'])}",
                      f"ahc_S_per_cm.xy within 3.5% of {r(CAL['ahc_S_per_cm']['reference_xy_300'])}; |yz|, |zx| <= 1.5"],
    'length_estimate': {'oracle_lab_evaluations': ORACLE.get('evaluations'), 'oracle_charge': ORACLE.get('charge_used'),
                        'budget': 28,
                        'catalog_units': 'archived items worth 32 units (4 DFT runs x 3, 6 Kubo folders x 2, 4 WannierTools runs x 2) '
                                         'plus any number of Wannier evaluations (1-2 units each)',
                        'tool_calls': '50-80', 'gpt6astra_cost_usd': '3-5',
                        'why': 'Two Wannier censuses (scan, minimisation, chirality), band-index mapping against DFT, a magnetic '
                               'symmetry analysis with antiunitary operations, discovery that the DFT+U interpolation fails off its '
                               'fitting grid and extraction of the node from the archived DFT lines, and a Kubo convergence and '
                               'symmetry check; each wrong fork moves a graded value outside its tolerance.',
                        'oracle_errors': ORACLE.get('errors')},
    'runtime_dependencies': ['numpy and scipy (runtime image); no extra layer'],
    'data_provenance': {'archive_url': 'https://archive.materialscloud.org/records/hkprp-xw234',
                        'doi': '10.24435/materialscloud:ph-3c', 'paper_doi': '10.1103/PhysRevResearch.6.013140',
                        'files_used': ['upload.zip', 'files_description.md'],
                        'md5': {'upload.zip': '27d1aa76a5cc384765108beb2ae05289',
                                'files_description.md': '37d5be45010d2ae698b72c26c1dcd574'},
                        'members': 'benchmarks/construction/inmnti2-weyl-census/members.md5', 'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': 'Adapted from the candidate design (see the construction README). The archived DFT+U Wannier model is built on a '
             '9x9x9 grid and breaks the magnetic symmetry by up to ~0.1 eV between grid points; graded node energy, separation '
             'and kz offset therefore come from the archived DFT lines, and counts, orbits and chirality from the Wannier model. '
             'The optical onset is not graded (the archived Kubo onset is smeared by the Gamma pocket and the Wannier node '
             'scatter, so it is not uniquely recoverable); it is a rubric item. Memorised paper values: E_W = -96 meV and '
             '0.37 1/A pass; the paper kz shift 0.0047 and "all nodes equivalent" fail.',
    'calibration': {
        'node_energy_meV': {'reference': r(CAL['node_energy_meV']['reference']),
                            'definition': CAL['node_energy_meV']['definition'],
                            'robustness': {k: r(v) for k, v in CAL['node_energy_meV']['robustness'].items()},
                            'nearest_ablation': {k: r(v) for k, v in CAL['node_energy_meV']['nearest_ablation'].items()},
                            'wannier_family_energies_meV': [r(x, 1) for x in CAL['node_energy_meV']['wannier_family_energies_meV']],
                            'tolerance': '+-3 meV'},
        'pair_separation_invA': {'reference': r(CAL['pair_separation_invA']['reference'], 5),
                                 'robustness': {k: r(v, 5) for k, v in CAL['pair_separation_invA']['robustness'].items()},
                                 'nearest_ablation': {k: r(v, 4) for k, v in CAL['pair_separation_invA']['nearest_ablation'].items()},
                                 'tolerance': '+-2%'},
        'kz_offset_plus_invA': {'reference': r(CAL['kz_offset_plus_invA']['reference'], 5),
                                'robustness': {k: r(v, 5) for k, v in CAL['kz_offset_plus_invA']['robustness'].items()},
                                'nearest_ablation': CAL['kz_offset_plus_invA']['nearest_ablation'],
                                'wannier_plus_family_kz': [r(x, 4) for x in CAL['kz_offset_plus_invA']['wannier_plus_family_kz']],
                                'not_ablation': CAL['kz_offset_plus_invA']['not_ablation'], 'tolerance': '+-0.002 1/A'},
        'pbe_node_energy_meV': {'reference': r(CAL['pbe_node_energy_meV']['reference']),
                                'spread': [r(x) for x in CAL['pbe_node_energy_meV']['spread']],
                                'nearest_ablation': {k: r(v) for k, v in CAL['pbe_node_energy_meV']['nearest_ablation'].items()},
                                'tolerance': '+-4 meV'},
        'ahc_S_per_cm': {'reference_xy': r(CAL['ahc_S_per_cm']['reference_xy_300'], 3),
                         'archived_by_folder': {k: {c: r(x, 2) for c, x in v.items()} for k, v in CAL['ahc_S_per_cm']['by_folder'].items()},
                         'robustness': {'mesh_210': r(CAL['ahc_S_per_cm']['by_folder']['210']['xy'], 2),
                                        'mesh_120_passes_by_2.8%': r(CAL['ahc_S_per_cm']['by_folder']['120']['xy'], 2)},
                         'nearest_ablation': {'mesh_30_fails_by_4.6%': r(CAL['ahc_S_per_cm']['by_folder']['30']['xy'], 2),
                                              'mesh_58': r(CAL['ahc_S_per_cm']['by_folder']['58']['xy'], 2),
                                              'folder_10': r(CAL['ahc_S_per_cm']['by_folder']['10']['xy'], 2),
                                              'raw_zx_300': r(CAL['ahc_S_per_cm']['by_folder']['300']['zx'], 2)},
                         'tolerance': CAL['ahc_S_per_cm']['tolerance']},
        'node_count': CAL['node_count'], 'orbit_sizes': CAL['orbit_sizes']}}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print('design.json written with', len(DPS), 'decision points')
