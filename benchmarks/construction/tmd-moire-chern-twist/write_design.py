"""Write benchmarks/tasks/tmd-moire-chern-twist/design.json from calibration.json (make_task.py) and, when present,
oracle_run.json (record_oracle.py).

    python3 -B benchmarks/construction/tmd-moire-chern-twist/write_design.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/tmd-moire-chern-twist'
cal = json.loads((HERE/'calibration.json').read_text())
T, A, R = cal['targets'], cal['ablations'], cal['robustness']
oracle = json.loads((HERE/'oracle_run.json').read_text()) if (HERE/'oracle_run.json').exists() else {}

DPS = [
    ('dp1', 'Which member of each near-degenerate DFT doublet carries the K-valley spin-up topology',
     'Use the spin-down (K\') entries, or both', 'All DFT Chern signs flipped ((+1,+1) read as (-1,-1))', 'spin_down_sector'),
    ('dp2', 'Spinful symmetry indicator: exp(2 pi i C/3) = -xi_gamma xi_kappa xi_kappa\' (C3z^3 = -1)',
     'Spinless formula without the minus sign', 'Chern numbers shifted mod 3 (e.g. (+1,+1) -> (-1,0))', 'indicator_without_minus_sign'),
    ('dp3', 'The kappa\' eigenvalue: buy it (equal to kappa, by C2y T) or assume it',
     'Assume xi_kappa\' = conj(xi_kappa) by time reversal', 'C fixed by gamma alone; 7 of 8 pairs wrong', 'kappa_prime_as_conjugate'),
    ('dp4', 'Pairing of the K-up / K\'-down doublets into moire bands',
     'Read band n as the n-th state of the energy-ordered list', 'Band 2 becomes the K\' partner of band 1; wrong pairs', 'state_index_as_band_index'),
    ('dp5', 'Real-space basis and layer assignment of the archived continuum maps',
     'Use the 120-deg plotting vectors of plot-sky.py, a left-handed basis, or Delta-top as the kappa_+ block',
     'Local-stacking set gets (-1,0) or (-1,+1); swapped layers flip both sets', 'continuum_plot_script_basis'),
    ('dp6', 'Which continuum parameter set describes the 3.89 deg DFT bands (literature prior vs data)',
     'Adopt the familiar local-stacking parameters', 'Wrong set (RMS 7.6 vs 1.6 meV); top-band Chern -1 instead of +1',
     'literature_local_stacking_set'),
    ('dp7', 'Transfer of a fixed continuum model to small angles vs buying the DFT symmetry data there',
     'Use the 3.89-deg-fitted model for 2.14 deg (or assume the model tracks DFT)', 'tMoTe2 2.14 deg reported (+1,-1) instead of (+1,+1)',
     'fixed_continuum_transferred_to_dft'),
    ('dp8', 'Location of XM vs MX in the Hartree maps (from the relaxed structures)',
     'Take (1/3,1/3) as MX', 'Every XM value has the wrong sign', 'xm_mx_swapped'),
    ('dp9', 'First-star expansion vs the full potential at small angles',
     'Reconstruct the layer potential from the first star only (or report the full drop as first-star)',
     'XM drop at 1.70 deg 0.34 meV instead of 8.9 meV', 'first_star_reconstruction_as_xm'),
    ('dp10', 'Signed Fourier data vs the archived magnitude-only Fourier file (budget: map vs cheap file)',
     'Use |Delta v_H(g)| with the sign of the large-angle map', 'First-star part at 1.70 / 1.25 deg with the wrong sign',
     'fourier_magnitudes_with_large_angle_sign'),
]

design = {
    'task_id': 'tmd-moire-chern-twist', 'catalog_id': 'tmd-moire-chern-twist', 'start': 'released',
    'source_group': 'zhang_2024_tmd_polarization_topology',
    'title': 'Twisted MoTe2 and WSe2: twist-angle evolution of moire Chern bands',
    'decision_points': [{'id': i, 'fork': f, 'wrong_choice': w, 'consequence': c, 'ablation': a} for i, f, w, c, a in DPS],
    'graded_fields': ['dft_chern_K_up: (C1, C2) of 8 DFT calculations (exact)',
                      'continuum_chern_3p89: (C1, C2) of sets P and Q (exact)',
                      'continuum_set_matching_dft (exact)', 'continuum_P_chern_2p14 (exact)',
                      'hartree_xm_meV at 3.15/1.70/1.54/1.25 deg (+-1.0 meV)',
                      'hartree_first_star_xm_meV at 1.70/1.25 deg (+- max(0.25 meV, 5%))'],
    'length_estimate': {'oracle_lab_evaluations': oracle.get('evaluations'), 'oracle_charge': oracle.get('charge_used'),
                        'budget': 27, 'catalog_total_charge': 56, 'tool_calls': '45-75',
                        'oracle_job_seconds': oracle.get('total_job_seconds'), 'gpt6astra_cost_usd': '3-4.5'},
    'runtime_dependencies': ['numpy (runtime image)'],
    'data_provenance': {'archive_url': 'https://zenodo.org/records/10888769', 'files_used': ['source-data.tar.gz'],
                        'md5': {'source-data.tar.gz': '038b59d589a8f254d76599e6b3e96f74'}, 'license': 'CC BY 4.0',
                        'paper_figure_transcription': 'C3z labels of Fig. 2c/2d (arXiv:2311.12776v3), see SOURCE.md'},
    'exceptions': {},
    'notes': ('Chern numbers of the DFT calculations exist only as figure labels; they were transcribed (C3z labels and '
              'printed C_K-up agree in every panel). Continuum targets are validated-evaluator values of the archived maps. '
              'Graded tWSe2 1.25 deg pair (-1,-1) appears only in the figure, not the text, which states (0,0) for the '
              'change between 1.70 and 1.54 deg.'),
    'calibration': {
        'dft_chern_K_up': {'reference': T['dft_chern_K_up'], 'nearest_ablation': {k: A[k] for k in (
            'dft_spin_down_sector', 'dft_no_minus_sign', 'dft_kappa_prime_conjugate', 'dft_state_index_as_band')},
                           'tolerance': 'exact'},
        'continuum_chern_3p89': {'reference': T['continuum_chern_3p89'],
                                 'robustness': 'identical for 61/91/169 plane waves per layer and 12^2/18^2/24^2 FHS meshes',
                                 'nearest_ablation': {k: A[k] for k in ('continuum_plot_script_basis',
                                                                        'continuum_left_handed_basis', 'continuum_layers_swapped')},
                                 'tolerance': 'exact'},
        'continuum_set_matching_dft': {'reference': T['continuum_set_matching_dft'],
                                       'rms_vs_dft_3p89_meV': cal['continuum_rms_vs_dft_3p89_meV'], 'tolerance': 'exact'},
        'continuum_P_chern_2p14': {'reference': T['continuum_P_chern_2p14'],
                                   'robustness': '91/169 plane waves, 18^2/24^2/30^2 meshes; min direct gap 1-2 = 0.86 meV',
                                   'nearest_ablation': {'equals_dft': A['continuum_2p14_equals_dft']}, 'tolerance': 'exact'},
        'hartree_xm_meV': {'reference': T['hartree_xm_meV'],
                           'robustness': {k: R[k] for k in ('hartree_xm_nearest_grid', 'hartree_xm_bilinear', 'hartree_xm_archived_linecut')},
                           'nearest_ablation': {k: A[k] for k in ('hartree_xm_mx_swapped', 'hartree_xm_first_star_only',
                                                                  'hartree_xm_no_sign_flip')},
                           'tolerance': '+-1.0 meV (estimator spread <= 0.035 meV; nearest wrong path 6.7 meV away at 1.54 deg)'},
        'hartree_first_star_xm_meV': {'reference': T['hartree_first_star_xm_meV'],
                                      'robustness': {'full_grid_fft': R['first_star_full_grid_fft']},
                                      'nearest_ablation': {k: A[k] for k in ('first_star_xm_mx_swapped', 'first_star_equals_full_xm',
                                                                             'first_star_from_magnitudes_3p15_sign')},
                                      'tolerance': '+- max(0.25 meV, 5%) (spread 0.055 meV; nearest wrong path 0.64 meV away at 1.70 deg)'}},
}
if oracle:
    design['length_estimate']['oracle_requests'] = {}
    for r in oracle['requests']:
        design['length_estimate']['oracle_requests'][r['kind']] = design['length_estimate']['oracle_requests'].get(r['kind'], 0) + r['charge']
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print('design.json written; oracle', bool(oracle))
