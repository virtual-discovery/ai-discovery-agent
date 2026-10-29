"""Write benchmarks/tasks/cu2wse4-orbital-edelstein/design.json from calibration.json (make_task.py) and
oracle_run.json (record_oracle.py, optional).

    python3 -B benchmarks/construction/cu2wse4-orbital-edelstein/write_design.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/cu2wse4-orbital-edelstein'
CAL = json.loads((HERE/'calibration.json').read_text())
ORACLE = json.loads((HERE/'oracle_run.json').read_text()) if (HERE/'oracle_run.json').exists() else {}
A, R, TG, PM = CAL['ablations'], CAL['robustness'], CAL['targets'], CAL['per_mesh']
f = lambda x, n=3: float(f'%.{n}g' % x)

DPS = [
    {'id': 'dp1', 'fork': 'Point group and axes: the staged CIF is written in P1; the D2d operations must be found from the '
                          'coordinates, with the twofold axes along a and b (I-42m setting), not on the diagonals',
     'wrong_choice': 'Assume the -4m2 setting (mirrors along a, b) or a chiral class',
     'consequence': f"allowed Edelstein components {A['rotated_axes_setting']['edelstein_allowed_components']} and nonlinear "
                    f"components {A['rotated_axes_setting']['nonlinear_hall_allowed_components']} instead of {TG['edelstein_allowed_components']} "
                    f"and {TG['nonlinear_hall_allowed_components']}", 'ablation': 'rotated_axes_setting'},
    {'id': 'dp2', 'fork': 'Orbital-moment operator: the wave-packet moment with the Wannier-centre position operator is the '
                          'Berry-phase moment; the service also offers the intra-atomic -mu_B<L>',
     'wrong_choice': 'Use the intra-atomic (atomic L) moment, or add it to the Berry-phase moment (both dev pilots did the latter: CB ratio -2.75)',
     'consequence': f"orbital/spin {f(A['intra_atomic_orbital_moment']['orbital_to_spin_ratio']['VBM-0.10'])} (VB) and "
                    f"{f(A['intra_atomic_orbital_moment']['orbital_to_spin_ratio']['CBM+0.10'])} (CB) instead of "
                    f"{f(TG['orbital_to_spin_ratio']['VBM-0.10'])} and {f(TG['orbital_to_spin_ratio']['CBM+0.10'])}: both conclusions reversed",
     'ablation': 'intra_atomic_orbital_moment'},
    {'id': 'dp3', 'fork': 'The no-SOC model has 54 spinless Wannier functions (header: collinear, spinsize 1): every state counts '
                          'twice for the crystal in alpha',
     'wrong_choice': 'Use the model sum as the crystal response',
     'consequence': f"SOC/no-SOC ratio {f(A['no_spin_factor_in_nosoc_alpha']['soc_to_nosoc_orbital_ratio'])} instead of {f(TG['soc_to_nosoc_orbital_ratio'])}",
     'ablation': 'no_spin_factor_in_nosoc_alpha'},
    {'id': 'dp4', 'fork': 'Rigid-band hole counting in the spinless model (cascade: mu(p) then alpha); the two errors partly cancel',
     'wrong_choice': 'Count one hole per band and k-point of the spinless model',
     'consequence': f"ratio {f(A['no_spin_factor_in_nosoc_holes']['soc_to_nosoc_orbital_ratio'])} (or "
                    f"{f(A['no_spin_factor_anywhere_nosoc']['soc_to_nosoc_orbital_ratio'])} with both spin factors omitted)",
     'ablation': 'no_spin_factor_in_nosoc_holes'},
    {'id': 'dp5', 'fork': 'Energy reference: the header Fermi levels of the two models lie mid-gap at different distances from '
                          'the band edges (VBM -0.599 vs -0.644 eV); the paper compares the models on that common axis',
     'wrong_choice': 'Compare the models at the same mu on the header scale instead of at equal hole density',
     'consequence': f"ratio {f(A['common_header_energy_axis']['soc_to_nosoc_orbital_ratio'])} ('agree well') instead of {f(TG['soc_to_nosoc_orbital_ratio'])}",
     'ablation': 'common_header_energy_axis'},
    {'id': 'dp6', 'fork': 'Mesh convergence: the orbital terms converge by N = 24-32, but the small valence-edge spin term keeps '
                          'drifting to the densest meshes; a mesh that looks converged for the orbital part is not converged for the ratio',
     'wrong_choice': 'Stop at the 32 mesh (1 unit) because M_orb is stable',
     'consequence': f"VB orbital/spin {f(A['coarse_mesh_32']['orbital_to_spin_ratio']['VBM-0.10'])} instead of {f(TG['orbital_to_spin_ratio']['VBM-0.10'])}",
     'ablation': 'coarse_mesh_32'},
    {'id': 'dp7', 'fork': 'Units and relaxation time: gamma = hbar/(2 tau), q = -|e|, sums per cell of the model divided by the '
                          'conventional-cell volume',
     'wrong_choice': 'tau = hbar/gamma, or a positive carrier charge',
     'consequence': f"M_x {f(A['tau_hbar_over_gamma']['orbital_magnetization_A_per_m'])} or {f(A['positive_charge_sign']['orbital_magnetization_A_per_m'])} "
                    f"instead of {f(TG['orbital_magnetization_A_per_m'])} A/m", 'ablation': 'tau_hbar_over_gamma'},
    {'id': 'dp8', 'fork': 'Valence versus conduction band: W-d spin-orbit coupling makes spin dominate at the conduction-band edge',
     'wrong_choice': 'Transfer the valence-band conclusion (orbital dominance) to the conduction band',
     'consequence': f"CB ratio {f(TG['orbital_to_spin_ratio']['VBM-0.10'])} instead of {f(TG['orbital_to_spin_ratio']['CBM+0.10'])}",
     'ablation': 'conduction_ratio_from_valence'},
    {'id': 'dp9', 'fork': 'Nonlinear Hall relation: for the Berry-curvature dipole D_xx = -D_yy, D_zz = 0 give sigma_zxy = -2 sigma_xyz',
     'wrong_choice': 'Assume full index-permutation symmetry of sigma_ijl (sigma_zxy = sigma_xyz)',
     'consequence': 'ratio +1 instead of -2', 'ablation': 'nlh_full_permutation_symmetry'},
    {'id': 'dp10', 'fork': 'Chirality indicator: tr G vanishes for the achiral, strongly gyrotropic class D2d',
     'wrong_choice': 'Treat the gyrotropic crystal as chiral (alpha_yy = alpha_xx)',
     'consequence': 'chirality indicator 2 and alpha_yy/alpha_xx = +1 instead of 0 and -1', 'ablation': 'chirality_from_chiral_class'},
]

design = {
    'task_id': 'cu2wse4-orbital-edelstein', 'catalog_id': 'cu2wse4-orbital-edelstein', 'start': 'released',
    'source_group': 'nakazawa_2026_cu2wse4_gyrotropic',
    'title': 'Cu2WSe4: orbital versus spin Edelstein response of a gyrotropic achiral crystal',
    'decision_points': DPS,
    'graded_fields': [
        'orbital_to_spin_ratio at VBM-0.10 and CBM+0.10 eV (SOC model), each within 25%',
        'orbital_magnetization_A_per_m (M_x, VBM-0.10 eV, E_x = 1e4 V/m, gamma = 30 meV) within 15%',
        'soc_to_nosoc_orbital_ratio at p = 5e20 cm^-3 within 10%',
        'edelstein_allowed_components == [xx, yy]; alpha_yy/alpha_xx = -1 and chirality indicator = 0 (+-0.05)',
        'nonlinear_hall_allowed_components == [xyz, xzy, yxz, yzx, zxy, zyx]; sigma_zxy/sigma_xyz = -2 (+-0.15)'],
    'length_estimate': {
        'oracle_lab_evaluations': ORACLE.get('evaluations'), 'oracle_charge': ORACLE.get('charge_used'), 'budget': 30,
        'catalog_units': 'fermi_surface meshes 16-128 cost 24 units (soc) + 14 units (nosoc) = 38, plus any number of bands / '
                         'band_properties requests (1 unit each)',
        'tool_calls': '45-75', 'gpt6astra_cost_usd': '2.5-4.5',
        'why': 'Symmetry analysis from a P1 CIF, band-edge location in two models, choice of orbital operator, a mesh study of a '
               'slowly converging spin term under a binding budget, rigid-band counting in a spinless model and a cascade into '
               'the SOC comparison, unit conversion of the Boltzmann response, and the BCD relation.',
        'oracle_errors': ORACLE.get('errors')},
    'runtime_dependencies': ['numpy and scipy (runtime image); no extra layer'],
    'data_provenance': {'archive_url': 'https://zenodo.org/records/20131063', 'doi': '10.5281/zenodo.20131063',
                        'paper': 'Phys. Rev. Materials 10, 064603 (2026); arXiv:2512.16387', 'files_used': ['CWS_PRM.zip'],
                        'md5': {'CWS_PRM.zip': '9a50ca688ae85f597776e7d6e5931765'},
                        'members': 'benchmarks/construction/cu2wse4-orbital-edelstein/members.md5', 'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': 'Graded values come from the validated evaluator, not directly from the archive: the archived orbital Edelstein '
             'curves equal -2 x the paper-convention (wave-packet) orbital response, while the archived spin curve matches the '
             'evaluator in absolute units (scale 1.002). Nonlinear Hall magnitudes are not graded (BCD at kT = 10 meV is '
             'dominated by near-degeneracy hot spots and the archive differs by a convention factor ~2); only exact symmetry '
             'relations are. The valence-edge orbital/spin ratio carries a mesh uncertainty (128 mesh -11.2, 1/N extrapolation '
             '-13.4), covered by the 25% window.',
    'calibration': {
        'orbital_to_spin_ratio': {'reference_mesh128': TG['orbital_to_spin_ratio'],
                                  'robustness': {'mesh_96': R['mesh_96']['orbital_to_spin_ratio'],
                                                 'extrapolated_VB': R['valence_spin_extrapolated_1_over_N']['orbital_to_spin_ratio']['VBM-0.10']},
                                  'per_mesh': {n: v['ratio'] for n, v in PM.items()},
                                  'nearest_ablation': {'coarse_mesh_32': A['coarse_mesh_32']['orbital_to_spin_ratio'],
                                                       'intra_atomic': A['intra_atomic_orbital_moment']['orbital_to_spin_ratio'],
                                                       'berry_plus_intra_atomic': A['berry_plus_intra_atomic_moment']['orbital_to_spin_ratio'],
                                                       'cb_from_vb': A['conduction_ratio_from_valence']['orbital_to_spin_ratio']['CBM+0.10']},
                                  'tolerance': '+-25% each'},
        'orbital_magnetization_A_per_m': {'reference_mesh128': TG['orbital_magnetization_A_per_m'],
                                          'robustness': {'mesh_96': R['mesh_96']['orbital_magnetization_A_per_m'],
                                                         'per_mesh': {n: v['mag'] for n, v in PM.items()}},
                                          'nearest_ablation': {'tau_hbar_over_gamma': A['tau_hbar_over_gamma']['orbital_magnetization_A_per_m'],
                                                               'positive_charge': A['positive_charge_sign']['orbital_magnetization_A_per_m'],
                                                               'intra_atomic': A['intra_atomic_orbital_moment']['orbital_magnetization_A_per_m']},
                                          'tolerance': '+-15%'},
        'soc_to_nosoc_orbital_ratio': {'reference_mesh128': TG['soc_to_nosoc_orbital_ratio'],
                                       'robustness': {'mesh_96': R['mesh_96']['soc_to_nosoc_orbital_ratio'],
                                                      'per_mesh': {n: v['soc_nosoc'] for n, v in PM.items()}},
                                       'nearest_ablation': {k: A[k]['soc_to_nosoc_orbital_ratio'] for k in
                                                            ('no_spin_factor_anywhere_nosoc', 'common_header_energy_axis',
                                                             'intra_atomic_orbital_moment', 'no_spin_factor_in_nosoc_holes',
                                                             'no_spin_factor_in_nosoc_alpha')},
                                       'tolerance': '+-10%'},
        'symmetry_fields': {'alpha_yy_over_alpha_xx': -1.0, 'chirality_indicator': 0.0, 'sigma_zxy_over_sigma_xyz': -2.0,
                            'raw_bcd_ratio_mesh128': R['raw_bcd_ratio_mesh128']['sigma_zxy_over_sigma_xyz'],
                            'tolerance': '+-0.05 (alpha), +-0.15 (sigma ratio); allowed components exact'}},
}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print('decision points', len(DPS))
