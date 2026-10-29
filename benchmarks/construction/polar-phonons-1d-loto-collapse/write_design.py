"""Write design.json of polar-phonons-1d-loto-collapse from calibration.json (and check.json when present).

    python3 -B write_design.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/polar-phonons-1d-loto-collapse'
C = json.loads((HERE/'calibration.json').read_text())
CHECK = json.loads((HERE/'check.json').read_text()) if (HERE/'check.json').exists() else None
E = C['ablation_errors']
R = C['robustness_errors']
O = C['oracle_errors']
T3 = C['three_d']


def worst(errs):
    return max(((tb, q, v) for tb in errs for q, v in errs[tb].items()), key=lambda x: abs(x[2]))


decision_points = [
    {'id': 'dp1', 'kind': 'competing explanations / boundary conditions',
     'fork': 'Which boundary conditions describe an isolated tube: the 3D-periodic runs (3D_code) or the 1D-cutoff runs (1D_code); and whether the ph.x Gamma frequency is the 3D q->0 LO',
     'wrong_choice': 'Use the 3D-periodic dispersion (finite LO-TO splitting at Gamma, bulk-like) for the isolated tubes; or report the analytic Gamma frequency as the 3D q->0 limit',
     'consequence': 'LO +95/+93 cm^-1 at 0.01|GZ| and +46/+40 at 0.05 (BN55/BN66); 3D q->0 value low by %.0f cm^-1' % (C['truth_3d'] - T3['abl_analytic_gamma_frequency']),
     'ablation': 'three_d_boundary_condition_dispersion'},
    {'id': 'dp2', 'kind': 'model class (dimensional fingerprint)',
     'fork': 'Asymptotic form of the dielectric shift: constant (3D), linear in q (2D) or q^2 ln(1/qt) (1D)',
     'wrong_choice': '2D-like linear collapse omega^2 = omega_G^2 + B q, B from the purchasable explicit 0.15|GZ| run of the hidden tube',
     'consequence': 'LO +3.6/+3.4 at 0.01 and -12/-16 at 0.05, -18/-21 at 0.075 cm^-1',
     'ablation': 'linear_2d_collapse'},
    {'id': 'dp3', 'kind': 'evidence to discover (defective archived reference)',
     'fork': 'Standard Fourier interpolation of the 1D force constants versus subtracting and re-adding the 1D long-range term; the archived BN(4,4) 1D_code interpolated dispersion is itself a plain interpolation (no long-range term), 33 cm^-1 below its own explicit 1D runs at 0.05|GZ|, visible only by cross-checking',
     'wrong_choice': 'Trust matdyn-like plain interpolation (or validate against the archived BN(4,4) interpolated dispersion)',
     'consequence': 'LO -5.6/-6.0 at 0.01, -40/-44 at 0.05, -46/-50 at 0.075 cm^-1',
     'ablation': 'standard_fourier_interpolation'},
    {'id': 'dp4', 'kind': 'model parameter / physical definition',
     'fork': 'Radius and dielectric constant of the cylinder model: the extent of the electronic polarization (about 1.4-1.9 times the atomic-shell radius, or the authors\' d/4) with eps1d = A(eps_zz - 1)/(pi t^2), versus the atomic-shell radius or the vacuum-diluted cell eps_zz used as the material constant',
     'wrong_choice': 'Atomic-shell radius (6.6 and 7.9 bohr) with the polarizability per length conserved; or eps1d = eps_zz of the cell (1.46) at the calibrated radius',
     'consequence': 'Shell radius: LO +6.6/+7.8 at 0.05 and +8.8/+10.2 at 0.075 cm^-1. Cell eps: +5.9/+5.3 for BN66 at 0.1/0.125 (thin margin; passes for BN55)',
     'ablation': 'atomic_shell_radius'},
    {'id': 'dp5', 'kind': 'regime judgement (kernel)',
     'fork': 'Full Bessel-function cylinder kernel versus its small-qt asymptotic form (the paper\'s q^2 ln q expansion) as the long-range term of the interpolation',
     'wrong_choice': 'Use the asymptotic form everywhere; it grows without bound at the grid q-points, so the subtracted long-range part is wrong',
     'consequence': 'LO +9.9/+18.8 at 0.075 and +11.3/+18.6 at 0.125 cm^-1',
     'ablation': 'asymptotic_kernel_in_interpolation'},
    {'id': 'dp6', 'kind': 'discriminating data under budget / regime of validity',
     'fork': 'Whether the Gamma run alone (analytic model with a constant short-range part, the paper\'s Eq. 2) suffices, or the 1D grid of each hidden tube must be bought to capture the mechanical dispersion',
     'wrong_choice': 'Gamma-only model D_Gamma + 1D long-range term at all q',
     'consequence': 'Fine below 0.05|GZ|, but +4.2/+5.9 at 0.075, +11.5 (BN66) at 0.1 and +14/+17 at 0.125 cm^-1',
     'ablation': 'gamma_only_long_range_model'},
    {'id': 'dp7', 'kind': 'conventions and units',
     'fork': 'q conversions (archived 2pi/alat with alat the lateral cell; |GZ| = pi/L) and the Rydberg e^2 = 2 of the long-range prefactor',
     'wrong_choice': 'Take q/|GZ| as a fraction of 2pi/L; or use e^2 = 1 with Rydberg force constants',
     'consequence': '+12/+13 at 0.01 and +54 at 0.05 (zone convention); -3.0 at 0.01 and -21/-22 at 0.05 cm^-1 (e^2)',
     'ablation': 'q_fraction_of_2pi_over_L'},
    {'id': 'dp8', 'kind': 'which quantities set the artefact',
     'fork': 'The 3D-periodic q->0 shift is set by the vacuum-diluted eps_zz of the supercell and the supercell volume, not by the cylinder model, the perpendicular component or no screening',
     'wrong_choice': 'Cylinder dielectric constant eps1d = A(eps_zz - 1)/(pi t^2) in the 3D term (or with the tube volume; or eps_xx; or eps = 1)',
     'consequence': '3D q->0 LO %.1f (cylinder eps), %.1f (cylinder eps and volume), %.1f (eps_xx), %.1f (unscreened) vs %.1f cm^-1' % (
         T3['abl_cylinder_dielectric'], T3['abl_cylinder_dielectric_and_volume'], T3['abl_perpendicular_dielectric'], T3['abl_unscreened'], C['truth_3d']),
     'ablation': 'cylinder_dielectric_for_3d_artefact'},
    {'id': 'dp9', 'kind': 'transfer across systems',
     'fork': 'Whether the polar shift curve is a property of the BN tube family or depends on each tube\'s mode charge per length and radius',
     'wrong_choice': 'Add the released BN(4,4) explicit shift curve to the hidden tube\'s Gamma LO',
     'consequence': 'LO -6.9/-11.5 at 0.05 and -7.8/-12.7 at 0.075 cm^-1 (scaling by the mode charge per length is worse: +7.5/+17 at 0.05)',
     'ablation': 'family_transfer_from_44'},
]

calibration = {}
for tb, field in (('BN55', 'BN55_LO_1D_cm-1'), ('BN66', 'BN66_LO_1D_cm-1')):
    for q, ref in C['truth'][tb].items():
        robust = {k: round(C['truth'][tb][q] + v[tb][q], 2) for k, v in R.items()}
        abls = sorted(((k, v[tb][q]) for k, v in E.items()), key=lambda kv: abs(kv[1]))
        calibration[f'{field}[{q}]'] = {
            'reference_archived_dfpt': ref, 'oracle': round(ref + O[tb][q], 2), 'oracle_error': O[tb][q],
            'robustness': robust, 'max_robustness_error': max(abs(v[tb][q]) for v in R.values()),
            'nearest_ablations': {k: round(ref + e, 2) for k, e in abls[:3]},
            'tolerance': f"+-{C['tolerances']['tolerance_cm-1'][q]} cm^-1"}
calibration['BN55_LO_3D_q0_cm-1'] = {
    'reference_archived_3d_matdyn_q0': C['truth_3d'], 'oracle_3D_code_gamma_plus_nonanalytic': T3['oracle_3D_code_gamma'],
    'robustness': {'1D_code_gamma_plus_nonanalytic': T3['robust_1D_code_gamma']},
    'nearest_ablations': {k: v for k, v in T3.items() if k.startswith('abl_')},
    'tolerance': f"+-{C['tolerances']['tolerance_3d_cm-1']} cm^-1"}

oracle_check = (CHECK or {}).get('summary', {}).get('oracle', {})
design = {
    'task_id': 'polar-phonons-1d-loto-collapse',
    'catalog_id': 'polar-phonons-1d-loto-collapse',
    'start': 'released',
    'source_group': 'rivano_2023_1d_polar_phonons',
    'title': '1D polar phonons: collapse of the LO-TO splitting in BN nanotubes',
    'decision_points': decision_points,
    'graded_fields': [
        'BN55_LO_1D_cm-1 at q/|GZ| = 0.01, 0.025, 0.05, 0.075, 0.125 (tolerances 2.0, 3.0, 4.5, 4.5, 4.5 cm^-1)',
        'BN66_LO_1D_cm-1 at q/|GZ| = 0.01, 0.025, 0.05, 0.075, 0.1, 0.125 (2.0, 3.0, 4.5, 4.5, 4.5, 4.5 cm^-1)',
        'BN55_LO_3D_q0_cm-1 within 3.0 cm^-1',
    ],
    'length_estimate': {
        'oracle_lab_evaluations': oracle_check.get('evaluations', 27),
        'oracle_charge': oracle_check.get('charge_used', 27),
        'budget': C['budget'],
        'catalog_units': C['catalog_units'],
        'oracle_container_seconds': oracle_check.get('container_seconds'),
        'tool_calls': '45-75',
        'wall_time_h': '1-2',
        'gpt6astra_cost_usd': '2.5-4.5',
        'why': 'Parsing raw QE dynamical-matrix files of 16-24-atom tubes, building a 1D long-range-corrected axial '
               'interpolation, choosing and validating the cylinder radius on a released tube, budgeted purchases '
               'across 12 folders, and the 3D artefact.',
        'harbor_oracle': {'campaign': 'build-polar-phonons-1d-loto-collapse', 'reward': 1.0, 'solver_charge': 27,
                          'total_runtime_s': 50},
        'pilot': 'pending (coordinator)'},
    'runtime_dependencies': ['numpy and scipy of the runtime image (lab only decompresses and serves archived files)'],
    'data_provenance': {
        'archive_url': 'https://archive.materialscloud.org/records/08gsg-gpd20',
        'doi': '10.24435/materialscloud:46-wj',
        'files_used': ['materialscloud_data_Rivano.zip', 'files_description.md'],
        'md5': {'materialscloud_data_Rivano.zip': 'f3ee3048a131eb831eecee865d51e705',
                'files_description.md': '0c399e8707d568e68028464c7cb683b0'},
        'member_sha256': 'benchmarks/construction/polar-phonons-1d-loto-collapse/members.json',
        'license': 'CC BY 4.0',
        'hidden_truth': 'Archived explicit 1D DFPT runs of BN(5,5) and BN(6,6) below 0.15|GZ| (longwave folders) and the '
                        'q->0 point of the archived BN(5,5) 3D-periodic matdyn.modes; stored only in environment/lab/private/'},
    'exceptions': {},
    'notes': ('Truth = the authors\' explicit DFPT runs. The authors\' own 1D-corrected interpolation (reproduced by the '
              'construction code to <0.15 cm^-1 for BN(5,5) and BN(6,6)) differs from them by up to 2.2 cm^-1, which sets '
              'the scale of the defensible-method spread. A fitted effective radius is not graded: the fit is shallow '
              '(t from 9 to 12.5 bohr for BN(5,5) passes all checks) and model dependent.'),
    'calibration': calibration,
}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print('decision points', len(decision_points), 'calibrated fields', len(calibration))
