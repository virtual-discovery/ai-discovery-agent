"""Write design.json of 2d-qsh-shc-plateau from calibration.json (make_task.py), the reference calculations and
check.json (gate). Run on the host: python3 -B write_design.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/2d-qsh-shc-plateau'
cal = json.loads((HERE/'calibration.json').read_text())
abl = json.loads((TASK/'solution/ablation_values.json').read_text())
check = json.loads((HERE/'check.json').read_text()) if (HERE/'check.json').exists() else {}
oracle_run = json.loads((HERE/'oracle_run.json').read_text()) if (HERE/'oracle_run.json').exists() else {}
P, Z, D = cal['plateau'], cal['z2'], cal['spin_direction']
A, R = abl['ablations'], abl['robustness']
f3 = lambda x: round(x, 3)


def pv(name, mat):
    return f3(A[name]['plateau_sigma_z'][mat])


dps = [
    {'id': 'dp1', 'fork': 'Where to put the chemical potential: the scf Fermi level (cold smearing 0.01 Ry) of the small-gap '
                          'models lies inside the valence band, and the .win fermi_energy is the valence-band maximum',
     'wrong_choice': 'Evaluate the SHC at the Fermi energy printed by aiida.out (or at the .win fermi_energy)',
     'consequence': f"CuLi2As {pv('plateau_at_scf_fermi', 'CuLi2As_P-6m2')} instead of {f3(P['CuLi2As_P-6m2'])}, "
                    f"Ti2I2N2 {pv('plateau_at_scf_fermi', 'Ti2I2N2_Pmm2')} instead of {f3(P['Ti2I2N2_Pmm2'])}, "
                    f"Li2Tl2 {pv('plateau_at_scf_fermi', 'Li2Tl2_P3m1')} instead of {f3(P['Li2Tl2_P3m1'])}",
     'ablation': 'plateau_at_scf_fermi'},
    {'id': 'dp2', 'fork': 'Gap of the Wannier model versus the DFT band edges written in aiida.win: they come from the coarse '
                          'scf grid; for Bi2In2 the Wannier valence-band maximum at a generic k lies above the .win mid-gap',
     'wrong_choice': 'Place mu at the mid-gap of the .win VBM/CBM comments',
     'consequence': f"Bi2In2 {pv('plateau_at_dft_mid_gap', 'Bi2In2_P4_nmm')} instead of {f3(P['Bi2In2_P4_nmm'])}",
     'ablation': 'plateau_at_dft_mid_gap'},
    {'id': 'dp3', 'fork': 'k-mesh convergence inside small gaps: the spin Berry curvature of Ti2I2N2 (22 meV gap) and Li2Tl2 '
                          '(31 meV) is concentrated in hot spots that coarse meshes alias; dense or adaptively refined meshes are needed',
     'wrong_choice': 'Accept a 120 x 120 (or 60 x 60) mesh without refinement or a convergence test',
     'consequence': f"Ti2I2N2 {pv('mesh_120_unrefined', 'Ti2I2N2_Pmm2')} (120^2) or {pv('mesh_60_unrefined', 'Ti2I2N2_Pmm2')} (60^2) "
                    f"instead of {f3(P['Ti2I2N2_Pmm2'])}; Li2Tl2 {pv('mesh_60_unrefined', 'Li2Tl2_P3m1')} at 60^2",
     'ablation': 'mesh_120_unrefined'},
    {'id': 'dp4', 'fork': 'Velocity-operator convention: Wannier functions at the lattice points (phase convention of the hr.dat) '
                          'versus at their centres; the lattice convention breaks the crystal symmetry of the Kubo response, '
                          'visible as a sigma^x_xy forbidden by the mirror of Bi2',
     'wrong_choice': 'Use the lattice convention',
     'consequence': f"Bi2 spin direction theta/phi {f3(A['lattice_position_convention']['spin_direction']['Bi2_P-3m1']['theta_deg'])}/"
                    f"{f3(A['lattice_position_convention']['spin_direction']['Bi2_P-3m1']['phi_deg'])} instead of "
                    f"{f3(D['Bi2_P-3m1']['theta_deg'])}/{f3(D['Bi2_P-3m1']['phi_deg'])} deg; Hg4O2 plateau "
                    f"{pv('lattice_position_convention', 'Hg4O2_P1')} instead of {f3(P['Hg4O2_P1'])}",
     'ablation': 'lattice_position_convention'},
    {'id': 'dp5', 'fork': 'Which calculation the archived spin Hall files describe: Top_materials holds relaxed structures '
                          '(other lattice constants, Hf2Br2 cell rotated), not the database models; the paper table matches them',
     'wrong_choice': 'Report the archived Top_materials plateaus (or the paper table) for the database models',
     'consequence': f"Bi2 {pv('relaxed_archive_values', 'Bi2_P-3m1')}, Hf2Br2 {pv('relaxed_archive_values', 'Hf2Br2_P-3m1')}, "
                    f"Ti2I2N2 {pv('relaxed_archive_values', 'Ti2I2N2_Pmm2')} instead of {f3(P['Bi2_P-3m1'])}, "
                    f"{f3(P['Hf2Br2_P-3m1'])}, {f3(P['Ti2I2N2_Pmm2'])}",
     'ablation': 'relaxed_archive_values'},
    {'id': 'dp6', 'fork': 'Topology versus SHC magnitude: a large non-quantized in-gap SHC does not imply Z2 = 1 (Bi2In2 is '
                          'trivial with a plateau near 0.8) and a plateau far from 1 does not imply a trivial insulator',
     'wrong_choice': 'Classify by the plateau magnitude (|sigma^z| >= 0.5 means QSHI)',
     'consequence': 'Bi2In2 classified Z2 = 1',
     'ablation': 'z2_from_plateau_magnitude'},
    {'id': 'dp7', 'fork': 'Resolution of the Wilson-loop flow: near the small direct gaps of Ti2I2N2 and Li2Tl2 the hybrid '
                          'Wannier centres wind quickly; a coarse sampling misses crossings',
     'wrong_choice': 'Trust a single 51 x 80 Wilson-loop calculation',
     'consequence': 'Z2 = 0 for ' + ', '.join(n for n, v in A['z2_from_coarse_wcc']['z2'].items() if v != Z[n]),
     'ablation': 'z2_from_coarse_wcc'},
    {'id': 'dp8', 'fork': 'Normalization: the service (like Wannier90) returns (hbar/e) S/cm with the vacuum-padded cell volume; '
                          'the sheet value needs each material\'s cell height c and the unit e/2pi',
     'wrong_choice': 'Normalize by the atomic thickness of the layer, or report in units of e/4pi',
     'consequence': 'plateaus off by c/t (factors 3-10) or by 2',
     'ablation': 'atomic_thickness_normalisation'},
    {'id': 'dp9', 'fork': 'Spin axis of the plateau: low symmetry allows in-plane spin components (Bi2: sigma^y_xy; Hg4O2, P1: '
                          'sigma^x_xy), so all three components must be computed',
     'wrong_choice': 'Compute only sigma^z_xy and take the spin along z',
     'consequence': f"theta 0 or 180 deg instead of {f3(D['Bi2_P-3m1']['theta_deg'])} (Bi2) and {f3(D['Hg4O2_P1']['theta_deg'])} (Hg4O2)",
     'ablation': 'sigma_z_only_spin_direction'},
]

calibration = {}
for n in P:
    near = {k: f3(v['plateau_sigma_z'][n]) for k, v in A.items() if 'plateau_sigma_z' in v and abs(v['plateau_sigma_z'][n] - P[n]) > 1e-9}
    nearest = dict(sorted(near.items(), key=lambda kv: abs(kv[1] - P[n]))[:3])
    calibration[f'plateau_sigma_z:{n}'] = {
        'reference': round(P[n], 4),
        'robustness': {k: round(v['plateau_sigma_z'][n], 4) for k, v in R.items() if 'plateau_sigma_z' in v},
        'nearest_ablations': nearest,
        'tolerance': f"+-max(3% of |ref|, 0.02) = +-{round(max(0.03 * abs(P[n]), 0.02), 4)}"}
for n in D:
    for ang, tol in (('theta_deg', 2.5), ('phi_deg', 5.0)):
        calibration[f'{ang}:{n}'] = {
            'reference': round(D[n][ang], 2),
            'robustness': {k: round(v['spin_direction'][n][ang], 2) for k, v in R.items() if 'spin_direction' in v},
            'nearest_ablations': {k: round(v['spin_direction'][n][ang], 2) for k, v in A.items() if 'spin_direction' in v},
            'tolerance': f'+-{tol} deg'}
calibration['z2'] = {'reference': Z, 'robustness': {k: v['z2'] for k, v in R.items() if 'z2' in v},
                     'nearest_ablations': {k: v['z2'] for k, v in A.items() if 'z2' in v}, 'tolerance': 'exact'}

summary = check.get('summary', {})
design = {
    'task_id': '2d-qsh-shc-plateau', 'catalog_id': '2d-qsh-shc-plateau', 'start': 'released',
    'source_group': 'zhou_2025_2d_shc_database', 'title': '2D quantum spin Hall insulators: spin Hall plateaus',
    'decision_points': dps,
    'graded_fields': ['z2 of 8 materials (exact)',
                      'plateau_sigma_z of 8 materials within max(3% of |ref|, 0.02) (hbar/e)(e^2/h)',
                      'spin_direction theta within 2.5 deg and phi within 5 deg for Bi2_P-3m1 and Hg4O2_P1'],
    'length_estimate': {'oracle_lab_evaluations': summary.get('oracle', {}).get('evaluations'),
                        'oracle_charge': summary.get('oracle', {}).get('charge_used'),
                        'budget': summary.get('oracle', {}).get('budget'),
                        'oracle_job_seconds': oracle_run.get('total_job_seconds'), 'oracle_gate_container_seconds_full_run': 642,
                        'tool_calls': '50-80', 'gpt6astra_cost_usd': '3-5',
                        'why': 'Eight models x (band edges, Wilson loop, converged spin Hall plateau) with resolution-dependent '
                               'costs; the Fermi-level traps need the staged outputs read against the Wannier band edges, the '
                               'position convention needs a symmetry test, and two canted plateaus need all three components.',
                        'oracle_errors': oracle_run.get('errors')},
    'runtime_dependencies': ['numpy and scipy (runtime image); no extra layer'],
    'data_provenance': {'archive_url': 'https://archive.materialscloud.org/records/rfcdy-wcy32',
                        'doi': '10.24435/materialscloud:yv-h7', 'paper_doi': '10.1038/s41699-025-00562-4',
                        'files_used': ['Data.7z', 'README.txt', 'files_description.md'],
                        'md5': {'Data.7z': '1d261e780f092b8c7c1f2e051855c798', 'README.txt': 'aedc173ccb0587b13fad1cbbcd5f58cc',
                                'files_description.md': '10c58cb5ed5d38de19f445f2fca21c19'},
                        'members': 'benchmarks/construction/2d-qsh-shc-plateau/members.md5', 'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': 'Adapted from the candidate design (see the construction README): graded values are validated-evaluator values '
             'of the archived database Wannier models (the evaluator reproduces the archived aiida_band.dat to <0.1 meV '
             'with the Wigner-Seitz image selection, and gives sigma^x_xy = 0.002 / sigma^y_xy = 0.352 for Bi2, where the '
             'archived relaxed-structure Wannier90 run gives 0.009 / 0.356). The archive has no hr.dat for the relaxed '
             'Top_materials structures and no .spn files; Y2C2I2 (a metal, P1 in the database, in a rotated cell) was '
             'replaced by the canted plateaus of Bi2 and Hg4O2. Edge-state counts are a rubric item.',
    'calibration': calibration,
}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print(json.dumps({k: design[k] for k in ('length_estimate',)}, indent=1))
