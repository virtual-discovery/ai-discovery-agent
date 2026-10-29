"""Write the task's design.json (decision points, graded fields, calibration) from calibration.json.

    python3 -B write_design.py          # after make_task.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/gan-irradiation-defect-levels'
C = json.loads((HERE/'calibration.json').read_text())
T = json.loads((TASK/'environment/lab/private/hidden.json').read_text())['thresholds']
DIEL = json.loads((HERE/'dielectric.json').read_text())

DPS = [
    {'id': 'dp1', 'kind': 'evidence in run metadata (k-point) and eigenvalues',
     'fork': 'Band edges: every defect cell and `perfect` use the single k-point (1/2,1/2,1/2), where the supercell '
             'eigenvalues are folded states (gap ~6.5 eV); the edges must come from the Gamma-only perfect cell with the '
             'same PAW set',
     'wrong_choice': 'Take the VBM/CBM from the HOMO/LUMO of `perfect` (R point), or the perfect-cell Fermi energy as VBM',
     'consequence': 'Gap 6.53 eV and every level +0.66 eV (R point); levels +0.33 eV (Fermi energy); pair energies at '
                    'E_CBM - 0.5 eV off by eV',
     'ablation': 'band_edges_from_R_point_supercell', 'other_ablations': ['vbm_from_perfect_cell_fermi_energy']},
    {'id': 'dp2', 'kind': 'metadata consistency, cascade into pair energies',
     'fork': 'Host total energy: the perfect cell computed like the defect cells (same k-point and PAW set) versus the '
             'Gamma-only perfect cell that supplies the band edges',
     'wrong_choice': 'Use the Gamma-only perfect cell (4.71 eV higher) as the host energy',
     'consequence': 'Levels unchanged but E_b -4.71 eV (divacancy apparently unbound) and both Frenkel energies -9.42 eV',
     'ablation': 'host_energy_from_gamma_only_supercell'},
    {'id': 'dp3', 'kind': 'model choice (screening of relaxed charged defects)',
     'fork': 'Dielectric screening for relaxed charged defects: static eps_0 versus electronic eps_inf, or no correction',
     'wrong_choice': 'eps_inf (corrections ~1.8x larger), or no finite-size correction',
     'consequence': 'eps_inf makes the 2+ states of V_N (0.12 eV window) and Ga_i (0.17 eV) stable, moves V_Ga (2-/3-) by '
                    '+0.46 eV, E_b by +0.55 eV and the Ga Frenkel energy by +0.92 eV; no correction compresses the V_Ga '
                    'ladder (levels within 0.17 eV) and removes the V_N neutral state',
     'ablation': 'electronic_dielectric_tensor', 'other_ablations': ['no_finite_size_correction']},
    {'id': 'dp4', 'kind': 'model choice, cascade',
     'fork': 'Potential alignment from the core potentials far from the defect (eFNV) versus the point-charge term alone',
     'wrong_choice': 'Makov-Payne-like point-charge correction without alignment',
     'consequence': 'V_Ga levels -0.10 to -0.13 eV, V_N (3+/+) -0.13 eV, Frenkel energies -0.41 (Ga) and +0.15 eV (N)',
     'ablation': 'point_charge_term_without_alignment'},
    {'id': 'dp5', 'kind': 'evidence in raw outputs (sign convention)',
     'fork': "Sign of VASP's 'average (electrostatic) potential at core': it is an electron potential energy, so the "
             'defect-induced change anticorrelates with the point-charge potential (r ~ -0.85 far from localized defects)',
     'wrong_choice': 'Read the core potentials as electrostatic potentials, flipping the sign of the DFT part of the alignment',
     'consequence': 'V_Ga ladder compressed (1.27-1.70 eV), V_N loses its neutral state, Ga_i (3+/+) +0.61 eV, E_b -0.59 eV, '
                    'Ga Frenkel energy -2.1 eV',
     'ablation': 'core_potentials_read_as_electrostatic_potentials'},
    {'id': 'dp6', 'kind': 'discriminating purchases under budget (site switching with charge)',
     'fork': 'Lowest-energy configuration per charge state: buy both interstitial sites for each charge',
     'wrong_choice': 'Use one archived site (i1) for every charge state',
     'consequence': 'Ga_i (3+/+) +0.52 eV plus a spurious (+/0) at 3.25 eV; no N_i level left in the graded window; Frenkel '
                    'energies +1.04 (Ga) and +1.11 eV (N)',
     'ablation': 'single_interstitial_site_i1'},
    {'id': 'dp7', 'kind': 'thermodynamic stability (negative U)',
     'fork': 'Levels between consecutive thermodynamically stable states versus every adjacent pair (q/q-1)',
     'wrong_choice': 'Report (3+/2+) and (2+/+) for V_N and Ga_i although 2+ is never stable',
     'consequence': 'Two levels in inverted order instead of the single (3+/+) level (V_N: U = -0.07 eV; Ga_i: U ~ 0, '
                    'with any 2+ window below the 0.1 eV merge width)',
     'ablation': 'adjacent_charge_levels_without_stability'},
    {'id': 'dp8', 'kind': 'metadata inconsistency (PAW set)',
     'fork': 'Which Gamma-point perfect run gives the band edges: WAV_ref (Gamma, on the exact host geometry) and '
             'kpt_gam/d10 use the Ga_d PAW set (864 electrons) with an incompatible eigenvalue scale; kpt_gam uses the '
             "defect cells' Ga set",
     'wrong_choice': 'Take the band edges from WAV_ref (its gap, 3.51 eV, looks reasonable)',
     'consequence': 'VBM 2.33 eV too high: most levels collapse below midgap or vanish; gap -0.05 eV',
     'ablation': 'band_edges_from_Ga_d_paw_run'},
]


def lv(r, d):
    return r['levels'][d]


def main():
    ref = C['reference']
    cal = {'band_gap_eV': {'reference': ref['gap'], 'robustness': {k: C[k]['gap'] for k in ROB},
                           'nearest_ablation': {'Ga_d (WAV_ref)': C['edges_Ga_d']['gap'], 'R_point': C['edges_R_point']['gap'],
                                                'fermi_energy_as_vbm': C['edges_efermi']['gap']},
                           'tolerance': '+-%g' % T['gap_tol_eV']}}
    for d in ['V_Ga', 'V_N', 'Ga_i', 'N_i', 'C_N']:
        cal['levels_' + d] = {'reference': lv(ref, d), 'robustness': {k: lv(C[k], d) for k in ROB},
                              'ablations': {k: lv(C[k], d) for k in ABL},
                              'tolerance': '+-%g eV per level, same charge pair; levels within %g eV of a band edge not graded'
                                           % (T['level_tol_eV'], T['level_window_eV'])}
    for key, name in [('Eb', 'divacancy_binding_energy_eV'), ('FP_N', 'frenkel_pair_energy_eV.N'),
                      ('FP_Ga', 'frenkel_pair_energy_eV.Ga')]:
        tol = {'Eb': T['binding_tol_eV'], 'FP_N': T['frenkel_tol_eV']['N'], 'FP_Ga': T['frenkel_tol_eV']['Ga']}[key]
        cal[name] = {'reference': ref[key], 'robustness': {k: C[k][key] for k in ROB},
                     'ablations': {k: C[k][key] for k in ABL}, 'tolerance': '+-%g' % tol}
    design = json.loads((TASK/'design.json').read_text())
    design.update({
        'decision_points': DPS,
        'graded_fields': [
            'band_gap_eV within +-%g eV' % T['gap_tol_eV'],
            'transition_levels of V_Ga, V_N, Ga_i, N_i, C_N: every reference level inside [%g, gap-%g] eV reported with the '
            'same (from, to) pair within +-%g eV; no other level inside [%g, gap-%g] eV'
            % (T['level_window_eV'], T['level_window_eV'], T['level_tol_eV'], T['level_window_eV'] + T['level_tol_eV'],
               T['level_window_eV'] + T['level_tol_eV']),
            'divacancy_binding_energy_eV within +-%g eV' % T['binding_tol_eV'],
            'frenkel_pair_energy_eV: N within +-%g eV, Ga within +-%g eV' % (T['frenkel_tol_eV']['N'], T['frenkel_tol_eV']['Ga']),
        ],
        'calibration': cal,
    })
    (TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
    print('design.json written:', len(DPS), 'decision points')


ROB = ['edges_aligned_all_sites', 'eps_0_isotropic', 'region_beyond_5.5A', 'region_beyond_4.0A']
ABL = ['edges_R_point', 'edges_efermi', 'host_gamma_only', 'eps_inf', 'no_correction', 'pc_only', 'alignment_wrong_sign',
       'single_site_i1', 'edges_Ga_d', 'single_site_i2', 'neutral_states', 'edges_aligned_N_sites']

if __name__ == '__main__':
    main()
