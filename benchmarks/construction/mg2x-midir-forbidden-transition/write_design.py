"""Write design.json of mg2x-midir-forbidden-transition from calibration.json (make_task.py) and check.json (gate).

    python3 -B benchmarks/construction/mg2x-midir-forbidden-transition/write_design.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/mg2x-midir-forbidden-transition'
C = json.loads((HERE/'calibration.json').read_text())
REF, ABL, ROB, X, TOL = C['reference'], C['ablations'], C['robustness'], C['extra'], C['tolerances']
r3 = lambda v: round(float(v), 4)
mev = lambda v: round(1e3 * float(v), 1)


def main():
    check = json.loads((HERE/'check.json').read_text())['summary'] if (HERE/'check.json').exists() else {}
    oracle = check.get('oracle') or {}
    s, sa = REF['x3_minus_x1_eV'], ABL['split']
    rt, ra, rr = REF['ir_raman_ratio_300K_over_77K'], ABL['ratio'], ROB['ratio']
    pt, pa, pr = REF['predicted_splitting_eV'], ABL['pred'], ROB['pred']
    design = {
        'task_id': 'mg2x-midir-forbidden-transition',
        'catalog_id': 'mg2x-midir-forbidden-transition',
        'start': 'released',
        'source_group': 'sugawara_2026_mg2x_qsgw',
        'title': 'Mg2X mid-IR band: forbidden X-point transition, QSGW levels and temperature',
        'decision_points': [
            {'id': 'dp1', 'fork': 'Identity of the two X states: dispersion (longitudinal mass ~0.54 for X1, ~1.03 for X3; transverse masses) versus energy order, given only that X1 is lower in Mg2Si',
             'wrong_choice': 'Label by energy order (the lower state is X1 everywhere)',
             'consequence': f"Mg2Sn splittings come out positive (+{s['Mg2Sn']['LDA'] * -1:.3f}, +{-s['Mg2Sn']['G0W0']:.3f}, +{-s['Mg2Sn']['QSGW']:.3f} eV) instead of negative at all three levels",
             'ablation': 'energy_order_labels'},
            {'id': 'dp2', 'fork': 'Which TOTE column is the defined quantity: the README labels column 7 "eQP (Z = 1)" but it equals e_in + Z (column 8 - e_in), the linearized Z-renormalized G0W0 energy; in QSGW iteration files column 6 is the QSGW level and the QP columns are not',
             'wrong_choice': 'Take the column described as eQP without Z (Z = 1), or the input eigenvalue column; read the QP column of the QSGW TOTE at the 100 K lattice',
             'consequence': f"G0W0 X3 - X1 {r3(sa['g0w0_z_equal_1_column']['Mg2Si']['G0W0'])} / {r3(sa['g0w0_z_equal_1_column']['Mg2Ge']['G0W0'])} / {r3(sa['g0w0_z_equal_1_column']['Mg2Sn']['G0W0'])} eV instead of {r3(s['Mg2Si']['G0W0'])} / {r3(s['Mg2Ge']['G0W0'])} / {r3(s['Mg2Sn']['G0W0'])}; QP column of the QSGW file gives {r3(pa['qsgw_column7_low_temperature']['Mg2Si_T0'])} eV at T0",
             'ablation': 'g0w0_z_equal_1_column'},
            {'id': 'dp3', 'fork': 'Ladder (vertex) correction: QSGW+ladder minus QSGW within the same code and basis, versus across codes or the QSGW80 scaling that mimics the vertex for the gap',
             'wrong_choice': 'Subtract the ecalj QSGW splitting from the Questaal ladder run, or take QSGW80 - QSGW',
             'consequence': f"{X['checks']['vertex']['cross_code_questaal_ladder_minus_ecalj_qsgw']:+.1f} or {X['checks']['vertex']['qsgw80_minus_qsgw']:+.1f} meV instead of {REF['vertex_correction_Mg2Si_meV']:+.1f} meV; QSGW80 taken as QSGW also fails the QSGW splittings",
             'ablation': 'vertex_cross_code'},
            {'id': 'dp4', 'fork': 'Lowest non-vanishing order of the X_lo -> X_hi matrix element: forbidden at X, |p|^2 = beta^2 q_perp^2 (from the PBE matrix elements: zero at X and along q_par, quadratic in q_perp)',
             'wrong_choice': 'Treat the transition as allowed (constant |p|^2), or take |p|^2 ~ |q|^2 including q_par, or average the full PBE |p|^2',
             'consequence': f"Q = 1 (constant), Mg2Si {r3(ra['isotropic_q2']['Mg2Si_1e16'])} / {r3(ra['isotropic_q2']['Mg2Si_3e18'])} and Mg2Sn {r3(ra['isotropic_q2']['Mg2Sn_1e16'])} (|q|^2), Mg2Si {r3(ra['full_pbe_matrix_element']['Mg2Si_1e16'])} (full PBE |p|^2) instead of {r3(rt['Mg2Si_1e16'])} / {r3(rt['Mg2Si_3e18'])} / {r3(rt['Mg2Sn_1e16'])}",
             'ablation': 'constant_matrix_element'},
            {'id': 'dp5', 'fork': 'Thermal average over the calculated dispersion of the whole transverse plane (box) versus a parabolic valley or the X-W line taken as isotropic; the lower state flattens along the [011] diagonal',
             'wrong_choice': 'Use the parabolic masses fitted to the scans, or the X-W scan as an isotropic transverse dispersion',
             'consequence': f"Mg2Si Q {r3(ra['parabolic']['Mg2Si_1e16'])} / {r3(ra['parabolic']['Mg2Si_3e18'])}, Mg2Sn {r3(ra['parabolic']['Mg2Sn_1e16'])} (parabolic); X-W line {r3(ra['xw_line_isotropic']['Mg2Si_1e16'])} / {r3(ra['xw_line_isotropic']['Mg2Sn_1e16'])}",
             'ablation': 'parabolic_valley'},
            {'id': 'dp6', 'fork': 'Degeneracy at 3e18 cm^-3 (77 K) and the counting of the box: Fermi-Dirac with 2 spins x 3 valleys x 8 octants, half weights on the symmetry planes',
             'wrong_choice': 'Boltzmann statistics at every density; forgetting the spin, the octants, or the half weights of the symmetry-plane points',
             'consequence': f"Mg2Si 3e18: {r3(ra['boltzmann']['Mg2Si_3e18'])} (Boltzmann), {r3(ra['spin_not_counted']['Mg2Si_3e18'])} (no spin), {r3(ra['one_octant_only']['Mg2Si_3e18'])} (one octant); no plane weights {r3(ra['no_symmetry_plane_weights']['Mg2Si_1e16'])} at 1e16; reference {r3(rt['Mg2Si_3e18'])}",
             'ablation': 'boltzmann_statistics'},
            {'id': 'dp7', 'fork': 'Whether the Mg2Si result transfers to Mg2Sn, whose lower state is X3 and whose small splitting makes the valley strongly non-parabolic (a second box purchase)',
             'wrong_choice': 'Reuse the Mg2Si Q for Mg2Sn, or occupy the X1 state of Mg2Sn',
             'consequence': f"Mg2Sn Q {r3(ra['mg2si_value_for_mg2sn']['Mg2Sn_1e16'])} or {r3(ra['x1_occupied_in_mg2sn']['Mg2Sn_1e16'])} instead of {r3(rt['Mg2Sn_1e16'])}",
             'ablation': 'mg2si_ratio_for_mg2sn'},
            {'id': 'dp8', 'fork': 'Lattice constant matched to the temperature: the production QSGW uses the room-temperature lattice; the T -> 0 prediction needs the 100 K lattice run (+23 meV)',
             'wrong_choice': 'Add the zero-point renormalization to the production (room-temperature lattice) splitting; or use the 100 K lattice at 300 K',
             'consequence': f"Mg2Si T0 {r3(pa['room_temperature_lattice_at_T0']['Mg2Si_T0'])} eV instead of {r3(pt['Mg2Si_T0'])}; 300 K {r3(pa['low_temperature_lattice_at_300K']['Mg2Si_300K'])} instead of {r3(pt['Mg2Si_300K'])}",
             'ablation': 'room_temperature_lattice_at_T0'},
            {'id': 'dp9', 'fork': 'Electron-phonon renormalization of X3 - X1 (positive: X1 drops more than X3), zero point at T -> 0 and the full 300 K value at 300 K',
             'wrong_choice': 'Ignore it, take the zero-point value at 300 K, or apply it with the sign of X1 - X3',
             'consequence': f"Mg2Si T0 {r3(pa['no_electron_phonon']['Mg2Si_T0'])}, 300 K {r3(pa['no_electron_phonon']['Mg2Si_300K'])} / {r3(pa['zero_point_only_at_300K']['Mg2Si_300K'])}, Mg2Ge 300 K {r3(pa['no_electron_phonon']['Mg2Ge_300K'])} / {r3(pa['zero_point_only_at_300K']['Mg2Ge_300K'])} eV instead of {r3(pt['Mg2Si_T0'])}, {r3(pt['Mg2Si_300K'])}, {r3(pt['Mg2Ge_300K'])}",
             'ablation': 'no_electron_phonon'},
        ],
        'graded_fields': [
            f"x3_minus_x1_eV: signed E(X3) - E(X1) at LDA, G0W0 and QSGW for Mg2Si, Mg2Ge, Mg2Sn (9 values), +-{TOL['x3_minus_x1_eV']} eV",
            f"vertex_correction_Mg2Si_meV, +-{TOL['vertex_meV']} meV",
            f"ir_raman_ratio_300K_over_77K: Mg2Si at 1e16 and 3e18 cm^-3, Mg2Sn at 1e16, within {100 * TOL['ratio_rel']:.0f}% (relative)",
            f"predicted_splitting_eV: Mg2Si T0 and 300 K, Mg2Ge 300 K, +-{TOL['predicted_eV']} eV",
        ],
        'length_estimate': {
            'oracle_lab_evaluations': oracle.get('evaluations', 19), 'oracle_charge': oracle.get('charge_used', 27), 'budget': 32,
            'catalog_units': 'about 114 (bands 6, x_scan 12, qsgw_box 5 x 5, momentum_matrix 27 for the full grid, tote 11, questaal 8, pseudopotential_g0w0 7, electron_phonon 18)',
            'minimum_charge_for_all_targets': 24,
            'tool_calls': '45-75', 'gpt6astra_cost_usd': '2.5-4.5',
            'why': 'Identifying X1/X3 from dispersions in three compounds and two levels, decoding the TOTE columns against the README, a same-code ladder comparison, deriving the selection rule from the PBE matrix elements, Fermi-Dirac averages over octant boxes for two compounds (two 5-unit purchases), and assembling temperature-consistent predictions from lattice and electron-phonon runs under a budget of about 28% of the catalog.',
        },
        'runtime_dependencies': ['numpy (runtime image); no extra layer'],
        'data_provenance': {
            'archive_url': 'https://zenodo.org/records/23319151', 'doi': '10.5281/zenodo.23319151',
            'files_used': ['20261012aa-paper1-zenodo-archive.zip', 'README.md', 'LICENSES.md', 'MD5SUMS.txt'],
            'md5': {'20261012aa-paper1-zenodo-archive.zip': '865ba37d739a8118a8a556cb95c643bb', 'README.md': '246ab9331d014b4665cb1629ca218586',
                    'LICENSES.md': 'aed582f0574ab7043a3627ea9ac59b82', 'MD5SUMS.txt': '91946a35a975e8bfcda1d31e4a6db173'},
            'members': 'benchmarks/construction/mg2x-midir-forbidden-transition/members.md5 (776 members)',
            'license': 'Mixed per file: computed outputs and docs CC BY 4.0 (served), digitized points CC0 1.0 (staged), programs MIT (ported), ecalj/Questaal inputs AGPL-3.0 (not used in the package)'},
        'exceptions': {},
        'notes': ('Adapted from the candidate design (see the construction README). The design staged the QSGW path bands, which shows the '
                  'QSGW splitting directly, and graded it; here every level is purchasable and the splittings are graded with their sign '
                  '(the Mg2Sn inversion), while the difficulty comes from evidence in the files: the garbled TOTE column label, the '
                  'same-code ladder comparison, the selection rule from the PBE matrix elements, the box versus parabolic valley, '
                  'degeneracy, the Mg2Sn non-parabolicity, and the temperature consistency of the lattice. The design\'s experimental '
                  'peak reading and absorption-peak targets were dropped (they reduce to the splitting within a few meV, and the reading '
                  'is estimator dependent); they remain in the rubric.'),
        'calibration': {
            'x3_minus_x1_eV': {
                'reference': s, 'robustness': 'scan, path and box agree at X to < 1e-5 eV; TOTE column 6 of the production QSGW and LDA within 0.5 meV',
                'nearest_ablation': {'G0W0 from the PBE start (Mg2Si, Mg2Ge)': [r3(sa['g0w0_from_pbe']['Mg2Si']['G0W0']), r3(sa['g0w0_from_pbe']['Mg2Ge']['G0W0'])],
                                     'QSGW80 as QSGW': {k: r3(v['QSGW']) for k, v in sa['qsgw80_as_qsgw'].items()},
                                     'energy-order labels (Mg2Sn)': {k: r3(v) for k, v in sa['energy_order_labels']['Mg2Sn'].items()},
                                     'G0W0 Z = 1 column': {k: r3(v['G0W0']) for k, v in sa['g0w0_z_equal_1_column'].items()}},
                'tolerance': f"+-{TOL['x3_minus_x1_eV']} eV"},
            'vertex_correction_Mg2Si_meV': {
                'reference': round(REF['vertex_correction_Mg2Si_meV'], 3),
                'robustness': {'ladder window Nc 8': round(X['checks']['vertex']['nc8'], 2), 'ladder window Nc 12': round(X['checks']['vertex']['nc12'], 2),
                               'larger basis': round(X['checks']['vertex']['basis2'], 2)},
                'nearest_ablation': {'cross code': round(ABL['vertex']['cross_code'], 1), 'QSGW80 - QSGW': round(ABL['vertex']['qsgw80_scaling'], 1),
                                     'gap change with ladder': round(ABL['vertex']['gap_change'], 1)},
                'tolerance': f"+-{TOL['vertex_meV']} meV"},
            'ir_raman_ratio_300K_over_77K': {
                'reference': {k: r3(v) for k, v in rt.items()},
                'robustness': {name: {k: r3(v) for k, v in vals.items()} for name, vals in rr.items()},
                'nearest_ablation': {name: {k: r3(v) for k, v in vals.items()} for name, vals in ra.items()},
                'tolerance': f"+-{100 * TOL['ratio_rel']:.0f}% relative"},
            'predicted_splitting_eV': {
                'reference': {k: r3(v) for k, v in pt.items()},
                'components_meV': {'lattice 300 K -> 100 K (Mg2Si QSGW)': round(X['lattice_shift_meV'], 1),
                                   'AHC dE (24^3, eta 10 meV)': {k: round(v, 2) for k, v in X['eph_dE_meV'].items()}},
                'robustness': {name: {k: r3(v) for k, v in vals.items()} for name, vals in pr.items()},
                'nearest_ablation': {name: {k: r3(v) for k, v in vals.items()} for name, vals in pa.items()},
                'tolerance': f"+-{TOL['predicted_eV']} eV"},
        },
    }
    (TASK/'design.json').write_text(json.dumps(design, indent=1, ensure_ascii=False) + '\n')
    print('wrote design.json')


if __name__ == '__main__':
    main()
