"""Write design.json of hydride-coulomb-mustar from calibration.json and oracle_run.json (run after the gate)."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
TASK = os.path.join(HERE, '../../tasks/hydride-coulomb-mustar')
C = json.load(open(os.path.join(HERE, 'calibration.json')))
O = json.load(open(os.path.join(HERE, 'oracle_run.json'))) if os.path.exists(os.path.join(HERE, 'oracle_run.json')) else {}
ref, ab, rob = C['reference'], C['ablations'], C['robustness']
r4 = lambda x: round(x, 4)
rel = lambda d: {k: r4(v / ref['tc_full_W_K'][k] - 1) for k, v in d.items()}

dps = [
    {'id': 'dp1', 'fork': 'Bare mu = N(E_F)W(E_F,E_F) versus the pseudopotential mu* at a 1 eV cutoff',
     'wrong_choice': 'Use mu itself as the constant-DOS mu*, or Morel-Anderson-rescale mu to the cutoff',
     'consequence': f"mu* = {ab['mu_star_bare_mu']:.3f} or {ab['mu_star_morel_anderson']:.3f} instead of {ref['mu_star_cdos_1eV_Mg2IrH6_20GPa']:.3f}",
     'ablation': 'mu_star_equals_bare_mu'},
    {'id': 'dp2', 'fork': 'Smeared versus dense-grid N(E_F) in mu (mu_results "mu" vs "mu rescaled")',
     'wrong_choice': 'Transfer the sigma = 0.3 eV smeared mu (0.33) to the withheld Mg2IrH6 pressures',
     'consequence': f"Mg2IrH6 Tc +{100 * rel(ab['tc_smeared_NF_mu_sigma0.3'])['Mg2IrH6_25GPa']:.0f}-{100 * rel(ab['tc_smeared_NF_mu_sigma0.3'])['Mg2IrH6_15GPa']:.0f} %",
     'ablation': 'mu_from_smeared_NF'},
    {'id': 'dp3', 'fork': 'Which smearing row of mu_results to trust',
     'wrong_choice': 'Take the smallest sigma (0.01 eV) as the least-smeared value',
     'consequence': f"W(0,0) blows up, mu ~0.70-0.74, Mg2IrH6 Tc {100 * rel(ab['tc_mu_rescaled_sigma0.01'])['Mg2IrH6_15GPa']:.0f} / {100 * rel(ab['tc_mu_rescaled_sigma0.01'])['Mg2IrH6_25GPa']:.0f} %",
     'ablation': 'mu_from_smallest_sigma'},
    {'id': 'dp4', 'fork': 'Constant versus energy-dependent DOS in the Eliashberg solver',
     'wrong_choice': 'Constant-DOS solver with the transferred mu',
     'consequence': 'Tc +41-45 % (Mg2IrH6), +6-8 % (H3S)', 'ablation': 'constant_dos_solver'},
    {'id': 'dp5', 'fork': 'Allen-Dynes versus Migdal-Eliashberg mu* (the archived Info.txt prints a muc_AD)',
     'wrong_choice': 'Allen-Dynes Tc with the archived muc_AD of neighbouring full-W runs',
     'consequence': 'Tc +30-32 % (Mg2IrH6), +3-6 % (H3S)', 'ablation': 'allen_dynes_with_archived_mucAD'},
    {'id': 'dp6', 'fork': 'Interpolate Tc(P) or compute Tc from each pressure\'s own alpha2F',
     'wrong_choice': 'Linear interpolation of archived full-W Tc between released pressures',
     'consequence': 'H3S 180 GPa -10 %, 190 GPa -5.3 % (the H3S Tc(P) peak near 180 GPa)', 'ablation': 'tc_interpolated_in_pressure'},
    {'id': 'dp7', 'fork': 'Hold the Coulomb parameter fixed (conventional mu* or the 0 GPa mu) or let mu track N(E_F)',
     'wrong_choice': 'mu* = 0.16 everywhere; or mu = 0.62 from 0 GPa for all Mg2IrH6 pressures',
     'consequence': 'mu* = 0.16: Mg2IrH6 Tc +46-47 %, H3S +7-9 %; 0 GPa mu: Mg2IrH6 -3.9 / -6.9 %',
     'ablation': 'conventional_mustar_0.16'},
    {'id': 'dp8', 'fork': 'Cell content of the H3S calculations (conventional bcc cell, two formula units)',
     'wrong_choice': 'Report the DOS file per cell (or spin-summed) as per formula unit per spin',
     'consequence': 'N(E_F) off by a factor 2', 'ablation': 'NF_per_conventional_cell'},
    {'id': 'dp9', 'fork': 'Which solver defines the equivalent mu*',
     'wrong_choice': 'Match the full-W Tc with the vDOS solver at 7 eV instead of the constant-DOS solver at 1 eV',
     'consequence': f"mu* = {ab['mu_star_vdos_7eV']:.3f}", 'ablation': 'mu_star_equivalent_in_vdos_solver'},
]
design = {
    'task_id': 'hydride-coulomb-mustar', 'catalog_id': 'hydride-coulomb-mustar', 'start': 'released',
    'source_group': 'bozier_2026_hydride_mustar', 'title': 'Hydrides: first-principles Coulomb repulsion and Tc',
    'decision_points': dps,
    'graded_fields': ['tc_full_W_K at Mg2IrH6 15/25 GPa and H3S 180/190 GPa, each within 3.5 % of the archived IsoME full-W Tc',
                      'mu_star_cdos_1eV_Mg2IrH6_20GPa within 0.025 of the validated-solver value',
                      'N_EF_H3S_180GPa_per_fu_per_spin within 3 %'],
    'length_estimate': {'oracle_lab_evaluations': O.get('evaluations'), 'oracle_charge': O.get('charge_used'),
                        'budget': 38, 'catalog_retrieval_units': 'about 115 retrievals plus unlimited solver runs (2 units each)',
                        'tool_calls': '45-75', 'gpt6astra_cost_usd': '2.5-4.5',
                        'solver_seconds': 'cdos <1 s, vdos 5-20 s per mu*, vdos_W 10-40 s (2 cores)'},
    'runtime_dependencies': ['numpy, scipy (runtime image)'],
    'data_provenance': {'archive_url': 'https://zenodo.org/records/23235290', 'doi': '10.5281/zenodo.23235290',
                        'files_used': ['mantel_data_files_v1.tar.gz (459 members, see members.md5)'],
                        'md5': {'mantel_data_files_v1.tar.gz': '31b5f98dafc279ce627f7a1473d08d3d'}, 'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': ('Graded Tc are the authors\' archived IsoME values (released start). The lab solver reproduces 15 archived '
              'full-W runs within 1.6 % (0.1-0.25 % at the withheld pressures that have W) and the archived fixed-mu* runs '
              'within 0.5 %. The equivalent mu* reference is a validated-solver value. Abstract numbers are not graded.'),
    'calibration': {
        'tc_full_W_K': {
            'reference_archived': ref['tc_full_W_K'],
            'oracle_route_rel': C['oracle_route_rel'],
            'robustness_rel': C['robustness_rel'],
            'nearest_ablations_rel': {k: C['ablations_rel'][k] for k in C['ablations_rel']},
            'tolerance': '|rel| <= 0.035'},
        'mu_star_cdos_1eV_Mg2IrH6_20GPa': {
            'reference': r4(ref['mu_star_cdos_1eV_Mg2IrH6_20GPa']),
            'robustness': {'matched_to_lab_full_W_Tc': r4(rob['mu_star_vs_lab_TcW'])},
            'nearest_ablation': {'bare_mu': r4(ab['mu_star_bare_mu']), 'vdos_7eV_equivalent': r4(ab['mu_star_vdos_7eV']),
                                 'morel_anderson': r4(ab['mu_star_morel_anderson'])},
            'tolerance': '+-0.025'},
        'N_EF_H3S_180GPa_per_fu_per_spin': {
            'reference': r4(ref['N_EF_H3S_180GPa_per_fu_per_spin']),
            'robustness': {'nearest_grid_point': r4(rob['nef_nearest_grid'])},
            'nearest_ablation': {'per_conventional_cell': r4(ab['nef_per_cell'])},
            'tolerance': '|rel| <= 0.03'},
    },
}
json.dump(design, open(os.path.join(TASK, 'design.json'), 'w'), indent=1)
print('wrote design.json')
