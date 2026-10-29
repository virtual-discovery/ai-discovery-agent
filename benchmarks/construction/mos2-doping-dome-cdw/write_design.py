"""Write design.json (decision points, calibration table, length estimate) from reference.json and oracle_run.json.

    python3 -B write_design.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/mos2-doping-dome-cdw'
REF = json.loads((HERE/'reference.json').read_text())
ORACLE = json.loads((HERE/'oracle_run.json').read_text()) if (HERE/'oracle_run.json').exists() else {}
VALID = json.loads((HERE/'solver_validation.json').read_text()) if (HERE/'solver_validation.json').exists() else {}
T, A, R = REF['targets'], REF['ablations'], REF['robustness']


def rel(v, ref):
    return round(v / ref - 1, 4)


def catalog_cost():
    c = json.loads((TASK/'environment/agent/catalog.json').read_text())
    k = c['purchase_costs']
    return (len(c['spectra_available']) * k['spectrum'] + len(c['dos_available']) * k['dos'] + len(c['pes_available']) * k['pes']
            + 2 * k['soft_mode'] + 2 * k['mo_displacement'] + 3 * k['archived_tc_table'] + k['moments_table'] + k['a1g_raman']
            + k['experiment_dome'])


def main():
    t1, tcdw, xc, tsoc = T['tc_1x1_K'], T['tc_cdw_high_K'], T['critical_doping_e_per_cell'], T['tc_1x1_soc_at_critical_doping_K']
    dps = [
        {'id': 'dp1', 'fork': 'How to estimate Tc at the withheld dopings, given that lambda_FA jumps from 0.77 (0.14) to 1.41 (0.16)',
         'wrong_choice': 'Interpolate neighbouring Tc linearly, or average neighbouring spectra without rescaling to the tabulated lambda',
         'consequence': 'Tc(0.15) %.1f%% / %.1f%% low (tolerance 3%%)' % (100 * rel(A['tc_linear_interpolation']['0.15'], t1['0.15']),
                                                                           100 * rel(A['tc_unscaled_average_spectrum']['0.15'], t1['0.15'])),
         'ablation': 'tc_linear_interpolation_across_lambda_jump', 'kind': 'regime found in the data; discriminating purchase (moments table)'},
        {'id': 'dp2', 'fork': 'Which spectra: adiabatic, NA (shifts) or FA (shifts + linewidths)',
         'wrong_choice': 'NA spectra (the paper\'s Fig. 2d level) or adiabatic spectra',
         'consequence': 'Tc +20 to +45%; SOC maximum +23%', 'ablation': 'nonadiabatic_shifts_only_NA', 'kind': 'convention (stated)'},
        {'id': 'dp3', 'fork': 'Constant vs energy-dependent DOS for the 1x1 H phase',
         'wrong_choice': 'Constant-DOS Eliashberg (cheaper to reason about; ignores the Q-valley edge and the shift of mu at low T)',
         'consequence': 'Tc -8.3%% / -6.2%%; SOC maximum %.1f%%' % (100 * rel(A['soc_constant_dos'], tsoc)),
         'ablation': 'constant_dos_for_1x1', 'kind': 'model class'},
        {'id': 'dp4', 'fork': 'Allen-Dynes / McMillan vs Eliashberg at lambda > 1.5',
         'wrong_choice': 'Allen-Dynes formula from lambda, omega_log, omega_2',
         'consequence': 'Tc(0.17) %.2f%%, SOC maximum %.2f%%, high-doping CDW +6.3%%/+10.7%%' % (
             100 * rel(A['tc_allen_dynes']['0.17'], t1['0.17']), 100 * rel(A['soc_allen_dynes'], tsoc)),
         'ablation': 'allen_dynes_formula', 'kind': 'model class'},
        {'id': 'dp5', 'fork': 'Which criterion marks the end of the 1x1 H ground state: PES energy crossing, 1x1 soft-mode instability or CDW lower stability edge',
         'wrong_choice': 'Imaginary Gamma* (M) soft mode of the 1x1 phase (x = %.4f) or the lower edge of the CDW branch (%.4f)' % (
             A['xc_soft_mode_instability'], A['xc_cdw_spinodal']),
         'consequence': 'x_c off by +0.014 / -0.004; cascades into the SOC maximum (+28%% / %.1f%%)' % (100 * rel(A['soc_at_cdw_spinodal'], tsoc)),
         'ablation': 'boundary_at_soft_mode_instability', 'kind': 'regime judgement; coexistence and hysteresis found in the data; cascade'},
        {'id': 'dp6', 'fork': 'Doping units: PES labelled by tot_charge of the 2x2 cell; soft-mode and hysteresis tables in 10^14 cm^-2',
         'wrong_choice': 'Report the 2x2-cell charge (or the areal density) as x_c',
         'consequence': 'x_c = %.3f or %.2f instead of %.4f' % (A['xc_per_2x2_cell'], A['xc_areal_density_1e14'], xc),
         'ablation': 'doping_per_2x2_cell', 'kind': 'convention found in metadata'},
        {'id': 'dp7', 'fork': 'Whether the ground-state maximum includes spin-orbit coupling (SOC spectra exist only at 0.18 and 0.19)',
         'wrong_choice': 'Use the scalar-relativistic curve', 'consequence': 'SOC maximum +28%',
         'ablation': 'soc_omitted', 'kind': 'data coverage'},
        {'id': 'dp8', 'fork': 'Dopings of the archive spectra 2x2cdw_dop1-3, which carry no doping; the paper caption quotes n1 = 3.0, n2 = 4.3, n3 = 4.6 x 10^14 cm^-2',
         'wrong_choice': 'Trust the caption, or assign dop1-3 to the first three rows of the authors\' table',
         'consequence': 'Tc(0.35) %.0f%% / %.0f%% off' % (100 * rel(A['cdw_caption_dopings']['0.35'], tcdw['0.35']),
                                                         100 * rel(A['cdw_table_order']['0.35'], tcdw['0.35'])),
         'ablation': 'cdw_dopings_from_caption', 'kind': 'label inconsistency found by a cross-check (Allen-Dynes of the spectra vs the archived table)'},
        {'id': 'dp9', 'fork': 'What the archived Tc tables are: their header says only "Tc"; they are Allen-Dynes values (and their 1x1 A/NA columns do not match the archived A/NA spectra)',
         'wrong_choice': 'Use the archived tables as Eliashberg Tc (interpolating the 1x1 table; reading the CDW table at 1.4/1.6)',
         'consequence': 'Tc(0.15) -12.3%, Tc(0.17) -5.1%; CDW +6.3%/+10.7%',
         'ablation': 'archived_tables_taken_as_eliashberg', 'kind': 'metadata inconsistency found by a cross-check'}]
    calibration = {
        'tc_1x1_K (relative, tolerance 3%)': {
            'reference': t1, 'robustness': {k: {x: rel(v, t1[x]) for x, v in R[k].items()} for k in
                                            ('tc_custom_rigid_band', 'tc_custom_rigid_band_upper_dos', 'tc_moments_no_shape', 'tc_moments_ratio',
                                             'tc_lambda_weighted_mix')},
            'gray_zone': {k: {x: rel(v, t1[x]) for x, v in REF['gray_zone'][k].items()} for k in REF['gray_zone']},
            'gray_zone_note': 'one neighbour\'s spectral shape rescaled to the tabulated lambda ignores the tabulated omega_log; '
                              'not defensible, not graded as a wrong path',
            'nearest_ablations': {k: {x: rel(v, t1[x]) for x, v in A[k].items()} for k in
                                  ('tc_allen_dynes', 'tc_constant_dos', 'tc_unscaled_average_spectrum', 'tc_linear_interpolation',
                                   'tc_archived_table', 'tc_mcmillan', 'tc_NA_spectra', 'tc_A_spectra')},
            'oracle': ORACLE.get('errors', {}).get('tc_1x1_rel'), 'tolerance': '|rel| <= 0.03'},
        'critical_doping_e_per_cell (absolute, tolerance 0.002)': {
            'reference': round(xc, 5), 'robustness': {'quadratic': round(R['xc_quadratic'] - xc, 5), 'well_minimum': round(R['xc_cdw_well_minimum'] - xc, 5)},
            'nearest_ablations': {'cdw_spinodal': round(A['xc_cdw_spinodal'] - xc, 4), 'soft_mode': round(A['xc_soft_mode_instability'] - xc, 4),
                                  'per_2x2_cell': round(A['xc_per_2x2_cell'], 4), 'areal_density': round(A['xc_areal_density_1e14'], 3)},
            'oracle': ORACLE.get('errors', {}).get('critical_doping_abs'), 'tolerance': '|abs| <= 0.002'},
        'tc_1x1_soc_at_critical_doping_K (relative, tolerance 3%)': {
            'reference': round(tsoc, 3), 'robustness': {k: rel(R[k], tsoc) for k in ('soc_at_xc_quadratic', 'soc_at_xc_well_minimum', 'soc_ratio_method')},
            'nearest_ablations': {k: rel(A[k], tsoc) for k in ('soc_allen_dynes', 'soc_constant_dos', 'soc_at_cdw_spinodal', 'soc_mcmillan',
                                                               'soc_NA_spectra', 'soc_omitted', 'soc_at_soft_mode_instability')},
            'oracle': ORACLE.get('errors', {}).get('tc_soc_rel'), 'tolerance': '|rel| <= 0.03'},
        'tc_cdw_high_K (relative, tolerance 4%)': {
            'reference': {k: round(v, 3) for k, v in tcdw.items()}, 'robustness': 'direct lab output once the labels are assigned (spread < 0.1%)',
            'nearest_ablations': {k: {x: rel(v, tcdw[x]) for x, v in A[k].items()} for k in
                                  ('cdw_mcmillan', 'cdw_allen_dynes', 'cdw_archived_table', 'cdw_caption_dopings', 'cdw_table_order')},
            'oracle': ORACLE.get('errors', {}).get('tc_cdw_rel'), 'tolerance': '|rel| <= 0.04'}}
    design = {
        'task_id': 'mos2-doping-dome-cdw', 'catalog_id': 'mos2-doping-dome-cdw', 'start': 'released',
        'source_group': 'girottoerhardt_2025_mos2_dome',
        'title': 'Electron-doped MoS2: superconducting dome and the 1x1 to 2x2 CDW transition',
        'decision_points': dps,
        'graded_fields': ['tc_1x1_K at x = 0.15 and 0.17 (within 3%)', 'critical_doping_e_per_cell (within 0.002)',
                          'tc_1x1_soc_at_critical_doping_K (within 3%)', 'tc_cdw_high_K at x = 0.35 and 0.40 (within 4%)'],
        'length_estimate': {'oracle_lab_evaluations': ORACLE.get('evaluations'), 'oracle_charge': ORACLE.get('charge_used'),
                            'budget': 34, 'catalog_cost_retrievals_only': catalog_cost(),
                            'catalog_cost_note': 'retrievals only; every Eliashberg run costs 2 more',
                            'tool_calls': '45-75', 'gpt6astra_cost_usd': '2.5-4.5',
                            'oracle_container_seconds': ORACLE.get('container_seconds')},
        'runtime_dependencies': ['numpy, scipy (runtime image); the lab runs eliashberg.py (numpy/scipy only)'],
        'data_provenance': {'archive_url': 'https://archive.materialscloud.org/records/51yyy-ksf72',
                            'doi': '10.24435/materialscloud:z1-aw', 'license': 'CC BY 4.0',
                            'files_used': ['fig2a-d.tar.gz', 'fig2e.tar.gz', 'fig3.tar.gz', 'fig4.tar.gz', 'fig5.tar.gz', 'input.tar.gz', 'README.md'],
                            'md5': {'fig2a-d.tar.gz': '4cdb4ac5e83aaa4e3067d6304cb5e97e', 'fig2e.tar.gz': '0514e9cf3bd294ef9f49e5d79912aefb',
                                    'fig3.tar.gz': 'dbf67323c8b26e984707a9ab512f979c', 'fig4.tar.gz': 'ad524c1f940d23e3280ac430a32fe286',
                                    'fig5.tar.gz': '25bf433a60498198d5aec2e36ba287bd', 'input.tar.gz': '6f541b537c66c7de1f3f92ef30242a17',
                                    'README.md': '3b22839753c99e088435274842d0449d'},
                            'not_used': 'fig6.tar.gz (elphmod model III, shared with the excluded Schobert 2024 study), figS11.tar.gz'},
        'solver_validation': {'n_values': VALID.get('n_values'), 'max_abs_relative_deviation': VALID.get('max_abs_relative_deviation'),
                              'against': 'archived fig2e/tc_*.dat (MM, AD, CE, CA, DE, DA at mu* = 0-0.20)'},
        'exceptions': {},
        'notes': 'Dome maximum and its doping (design target 3) are not graded: the archived 2x2 CDW Tc in the coexistence window are '
                 'Allen-Dynes values without spectra, slightly above the 1x1 values, and the paper\'s headline 34.6 K (1x1 at 0.20, '
                 'Allen-Dynes) would pass any such check by recall. The graded SOC maximum of the 1x1 ground-state branch replaces it.',
        'calibration': calibration}
    (TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
    print('design.json written:', len(dps), 'decision points')


if __name__ == '__main__':
    main()
