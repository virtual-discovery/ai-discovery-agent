"""Write design.json for gan-hole-piezomobility, including the calibration table computed from the staged files.

    python3 -B write_design.py      # after make_task.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/gan-hole-piezomobility'
H = json.loads((TASK/'environment/lab/private/hidden.json').read_text())
V = json.loads((TASK/'solution/ablation_values.json').read_text())
T, G, TH = H['targets'], H['graded'], H['thresholds']
ORACLE = json.loads((HERE/'oracle_run.json').read_text()) if (HERE/'oracle_run.json').exists() else {}


def rel(pred, lab, k):
    return round(pred[lab][k] / T[lab][k] - 1, 4)


def field_row(lab, k):
    row = {'reference_value': round(T[lab][k], 2), 'oracle_relative_error': rel(V['reference'], lab, k)}
    if lab in TH['window'] and k in TH['window'][lab]:
        row = {'reference_value': round(T[lab][k], 2), 'oracle_value': round(V['reference'][lab][k], 1),
               'robustness_values': {n: round(v[lab][k], 1) for n, v in V['robustness'].items()
                                     if isinstance(v, dict) and lab in v and round(v[lab][k], 3) != round(V['reference'][lab][k], 3)},
               'nearest_ablation_values': {n: round(v[lab][k], 1) for n, v in V['ablations'].items()
                                           if isinstance(v, dict) and lab in v and round(v[lab][k], 3) != round(V['reference'][lab][k], 3)},
               'tolerance': f"window {TH['window'][lab][k]}"}
        return row
    row['robustness_relative_errors'] = {n: rel(v, lab, k) for n, v in V['robustness'].items()
                                         if isinstance(v, dict) and lab in v and rel(v, lab, k) != rel(V['reference'], lab, k)}
    row['ablation_relative_errors'] = {n: rel(v, lab, k) for n, v in V['ablations'].items()
                                       if isinstance(v, dict) and lab in v and rel(v, lab, k) != rel(V['reference'], lab, k)}
    row['tolerance'] = f"+/-{int(100 * (TH['diag_rel'] if k in ('xx', 'yy', 'zz') else TH['offdiag_rel']))}% (relative, signed)"
    return row


calibration = {f'{lab}:{k}': field_row(lab, k) for lab, keys in G.items() for k in keys}
bz, bi = H['best_mu_zz_2pct'], H['best_inplane_2pct']
calibration['best_mu_zz_2pct'] = {
    'reference': {'state': bz['state'], 'strain_voigt': bz['strain_voigt'], 'mu_zz': round(bz['mu_zz'], 2)},
    'runner_up': {k: round(v, 2) for k, v in sorted(H['two_pct_ranking']['mu_zz'].items(), key=lambda x: -x[1])[1:4]},
    'nearest_ablation': {'serta_survey': V['ablations']['best_mu_zz_from_serta']},
    'tolerance': 'strain: dominant component e3 = -0.02 +/- 0.002, others |e| <= 0.002; mu_zz +/-5%'}
calibration['best_inplane_2pct'] = {
    'reference': {'states': bi['states'], 'strain_voigt': 'e6 = +/-0.04', 'mu_ref': round(bi['mu_ref'], 2),
                  'mu_values': {k: round(v, 2) for k, v in bi['mu_values'].items()}},
    'runner_up': {k: round(v, 2) for k, v in sorted(H['two_pct_ranking']['inplane_max'].items(), key=lambda x: -x[1])[2:4]},
    'nearest_ablation': {'diagonal_only': V['ablations']['best_inplane_diagonal_only'],
                         'label_read_as_voigt': {'strain_voigt': V['ablations']['best_inplane_label_as_voigt']['strain_voigt']}},
    'tolerance': 'strain: dominant component |e6| = 0.04 +/- 0.002 (either sign), others |e| <= 0.002; mu +/-8%'}

design = {
    'task_id': 'gan-hole-piezomobility', 'catalog_id': 'gan-hole-piezomobility', 'start': 'released',
    'source_group': 'chen_2025_gan_piezomobility',
    'title': 'Wurtzite GaN: strain engineering of the hole mobility tensor',
    'decision_points': [
        {'id': 'dp1', 'fork': 'How to obtain a withheld shear state (expand_3%_e4, compress_3%_e6)',
         'wrong_choice': 'Extrapolate the same-sign 1%/2% series (or a linear tensor) instead of using the exact mirror '
                         'partner; or copy the partner without changing the sign of the strained off-diagonal element',
         'consequence': 'expand_3%_e4 mu_yy +19%, mu_zz +14%, mu_yz +21%; unflipped off-diagonals have the wrong sign',
         'ablation': 'W4_same_sign_series', 'difficulty_type': 'symmetry that must be recognised and verified; '
         'the exact partner is purchasable but not labelled as such'},
        {'id': 'dp2', 'fork': 'Whether e2 can be obtained from e1 by swapping x and y',
         'wrong_choice': 'Use the transverse-isotropy (linear-order) swap of the purchasable expand_3%_e1 state',
         'consequence': 'expand_3%_e2 mu_xx -18%, mu_yy +17%: uniaxial strain along a1 and perpendicular to it are '
                        'inequivalent at finite strain (no 90-degree rotation in 6mm)',
         'ablation': 'W3_swap_from_e1', 'difficulty_type': 'textbook symmetry that fails in the nonlinear regime; '
         'discovered only by buying both e1 and e2 at one strain'},
        {'id': 'dp3', 'fork': 'Model class: linear piezomobility tensor vs nonlinear, one-sided response',
         'wrong_choice': 'Fit a linear (symmetry-constrained) piezomobility tensor to the +-1% states and apply it at 3%',
         'consequence': 'All targets fail (expand_3%_e3 mu_zz -48%, expand_3%_e2 mu_yy negative, shear diagonals '
                        'unchanged): HH/LH are quasi-degenerate, so the response is non-analytic at zero strain',
         'ablation': 'linear_piezomobility_from_1pct', 'difficulty_type': 'regime judgement'},
        {'id': 'dp4', 'fork': 'Which regime compress_3%_e3 is in (c-axis band crossing)',
         'wrong_choice': 'Extrapolate the modest 0->1% c-axis response',
         'consequence': 'mu_zz ~94 cm2/Vs instead of ~214: the crystal-field level crosses to the top between 1% and '
                        '2%, visible in the purchasable GW bands and in the 2% mobility',
         'ablation': 'W1_precrossing_extrapolation', 'difficulty_type': 'regime judgement; discriminating purchase'},
        {'id': 'dp5', 'fork': 'SERTA (1 unit) as a cheap survey vs IBTE (3 units)',
         'wrong_choice': 'Rank the 2% states by SERTA mu_zz and rescale by the unstrained IBTE/SERTA ratio',
         'consequence': 'SERTA ranks compress_2%_e4 (71.0) above compress_2%_e3 (64.1); IBTE reverses the order '
                        '(148.0 vs 163.9) because the IBTE/SERTA ratio depends on band character (2.1 vs 2.6)',
         'ablation': 'best_mu_zz_from_serta', 'difficulty_type': 'budget allocation; consequence visible only after '
         'buying IBTE for competing states'},
        {'id': 'dp6', 'fork': 'Diagonal elements vs principal values for the in-plane maximum',
         'wrong_choice': 'Compare mu_xx and mu_yy only',
         'consequence': 'Picks compress_2%_e1 (158 cm2/Vs) instead of in-plane shear e6 (~205 cm2/Vs along ~52 deg)',
         'ablation': 'best_inplane_diagonal_only', 'difficulty_type': 'off-diagonal elements / principal axes'},
        {'id': 'dp7', 'fork': 'Meaning of the archive shear labels',
         'wrong_choice': "Read 'compress_2%_e6' as Voigt e6 = -0.02",
         'consequence': 'The archive applies eps_xy = -0.02 (visible in the cell), i.e. Voigt e6 = -0.04; the reported '
                        'strain is off by a factor 2',
         'ablation': 'best_inplane_label_as_voigt', 'difficulty_type': 'convention discovered in metadata (cells)'},
        {'id': 'dp8', 'fork': 'Whether IBTE runs flagged "reached the maximum but did not converge" are usable',
         'wrong_choice': 'Distrust the flagged runs (15 archived, 12 purchasable: e3 expansion and most e4/e5 states) '
                         'and substitute SERTA scaled by the unstrained IBTE/SERTA ratio',
         'consequence': 'expand_3%_e3 mu_xx +17%, expand_3%_e4 mu_yy +47%; the flagged runs change by <= 3e-4 '
                        'cm2/Vs per iteration at the end and drift by < 0.06% (300 K diagonal) after iteration 100',
         'ablation': 'flagged_runs_replaced_by_scaled_serta', 'difficulty_type': 'data-quality judgement from raw solver logs'},
    ],
    'graded_fields': [
        'mobility_predictions_cm2_per_Vs: compress_3%_e3 mu_zz in [150, 280]',
        'expand_3%_e3 xx, yy, zz; expand_3%_e2 xx, yy; expand_3%_e4 xx, yy, zz; compress_3%_e6 xx, yy, zz within +/-10%',
        'expand_3%_e4 yz and compress_3%_e6 xy within +/-15% (signed)',
        'best_mu_zz_2pct: Voigt strain e3 = -0.02 and mu_zz within 5%',
        'best_inplane_2pct: Voigt strain e6 = +/-0.04 and in-plane principal mobility within 8%',
    ],
    'length_estimate': {
        'oracle_lab_evaluations': ORACLE.get('evaluations', 18), 'oracle_charge': ORACLE.get('charge_used', 44),
        'budget': 50, 'catalog_charge_all_bands_and_ibte': 37 + 3 * 32, 'tool_calls': '40-70',
        'gpt6astra_cost_usd': '2-4', 'oracle_runs': {k: ORACLE[k] for k in ('gate', 'harbor') if k in ORACLE},
        'pilot': 'pending (coordinator)'},
    'runtime_dependencies': ['numpy (runtime image); the lab is a lookup of archived values'],
    'data_provenance': {'archive_url': 'https://archive.materialscloud.org/records/qfdfa-cs517',
                        'doi': '10.24435/materialscloud:zy-qw',
                        'files_used': ['strain_GaN.tar.gz (scf.in, epw2.in, gan.eig, epw2.out of 37 states; see '
                                       'members.md5)', 'README.txt'],
                        'md5': {'strain_GaN.tar.gz': '9a8747242f4759d235c5552aaea0ba52',
                                'README.txt': 'e72a76b7c622cd39d0d3a91ea75811e6'},
                        'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': ('Adapted from the design: the linear piezomobility coefficients are not graded (the response is '
              'non-analytic at zero strain, so a linear coefficient is not uniquely defined by the data); the withheld '
              'compress_2%_e2 is replaced by expand_3%_e2 (compress_2%_e2 is not recoverable within 10% by any '
              'defensible route); compress_3%_e3 is graded with a regime window because its exact value is not '
              'recoverable from the purchasable data; the GW-vs-DFT fork is dropped (only GW energies are archived) and '
              'replaced by the shear-label convention, the SERTA survey trap and the unconverged-flag judgement.'),
    'calibration': calibration,
}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print(json.dumps(calibration, indent=1)[:3000])
