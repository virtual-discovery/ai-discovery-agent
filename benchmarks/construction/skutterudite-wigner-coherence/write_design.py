"""Write design.json (decision points, calibration table, length estimate) from ablation_values.json.

    python3 -B write_design.py      (after make_task.py and the gate; reads check.json if present)
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/skutterudite-wigner-coherence'
V = json.loads((TASK/'solution/ablation_values.json').read_text())
H = json.loads((TASK/'environment/lab/private/hidden.json').read_text())
CHECK = json.loads((HERE/'check.json').read_text()) if (HERE/'check.json').exists() else {}
REF, A, R, T = V['reference'], V['ablations'], V['robustness'], H['thresholds']
FIELDS = {'YbIr4Sb12 kappa_P 300 K': ('YbIr4Sb12', 'kappa_P_W_mK', '300'), 'YbIr4Sb12 kappa_P 700 K': ('YbIr4Sb12', 'kappa_P_W_mK', '700'),
          'YbIr4Sb12 kappa_C 100 K': ('YbIr4Sb12', 'kappa_C_W_mK', '100'), 'YbIr4Sb12 kappa_C 300 K': ('YbIr4Sb12', 'kappa_C_W_mK', '300'),
          'YbIr4Sb12 kappa_C 700 K': ('YbIr4Sb12', 'kappa_C_W_mK', '700'), 'YbIr4Sb12 T_half (K)': ('YbIr4Sb12', 'T_half_K'),
          'YbIr4Sb12 rattler frequency (THz)': ('YbIr4Sb12', 'rattler_frequency_THz'),
          'IrSb3 kappa_P 50 K': ('IrSb3', 'kappa_P_50K_W_mK'), 'IrSb3 kappa_TOT 300 K': ('IrSb3', 'kappa_TOT_300K_W_mK')}
TOL = {'YbIr4Sb12 kappa_P 300 K': 'within %d%%' % round(100 * T['kappa_P_rel']), 'YbIr4Sb12 kappa_P 700 K': 'within %d%%' % round(100 * T['kappa_P_rel']),
       'YbIr4Sb12 kappa_C 100 K': 'within %d%%' % round(100 * T['kappa_C_rel']), 'YbIr4Sb12 kappa_C 300 K': 'within %d%%' % round(100 * T['kappa_C_rel']),
       'YbIr4Sb12 kappa_C 700 K': 'within %d%%' % round(100 * T['kappa_C_rel']), 'YbIr4Sb12 T_half (K)': '+-%d K' % T['T_half_K'],
       'YbIr4Sb12 rattler frequency (THz)': 'within %d%%' % round(100 * T['rattler_rel']),
       'IrSb3 kappa_P 50 K': 'within %d%%' % round(100 * T['irsb3_kappa_P_50K_rel']), 'IrSb3 kappa_TOT 300 K': 'within %d%%' % round(100 * T['irsb3_kappa_TOT_300K_rel'])}


def get(d, path):
    for k in path:
        d = d[k]
    return d


def fmt(x, ref, field):
    return round(x, 1) if 'T_half' in field else round(x, 4)


def dev(x, ref, field):
    return ('%+.1f K' % (x - ref)) if 'T_half' in field else ('%+.1f%%' % (100 * (x / ref - 1)))


calibration = {}
for field, path in FIELDS.items():
    ref = get(REF, path)
    abl = {k: '%s (%s)' % (fmt(get(v, path), ref, field), dev(get(v, path), ref, field)) for k, v in A.items()
           if abs(get(v, path) / ref - 1) > 1e-3}
    rob = {k: '%s (%s)' % (fmt(get(v, path), ref, field), dev(get(v, path), ref, field)) for k, v in R.items()}
    calibration[field] = {'reference': fmt(ref, ref, field), 'robustness': rob, 'ablations_that_move_it': abl, 'tolerance': TOL[field]}

oracle = CHECK.get('summary', {}).get('oracle', {})
design = {
    'task_id': 'skutterudite-wigner-coherence', 'catalog_id': 'skutterudite-wigner-coherence', 'start': 'released',
    'source_group': 'dilucente_2024_skutterudite_wigner',
    'title': 'Yb-filled skutterudites: Boltzmann-to-Wigner crossover',
    'decision_points': [
        {'id': 'dp1', 'fork': 'Whether the archived isotope rates of the Yb-filled runs are physical',
         'wrong_choice': 'Use the archived gamma_isotope as is (it reproduces the archived YbCo4Sb12 table exactly)',
         'consequence': '21 entries of 1.7e3-2.0e7 THz at the q-point nearest Gamma kill the long-wavelength acoustic modes: '
                        'kappa_P of YbIr4Sb12 11% low, T_half 58 K low',
         'evidence': 'gamma_isotope >> frequency only at q-index 1; all other entries < 0.01 THz', 'ablation': 'archived_isotope_rates_as_is'},
        {'id': 'dp2', 'fork': 'Whether isotope scattering must be added to the archived (ph-ph) gamma',
         'wrong_choice': 'Use gamma alone (or conclude from the byte-identical 1.19 "isotope" files that isotopes do nothing)',
         'consequence': 'Passes for YbIr4Sb12 (0.6%) but IrSb3 kappa_P at 50 K 6.3% high',
         'evidence': 'FeSb3/YbCo4Sb12 logs reproduce only with gamma + gamma_isotope; the 1.19 isotope kappa files equal the no-isotope ones while their logs differ',
         'ablation': 'isotope_scattering_omitted'},
        {'id': 'dp3', 'fork': 'Which archived CoSb3 Wigner table validates the lifetime convention',
         'wrong_choice': 'Calibrate on tc_data.dat (kP twice the phono3py log) and take tau = 1/(2 pi gamma)',
         'consequence': 'All kappa_P doubled, Lorentzian widths halved (kappa_C 10-20% low), T_half beyond 800 K',
         'evidence': 'phono3py log in the same folder; the script that wrote tc_data.dat uses a 0.5 populations prefactor and omega_1 for omega_2',
         'ablation': 'lifetime_convention_from_tc_data'},
        {'id': 'dp4', 'fork': 'Where velocity-operator elements within degenerate subspaces belong',
         'wrong_choice': 'Count every off-diagonal element as a coherence',
         'consequence': 'kappa_C of YbIr4Sb12 at 100 K 10% high (2.5% at 300 K)',
         'evidence': 'archived kappa_C of the released runs reproduces only with degenerate pairs excluded',
         'ablation': 'degenerate_pairs_counted_as_coherences'},
        {'id': 'dp5', 'fork': 'Band-overlap function and width in the coherence term',
         'wrong_choice': 'Gaussian overlap, or a Lorentzian with the phono3py gamma as full width',
         'consequence': 'kappa_C 15-28% (Gaussian) or 16-20% (half width) low; T_half 47-66 K high',
         'evidence': 'archived kappa_C of released runs', 'ablation': 'gaussian_band_overlap'},
        {'id': 'dp6', 'fork': 'How to repair the corrupted isotope entries',
         'wrong_choice': 'Discard the affected q-point from the Brillouin-zone sum',
         'consequence': 'kappa_P of YbIr4Sb12 15% low, T_half 86 K low (the q-point carries the long-wavelength acoustic heat)',
         'evidence': 'the ph-ph linewidths and velocities at that q-point are normal; only gamma_isotope failed',
         'ablation': 'corrupted_qpoint_discarded'},
        {'id': 'dp7', 'fork': 'Units of the archived projected-DOS frequency axis',
         'wrong_choice': 'Read the projected_dos.dat axis as THz',
         'consequence': 'Rattler frequency 6.75 THz instead of 0.97 THz (factor 6.97); the rattler band appears in the optical range',
         'evidence': 'PDOS extends to ~45 on its axis while modes stop at 6.4 THz; the PDOS run used phonopy --qe (factor 108.97) '
                     'on eV/A^2 force constants (phono3py factor 15.633); the Yb MSD Einstein frequency is 0.98 THz',
         'ablation': 'pdos_axis_taken_as_THz'},
    ],
    'graded_fields': ['%s (%s)' % (k, TOL[k]) for k in FIELDS],
    'length_estimate': {
        'oracle_lab_evaluations': oracle.get('evaluations'), 'oracle_charge': oracle.get('charge_used'), 'budget': oracle.get('budget'),
        'catalog_units': '10 runs x 16 temperatures of mode_data (160) + 6 velocity operators (24) + 11 archived files + 7 harmonic = 202',
        'tool_calls': '45-75', 'gpt6astra_cost_usd': '2.5-4.5',
        'why': 'Parsing large per-mode JSON, implementing the Wigner populations/coherences formula with degenerate handling, '
               'validating it against released logs, resolving conflicting archived tables, detecting corrupted isotope rates, '
               'and locating the crossover from several temperatures under a budget of about 15% of the catalog.',
        'pilot': 'pending (coordinator)'},
    'runtime_dependencies': ['numpy (runtime image)'],
    'data_provenance': {'archive_url': 'https://archive.materialscloud.org/records/5jsfj-cdg70', 'doi': '10.24435/materialscloud:xn-qr',
                        'files_used': 'materials_cloud_archive_submission.zip, 87 members read by HTTP range (list and sha256 in construction/members.json)',
                        'md5': {'materials_cloud_archive_submission.zip': '9dc2bd3b8659022e5703862f82adaa40'}, 'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': ('Graded values come from the archived inputs of the withheld runs with the authors\' formula (reproduced to 1e-6 on every '
              'released 2.1.0 run); the only change is that corrupted gamma_isotope entries (> 1 THz, all at the q-point nearest Gamma) '
              'are set to zero. Any repair that keeps the q-point (zero, band median, omega^4 fit) agrees within 0.4%. YbFe4Sb12 is '
              'not graded: its two archived fc3 variants differ by 2.4x in kappa_P and neither uses the paper\'s DFPT-fc2 + hiPhive-fc3 '
              'combination. The design\'s velocity-operator transfer to the hidden compound was dropped (not uniquely recoverable).'),
    'calibration': calibration,
}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print(json.dumps(calibration, indent=1))
