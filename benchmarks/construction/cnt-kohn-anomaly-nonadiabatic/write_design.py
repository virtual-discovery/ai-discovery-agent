"""Write design.json for cnt-kohn-anomaly-nonadiabatic from calibration.json (and check.json when present).

    python3 -B write_design.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/cnt-kohn-anomaly-nonadiabatic'
CAL = json.loads((HERE/'calibration.json').read_text())
CHECK = json.loads((HERE/'check.json').read_text()) if (HERE/'check.json').exists() else None
REF, ROB, ABL, FACTS = CAL['reference'], CAL['robustness'], CAL['ablations'], CAL['facts']


def r(x, n=4):
    return None if x is None else float('%.*g' % (n, x))


DPS = [
    {'id': 'dp1', 'fork': 'Does the archived DFPT frequency omega_Q already contain an electronic self-energy?',
     'wrong_choice': 'Take omega_b = omega_Q and add Pi(T) on top (double counting)',
     'consequence': 'nu = 20 softens to ~100 meV and the RBM is unstable at every archived temperature',
     'ablation': 'no_dfpt_subtraction', 'kind': 'physics, supported by the archived DFPT smearing scan'},
    {'id': 'dp2', 'fork': 'Which static self-energy is contained in omega_Q: the one at the DFPT smearing (0.01 Ry = 1578.87 K), not at the physical T or at 0.01 eV',
     'wrong_choice': 'Subtract the static self-energy at the same temperature T (or read degauss in eV)',
     'consequence': 'nu = 20 hardens above omega_Q (Omega ~ 145 meV) and the RBM never becomes unstable',
     'ablation': 'same_temperature_reference', 'kind': 'metadata (degauss units, smearing) and released 1578.87 K table'},
    {'id': 'dp3', 'fork': 'Static (adiabatic) versus dynamic (nonadiabatic) self-energy for the mode-20 frequency',
     'wrong_choice': 'Use the adiabatic static result for the nonadiabatic fields (or the converse)',
     'consequence': 'Omega(100 K) off by ~19 meV, Omega(300 K) by ~7 meV',
     'ablation': 'adiabatic_for_nonadiabatic', 'kind': 'model class'},
    {'id': 'dp4', 'fork': 'Column identity of the archived self-energy tables: the header names 10 columns in an order that does not match the 8 data columns',
     'wrong_choice': "Follow the header: take column 8 ('Real Sigma(w=0)') as the static self-energy",
     'consequence': 'The subtracted static part is ~0 (column 8 is Im Pi, which vanishes at w = 0): equivalent to double counting',
     'ablation': 'header_column_order', 'kind': 'archive defect visible only in the raw data'},
    {'id': 'dp5', 'fork': 'Scattering processes: k -> k+Q alone covers the zone; the authors\' coupling-plot script sums k+q and k-q',
     'wrong_choice': 'Sum both processes',
     'consequence': 'Self-energy doubled: Omega(300 K) 7 meV low, linewidth doubled, RBM T* ~ 910 K',
     'ablation': 'both_processes_summed', 'kind': 'bookkeeping; discriminated only by validation against an archived table'},
    {'id': 'dp6', 'fork': 'Band pairs in the energy-ordered Wannier basis: the Dirac crossing at E_F makes the 2k_F backscattering the 1<->2 pair',
     'wrong_choice': 'Keep only the intraband pairs (1,1) and (2,2)',
     'consequence': 'No Kohn anomaly: Omega ~ omega_Q, zero linewidth, no RBM instability',
     'ablation': 'sorted_intraband_only', 'kind': 'bookkeeping found in the coupling data'},
    {'id': 'dp7', 'fork': 'Convergence of the dynamic self-energy at low T (eta = k_B T/100): validation at high T does not transfer',
     'wrong_choice': 'Coarse sampling (nk = 1e5) at 100 K, or the raw 2000-point archived grid',
     'consequence': 'Linewidth at 100 K 10% low (1e5); raw grid: Omega 6-7 meV high and the linewidth almost zero',
     'ablation': 'coarse_mesh_1e5', 'kind': 'convergence regime; needs a low-temperature validation (released 50 K table)'},
    {'id': 'dp8', 'fork': 'Reference frequency and smearing scheme: omega_Q (cold-smearing DFPT behind the coupling data) versus the Fermi-Dirac DFPT smearing scan',
     'wrong_choice': 'Anchor on the Fermi-Dirac 10 mRy frequency, or use the scan with sigma as temperature for T*',
     'consequence': 'RBM T* ~ 330 K or ~ 204 K instead of ~ 527 K; nu = 20 frequencies 8 meV high',
     'ablation': 'fd_scan_reference_frequency', 'kind': 'inconsistent archived calculations (smearing scheme)'},
    {'id': 'dp9', 'fork': 'Mass parameter of the archived threshold program: the variable called mass is a = 1/m*',
     'wrong_choice': 'Pass m* = 2 as the program parameter',
     'consequence': 'g0 doubled (121 and 305 meV)', 'ablation': 'inverted_mass_parameter',
     'kind': 'archive naming inconsistency (main() mass_params, function doc a = 1/m*, fig1.py reversed labels)'},
    {'id': 'dp10', 'fork': 'Omega(T) of nu = 20 is non-monotonic (minimum near omega_Q/4.6 k_B): withheld temperatures must be computed',
     'wrong_choice': 'Interpolate the released tables (50 K and 700 K)',
     'consequence': 'Omega(300 K) 1.9 meV high; RBM T* ~ 370 K', 'ablation': 'interpolate_released_tables',
     'kind': 'regime judgement under a binding budget'},
]

def field_values(group, path):
    out = {}
    for name, d in group.items():
        v = d
        for p in path:
            v = v.get(p) if isinstance(v, dict) else None
            if v is None:
                break
        if v is not None or (path[0] in d and len(path) == 1):
            out[name] = r(v)
    return out


def main():
    tol = {'omega_nonadiabatic_meV': '+-1.2 meV', 'linewidth_meV': '+-6%', 'omega_adiabatic_100K_meV': '+-1.5 meV',
           'rbm_instability_temperature_K': '+-8%', 'threshold_coupling_meV': '+-6%'}
    paths = [('omega_nonadiabatic_meV', '100'), ('omega_nonadiabatic_meV', '300'), ('linewidth_meV', '100'),
             ('linewidth_meV', '300'), ('omega_adiabatic_100K_meV',), ('rbm_instability_temperature_K',),
             ('threshold_coupling_meV', '1'), ('threshold_coupling_meV', '20')]
    calibration = {}
    for p in paths:
        ref = REF
        for q in p:
            ref = ref[q]
        calibration['/'.join(p)] = {'reference': r(ref, 6), 'robustness': field_values(ROB, p),
                                    'ablations': field_values(ABL, p), 'tolerance': tol[p[0]]}
    oracle = CHECK['summary']['oracle'] if CHECK else {}
    design = {
        'task_id': 'cnt-kohn-anomaly-nonadiabatic', 'catalog_id': 'cnt-kohn-anomaly-nonadiabatic', 'start': 'released',
        'source_group': 'marazzi_2025_kohn_1d',
        'title': 'Kohn anomaly in a 4 A carbon nanotube: nonadiabatic softening and the RBM instability',
        'decision_points': DPS,
        'graded_fields': ['omega_nonadiabatic_meV {100, 300} (nu = 20) within 1.2 meV',
                          'linewidth_meV {100, 300} (nu = 20) within 6%',
                          'omega_adiabatic_100K_meV (nu = 20) within 1.5 meV',
                          'rbm_instability_temperature_K within 8% (null fails)',
                          'threshold_coupling_meV {1, 20} within 6%'],
        'length_estimate': {'oracle_lab_evaluations': oracle.get('evaluations'), 'oracle_charge': oracle.get('charge_used'),
                            'budget': oracle.get('budget'), 'tool_calls': '45-75', 'gpt6astra_cost_usd': '2-4',
                            'why': ('validation of the evaluator at low and high T, self-energies at two withheld '
                                    'temperatures, a temperature bracket for T*, the reference smearing, the threshold '
                                    'program, and archive inspection')},
        'runtime_dependencies': ['numpy and scipy (runtime image); gcc (runtime image) only if the agent compiles the archived C model'],
        'data_provenance': {'archive_url': 'https://archive.materialscloud.org/records/vs5fd-dbv98',
                            'doi': '10.24435/materialscloud:yx-x1', 'license': 'CC BY 4.0',
                            'files_used': ['README.md', 'dft_inputs.zip', 'Fig1.zip', 'Fig2.zip'],
                            'md5': {'README.md': '1c6dbfd0eb10564470d215e3eca226fe', 'dft_inputs.zip': 'c85fa7a8af0ffa84597e7fe76e1749b5',
                                    'Fig1.zip': '4c1760d1aaf80add6900cd65233ca571', 'Fig2.zip': 'b17a42c46d5a7856356d71eab691875c'}},
        'exceptions': {},
        'notes': ('Truth: Eq. (6)-(7) complex pole on the archived self-energy tables at the withheld temperatures, with the '
                  'static part of the released 1578.87 K table subtracted (the authors\' protocol, fig2_d_f_h.py); T* from '
                  'the Omega -> 0 limit of the pole across the archived 500-600 K tables; g0 from the archived '
                  'coupling_calculator.c rebuilt with a parametrised main(). The lab evaluator reproduces the archived '
                  f"tables to {FACTS['evaluator_max_deviation_meV']:.2f} meV. Paper-text thresholds (54.5 and 129.0 meV) "
                  'disagree with the archived program and with Eq. (9) and fail the threshold check.'),
        'calibration': calibration,
    }
    if CHECK:
        design['gate'] = {'passed': CHECK['passed'], 'checked_utc': CHECK['checked_utc']}
    (TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
    print(json.dumps(calibration, indent=1))


if __name__ == '__main__':
    main()
