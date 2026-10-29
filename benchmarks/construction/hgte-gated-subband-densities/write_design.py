"""Write the task's design.json from calibration.json (and oracle_run.json when present).

    python3 -B write_design.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/hgte-gated-subband-densities'
calib = json.loads((HERE/'calibration.json').read_text())
oracle = json.loads((HERE/'oracle_run.json').read_text()) if (HERE/'oracle_run.json').exists() else {}

design = {
    'task_id': 'hgte-gated-subband-densities',
    'catalog_id': 'hgte-gated-subband-densities',
    'start': 'released',
    'source_group': 'hofer_2026_hgte_selfconsistent',
    'title': 'Gated HgTe TI films: SdH subband densities vs self-consistent k.p',
    'decision_points': [
        {'id': 'dp1', 'fork': 'SdH frequency-to-density conversion for bulk subbands and hole pockets',
         'wrong_choice': 'n = 2 e f / h (spin-degenerate) for the bulk subband and the hole pocket',
         'consequence': 'first bulk subband 2.45 instead of 1.22, hole pocket 4.98 instead of 2.49; charge balance violated',
         'ablation': 'bulk_spin_degenerate'},
        {'id': 'dp2', 'fork': 'Which SdH feature is the majority hole pocket at -0.70 V',
         'wrong_choice': 'Strongest peak above 1e11 cm^-2 at that gate voltage (a harmonic / sum frequency), or the peak nearest the model',
         'consequence': '1.17 or 2.15 instead of 2.49; only ridge tracking across gate voltages finds it',
         'ablation': 'hole_strongest_peak'},
        {'id': 'dp3', 'fork': 'Top vs bottom surface on the hole side',
         'wrong_choice': 'Top (gate-side) surface = the larger surface density, as on the electron side; kdotpy band 1 = top everywhere',
         'consequence': 'top/bottom swapped at -0.45 V (0.63 vs 0.45) and model top surface 0.68 instead of 0.48',
         'ablation': 'top_is_larger_density'},
        {'id': 'dp4', 'fork': 'Multiplicity of the majority hole pocket',
         'wrong_choice': 'Four equivalent side maxima (cubic symmetry) or a single pocket',
         'consequence': 'charge balance p_total = 4.6 = 1.85 x 2.49 fails; BIA splits the side maxima into 2 + 2',
         'ablation': 'four_equivalent_pockets'},
        {'id': 'dp5', 'fork': 'Fermi level at which the wedge (pocket) runs are read',
         'wrong_choice': "The wedge run's own E_F(n) from its .kpout",
         'consequence': 'majority pocket at -0.60 V 1.26 instead of 1.60 (the wedge E_F assumes an isotropic extrapolation)',
         'ablation': 'wedge_own_fermi_level'},
        {'id': 'dp6', 'fork': 'Normalisation of the wedge IDOS / pocket counting in the main run',
         'wrong_choice': 'Main-run band -1 divided by 2 or 4, or the wedge IDOS without dividing by 360/width',
         'consequence': '1.86, 0.93 or 4.81 instead of 1.60',
         'ablation': 'band_minus1_split_in_two'},
        {'id': 'dp7', 'fork': 'Self-consistent runs at each density vs rigid-band filling of one run',
         'wrong_choice': 'Fill the n = 0 run rigidly to n_tot(+0.20 V)',
         'consequence': 'model bottom surface 1.59 instead of 0.97: no screening of the bottom surface',
         'ablation': 'rigid_band_from_neutral_run'},
        {'id': 'dp8', 'fork': 'Where bulk hole pockets start to fill',
         'wrong_choice': 'At charge neutrality of the gate line, or at the Hall-slope sign change',
         'consequence': '-0.33 / -0.35 V instead of -0.21 V: holes coexist with surface electrons on the n side',
         'ablation': 'holes_appear_at_neutrality'},
        {'id': 'dp9', 'fork': 'Top surface vs first bulk subband on the electron side',
         'wrong_choice': 'Take the steep ridge (2.8 at +0.20 V, 3.2 at +0.30 V) as the bulk subband and a flat ~1.4-1.7 feature '
                         '(2x bottom-surface harmonic or a sum frequency) as the top surface; the Hall sum rule is still met',
         'consequence': 'top surface 1.42/1.72 instead of 2.81, first bulk subband 3.25 instead of 1.22. Ruled out by the data: '
                        'the steep ridge is already at 1.6-2.0 at Vg = -0.05..0 V, while a newly filled subband starts at zero '
                        '(the 0.34 -> 1.74 ridge from +0.05 V extrapolates to zero near -0.07 V); the ~1.5 feature is flat at '
                        '2x the bottom surface; the model band continuous across neutrality is band 1 (2.76), bands 3/4 start from zero',
         'ablation': 'steep_ridge_as_bulk_harmonic_as_top'},
    ],
    'graded_fields': [
        'experiment_1e11_cm2: six 70 nm SdH densities, each within max(6 %, 0.07)',
        'model_1e11_cm2: three densities of the released 70 nm runs, each within max(6 %, 0.04)',
        'hole_onset_gate_V within 0.06 V',
        'majority_pocket_multiplicity exact (2)',
    ],
    'length_estimate': {
        'oracle_lab_evaluations': oracle.get('evaluations'),
        'oracle_charge': oracle.get('charge_used'),
        'budget': 25,
        'catalog_units': 162,
        'tool_calls': '50-80',
        'gpt6astra_cost_usd': '3-5',
        'why': 'Raw-data SdH pipeline (sigma_xx, 1/B resampling, FFT), assignment of 6 ridges using purchased model runs, '
               'ridge tracking of the hole pocket over 5 gate blocks, kdotpy output parsing incl. wedge runs, gate calibration; '
               'budget allows ~15 % of the catalog.',
    },
    'runtime_dependencies': ['numpy, scipy (runtime image)'],
    'data_provenance': {
        'archive_url': 'https://zenodo.org/records/18709991 ; https://zenodo.org/records/18429111',
        'files_used': ['data_repository.zip (members, range reads)', 'data_and_code.zip (members, range reads)'],
        'md5': {'data_repository.zip': 'a948d21c70f3f911d3bbc8554acccbce', 'data_and_code.zip': 'e777270ff41037b2072b819a8d4b9961',
                'members': 'benchmarks/construction/hgte-gated-subband-densities/members.md5'},
        'license': 'CC BY 4.0 (both records)',
    },
    'exceptions': {},
    'notes': ('Adapted design: no kdotpy runs on request (GPU/CPU-days; the archive holds one parameter set), so the '
              'dielectric/VBO calibration fork is a rubric item. The graded experimental values are measured 70 nm SdH '
              'densities that the agent derives from the raw record; model values come from the archived runs. '
              '45 nm (CompMat-Bench) is not used.'),
    'calibration': {k: v for k, v in calib.items() if not k.startswith('_')},
}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print('design.json written; oracle', oracle.get('evaluations'))
