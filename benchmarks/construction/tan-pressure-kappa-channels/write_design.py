"""Write the task's design.json from calibration.json and (if present) check.json. python3 -B write_design.py"""
import json
from pathlib import Path

HERE = Path(__file__).parent
TASK = HERE.parents[2]/'benchmarks/tasks/tan-pressure-kappa-channels'
cal = json.loads((HERE/'calibration.json').read_text())
check = json.loads((HERE/'check.json').read_text()) if (HERE/'check.json').exists() else {}
oracle = check.get('summary', {}).get('oracle', {})

DP = [
    ('dp1', 'RTA vs iterative BTE (strong Normal processes)', 'Report the RTA conductivity',
     'kappa about 45% too low at every target', 'abl_RTA'),
    ('dp2', 'Three-phonon only vs adding four-phonon scattering (large a-o gap and acoustic bunching freeze out '
            'three-phonon channels)', 'Omit four-phonon rates (they are withheld at 40/80 GPa)',
     'kappa 3-4x too high; with scalebroad 0.1 the iteration does not converge', 'abl_no_4ph'),
    ('dp3', 'Insulator vs semimetal: include phonon-electron scattering at the withheld 80 GPa',
     'Drop phonon-electron scattering where its table is withheld', 'kappa_80 +31%, resistivity share 0',
     'abl_no_phel'),
    ('dp4', 'Units of BTE.sr_phel_metal (column 2 is 1/tau in ps^-1; column 1 is omega in rad/ps)',
     'Treat the rate column as angular (multiply by 2 pi) when merging with the lab rates',
     'kappa_80 -19%, share +0.14', 'abl_phel_times_2pi'),
    ('dp5', 'Row order of the archived phonon-electron tables (ShengBTE irreducible list, third grid index fastest),'
            ' discoverable only by comparing the archived frequency column with the lab phonons',
     'Map rows onto the catalog irreducible list in its own order', 'rates on the wrong modes: kappa -20%, share +0.17',
     'abl_phel_row_order'),
    ('dp6', 'Calculation settings: the archived CONTROL fixes 16^3 and scalebroad 0.1; the three-phonon conductivity '
            'is strongly broadening-dependent', 'Use the ShengBTE default scalebroad = 1.0',
     'kappa -14 to -17%', 'abl_scalebroad_1'),
    ('dp7', 'Transfer to withheld pressures: transfer rates and solve one BTE vs interpolate conductivities '
            '(Matthiessen-like combination of separately solved results)',
     'Interpolate released all-channel kappa in pressure (and 1/T for 600 K)',
     'kappa_40 -7%, kappa_60(600 K) +51%', 'abl_kappa_level_interpolation'),
    ('dp8', 'Temperature transfer of four-phonon rates to 600 K (60 GPa withheld)',
     'Reuse 300 K rates, or scale classically by (T/300)^2', 'kappa_60(600 K) +88% or +21%',
     'abl_4ph_not_rescaled_to_600K'),
    ('dp9', 'Anisotropy (hexagonal): xx vs zz components', 'Report the trace average',
     'kappa_zz +15%', 'abl_trace_average'),
]
design = {
    'task_id': 'tan-pressure-kappa-channels', 'catalog_id': 'tan-pressure-kappa-channels', 'start': 'released',
    'source_group': 'kundu_2024_tan_kappa',
    'title': 'theta-TaN under pressure: competing phonon scattering channels',
    'decision_points': [{'id': i, 'fork': f, 'wrong_choice': w, 'consequence': c, 'ablation': a}
                        for i, f, w, c, a in DP],
    'extra_ablations': {'dp8': 'abl_4ph_classical_T2_scaling'},
    'graded_fields': ['kappa_xx_80GPa_300K_W_mK (+-10%)', 'kappa_zz_80GPa_300K_W_mK (+-10%)',
                      'kappa_xx_40GPa_300K_W_mK (+-5%)', 'kappa_xx_60GPa_600K_W_mK (+-4%)',
                      'phel_resistivity_fraction_80GPa_300K (+-0.06 absolute)'],
    'length_estimate': {'oracle_lab_evaluations': oracle.get('evaluations'), 'oracle_charge': oracle.get('charge_used'),
                        'budget': 30, 'catalog_value_units': 'unbounded (any pressure/mesh/T/broadening/channel '
                        'combination); fixed items alone: 6 pressures x 4 files + 6x7 harmonic + many rates',
                        'tool_calls': '45-75', 'oracle_container_seconds': oracle.get('container_seconds'),
                        'gpt6astra_cost_usd': '3-5'},
    'runtime_dependencies': ['numpy, scipy (runtime image); no extra packages'],
    'data_provenance': {'archive_url': 'https://zenodo.org/records/10564835', 'doi': '10.5281/zenodo.10564835',
                        'files_used': ['IFCs-data-pressure-TaN.tar.gz'],
                        'md5': {'IFCs-data-pressure-TaN.tar.gz': '0e91973a88688ca03ac776d379486f9c'},
                        'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': 'References are new calculations from the archived inputs with the lab code path (tanbte.py); '
             'four-phonon rates by Monte Carlo sampling (8192 pairs per irreducible q; two seeds at 0 and 80 GPa '
             'agree to 0.1% in kappa). The paper\'s kappa(P) figure and the ~60 GPa maximum are not graded. The '
             'dominant defensible spread is the pressure interpolation of the withheld 80 GPa rates (+-6%).',
    'calibration': cal,
}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print('design.json written', oracle)
