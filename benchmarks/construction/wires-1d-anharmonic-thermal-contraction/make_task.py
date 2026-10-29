"""Write the lab data, catalog, private targets and ablation values for wires-1d-anharmonic-thermal-contraction.

    python3 -B make_task.py [POOL]   # default POOL: ~/.cache/ai-theorist/wires-1d-anharmonic-thermal-contraction/pool.json
"""
import json
import sys
from pathlib import Path
import model as M
import calib

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/wires-1d-anharmonic-thermal-contraction'
POOL = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/wires-1d-anharmonic-thermal-contraction/pool.json'
BUDGET = 36
WITHHELD_HESSIANS = [('TaSe3', '650K/10000'), ('TaSe3', '650K/100000'), ('TaSe3', '700K/10000'), ('TaSe3', '700K/100000')]
INFO = {
    'CuC2': {'atoms_per_unit_cell': 3, 'supercell': [1, 1, 6], 'q_star_files': 4,
             'description': 'straight Cu-C-C chain, metallic; DFT cell a = b = 27.0 A (vacuum), c = 4.887 A'},
    'TaSe3': {'atoms_per_unit_cell': 4, 'supercell': [1, 1, 6], 'q_star_files': 4,
              'description': 'TaSe3 single chain (Ta + 3 Se), metallic; DFT cell 28.22 x 27.66 A (vacuum), c = 3.208 A'},
    'AuSe2': {'atoms_per_unit_cell': 6, 'supercell': [1, 1, 10], 'q_star_files': 6,
              'description': 'AuSe2 in its insulating double-cell phase (6 atoms); DFT cell 25.37 x 25.37 A (vacuum), c = 7.301 A'},
}
THRESH = {'k_rel': 0.10, 'alpha_rel': 0.15, 'cv_cuc2_rel': 0.035, 'cv_tase3_rel': 0.015, 'T_inst_abs_K': 75.0}


def main():
    pool = M.load_pool(POOL)
    cal = calib.compute(pool)
    ref, rob, abl = cal['ref'], cal['rob'], cal['abl']
    # ---------------- lab data (served) ----------------
    data = {'runs': {}, 'dfpt': {}, 'dft_files': {}, 'strained': None}
    catalog_runs = {}
    for wire, w in pool['wires'].items():
        data['runs'][wire], catalog_runs[wire] = {}, []
        for key, r in w['runs'].items():
            hess = r['hessian']
            withheld = (wire, key) in WITHHELD_HESSIANS
            rec = {'temperature_K': r['T'], 'configurations': r['configurations'],
                   'log': ({'archive_file': r['log']['file'], 'section_header': r['log']['header'],
                            'populations_in_section': r['log']['populations'],
                            'final_population_excerpt': r['log']['final_excerpt']} if 'log' in r else None),
                   'auxiliary': r['auxiliary'], 'hessian': None if withheld else hess}
            data['runs'][wire][key] = rec
            catalog_runs[wire].append({'temperature_K': r['T'], 'configurations': r['configurations'],
                                       'sscha_log': rec['log'] is not None, 'auxiliary': r['auxiliary'] is not None,
                                       'hessian': 'withheld by this task' if withheld else hess is not None})
        data['dfpt'][wire] = w['dfpt']
        data['dft_files'][wire] = w['dft_files']
    data['strained'] = pool['strained_TaSe3']
    lab = TASK/'environment/lab'
    (lab/'data.json').write_text(json.dumps(data))
    # ---------------- catalog (agent-visible) ----------------
    catalog = {
        'source': 'Author-archived SSCHA and DFT data for exfoliable 1D wires (Cignarella, Bastonero, Monacelli and Marzari, '
                  'Materials Cloud Archive 2025, DOI 10.24435/materialscloud:fj-20). Purchases return archived files or '
                  'verbatim excerpts; nothing is recomputed.',
        'protocol': {
            'sscha': 'Stochastic self-consistent harmonic approximation (python-sscha) with NequIP machine-learned forces, '
                     'at fixed cell (the exfoliated DFT geometry, same cell at every temperature). Each run is labelled by temperature and by the number '
                     'of configurations per population (10000 or 100000).',
            'auxiliary_matrices': 'final SSCHA auxiliary dynamical matrices (archive names final_sscha_*), one file per q-star of the supercell grid',
            'hessian_matrices': 'Hessian of the SSCHA positional free energy (archive names hessian_*), one file per q-star',
            'dfpt': 'zero-temperature harmonic DFPT (PBE, Quantum ESPRESSO 7.2) on the same q grid: "none" = ph.x output '
                    'files as written (incl. the dyn0 q-list), "1d" = the archived matrices with the authors\' 1D acoustic sum rule applied',
            'strained_calculation': 'the only strained calculation in the archive: TaSe3 pw.x + ph.x (1 x 1 x 10 q grid) at one '
                                    'elongated cell, from the authors\' QHA AiiDA archive',
        },
        'file_formats': {
            'dynamical_matrices': 'Quantum ESPRESSO / CellConstructor text: header (ntyp, nat, ibrav, celldm(1) = alat in bohr), '
                                  'cell in units of alat, species masses in Rydberg atomic units (amu x 911.444), positions in '
                                  'units of alat; then for each q (2 pi / alat, Cartesian) the 3 x 3 blocks of the force-constant '
                                  'matrix Phi(q) in Ry/bohr^2 (real, imaginary pairs), Phi(q) = sum_R Phi(0, R) exp(-i q.R).',
            'sscha_log': 'verbatim excerpt of the minimization log: the last minimization step and the final RESULTS block '
                         '(free energy per supercell in meV with stochastic error, gradients, Kong-Liu effective sample size, '
                         'centroid forces, stress tensor and ab initio average stress in GPa) of the requested run. '
                         'Configuration-progress and debug lines are removed.',
        },
        'wires': INFO,
        'sscha_runs': catalog_runs,
        'dfpt_sets': {w: ['none', '1d'] for w in pool['wires']},
        'dft_files': {w: sorted(v) for w, v in data['dft_files'].items()},
        'requests': {
            'sscha_log': {'kind': 'sscha_log', 'wire': '<CuC2|TaSe3|AuSe2>', 'temperature_K': 300, 'configurations': 10000, 'cost': 1},
            'sscha_matrices': {'kind': 'sscha_matrices', 'wire': '<wire>', 'temperature_K': 300, 'configurations': 10000,
                               'matrix': '<auxiliary|hessian>', 'cost': 1},
            'dfpt': {'kind': 'dfpt', 'wire': '<wire>', 'asr': '<none|1d>', 'cost': 1},
            'dft_file': {'kind': 'dft_file', 'wire': '<wire>', 'file': '<name from dft_files>', 'cost': 1},
            'strained_calculation': {'kind': 'strained_calculation', 'cost': 1},
        },
        'budget': BUDGET,
    }
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    # ---------------- private ----------------
    hidden = {'targets': {'k_eV': ref['k'], 'alpha_1e6_per_K': ref['alpha'], 'cv_cuc2': ref['cv_cuc2'],
                          'cv_tase3': ref['cv_tase3'], 'T_inst_K': ref['T_inst']},
              'thresholds': THRESH,
              'notes': {'T_inst_100k_crossing': cal['info']['T_inst_100k'], 'T_inst_10k_crossing': cal['info']['T_inst_10k']}}
    (lab/'private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    (lab/'private/hidden_values.json').write_text(json.dumps(
        [ref['k'], ref['alpha'], ref['cv_cuc2'], ref['cv_tase3'], ref['T_inst'],
         cal['info']['T_inst_100k'], cal['info']['T_inst_10k']]) + '\n')
    (lab/'public.json').write_text(json.dumps({'task_id': TASK.name, 'solver_budget': BUDGET, 'reserved_validation_charge': 0,
                                               'requests': catalog['requests'], 'catalog': '/app/catalog.json'}, indent=1) + '\n')
    (TASK/'tests/public.json').write_text(json.dumps({'task_id': TASK.name, 'solver_budget': BUDGET, 'min_report_chars': 400},
                                                     indent=1) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob, 'reference': ref},
                                                                 indent=1) + '\n')
    (HERE/'calibration.json').write_text(json.dumps(cal, indent=1, default=float) + '\n')
    print(json.dumps({'ref': ref, 'thresholds': THRESH}, indent=1))
    print('lab data MB', round((lab/'data.json').stat().st_size / 2**20, 2))


if __name__ == '__main__':
    main()
