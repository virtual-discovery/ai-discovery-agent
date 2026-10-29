"""DeePMD reference calculations through the lab's own code path (task_lab.run_structure), runtime image.

    python3 -B reference.py      -> ref_dp.json

1. For every near-equilibrium frame of the metastable phases (ST12, BC8, R8) and of CD and bSn (see
   common.clean_frames), the DP energy of the archived geometry and after a fixed-cell ionic relaxation
   ("relax_ions"). Their difference is the internal-relaxation correction the archived DFT data lack.
2. DP cell + ion relaxations ("relax_cell") of every phase at -2 ... 12 GPa (for the surrogate-only ablations).
"""
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
LAB = HERE.parents[1]/'tasks/ge-pressure-polymorph-pathways/environment/lab'
sys.path.insert(0, str(LAB))
sys.path.insert(0, str(HERE))
import task_lab  # noqa: E402
from common import POOL, PHASES, clean_frames  # noqa: E402

START = {'CD': 'training_set/CD/config_98', 'bSn': 'training_set/bSn/config_61', 'ST12': 'training_set/ST12/config_1',
         'BC8': 'training_set/BC8/config_61', 'R8': 'training_set/R8/config_153'}
PRESSURES = [-2, 0, 2, 4, 6, 8, 10, 12]


def main():
    t0 = time.time()
    corr = {}
    for g in PHASES:
        for o in clean_frames(g):
            fr = POOL[o['key']]['frames'][o['j']]
            r = task_lab.run_structure({'cell': fr['cell'], 'positions': fr['positions'], 'mode': 'relax_ions'})
            corr[o['tag']] = {'dp_unrelaxed_eV': r['energy_initial_eV'], 'dp_relaxed_eV': r['energy_eV'], 'natoms': o['n'],
                              'max_force_after_eV_A': r['max_force_eV_A']}
        print(g, 'corrections done', round(time.time() - t0), 's', flush=True)
    curves = {}
    for g, key in START.items():
        fr = POOL[key]['frames'][-1]
        cell, pos = fr['cell'], fr['positions']
        pts = []
        for p in PRESSURES:
            r = task_lab.run_structure({'cell': cell, 'positions': pos, 'mode': 'relax_cell', 'pressure_GPa': p})
            pts.append({'P_GPa': p, 'H_eV_atom': r['enthalpy_eV'] / r['natoms'], 'V_A3_atom': r['volume_A3'] / r['natoms'],
                        'E_eV_atom': r['energy_eV'] / r['natoms'], 'max_force': r['max_force_eV_A'], 'stress_kB': r['stress_kB']})
            cell, pos = r['cell_A'], r['positions_A']
        curves[g] = pts
        print(g, 'curve done', round(time.time() - t0), 's', flush=True)
    (HERE/'ref_dp.json').write_text(json.dumps({'relax_corrections': corr, 'dp_relaxed_curves': curves}, indent=0) + '\n')


if __name__ == '__main__':
    main()
