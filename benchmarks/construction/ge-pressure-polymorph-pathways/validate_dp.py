"""Validate the lab's numpy DeePMD evaluator (environment/lab/dp_numpy.py) against the released DFT data.

    python3 -B validate_dp.py     (runtime image; writes dp_validation.json)

Checks: (1) analytic forces and virial against central finite differences; (2) energies, forces and pressures
against the archived VASP results for a random sample of configurations of every folder.
"""
import gzip
import json
import random
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
LAB = HERE.parents[1]/'tasks/ge-pressure-polymorph-pathways/environment/lab'
sys.path.insert(0, str(LAB))
from dp_numpy import DP  # noqa: E402

POOL = json.load(gzip.open(Path.home()/'.cache/ai-theorist/ge-pressure-polymorph-pathways/pool.json.gz'))


def main():
    dp = DP(LAB/'potential')
    out = {'finite_difference': [], 'by_folder': {}}
    for key in ['training_set/BC8/config_10', 'training_set/ST12/config_174', 'training_set/bSn/config_7']:
        f = POOL[key]['frames'][-1]
        cell, pos = np.array(f['cell']), np.array(f['positions'])
        r = dp.compute(cell, pos)
        h = 1e-5
        fd = np.zeros((3, 3))
        for a in range(3):
            for c in range(3):
                p = pos.copy(); p[a, c] += h; ep = dp.energy(cell, p); p[a, c] -= 2 * h; em = dp.energy(cell, p)
                fd[a, c] = -(ep - em) / 2 / h
        vd = np.zeros((3, 3))
        for a in range(3):
            for b in range(3):
                e = np.zeros((3, 3)); e[a, b] = h
                M = np.eye(3) + e; ep = dp.energy(cell @ M.T, pos @ M.T)
                M = np.eye(3) - e; em = dp.energy(cell @ M.T, pos @ M.T)
                vd[a, b] = (ep - em) / 2 / h
        out['finite_difference'].append({'config': key, 'max_force_error_eV_A': float(np.abs(r['forces'][:3] - fd).max()),
                                         'max_virial_error_eV': float(np.abs(r['dE_dstrain'] - vd).max())})
    random.seed(1)
    keys = sorted(k for k, v in POOL.items() if v['frames'] and v['settings']['NIONS'] <= 70)
    for key in random.sample(keys, 300):
        v = POOL[key]
        f = v['frames'][-1]
        n = v['settings']['NIONS']
        r = dp.compute(f['cell'], f['positions'])
        g = 'validation_set' if key.startswith('validation') else key.split('/')[1]
        row = out['by_folder'].setdefault(g, {'dE_meV': [], 'dF': [], 'dP_GPa': []})
        row['dE_meV'].append((r['energy'] - f['energy']) / n * 1000)
        row['dF'] += list((r['forces'] - np.array(f['forces'])).ravel())
        row['dP_GPa'].append(r['pressure_GPa'] - float(np.mean(f['stress_kB'][:3])) / 10)
    summary = {}
    allE = []
    for g, row in sorted(out['by_folder'].items()):
        e, df, dpp = np.array(row['dE_meV']), np.array(row['dF']), np.array(row['dP_GPa'])
        allE += list(e)
        summary[g] = {'n': len(e), 'energy_mean_meV': round(float(e.mean()), 2), 'energy_rmse_meV': round(float(np.sqrt((e ** 2).mean())), 2),
                      'force_rmse_eV_A': round(float(np.sqrt((df ** 2).mean())), 3), 'pressure_rmse_GPa': round(float(np.sqrt((dpp ** 2).mean())), 2)}
        print(g, summary[g])
    summary['ALL'] = {'n': len(allE), 'energy_rmse_meV': round(float(np.sqrt(np.mean(np.square(allE)))), 2)}
    print('ALL', summary['ALL'], out['finite_difference'])
    (HERE/'dp_validation.json').write_text(json.dumps({'finite_difference': out['finite_difference'], 'vs_dft': summary}, indent=1) + '\n')


if __name__ == '__main__':
    main()
