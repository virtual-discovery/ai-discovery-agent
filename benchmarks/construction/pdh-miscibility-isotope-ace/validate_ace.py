"""Validate ace_numpy.py against python-ace (pyace) on released test-set structures and perturbed cells.

    $CACHE/venv/bin/python -B validate_ace.py $CACHE   (build host, pyace installed in $CACHE/venv)
"""
import sys, time, json
import numpy as np
import pandas as pd
from pyace import PyACECalculator
sys.path.insert(0, '.')
from ace_numpy import ACE

C = sys.argv[1]
calc = PyACECalculator(C + '/potential_files/ace/PdH_PBE_v03a.yaml')
ace = ACE(C + '/work/PdH_ace_ctilde.json')
df = pd.read_pickle(C + '/potential_files/ace/df_test.pckl.gzip', compression='gzip')
rng = np.random.default_rng(0)
rows = df.sample(60, random_state=1)
emax = fmax = smax = 0
for at in rows.ase_atoms:
    at = at.copy()
    at.calc = calc
    e_ref, f_ref, s_ref = at.get_potential_energy(), at.get_forces(), at.get_stress(voigt=False)
    t = time.time()
    out = ace.compute(at.cell[:], at.get_chemical_symbols(), at.positions)
    dt = time.time() - t
    emax = max(emax, abs(out['energy'] - e_ref) / len(at))
    fmax = max(fmax, np.abs(out['forces'] - f_ref).max())
    smax = max(smax, np.abs(out['stress'] - s_ref).max())
print('max |dE|/atom', emax, 'max |dF|', fmax, 'max |dstress|', smax)
# close-contact structure exercising the ZBL switching region
from ase import Atoms
at = Atoms('Pd2H2', positions=[[0, 0, 0], [1.25, 0.2, 0.1], [2.5, 2.0, 0.3], [2.5, 2.5, 0.1]], cell=[6, 6, 6], pbc=True)
at.calc = calc
out = ace.compute(at.cell[:], at.get_chemical_symbols(), at.positions)
print('zbl case', at.get_potential_energy(), out['energy'], np.abs(out['forces'] - at.get_forces()).max())
# finite-difference force check of the numpy evaluator
x = at.positions.copy(); h = 1e-5; fd = np.zeros_like(x)
for i in range(len(x)):
    for k in range(3):
        xp = x.copy(); xp[i, k] += h; xm = x.copy(); xm[i, k] -= h
        fd[i, k] = -(ace.compute(at.cell[:], at.get_chemical_symbols(), xp, forces=False)['energy']
                     - ace.compute(at.cell[:], at.get_chemical_symbols(), xm, forces=False)['energy']) / (2 * h)
print('fd force check', np.abs(fd - out['forces']).max())
for pos in ([[0, 0, 0], [0.9, 0.2, 0.1], [2.5, 2.0, 0.3], [2.5, 2.5, 0.1]], [[0, 0, 0], [0.98, 0.0, 0.0], [2.5, 2.0, 0.3], [2.5, 2.46, 0.1]]):
    at = Atoms('PdH3', positions=pos, cell=[6, 6, 6], pbc=True)
    at.calc = calc
    out = ace.compute(at.cell[:], at.get_chemical_symbols(), at.positions)
    print('zbl region', at.get_potential_energy(), out['energy'], np.abs(out['forces'] - at.get_forces()).max())
# timing on a 4x4x4 conventional fcc Pd cell (256 Pd) with 8 H
from ase.build import bulk
big = bulk('Pd', 'fcc', a=3.94, cubic=True).repeat((4, 4, 4))
for k in range(8):
    big.append('H'); big.positions[-1] = [1.97 + 3.94 * (k % 4), 3.94 * (k // 4) * 2, 0]
t = time.time(); out = ace.compute(big.cell[:], big.get_chemical_symbols(), big.positions); print('264 atoms', time.time() - t, 's')
big.calc = calc; print('E', out['energy'], big.get_potential_energy())
