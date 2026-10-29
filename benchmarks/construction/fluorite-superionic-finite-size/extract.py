"""Build the run pool for fluorite-superionic-finite-size from the files fetched by fetch.py.

    python3 -B extract.py      # -> ~/.cache/ai-theorist/fluorite-superionic-finite-size/pool.npz + pool_meta.json

One record per released NVT run (material, supercell L, target T): the full msd_CM_B4.dat table, the thermo
columns of nvt_log (Time, TotEng, PotEng, Temp, Press), and the LAMMPS input deck with its timestep and dump
interval. PbF2/NVT_L_T_6p056_DW duplicates the main PbF2 set byte for byte at L = 2-4 (checked here), so only
its L = 1 and L = 5-7 runs are added.
"""
import json
import os
import re
from pathlib import Path
import numpy as np

CACHE = Path(os.path.expanduser('~/.cache/ai-theorist/fluorite-superionic-finite-size'))
RAW = CACHE/'raw'
SETS = [('CaF2', 'CaF2/NVT_L_T', 'CaF2.in', 12), ('PbF2', 'PbF2/NVT_L_T_6p056', 'PbF2.in', 12),
        ('PbF2', 'PbF2/NVT_L_T_6p056_DW', 'PbF2.in', 12), ('UO2', 'UO2/NVT_L_T', 'UO2.in', 12),
        ('AgI', 'alpha-AgI/NVT_L_T_5.37ang', 'AgI.in', 4)]


def deck_params(text):
    dt = float(re.search(r'variable\s+dt\s+equal\s+([\d.]+)', text).group(1))
    every = int(re.search(r'variable\s+traj_print_interval\s+equal\s+(\d+)', text).group(1))
    return dt, every


def main():
    meta, arrays = [], {}
    for material, setdir, deck, atoms_per_cell in SETS:
        for run in sorted((RAW/setdir).glob('*x*x*/T_*_K')):
            L, T = int(run.parent.name.split('x')[0]), int(run.name.split('_')[1])
            if setdir.endswith('_DW'):
                main_run = RAW/'PbF2/NVT_L_T_6p056'/run.parent.name/run.name
                if main_run.exists():
                    assert (main_run/'msd_CM_B4.dat').read_bytes() == (run/'msd_CM_B4.dat').read_bytes()
                    continue
            text = (run/deck).read_text()
            dt, every = deck_params(text)
            msd = np.loadtxt(run/'msd_CM_B4.dat')
            log = np.load(run/'nvt_log.npy')
            key = f'{material}-{L}-{T}'
            arrays[key + '/msd'] = msd
            arrays[key + '/log'] = log
            meta.append({'key': key, 'material': material, 'L': L, 'T': T, 'natoms': atoms_per_cell * L ** 3,
                         'set': setdir, 'timestep_ps': dt, 'dump_every': every, 'frame_ps': dt * every,
                         'msd_rows': len(msd), 'log_rows': len(log), 'deck': text})
    np.savez_compressed(CACHE/'pool.npz', **arrays)
    (CACHE/'pool_meta.json').write_text(json.dumps(meta))
    print(len(meta), 'runs ->', CACHE/'pool.npz')


if __name__ == '__main__':
    main()
