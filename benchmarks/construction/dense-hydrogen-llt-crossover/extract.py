"""Build the compact pool from the fetched archive members (run inside the runtime image).

    python3 -B extract.py      # reads ~/.cache/ai-theorist/dense-hydrogen-llt-crossover/raw, writes pool_thermo.npz,
                               # pool_pimd.npz and pool_meta.json in the same cache folder

Classical runs: every thermo block of log.lammps and slurm-*.out is parsed; blocks are stitched by step (a later job
overrides an earlier one at the same step; log.lammps repeats the last job). The served series is every 5th thermo
row (every 500 steps = 0.1 ps) of Temp, PotEng, TotEng, Press. Energies are stored as float32 offsets from a per-run
reference (E_ref, eV) that the lab adds back. PIMD: every 20th row of the archived thermodata (text despite the .npy
suffix) for the columns nmove, H, E_pot, Temp(K), quantum_kin_virial, pressure(GPa).
"""
import glob
import hashlib
import json
import os
import re
from pathlib import Path
import numpy as np
import lmplog

CACHE = Path(os.environ.get('HCACHE', os.path.expanduser('~/.cache/ai-theorist/dense-hydrogen-llt-crossover')))
RAW = CACHE/'raw'
STRIDE = 5
PIMD_STRIDE = 20
PIMD_COLS = [1, 2, 4, 7, 8, 10]   # nmove, H, E_pot, Temp(K), quantum_kin_virial, pressure (GPa)


def jobs_of(rsd):
    out = []
    for p in sorted(rsd.glob('slurm-*.out'), key=lambda p: int(re.search(r'(\d+)', p.name).group(1))):
        text = p.read_text(errors='replace')
        start = 'read_data init.lmp' if re.search(r'^\s*read_data', text, re.M) else None
        m = re.search(r'^\s*read_restart\s+(?!CPU)(\S+)', text, re.M)
        if m:
            start = 'read_restart ' + m.group(1)
        elif re.search(r'read_restart CPU', text):
            start = 'read_restart (restart file of the previous job)'
        steps = [rows[:, 0] for head, rows in lmplog.blocks(text) if len(rows)]
        err = re.search(r'^ERROR.*$', text, re.M)
        job = {'job': p.stem, 'start': start}
        if steps:
            job.update(first_step=int(steps[0][0]), last_step=int(steps[-1][-1]))
        if err:
            job['error'] = err.group(0)[:160]
        out.append(job)
    return out


def main():
    arrays, meta = {}, {'runs': {}, 'pimd': {}}
    for tard in sorted(RAW.glob('N_*_T*K.tar')):
        n, t = re.match(r'N_(\d+)_T(\d+)K', tard.name).groups()
        n, t = int(n), int(t)
        for rsd in sorted((tard/f'T{t}K').glob('rs*')):
            rs = rsd.name[2:]
            files = [p for p in rsd.iterdir() if p.name == 'log.lammps' or p.name.startswith('slurm-')]
            rows = {}
            for job, seg in lmplog.run_series(files):
                for r in seg:
                    rows[int(r[0])] = r
            steps = np.array(sorted(rows))
            X = np.array([rows[s] for s in steps])
            assert np.all(np.diff(steps) == 100), (n, t, rs)
            keep = (steps % (100 * STRIDE)) == 0
            X = X[keep]
            vol = float(X[0, 6])
            assert np.allclose(X[:, 6], vol)
            e_ref = float(np.round(X[:, 4].mean(), 1))
            key = f'{n}/{t}/{rs}'
            arrays[key] = np.c_[X[:, 1], X[:, 2] - e_ref, X[:, 4] - e_ref, X[:, 5]].astype(np.float32)
            meta['runs'][key] = {'N': n, 'T': t, 'rs': rs, 'first_step': int(steps[keep][0]),
                                 'last_step': int(steps[keep][-1]), 'stride_steps': 100 * STRIDE,
                                 'E_ref_eV': e_ref, 'volume_A3': vol, 'length_ps': round(float(steps[-1]) * 2e-4, 2),
                                 'input_deck': (rsd/'input.lammps').read_text(), 'jobs': jobs_of(rsd)}
    np.savez_compressed(CACHE/'pool_thermo.npz', **arrays)
    pimd, digests = {}, {}
    for f in sorted(glob.glob(str(RAW/'PIMD/*K/rs*/thermodata.npy'))):
        temp, d = f.split('/')[-3], f.split('/')[-2]
        rs, rep = d[2:].rsplit('.', 1)
        blob = Path(f).read_bytes()
        h = hashlib.sha256(blob).hexdigest()
        key = f'{temp}/{rs}/{rep}'
        meta['pimd'][key] = {'T_label': temp, 'rs': rs, 'replica': int(rep), 'sha256': h, 'rows_archived': None}
        if h not in digests:
            a = np.loadtxt(f)
            meta['pimd'][key]['rows_archived'] = len(a)
            digests[h] = key
            pimd[key.replace('/', '|')] = a[::PIMD_STRIDE][:, PIMD_COLS].astype(np.float32)
        meta['pimd'][key]['data_key'] = digests[h].replace('/', '|')
    for key, m in meta['pimd'].items():
        m['rows_archived'] = meta['pimd'][m['data_key'].replace('|', '/')]['rows_archived']
    np.savez_compressed(CACHE/'pool_pimd.npz', **pimd)
    (CACHE/'pool_meta.json').write_text(json.dumps(meta, indent=0))
    dup = sum(1 for k, m in meta['pimd'].items() if m['data_key'].replace('|', '/') != k)
    print(len(arrays), 'classical runs;', len(meta['pimd']), 'PIMD files,', dup, 'byte-identical to another label')


if __name__ == '__main__':
    main()
