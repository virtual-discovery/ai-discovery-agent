"""Extract the DFT validation subset served by the lab from the released training/test pickles.

    $CACHE/venv/bin/python -B extract_dft.py $CACHE dft_subset.json

Energies are the authors' `energy_corrected` (VASP PBE total energy minus the isolated-atom reference
energies used in the fit: Pd -1.47571017 eV, H -0.05660531 eV), i.e. the same energy zero as the ACE.
Forces in eV/A; stress (when archived) in eV/A^3 with the ASE sign convention.
"""
import json
import sys
import numpy as np
import pandas as pd

C, OUT = sys.argv[1], sys.argv[2]
frames = []
for split, f in (('train', 'df_train.pckl.gzip'), ('test', 'df_test.pckl.gzip')):
    d = pd.read_pickle(f'{C}/potential_files/ace/{f}', compression='gzip')
    d['split'] = split
    frames.append(d)
df = pd.concat(frames)
rows = []


def add(r, label):
    a = r.ase_atoms
    st = r.results.get('stress') if isinstance(r.results, dict) else None
    rows.append({'label': label, 'split': r.split, 'formula': a.get_chemical_formula(),
                 'cell': np.round(np.array(a.cell[:]), 8).tolist(), 'symbols': a.get_chemical_symbols(),
                 'positions': np.round(a.positions, 8).tolist(), 'energy_eV': float(r.energy_corrected),
                 'forces_eV_A': np.round(np.array(r.forces, float), 6).tolist(),
                 'stress_eV_A3': None if st is None else np.round(np.array(st, float), 8).tolist()})


for _, r in df.iterrows():
    n = r['name']
    if 'pdhX/' in n:
        tag = n.split('pdhX/')[1].split('/cfg')[0].split('/')[-1]          # e.g. PdH_0.5-octa
        kind = n.split('/cfg/')[1].split('/')[1]
        if kind in ('relax', 'static'):
            add(r, f'PdHx-16Pd {tag} {kind}')
    elif 'refs/pd_fcc' in n:
        add(r, 'Pd fcc ' + n.split('/cfg/Pd/')[1].rsplit('/data.json', 1)[0])
    elif 'refs/pdH' in n:
        add(r, 'PdH rocksalt ' + n.split('/cfg/PdH/')[1].rsplit('/data.json', 1)[0])
    elif r.n_Pd == 0 and r.n_H == 2 and r.ase_atoms.get_volume() > 100:
        d = r.ase_atoms.get_distance(0, 1, mic=True)
        if 0.6 <= d <= 1.0:
            add(r, f'H2 molecule d={d:.4f}')
    elif 'H_vac/' in n and '/static/' in n:
        add(r, 'Pd31 vacancy cell ' + n.split('H_vac/')[1].split('/cfg')[0].split('/')[-1])
rows.sort(key=lambda x: x['label'])
for i, x in enumerate(rows):
    x['id'] = f'DFT-{i + 1:03d}'
json.dump(rows, open(OUT, 'w'))
import collections
print(len(rows), collections.Counter(x['label'].split()[0] for x in rows))
