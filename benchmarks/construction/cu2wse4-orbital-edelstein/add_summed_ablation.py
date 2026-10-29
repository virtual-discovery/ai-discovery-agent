"""Add the 'Berry-phase + intra-atomic moment' ablation to solution/ablation_values.json and calibration.json from the
lab's own compact caches (identical to make_task.py's lab_cache), without the raw construction caches.

    bash run_in_image.sh add_summed_ablation.py
"""
import json
from pathlib import Path
import numpy as np
import make_task as M

LABDATA = M.LAB/'data'


def lab_cache(name, n):
    z = np.load(LABDATA/f'mesh_{name}_{n}.npz')
    c = {k: z[k].astype(float) for k in z.files}
    c['centres'] = c['e0'] + (np.arange(len(c['count'])) + 0.5) * M.T.EBIN
    return c


def fermi(name):
    return float(np.load(LABDATA/f'model_{name}.npz')['fermi_Ha']) * M.T.HARTREE


M.lab_cache, M.fermi = lab_cache, fermi
hidden = json.loads((M.LAB/'private/hidden.json').read_text())
edges = hidden['reference']['edges_eV']
ref = M.values(M.REF_MESH, M.REF_MESH, edges)
assert abs(ref['ratio']['CBM+0.10'] / hidden['targets']['orbital_to_spin_ratio']['CBM+0.10'] - 1) < 1e-9
s = M.values(M.REF_MESH, M.REF_MESH, edges, moment='orb+atomic')
entry = {'orbital_to_spin_ratio': s['ratio'], 'orbital_magnetization_A_per_m': s['mag'], 'soc_to_nosoc_orbital_ratio': s['soc_nosoc']}
for path in (M.TASK/'solution/ablation_values.json', M.HERE/'calibration.json'):
    d = json.loads(path.read_text())
    d['ablations']['berry_plus_intra_atomic_moment'] = entry
    path.write_text(json.dumps(d, indent=1) + '\n')
print(json.dumps(entry, indent=1))
