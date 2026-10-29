"""Extract the Fe configuration pool from the Materials Cloud archive 3p6fv-8hd59 (Fe.tar.gz, md5 3e76b9b6eaeba27b7873be4327806477).

Run from the repository root after downloading the archive into downloads/heis-configdesign/:
    python3 benchmarks/construction/heis-configdesign-fe/extract.py
Writes pool.json: 101 spin spirals and 89 random noncollinear states with DFT energies
(meV/atom relative to the ferromagnet, as in the authors' CSVs) and DFT local moments.
"""
import csv, glob, json, re
import numpy as np
from pathlib import Path

SRC = Path('downloads/heis-configdesign/Fe')
OUT = Path(__file__).with_name('pool.json')
E_FM = -8274.6594  # meV/atom; the authors' CSV energies are OUTCAR energy(sigma->0)/atom minus this


def outcar(path):
    text = path.read_text()
    nat = int(re.findall(r'NIONS\s*=\s*(\d+)', text)[0])
    moments = []
    for c in 'xyz':
        seg = text[text.rfind(f'magnetization ({c})'):].split('\n')[4:4 + nat]
        moments.append([float(line.split()[-1]) for line in seg])
    energy = float(re.findall(r'energy\(sigma->0\) =\s+([-\d.]+)', text)[-1]) / nat * 1000 - E_FM
    return nat, np.array(moments).T, energy


def poscar(path):
    lines = path.read_text().split('\n')
    scale = float(lines[1])
    cell = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]) * scale
    nat = int(lines[6].split()[0])
    frac = np.array([[float(x) for x in lines[8 + i].split()[:3]] for i in range(nat)])
    return cell, frac


spirals = []
for row in csv.DictReader(open(SRC/'spin_spiral/spin-spiral-energies.csv')):
    q = [float(row[k]) for k in ('qx', 'qy', 'qz')]
    spirals.append({'q_reciprocal_lattice': q, 'cone_angle_deg': float(row['theta_z_1']),
                    'moment_muB': float(row['magmom']), 'energy_meV_per_atom': float(row['energy'])})
randoms = []
for d in sorted(glob.glob(str(SRC/'random_states/0*-random-*/random-*'))):
    d = Path(d)
    nat, m, energy = outcar(d/'OUTCAR')
    cell, frac = poscar(d/'POSCAR')
    randoms.append({'supercell': d.parent.name.split('-')[-1], 'label': d.parent.name.split('-')[-1] + '-' + d.name.split('-')[1],
                    'cell_A': cell.round(10).tolist(), 'fractional_positions': frac.round(12).tolist(),
                    'moments_muB': m.round(3).tolist(), 'energy_meV_per_atom': round(energy, 3)})
csv_e = sorted(float(r[3]) for r in list(csv.reader(open(SRC/'random_states/random-states-energies.csv')))[1:])
assert np.abs(np.sort([r['energy_meV_per_atom'] for r in randoms]) - csv_e).max() < 2e-3
OUT.write_text(json.dumps({'source': 'Materials Cloud 3p6fv-8hd59 Fe.tar.gz', 'lattice_A': 2.87,
                           'primitive_cell_A': [[-1.435, 1.435, 1.435], [1.435, -1.435, 1.435], [1.435, 1.435, -1.435]],
                           'spirals': spirals, 'random_states': randoms}))
print(len(spirals), 'spirals;', len(randoms), 'random states ->', OUT)
