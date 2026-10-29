"""Read the fetched archive members into DATA/pool.json.gz (every number used by make_task.py, unmodified).

    python3 -B extract.py [DATA]     (DATA defaults to /data, the cache mounted by run_in_image.sh)
"""
import glob
import gzip
import json
import os
import sys
import numpy as np

DATA = sys.argv[1] if len(sys.argv) > 1 else '/data'
IN = os.path.join(DATA, 'mantel_data_files_v1/INPUTS')
FIG = os.path.join(DATA, 'mantel_data_files_v1/FIGURES')

SYSTEMS = {
    'Mg2IrH6': {0: 'hydrides/Mg2IrH6', **{p: f'Mg2IrH6_pressure_dependence/Mg2IrH6-{p}GPa' for p in (5, 10, 15, 20, 25, 30)}},
    'H3S': {p: f'H3S_pressure_dependence/H3S-{p}GPa' for p in (140, 160, 170, 180, 190, 200, 210, 220, 230, 240)},
    'CaH6': {160: 'hydrides/CaH6'},
    'LaH10': {250: 'hydrides/LaH10'},
    'YSbH6': {50: 'hydrides/YSbH6'},
}


def text(path):
    with open(path, errors='replace') as f:
        return f.read()


def read_dos(path):
    head = open(path).readline()
    ef = float(head.split('=')[-1].split()[0])
    d = np.loadtxt(path)
    return {'header': head.strip(), 'EF_eV': ef, 'E_eV': d[:, 0].tolist(), 'dos': d[:, 1].tolist(),
            'int_dos': d[:, 2].tolist()}


def read_a2f(path):
    lines = open(path).read().splitlines()
    head = lines[0]
    d = np.loadtxt(path)
    return {'header': head, 'table': d.tolist()}


def read_W(folder):
    raw = np.load(os.path.join(folder, 'W_ee_raw.npy'))
    dr = np.load(os.path.join(folder, 'dos_raw.npy'))
    return {'E_eV': dr[:, 0].tolist(), 'dos': dr[:, 1].tolist(), 'NNW': raw.tolist()}


def read_cell(folder, seed):
    for name in (f'{seed}.scf.in', 'scf.in'):
        p = os.path.join(folder, name)
        if os.path.exists(p):
            return text(p)
    return None


pool = {}
for system, pressures in SYSTEMS.items():
    for p, rel in pressures.items():
        folder = os.path.join(IN, rel)
        seed = system
        rec = {'folder': rel}
        iso = os.path.join(folder, 'IsoME')
        if os.path.exists(os.path.join(iso, 'rescaled-alpha2F.dat')):
            rec['a2f'] = read_a2f(os.path.join(iso, 'rescaled-alpha2F.dat'))
            rec['dos_file'] = read_dos(os.path.join(iso, f'{seed}.dos.dat'))
            rec['runs'] = {}
            for s in sorted(glob.glob(os.path.join(iso, '*/Summary.dat'))):
                run = os.path.basename(os.path.dirname(s))
                rec['runs'][run] = {'Summary.dat': text(s), 'Info.txt': text(os.path.join(os.path.dirname(s), 'Info.txt'))}
            rec['scripts'] = {os.path.basename(j): text(j) for j in sorted(glob.glob(os.path.join(iso, '*.jl')))}
            if os.path.exists(os.path.join(iso, 'README')):
                rec['isome_readme'] = text(os.path.join(iso, 'README'))
        elif os.path.exists(os.path.join(folder, f'DOS/{seed}.dos.dat')):
            rec['dos_file'] = read_dos(os.path.join(folder, f'DOS/{seed}.dos.dat'))
        if os.path.exists(os.path.join(folder, 'W_ee_raw.npy')):
            w = read_W(folder)
            if system not in ('Mg2IrH6', 'H3S'):     # 2000x2000 grids: keep only the diagonal
                nnw = np.array(w.pop('NNW'))
                w['NNW_diagonal'] = np.diag(nnw).tolist()
            rec['W'] = w
        if os.path.exists(os.path.join(folder, 'mu_results.dat')):
            rec['mu_results'] = text(os.path.join(folder, 'mu_results.dat'))
        logs = glob.glob(os.path.join(folder, '*_mantel.log'))
        if logs:
            rec['mantel_log_header'] = '\n'.join(l for l in text(logs[0]).splitlines()[:40] if 'Calculating' not in l)
        rec['scf_in'] = read_cell(folder, seed)
        relax = glob.glob(os.path.join(folder, 'relax/*.relax.in'))
        if relax:
            rec['relax_in'] = text(relax[0])
        pool[f'{system}/{p}'] = rec
        print(system, p, sorted(rec))

figs = {}
for f in sorted(glob.glob(os.path.join(FIG, 'pressure_dependence/*/*.txt')) + glob.glob(os.path.join(FIG, 'pressure_dependence/H3S/anharmonicity/*.txt'))):
    if not os.path.basename(f).startswith('._'):
        figs[os.path.relpath(f, FIG)] = text(f)
pool['_figures'] = figs
with gzip.open(os.path.join(DATA, 'pool.json.gz'), 'wt') as f:
    json.dump(pool, f)
print('wrote', os.path.join(DATA, 'pool.json.gz'))
