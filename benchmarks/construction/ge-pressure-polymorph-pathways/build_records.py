"""Write the lab's records.json.gz and the agent's catalog.json from pool.json.gz.

    python3 -B build_records.py

Catalog = every released training configuration of the five phase folders CD, bSn, ST12, BC8 and R8.
Each record holds what the OUTCAR shows: INCAR-type settings, k mesh, whether the relaxation reached the
required accuracy, and every ionic step (energy(sigma->0), free energy TOTEN, cell, volume, stress in kB with
VASP sign and order, Cartesian positions and forces, electronic steps). Nothing is renormalized.
"""
import gzip
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/ge-pressure-polymorph-pathways'
POOL = json.load(gzip.open(Path.home()/'.cache/ai-theorist/ge-pressure-polymorph-pathways/pool.json.gz'))
PHASES = {'CD': 'cubic diamond (Fd-3m), Ge-I', 'bSn': 'beta-tin (I4_1/amd), Ge-II',
          'ST12': 'simple tetragonal, 12 atoms per cell (P4_32_12), Ge-III',
          'BC8': 'body-centred cubic, 8 atoms per primitive cell (Ia-3), Ge-IV',
          'R8': 'rhombohedral, 8 atoms per primitive cell (R-3)'}
SETTING_KEYS = {'ENCUT': 'ENCUT_eV', 'PREC': 'PREC', 'ISMEAR': 'ISMEAR', 'SIGMA': 'SIGMA_eV', 'kmesh': 'kpoint_mesh',
                'PSTRESS': 'PSTRESS_kB', 'ISIF': 'ISIF', 'IBRION': 'IBRION', 'NSW': 'NSW', 'EDIFF': 'EDIFF_eV',
                'EDIFFG': 'EDIFFG', 'NELM': 'NELM', 'ISPIN': 'ISPIN', 'POTCAR': 'POTCAR'}


def record_id(folder):
    _, phase, cfg = folder.split('/')
    return '%s-%03d' % (phase, int(cfg.split('_')[1]))


def main():
    records, catalog = {}, []
    for folder in sorted(POOL, key=lambda k: (k.split('/')[1], int(k.split('_')[-1]) if k.split('_')[-1].isdigit() else 0)):
        parts = folder.split('/')
        if len(parts) != 3 or parts[0] != 'training_set' or parts[1] not in PHASES:
            continue
        v = POOL[folder]
        if not v['frames']:
            continue
        s = v['settings']
        steps = []
        for k, f in enumerate(v['frames']):
            steps.append({'step': k + 1, 'energy_sigma0_eV': f['energy'], 'free_energy_TOTEN_eV': f['toten'],
                          'volume_A3': f['volume'], 'cell_A': f['cell'], 'stress_kB': f['stress_kB'],
                          'positions_A': f['positions'], 'forces_eV_A': f['forces'], 'electronic_steps': f['n_electronic']})
        rid = record_id(folder)
        records[rid] = {'phase': parts[1], 'archive_folder': folder, 'natoms': s['NIONS'],
                        'settings': {SETTING_KEYS[k]: s.get(k) for k in SETTING_KEYS},
                        'reached_required_accuracy': v['reached_accuracy'], 'ionic_steps': steps,
                        'units': {'energy': 'eV per cell', 'volume': 'A^3 per cell', 'stress_kB': 'kB, order XX YY ZZ XY YZ ZX, positive = compressive',
                                  'positions': 'Cartesian A', 'forces': 'eV/A'}}
        f = v['frames'][-1]
        catalog.append({'id': rid, 'phase': parts[1], 'natoms': s['NIONS'], 'ionic_steps': len(v['frames']),
                        'volume_per_atom_A3': round(f['volume'] / s['NIONS'], 2)})
    with gzip.open(TASK/'environment/lab/records.json.gz', 'wt') as fh:
        json.dump(records, fh, separators=(',', ':'))
    cat = {'material': 'Ge, PBE (VASP 6, PAW Ge_d, ENCUT 500 eV), training database of the DeePMD potential of '
                       'Fantasia et al. (2024), Materials Cloud archive 238ad-1bg44',
           'phases': PHASES,
           'records_note': 'One record = one archived VASP calculation (a single point or a relaxation). '
                           'volume_per_atom_A3 refers to the last ionic step.',
           'purchase_costs': {'record': 1, 'dp': 1}, 'budget': 34, 'records': catalog}
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(cat, indent=0) + '\n')
    print(len(records), 'records;', {p: sum(r['phase'] == p for r in catalog) for p in PHASES})


if __name__ == '__main__':
    main()
