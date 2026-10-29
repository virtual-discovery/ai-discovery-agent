"""Extract the eight triangular dihalides from the AMaRaNTA AiiDA provenance databases.

    python3 -B extract.py [CACHE_DIR]      (run fetch.py first)

Writes pool.json next to this script. For each compound (the PK listed in the archive's
other_data/PK_list.txt) it stores the unit cell and sites, the authors' J_dict, the three
pairs ("couples") AMaRaNTA used, and every four-state calculation: parameter label, input
supercell (cell and sites), the constrained moment directions of the four states, the four
VASP total energies, the returned parameter and the convergence flags.

Asserts: the AiiDA J_dict equals the archive's *J_dict.txt to the printed precision, and the
four-state combination (E1 - E2 - E3 + E4)/4 of the stored energies equals each returned
parameter.
"""
import json
import re
import sqlite3
import sys
from pathlib import Path

CACHE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/amaranta-2d-triangular-helimagnets'
OUT = Path(__file__).with_name('pool.json')
COMPOUNDS = {'NiCl2': ('group10', 11815), 'NiBr2': ('group2', 12821), 'VCl2': ('group2', 16233), 'VBr2': ('group2', 17997),
             'VI2': ('group2', 14160), 'MnCl2': ('group3', 4157), 'MnBr2': ('group3', 3579), 'MnI2': ('group3', 2705)}
TXT_KEYS = {'Jxx': 'J_xx', 'Jyy': 'J_yy', 'Jzz': 'J_zz', 'Jxy': 'J_xy', 'Jyx': 'J_yx', 'Jyz': 'J_yz', 'Jzy': 'J_zy',
            'Jxz': 'J_xz', 'Jzx': 'J_zx', 'J2': 'J_scalar_2NN', 'J3': 'J_scalar_3NN', 'A': 'A_scalar'}


class DB:
    def __init__(self, path):
        self.c = sqlite3.connect(path)

    def node(self, i):
        r = self.c.execute('select id, process_type, attributes, ctime, uuid from db_dbnode where id=?', (i,)).fetchone()
        return {'id': r[0], 'ptype': r[1] or '', 'attrs': json.loads(r[2]) if r[2] else {}, 'ctime': r[3], 'uuid': r[4]}

    def inputs(self, i):
        return [(self.node(j), label, kind) for j, label, kind in
                self.c.execute('select input_id, label, type from db_dblink where output_id=?', (i,))]

    def outputs(self, i):
        return [(self.node(j), label, kind) for j, label, kind in
                self.c.execute('select output_id, label, type from db_dblink where input_id=?', (i,))]


def structure(attrs):
    return {'cell': attrs['cell'], 'sites': [{'species': s['kind_name'], 'position': s['position']} for s in attrs['sites']]}


def vasp_run(db, wc):
    run = {'pk': wc['id'], 'exit_status': wc['attrs'].get('exit_status')}
    for n, label, kind in db.outputs(wc['id']):
        if label == 'misc':
            rs = n['attrs'].get('run_status', {})
            run.update(energy_eV=n['attrs']['total_energies']['energy_extrapolated'], electronic_converged=rs.get('electronic_converged'),
                       nelm_breach=rs.get('contains_nelm_breach'), scf_steps=(rs.get('last_iteration_index') or [None, None])[1])
    for n, label, kind in db.inputs(wc['id']):
        if label == 'kpoints':
            run['kmesh'] = n['attrs'].get('mesh')
        if label == 'parameters':
            incar = n['attrs'].get('incar', {})
            run['incar'] = {k: incar.get(k) for k in ('ENCUT', 'LSORBIT', 'LNONCOLLINEAR', 'I_CONSTRAINED_M', 'LAMBDA', 'LDAUU', 'LDAUJ', 'ISPIN', 'EDIFF')}
    return run


def extract(db, pk):
    wc = db.node(pk)
    ins = {label: n for n, label, kind in db.inputs(pk)}
    outs, subs = {}, []
    for n, label, kind in db.outputs(pk):
        if kind == 'return':
            outs[label] = n
        elif kind == 'call_work' and 'SubWorkChain' in n['ptype']:
            subs.append(n)
    rec = {'pk': pk, 'uuid': wc['uuid'], 'process_state': wc['attrs'].get('process_state'), 'exit_status': wc['attrs'].get('exit_status'),
           'structure': structure(ins['structure']['attrs']), 'J_dict': outs['J_dict']['attrs'], 'couples': outs['couples']['attrs'],
           'spin_label': outs['spin']['attrs']['value'], 'units': outs['units']['attrs'], 'four_state_sets': []}
    for s in sorted(subs, key=lambda n: n['ctime']):
        sins = {label: n for n, label, kind in db.inputs(s['id'])}
        souts = db.outputs(s['id'])
        fe = [n for n, label, kind in souts if 'find_energy' in n['ptype']][0]
        energies = [n for n, label, kind in db.inputs(fe['id']) if label == 'list_of_energies'][0]['attrs']['list']
        runs = [vasp_run(db, n) for n, label, kind in sorted(souts, key=lambda x: x[0]['ctime'])
                if kind == 'call_work' and 'vasp' in n['ptype']]
        assert sorted(r['energy_eV'] for r in runs) == sorted(energies)
        rec['four_state_sets'].append({
            'label': sins['which_label']['attrs']['value'], 'soc': sins['which_incar']['attrs']['value'] == 2,
            'supercell': structure(sins['structure_sub']['attrs']), 'kmesh': sins['kpoints_sub']['attrs'].get('mesh'),
            'moment_directions': sins['magmom_sub']['attrs']['list'], 'energies_eV': energies,
            'parameter_meV': [n for n, label, kind in souts if label == 'J_parameter'][0]['attrs']['value'], 'vasp_runs': runs})
    return rec


def main():
    pool = {}
    for name, (group, pk) in COMPOUNDS.items():
        rec = extract(DB(CACHE/f'{group}_db.sqlite3'), pk)
        text = (CACHE/'output_readble'/f'{name}.J_dict.txt').read_text()
        assert f'(PK {pk})' in text
        for key, value in re.findall(r'^\s+(\w+)\s+([-+]\d+\.\d+) meV', text, re.M):
            assert abs(rec['J_dict'][TXT_KEYS[key]] - float(value)) < 6e-4, (name, key)
        for fs in rec['four_state_sets']:
            e = fs['energies_eV']
            assert abs((e[0] - e[1] - e[2] + e[3]) / 4 * 1000 - fs['parameter_meV']) < 1e-6
            assert all(r['electronic_converged'] and not r['nelm_breach'] for r in fs['vasp_runs'])
        pool[name] = rec
        print(name, pk, len(rec['four_state_sets']), 'four-state sets')
    OUT.write_text(json.dumps(pool, indent=None, separators=(',', ':')) + '\n')
    print('wrote', OUT, OUT.stat().st_size // 1024, 'kB')


if __name__ == '__main__':
    main()
