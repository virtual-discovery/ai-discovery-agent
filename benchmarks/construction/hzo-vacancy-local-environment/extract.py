"""Extract the HZO oxygen-vacancy pool from the Apollo archive (DOI 10.17863/CAM.121316, Data.zip).

Run fetch.py first, then:
    python3 benchmarks/construction/hzo-vacancy-local-environment/extract.py [--raw DIR]
Writes pool.json next to this file:
  - the three composition supercells (CIF text, verbatim) and the 3x3x3 CASTEP supercell CIF;
  - VASP-PBE neutral-vacancy total energies for all 64 O sites of each composition and the pristine-cell
    energies (Data_Fig2b_g.csv, rows at E_fermi = 0; the neutral columns do not depend on E_fermi);
  - the authors' x = 0.50 site table (FIGS12.csv), CASTEP multi-functional energies (FIGS3.csv),
    CASTEP-PBE 3x3x3 energies (FIGS4.txt) and the environment frequencies of FIGS5.txt.
Every derived relation used later is asserted against the archive here.
"""
import argparse
import collections
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT = Path.home()/'.cache/ai-theorist/hzo-vacancy-local-environment/raw/Apollo deposit'
CIFS = {'0.25': 'Main text/Structures/Fig2b_e_0.25_HZO.cif', '0.50': 'Main text/Structures/Fig2c_f_0.50_HZO.cif',
        '0.75': 'Main text/Structures/FIG2d_g_0.75_HZO.cif'}
MU = {'poor': -10.09, 'mid': -7.33, 'high': -4.57}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw', type=Path, default=DEFAULT)
    a = p.parse_args()
    raw = a.raw
    pool = {'source': 'Apollo 10.17863/CAM.121316 Data.zip (md5 70721b2c03a34c84e9a2deab61160afb)',
            'mu_O_eV': MU, 'cif': {x: (raw/f).read_text() for x, f in CIFS.items()},
            'cif_333': (raw/'SM material/FIGS4 3x3x3_supercell.cif').read_text()}
    # VASP PBE neutral vacancies, three compositions
    vac, bulk, seen = collections.defaultdict(dict), {}, collections.Counter()
    for r in csv.DictReader(open(raw/'Main text/Fig_2/Data_Fig2b_g.csv')):
        x = {'OUTCAR_0.25': '0.25', 'OUTCAR_0.5': '0.50', 'OUTCAR_0.75': '0.75'}[r['Concentration']]
        site, mu = int(r['vacancy']), r['Oxygen chem pot name']
        e0, eb = float(r['Neutral Defect energy']), float(r['Bulk_energy'])
        assert float(r['Oxygen chem pot']) == MU[mu]
        assert abs(e0 - eb + MU[mu] - float(r['Neutral Formation energy'])) < 1e-9
        if x in bulk:
            assert bulk[x] == eb
        bulk[x] = eb
        if site in vac[x]:
            assert vac[x][site] == e0
        vac[x][site] = e0
        seen[(x, site, mu)] += 1
    assert len(seen) == 3 * 64 * 3 and set(seen.values()) == {541}
    pool['vasp'] = {'bulk_eV': bulk, 'vacancy_eV': {x: {str(s): vac[x][s] for s in sorted(vac[x])} for x in sorted(vac)}}
    # Authors' x = 0.50 site table: E_vac is the O-rich neutral formation energy of the Fig2c_f cell
    table = []
    for r in csv.DictReader(open(raw/'SM material/FIGS12.csv')):
        s = int(r['O_index'])
        assert abs(float(r['E_vac']) - (vac['0.50'][s] - bulk['0.50'] + MU['high'])) < 1e-9
        table.append(r)
    assert len(table) == 64
    pool['table_x050'] = table
    # CASTEP multi-functional energies (each functional has its own site labels in the release)
    castep = collections.defaultdict(list)
    for r in csv.DictReader(open(raw/'SM material/FIGS3.csv')):
        castep[r['Source']].append({'site': int(r['Oxygen site']), 'final_enthalpy_eV': float(r['Final Enthalpy']),
                                    'frac_coord': [float(v) for v in r['frac_coord'].strip('[]').split()],
                                    'n_Hf': int(r['Hf Neighbors']), 'n_Zr': int(r['Zr Neighbors']),
                                    'vacancy_type': r['vacancy type'], 'local_environment': r['local environment']})
    assert sorted(castep) == ['B3LYP', 'LDA', 'PBE', 'PBEsol', 'RPBE', 'RSCAN', 'WC'] and all(len(v) == 64 for v in castep.values())
    pool['castep'] = dict(castep)
    rows = [l.split() for l in open(raw/'SM material/FIGS4.txt') if l.strip() and l.split()[0].isdigit()]
    pool['castep_333'] = {s: float(e) for s, e in rows}
    assert len(pool['castep_333']) == 45
    freq = list(csv.DictReader(open(raw/'SM material/FIGS5.txt')))
    pool['figs5'] = {r['site_type']: [float(r['frequency_25']), float(r['frequency_50']), float(r['frequency_75'])] for r in freq}
    (HERE/'pool.json').write_text(json.dumps(pool) + '\n')
    print('pool.json:', {x: len(v) for x, v in pool['vasp']['vacancy_eV'].items()}, 'castep', len(castep), 'x 64; 333:', 45)


if __name__ == '__main__':
    main()
