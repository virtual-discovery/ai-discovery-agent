"""Read the unpacked archive (fetch.sh) into pool.json: every number the task uses, unmodified.

    python3 -B extract.py [DEST]        # DEST defaults to ~/.cache/ai-theorist/mos2-doping-dome-cdw
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEST = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/mos2-doping-dome-cdw'
X = DEST/'x'


def rows(path, skip_hash=True):
    out = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or (skip_hash and line.startswith('#')):
            continue
        out.append([float(t) for t in line.split()])
    return out


def header(path):
    return [l.strip() for l in path.read_text().splitlines() if l.strip().startswith('#')]


def main():
    pool = {'spectra': {}, 'dos': {}, 'source': 'Materials Cloud record 51yyy-ksf72 (CC BY 4.0)'}
    for d in sorted((X/'fig2e/data').iterdir()):
        m = re.match(r'dop_0(\d\d)(_soc)?$', d.name)
        x = '0.' + m.group(1)
        phase = '1x1H-SOC' if m.group(2) else '1x1H'
        for fl in ('A', 'NA', 'FA'):
            f = d/f'800K_a2f_lam_{fl}.txt'
            pool['spectra'][f'{phase}/{x}/{fl}'] = {'file': str(f.relative_to(X)), 'header': header(f), 'rows': rows(f)}
        dos, dos2 = rows(d/'dos.dat'), rows(d/'dos2.dat')
        assert [r[0] for r in dos] == [r[0] for r in dos2]
        pool['dos'][f'{phase}/{x}'] = {'files': [str((d/'dos.dat').relative_to(X)), str((d/'dos2.dat').relative_to(X))],
                                       'energy_eV': [r[0] for r in dos], 'total': [r[1] for r in dos],
                                       'K_dz2': [r[1] for r in dos2], 'Q_rest': [r[2] for r in dos2]}
    for f in sorted((X/'fig4').glob('*_a2f_lam_*.txt')):
        m = re.match(r'(1tprime|2x2cdw_(dop\d))_a2f_lam_(A|NA|FA)\.txt', f.name)
        key = ('1Tprime/archived' if m.group(1) == '1tprime' else f'2x2CDW-high/{m.group(2)}') + '/' + m.group(3)
        pool['spectra'][key] = {'file': str(f.relative_to(X)), 'header': header(f), 'rows': rows(f)}
    pool['moments'] = {'lambda': {'header': header(X/'fig2a-d/lambda_2d.txt'), 'rows': rows(X/'fig2a-d/lambda_2d.txt')},
                       'omega_log_meV': {'header': header(X/'fig2a-d/omegalog_2d.txt'), 'rows': rows(X/'fig2a-d/omegalog_2d.txt')}}
    pool['tc_tables_fig2e'] = {}
    for f in sorted((X/'fig2e').glob('tc_*.dat')):
        lines = f.read_text().splitlines()
        pool['tc_tables_fig2e'][f.stem] = {'columns': lines[0].split(), 'rows': [[float(t) for t in l.split()] for l in lines[1:] if l.strip()]}
    pool['fig5'] = {}
    for name in ('tcs_800K.txt', 'tcs_CDW_800K.txt', 'tcs_1tprime.txt', 'a1g.txt', 'exp_a1g.txt', 'exp_dome_1'):
        f = X/'fig5'/name
        pool['fig5'][name] = {'header': header(f), 'rows': rows(f)}
    pool['fig3'] = {}
    for name in ('1x1_w_Gamma.txt', '2x2_w_Gamma.txt', 'Mo_hysteresis_1x1.txt', 'Mo_hysteresis_cdw.txt'):
        f = X/'fig3'/name
        pool['fig3'][name] = {'header': header(f), 'rows': rows(f)}
    pool['pes'] = {}
    for f in sorted((X/'fig3').glob('pes_*')):
        pool['pes'][f.name] = {'energies_Ry': [float(t) for t in f.read_text().split()],
                               'raw_lines': f.read_text().splitlines()}
    pool['qe_inputs'] = {}
    for f in sorted((X/'input').glob('*/sc.in')):
        text = f.read_text()
        system = text[text.index('&system'):text.index('/', text.index('&system'))]
        pool['qe_inputs'][f.parent.name] = {k: v.strip().strip(',') for k, v in
                                            re.findall(r'^\s*([A-Za-z_()0-9]+)\s*=\s*(.+?)\s*$', system, re.M)}
    (HERE/'pool.json').write_text(json.dumps(pool, separators=(',', ':')) + '\n')
    print('wrote pool.json:', len(pool['spectra']), 'spectra,', len(pool['dos']), 'DOS,', len(pool['pes']), 'PES')


if __name__ == '__main__':
    main()
