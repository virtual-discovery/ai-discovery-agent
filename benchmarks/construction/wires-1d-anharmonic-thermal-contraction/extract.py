"""Extract the archive into pool.json (cache): SSCHA logs (final-population excerpts), dynamical matrices,
DFPT initial guesses, DFT files and the archived strained TaSe3 calculation (from QHA/qha_1D.aiida).

    python3 -B extract.py [CACHE]      # default CACHE: ~/.cache/ai-theorist/wires-1d-anharmonic-thermal-contraction
"""
import io
import json
import re
import sqlite3
import sys
import tempfile
import zipfile
from pathlib import Path
import parselog

CACHE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/wires-1d-anharmonic-thermal-contraction'
ROOT = CACHE/'raw/Extreme-anharmonicity-thermal-contration-1D-wires'
SS = ROOT/'SSCHA'
WIRES = {'CuC2': {'prefix': 'c2cu', 'asr_prefix': 'cuc2_second_asr_', 'nq': 4, 'files': ['scf.in', 'scf.out', 'ph.in', 'ph.out']},
         'TaSe3': {'prefix': 'tase3', 'asr_prefix': 'tase3_second_asr_', 'nq': 4, 'files': ['scf.in', 'scf.out']},
         'AuSe2': {'prefix': 'ause2', 'asr_prefix': 'ause2_second_asr_', 'nq': 6, 'files': ['scf.in', 'ph.in']}}


def dyn_files(wire, T, nconf, kind):
    d = SS/wire/f'T{T}'
    nq = WIRES[wire]['nq']
    if wire == 'AuSe2':
        names = [f'final_sscha_T{T}_{nconf}conf{i}' for i in range(1, nq + 1)] if kind == 'auxiliary' else []
    elif kind == 'auxiliary':
        names = [f'final_sscha_{nconf}conf_T{T}_{i}' for i in range(1, nq + 1)]
    else:
        names = [f'hessian_{nconf}_T{T}_{i}' for i in range(1, nq + 1)]
    if not names or not all((d/n).exists() for n in names):
        return None
    return {n: (d/n).read_text() for n in names}


def log_excerpt(path, T):
    """Verbatim final-population excerpt of one temperature section (configuration-progress lines removed)."""
    lines = [l.rstrip('\n') for l in open(path, errors='ignore') if not l.startswith('Computing configuration')
             and not l.startswith('ENSEMBLE ALL')]
    starts = [i for i, l in enumerate(lines) if re.match(r'^\s+T = (\d+)\s*$', l) and i > 0 and lines[i - 1].startswith('=====')]
    for k, i in enumerate(starts):
        if int(lines[i].split('=')[1]) != T:
            continue
        end = starts[k + 1] - 1 if k + 1 < len(starts) else len(lines)
        sec = lines[i - 1:end]
        res = [j for j, l in enumerate(sec) if '*   RESULTS   *' in l]
        if not res:
            return None
        last = res[-1]
        step = max(j for j, l in enumerate(sec[:last]) if 'NEW MINIMIZATION STEP' in l)
        stop = last
        for j in range(last, len(sec)):
            stop = j
            if 'Ab initio average stress' in sec[j]:
                stop = j + 3
                break
        head = [l for l in sec[:6] if l.strip()]
        conf = [l for l in sec if l.startswith('Force computed shape')][:1]
        crit = [j for j in range(step, last) if 'Running =  False' in sec[j]]
        body = sec[step:crit[0] + 1] if crit else sec[step:step + 25]
        tail = sec[last - 2:stop + 1]
        return {'header': head + conf, 'populations': len(res),
                'final_excerpt': '\n'.join(body + ['[... debug output omitted ...]'] + tail)}
    return None


def main():
    pool = {'wires': {}, 'strained_TaSe3': None}
    for wire, info in WIRES.items():
        w = {'runs': {}, 'dfpt': {}, 'dft_files': {}}
        ig = SS/wire/'initial_guess'
        w['dfpt']['none'] = {f'{info["prefix"]}.dyn{i}': (ig/f'{info["prefix"]}.dyn{i}').read_text() for i in range(0, info['nq'] + 1)}
        w['dfpt']['1d'] = {f'{info["asr_prefix"]}{i}': (ig/f'{info["asr_prefix"]}{i}').read_text() for i in range(1, info['nq'] + 1)}
        for f in info['files']:
            w['dft_files'][f] = (ig/f).read_text()
        logs = {}
        for path in sorted((SS/wire).glob('log*')):
            for sec in parselog.parse_log(path):
                if not sec['results'] or 'stress_GPa' not in sec['results'][-1]:
                    continue
                nconf = sec['n_conf'] or (100000 if '100000' in path.name else 10000)
                key = (sec['T'], nconf)
                if key not in logs:
                    logs[key] = (path, sec)
        temps = sorted(int(p.name[1:]) for p in (SS/wire).glob('T*') if p.is_dir())
        keys = set(logs) | {(T, n) for T in temps for n in (10000, 100000) if dyn_files(wire, T, n, 'auxiliary')}
        for T, n in sorted(keys):
            r = {'T': T, 'configurations': n}
            if (T, n) in logs:
                path, sec = logs[(T, n)]
                ex = log_excerpt(path, T)
                r['log'] = {'file': f'SSCHA/{wire}/{path.name}', **ex, 'parsed': sec['results'][-1]}
            r['auxiliary'] = dyn_files(wire, T, n, 'auxiliary')
            r['hessian'] = dyn_files(wire, T, n, 'hessian')
            w['runs'][f'{T}K/{n}'] = r
        pool['wires'][wire] = w
    # Strained TaSe3 (QHA archive): one pw.x + ph.x calculation at a single strain.
    with zipfile.ZipFile(ROOT/'QHA/qha_1D.aiida') as z, tempfile.TemporaryDirectory() as tmp:
        (Path(tmp)/'db.sqlite3').write_bytes(z.read('db.sqlite3'))
        c = sqlite3.connect(str(Path(tmp)/'db.sqlite3'))
        node = {i: (t, json.loads(a), json.loads(rm)) for i, t, a, rm in
                c.execute('select id, node_type, attributes, repository_metadata from db_dbnode')}
        repo = lambda k: z.read('repo/' + k).decode()
        pw_calc = next(v for k, v in node.items() if v[0].endswith('CalcJobNode.') and 'pw' in v[1]['parser_name'])
        ph_calc = next(v for k, v in node.items() if v[0].endswith('CalcJobNode.') and v[1]['parser_name'].endswith('.ph'))
        pw_folder = next(v for k, v in node.items() if v[0].endswith('FolderData.') and 'data-file-schema.xml' in v[2].get('o', {}))
        ph_folder = next(v for k, v in node.items() if v[0].endswith('FolderData.') and 'DYN_MAT' in v[2].get('o', {}))
        pool['strained_TaSe3'] = {
            'pw_input': repo(pw_calc[2]['o']['aiida.in']['k']),
            'pw_output': repo(pw_folder[2]['o']['aiida.out']['k']),
            'ph_input': repo(ph_calc[2]['o']['aiida.in']['k']),
            'dyn_mat': {n: repo(v['k']) for n, v in ph_folder[2]['o']['DYN_MAT']['o'].items()},
            'structure_label': next(v for k, v in node.items() if v[0].endswith('StructureData.'))[0] and
            next(n[1].get('cell') for k, n in node.items() if n[0].endswith('StructureData.'))}
        labels = [r for r in c.execute("select label from db_dbnode where node_type like '%StructureData%'")]
        pool['strained_TaSe3']['structure_label'] = labels[0][0]
    out = CACHE/'pool.json'
    out.write_text(json.dumps(pool))
    for wire, w in pool['wires'].items():
        print(wire, len(w['runs']), 'runs;', sum(r['hessian'] is not None for r in w['runs'].values()), 'hessians;',
              sum('log' in r for r in w['runs'].values()), 'logs')
    print('pool.json', round(out.stat().st_size / 2**20, 1), 'MB')


if __name__ == '__main__':
    main()
