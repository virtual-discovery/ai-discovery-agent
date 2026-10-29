"""Parse every OUTCAR (+CONTCAR) of the archive zip into pool.json.gz (read in place, no unpacking).

    python3 -B extract.py [zip]   (run inside the runtime image)

Per configuration: INCAR-type settings, k mesh, and every ionic frame (energy(sigma->0), free energy TOTEN,
stress in kB (VASP order XX YY ZZ XY YZ ZX, positive = compressive), cell, volume, Cartesian positions, max force,
number of electronic steps, Cartesian forces), plus the CONTCAR text.
"""
import gzip
import json
import re
import sys
import zipfile
from pathlib import Path

ZIP = Path(sys.argv[1] if len(sys.argv) > 1 else
           Path.home()/'.cache/ai-theorist/ge-pressure-polymorph-pathways/allfiles.zip')
OUT = Path.home()/'.cache/ai-theorist/ge-pressure-polymorph-pathways/pool.json.gz'

PARAMS = {
    'ENCUT': r'ENCUT\s*=\s*([-\d.Ee+]+)\s*eV', 'ISMEAR': r'ISMEAR\s*=\s*(-?\d+);', 'SIGMA': r'SIGMA\s*=\s*([-\d.Ee+]+)',
    'PSTRESS': r'PSTRESS=\s*([-\d.Ee+]+)', 'ISIF': r'ISIF\s*=\s*(-?\d+)', 'IBRION': r'IBRION\s*=\s*(-?\d+)',
    'NSW': r'NSW\s*=\s*(\d+)', 'NELM': r'NELM\s*=\s*(\d+)', 'EDIFF': r'EDIFF\s*=\s*([-\d.Ee+]+)',
    'EDIFFG': r'EDIFFG\s*=\s*([-\d.Ee+]+)', 'PREC': r'\bPREC\s*=\s*([A-Za-z]+)', 'NIONS': r'NIONS\s*=\s*(\d+)',
    'NKPTS': r'NKPTS\s*=\s*(\d+)', 'ISPIN': r'ISPIN\s*=\s*(\d+)', 'POTIM': r'POTIM\s*=\s*([-\d.Ee+]+)',
    'IMAGES': r'IMAGES\s*=\s*(\d+)', 'SPRING': r'SPRING\s*=\s*([-\d.Ee+]+)', 'TEBEG': r'TEBEG\s*=\s*([-\d.Ee+]+)',
}
FLOAT = re.compile(r'-?\d+\.\d*(?:[Ee][-+]?\d+)?')


def num(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s


def parse(text):
    lines = text.splitlines()
    head = text[:60000]
    settings = {}
    for key, pat in PARAMS.items():
        m = re.search(pat, head)
        if m:
            settings[key] = num(m.group(1))
    m = re.search(r'generate k-points for:\s+(\d+)\s+(\d+)\s+(\d+)', head)
    settings['kmesh'] = [int(m.group(i)) for i in (1, 2, 3)] if m else None
    m = re.search(r'POTCAR:\s+(.+)', head)
    settings['POTCAR'] = m.group(1).strip() if m else None
    m = re.search(r' vasp\.([\d.]+)', head)
    settings['vasp'] = m.group(1) if m else None
    n = settings.get('NIONS', 0)
    frames, cur, nelec = [], {}, 0
    i = 0
    while i < len(lines):
        L = lines[i]
        if L.startswith('----') and 'Iteration' in L:
            nelec += 1
        elif L.strip().startswith('in kB'):
            cur['stress_kB'] = [float(x) for x in FLOAT.findall(L)[:6]]
        elif 'direct lattice vectors' in L and 'cell' not in cur and 'stress_kB' in cur:
            cur['cell'] = [[float(x) for x in FLOAT.findall(lines[i + k])[:3]] for k in (1, 2, 3)]
        elif L.strip().startswith('volume of cell') and 'stress_kB' in cur:
            cur['volume'] = float(L.split()[-1])
        elif L.startswith(' POSITION') and 'TOTAL-FORCE' in L:
            block = [FLOAT.findall(lines[i + 2 + k]) for k in range(n)]
            cur['positions'] = [[float(x) for x in b[:3]] for b in block]
            cur['forces'] = [[float(x) for x in b[3:6]] for b in block]
            cur['max_force'] = max((float(b[3]) ** 2 + float(b[4]) ** 2 + float(b[5]) ** 2) ** 0.5 for b in block)
            i += 2 + n
        elif 'free  energy   TOTEN' in L:
            cur['toten'] = float(L.split()[-2])
        elif 'energy(sigma->0)' in L and 'toten' in cur:
            cur['energy'] = float(L.split()[-1])
            cur['n_electronic'] = nelec
            frames.append(cur)
            cur, nelec = {}, 0
        i += 1
    return {'settings': settings, 'frames': frames,
            'reached_accuracy': 'reached required accuracy' in text,
            'terminated': 'General timing and accounting' in text, 'n_lines': len(lines)}


def main():
    pool = {}
    with zipfile.ZipFile(ZIP) as z:
        names = [x for x in z.namelist() if x.endswith('/OUTCAR')]
        for k, name in enumerate(sorted(names)):
            folder = name[len('allfiles/'):-len('/OUTCAR')]
            rec = parse(z.read(name).decode('latin-1'))
            cname = name[:-6] + 'CONTCAR'
            try:
                rec['contcar'] = z.read(cname).decode('latin-1')
            except KeyError:
                rec['contcar'] = None
            pool[folder] = rec
            if k % 500 == 0:
                print(k, folder, len(rec['frames']), flush=True)
    with gzip.open(OUT, 'wt') as f:
        json.dump(pool, f)
    print('configs', len(pool))


if __name__ == '__main__':
    main()
