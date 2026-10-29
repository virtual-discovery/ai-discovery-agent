"""Extract the MnF2 spin-spiral pool from the Materials Cloud archive 3p6fv-8hd59 (MnF2.tar.gz).

    python3 benchmarks/construction/heis-configdesign-mnf2-altermagnet/extract.py [CACHE_DIR]

Reads the authors' spin-spiral CSV and, for every spiral, its VASP OUTCAR (final energy(sigma->0) of the
6-atom cell; the plane-wave cutoff VASP asks for to represent the spiral), OSZICAR (final in-cell
constrained moments M_int of the two Mn), INCAR (MAGMOM, QSPIRAL, cutoffs) and POSCAR. Writes pool.json.

Checks performed here:
  * every CSV energy equals 1000*(E_cell(q) - E_cell(q=0)) of the matching OUTCAR to within 1.5e-3 meV,
    i.e. the CSV 'meV/atom' energies are in fact meV per 6-atom primitive cell (2 Mn);
  * the CSV moments equal the OSZICAR M_int of the final step.
"""
import csv
import json
import re
import sys
from pathlib import Path
import numpy as np

CACHE = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/heis-configdesign-mnf2-altermagnet')
SRC = CACHE/'MnF2/spin-spiral'
OUT = Path(__file__).with_name('pool.json')
DIRS = {'01-001': '[001]', '03-111': '[111]', '05-113': '[113]', '07-121': '[121]', '09-214': '[214]'}


def incar(path):
    out = {}
    for line in path.read_text().splitlines():
        if '=' in line:
            k, v = line.split('=', 1)
            out[k.strip()] = v.strip()
    return out


def outcar(path):
    text = path.read_text()
    energy = float(re.findall(r'energy\(sigma->0\) =\s+([-\d.]+)', text)[-1])
    need = re.findall(r'choose ENMAX >\s*([\d.e+]+)\s*eV', text.replace('\n', ' '))
    return energy, float(need[-1]) if need else None


def outcar_warnings(path):
    """Verbatim WARNING / ADVICE texts printed in the OUTCAR header (before the first ionic step)."""
    lines = path.read_text().split('\n')
    out, i = [], 0
    stop = next(k for k, l in enumerate(lines) if 'POSCAR:' in l)
    while i < stop:
        l = lines[i]
        if 'ADVICE to this user' in l:
            body = []
            i += 1
            while i < stop and not lines[i].strip().startswith('-----'):
                t = lines[i].strip().strip('|').strip()
                if t:
                    body.append(t)
                i += 1
            out.append('ADVICE: ' + ' '.join(body))
        elif l.strip().startswith('WARNING:'):
            body = [l.strip()]
            i += 1
            while i < stop and lines[i].strip() and not lines[i].strip().startswith(('WARNING', 'Optimization', 'LDA', 'POSCAR')):
                body.append(lines[i].strip())
                i += 1
            out.append(' '.join(body))
            continue
        i += 1
    return out


def scf_steps(path):
    return sum(1 for l in path.read_text().splitlines() if l.startswith(('DAV:', 'RMM:')))


def oszicar_moments(path, natoms=2):
    lines = path.read_text().splitlines()
    starts = [i for i, l in enumerate(lines) if l.split()[:3] == ['ion', 'MW_int', 'M_int']]
    block = lines[starts[-1] + 1: starts[-1] + 1 + natoms]
    return [[float(x) for x in l.split()[4:7]] for l in block]


def poscar(path):
    lines = path.read_text().split('\n')
    scale = float(lines[1])
    cell = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]) * scale
    counts = [int(x) for x in lines[6].split()]
    frac = np.array([[float(x) for x in lines[8 + i].split()[:3]] for i in range(sum(counts))])
    return cell, frac


def main():
    rows = []
    with open(SRC/'spin-spiral-energies.csv') as fh:
        reader = csv.reader(fh)
        next(reader)
        for r in reader:
            v = [float(x) for x in r]
            rows.append({'q': v[0:3], 'energy_csv': v[3], 'm1': v[4:7], 'mag1': v[7], 'theta1': v[8],
                         'm2': v[9:12], 'mag2': v[12], 'theta2': v[13]})
    calcs = []
    for d, label in DIRS.items():
        for qdir in sorted((SRC/d).glob('q-*')):
            inc = incar(qdir/'INCAR')
            q = [float(x) for x in inc['QSPIRAL'].split()]
            energy, need = outcar(qdir/'OUTCAR')
            calcs.append({'direction': label, 't': float(qdir.name[2:]), 'q': q, 'E_eV': energy, 'enmax_needed_eV': need,
                          'incar': inc, 'warnings': outcar_warnings(qdir/'OUTCAR'), 'scf_steps': scf_steps(qdir/'OSZICAR'),
                          'moments_all': oszicar_moments(qdir/'OSZICAR', natoms=6),
                          'M_int': oszicar_moments(qdir/'OSZICAR'), 'magmom': [float(x) for x in inc['MAGMOM'].split()[:6]],
                          'encut': float(inc['ENCUT']), 'enini': float(inc['ENINI']),
                          'kmesh': (qdir/'KPOINTS').read_text().split('\n')[3].split(), 'path': str(qdir.relative_to(CACHE))})
    e0 = [c['E_eV'] for c in calcs if np.allclose(c['q'], 0)]
    assert len(e0) == 2 and abs(e0[0] - e0[1]) < 1e-5, e0
    cell, frac = poscar(SRC/'01-001/q-0.00/POSCAR')
    spirals = []
    for r in rows:
        match = [c for c in calcs if np.allclose(c['q'], r['q'], atol=1e-6)]
        assert match, r['q']
        c = match[0]
        if np.allclose(r['q'], 0):
            c = [m for m in match if m['direction'] == '[001]'][0]
        de = 1000 * (c['E_eV'] - e0[0])
        assert abs(de - r['energy_csv']) < 1.5e-3, (r['q'], de, r['energy_csv'])
        assert np.allclose(c['M_int'][0], r['m1'], atol=1.5e-3) and np.allclose(c['M_int'][1], r['m2'], atol=1.5e-3), r['q']
        assert abs(r['theta1'] - 90) < 1e-6 and abs(r['theta2'] - 90) < 1e-6
        spirals.append({'direction': 'gamma' if np.allclose(r['q'], 0) else c['direction'], 't': c['t'],
                        'q_reciprocal_lattice': r['q'], 'total_energy_eV': c['E_eV'],
                        'energy_meV_per_cell': round(de, 4), 'energy_csv_meV': r['energy_csv'],
                        'mn_moments_muB': [r['m1'], r['m2']], 'moment_magnitudes_muB': [r['mag1'], r['mag2']],
                        'enmax_needed_eV': c['enmax_needed_eV'], 'encut_eV': c['encut'], 'enini_eV': c['enini'],
                        'kmesh': c['kmesh'], 'magmom_in': c['magmom'], 'incar': c['incar'], 'outcar_warnings': c['warnings'],
                        'scf_steps': c['scf_steps'], 'moments_all_atoms_muB': c['moments_all'], 'source_path': c['path']})
    OUT.write_text(json.dumps({'source': 'Materials Cloud 3p6fv-8hd59 MnF2.tar.gz (md5 4311cefc120956738f6d26440354f3fd)',
                               'cell_A': cell.round(10).tolist(), 'species': ['Mn', 'Mn', 'F', 'F', 'F', 'F'],
                               'fractional_positions': frac.round(12).tolist(), 'spirals': spirals}, indent=1))
    over = [(s['direction'], s['t'], s['enmax_needed_eV']) for s in spirals if s['enmax_needed_eV'] and s['enmax_needed_eV'] > s['encut_eV']]
    print(len(spirals), 'spirals ->', OUT, '; basis-truncated (ENMAX needed > ENCUT):', over)


if __name__ == '__main__':
    main()
