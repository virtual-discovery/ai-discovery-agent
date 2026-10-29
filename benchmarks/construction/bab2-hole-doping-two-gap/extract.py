"""Extract the BaB2 archive content used by the task into pool.json (run after fetch.py).

    python3 -B benchmarks/construction/bab2-hole-doping-two-gap/extract.py [DEST]

DEST is the fetch directory (default ~/.cache/ai-theorist/bab2-hole-doping-two-gap). Everything in
pool.json is archived content: input files verbatim, selected output lines verbatim, and numbers parsed
from the archived files without transformation. Folders keep the archive's names; the doping of each
folder is read from its scf.in (tot_charge) and asserted against the other headers.
"""
import json
import re
import sys
from pathlib import Path

DEST = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/bab2-hole-doping-two-gap')
ROOT = DEST/'x/bab2-analysis'
OUT = Path(__file__).resolve().with_name('pool.json')
FREQ = re.compile(r'freq \(\s*\d+\)\s*=\s*[-\d.]+ \[THz\] =\s*([-\d.]+) \[cm-1\]')
QLINE = re.compile(r'q = \(\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\)')


def text(path):
    return Path(path).read_text(errors='replace')


def lines_matching(path, patterns):
    return [line.rstrip() for line in text(path).splitlines() if any(p in line for p in patterns)]


def last_cell(relax_out):
    t = text(relax_out)
    block = t[t.rfind('CELL_PARAMETERS'):].splitlines()
    cell = [[float(v) for v in block[i].split()[:3]] for i in (1, 2, 3)]
    pos = t[t.rfind('ATOMIC_POSITIONS'):].splitlines()
    atoms = []
    for line in pos[1:]:
        p = line.split()
        if len(p) < 4 or p[0] == 'End':
            break
        atoms.append([p[0]] + [float(v) for v in p[1:4]])
    return {'header': block[0].strip(), 'cell': cell, 'positions_header': pos[0].strip(), 'positions': atoms}


def dyn(path):
    """q-point (as printed, cartesian, units 2 pi/alat) and harmonic frequencies (cm-1) of a dynamical-matrix file."""
    t = text(path)
    diag = t[t.find('Diagonalizing the dynamical matrix'):]
    q = [float(v) for v in QLINE.search(diag).groups()]
    freqs = [float(v) for v in FREQ.findall(diag)]
    return q, freqs


def dyn_header(path):
    lines = text(path).splitlines()
    return '\n'.join(lines[:lines.index(next(l for l in lines if 'Dynamical  Matrix in cartesian axes' in l))]).rstrip()


def a2f(path):
    lines = text(path).splitlines()
    rows, footer = [], []
    for line in lines[1:]:
        p = line.split()
        if len(p) == 21 and not rows or (rows and len(p) == 21 and not footer):
            rows.append([float(v) for v in p])
        else:
            footer.append(line.rstrip())
    return {'header': lines[0].rstrip(), 'rows': rows, 'footer': footer}


def iso_log(path):
    t = text(path)
    start = t.find('     ===================================================================\n     Solve isotropic')
    body = t[start:]
    end = body.find('     Unfolding on the coarse grid')
    body = body if end < 0 else body[:end]
    return '\n'.join(line.rstrip() for line in body.splitlines() if line.strip())


def aniso_log(path):
    """Per-temperature log lines of the archived anisotropic run (imaginary axis then Pade continuation)."""
    keys = ('temp(', 'Total number of frequency points', 'Convergence was', 'Chemical potential', 'Free energy',
            'Min. / Max. values of superconducting gap', 'Pade approximant', 'WARNING', 'tempsmax',
            'Electron-phonon coupling strength', 'Allen-Dynes', 'w_log', 'BCS superconducting gap', 'machine learning')
    out, current = {}, None
    head = []
    for line in text(path).splitlines():
        m = re.match(r'\s*temp\(\s*\d+\)\s*=\s*([\d.]+) K', line)
        if m:
            current = float(m.group(1))
            out[current] = []
        if any(k in line for k in keys):
            (out[current] if current is not None else head).append(line.rstrip())
    return head, out


def gap_file(path):
    lines = text(path).splitlines()
    return {'header': lines[0].rstrip(), 'rows': [[float(v) for v in l.split()] for l in lines[1:] if l.strip()]}


def run(folder):
    d = ROOT/'sensitivity'/folder
    scf_in = text(d/'epw/scf.in')
    x = float(re.search(r'tot_charge\s*=\s*([\d.]+)', scf_in).group(1))
    assert x == float(re.search(r'tot_charge\s*=\s*([\d.]+)', text(d/'ph/scf.in')).group(1))
    assert x == float(re.search(r'tot_charge\s*=\s*([\d.]+)', text(d/'relax/relax.in')).group(1))
    ef_iso = float(re.search(r'fermi_energy\s*=\s*([\d.]+)', text(d/'epw/iso.in')).group(1))
    ef_scf = float(re.findall(r'the Fermi energy is\s+([\d.]+)', text(d/'epw/scf.out'))[-1])
    assert abs(ef_iso - ef_scf) < 1e-4, (folder, ef_iso, ef_scf)
    nelec = float(re.search(r'number of electrons\s*=\s*([\d.]+)', text(d/'epw/scf.out')).group(1))
    assert abs(nelec - (16 - x)) < 0.0051, (folder, nelec, x)
    dfpt = []
    for n in range(1, 22):
        q, freqs = dyn(d/f'ph/bab2.dyn{n}')
        dfpt.append({'file': f'bab2.dyn{n}', 'q_cartesian_2pi_over_alat': q, 'frequencies_cm1': freqs})
    head, per_t = aniso_log(d/'epw/epw_int.out')
    aniso = {}
    for path in sorted((d/'epw').glob('bab2.imag_aniso_gap0_*')):
        T = float(path.name.split('_')[-1])
        aniso[f'{T:g}'] = {'imaginary': gap_file(path), 'real': gap_file(d/'epw'/path.name.replace('imag', 'pade')),
                           'log': per_t[T]}
    return {'folder': folder, 'tot_charge': x, 'scf_in': scf_in.rstrip(), 'ph_in': text(d/'ph/ph.in').rstrip(),
            'scf_out_lines': lines_matching(d/'epw/scf.out', ['number of electrons', 'the Fermi energy is',
                                                              'Exchange-correlation', 'number of k points=',
                                                              'kinetic-energy cutoff', 'unit-cell volume']),
            'relaxed_structure': last_cell(d/'relax/relax.out'), 'dfpt_header': dyn_header(d/'ph/bab2.dyn1'),
            'dfpt_grid': text(d/'ph/bab2.dyn0').rstrip(), 'dfpt': dfpt, 'a2f': a2f(d/'epw/bab2.a2f'),
            'iso_in': text(d/'epw/iso.in').rstrip(), 'iso_log': iso_log(d/'epw/iso.out'),
            'aniso_in': text(d/'epw/epw_int.in').rstrip(), 'aniso_log_head': head, 'aniso': aniso}


def supercell(name, scf_in, scf_out, relax=None):
    d = ROOT/name
    out = {'run': name, 'scf_in': text(d/scf_in).rstrip(), 'ph_in': text(d/'ph.in').rstrip(),
           'scf_out_lines': lines_matching(d/scf_out, ['number of electrons', 'the Fermi energy is',
                                                       'Exchange-correlation', 'number of atoms/cell',
                                                       'kinetic-energy cutoff', 'unit-cell volume']),
           'dfpt_grid': text(d/'bab2.dyn0').rstrip(), 'dfpt_header': dyn_header(d/'bab2.dyn1'), 'dfpt': []}
    for n in (1, 2):
        q, freqs = dyn(d/f'bab2.dyn{n}')
        out['dfpt'].append({'file': f'bab2.dyn{n}', 'q_cartesian_2pi_over_alat': q, 'frequencies_cm1': freqs})
    if relax:
        out['relax_in'] = text(d/relax).rstrip()
        out['relaxed_structure'] = last_cell(d/'relax.out')
    return out


def main():
    readme = text(DEST/'README.md')
    pool = {'source': 'Materials Cloud Archive 0kqen-w8h90, bab2-analysis.zip (md5 a63a59c6c913b4cab7d5050eb0f1710c)',
            'readme_bab2': readme[readme.find('`bab2-analysis.zip`'):readme.find('## Fig. 1 tables')].strip(),
            'runs': {f: run(f) for f in '0123456'},
            'supercells': {'331': supercell('331', 'scf_331.in', 'scf_331.out'),
                           '331_cs': supercell('331_cs', 'scf.in', 'scf.out', relax='relax.in')}}
    OUT.write_text(json.dumps(pool, separators=(',', ':')) + '\n')
    print({f: r['tot_charge'] for f, r in pool['runs'].items()}, f'-> {OUT} ({OUT.stat().st_size / 1e6:.2f} MB)')


if __name__ == '__main__':
    main()
