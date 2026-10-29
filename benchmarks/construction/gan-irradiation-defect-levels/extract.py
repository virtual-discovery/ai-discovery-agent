"""Extract every calculation of the GaN defect archive into pool.json.

Source: Zenodo record 21893631 (gan_defects.zip, md5 85d4a96b1bf0bcf032c89bf92dbb021f), unpacked by fetch.py into
~/.cache/ai-theorist/gan-irradiation-defect-levels/raw. Run in the runtime image:

    python3 -B extract.py [RAW_DIR]

For each calculation folder it records the INCAR tags printed in the OUTCAR, the POTCAR titles, the k-point,
NELECT/NUPDOWN, the ionic-step count and convergence flags, the final energies, E-fermi, the final cell and
positions (CONTCAR), the final Kohn-Sham eigenvalues and occupations of both spin channels, and the final
"average (electrostatic) potential at core" of every atom. No POTCAR content is copied.
"""
import json
import re
import sys
from pathlib import Path

DEFAULT = Path.home()/'.cache/ai-theorist/gan-irradiation-defect-levels/raw'
OUT = Path(__file__).with_name('pool.json')


def incar_tags(text):
    block = text[text.index(' INCAR:'):text.index(' POTCAR:', text.index(' INCAR:'))]
    tags = {}
    for line in block.splitlines()[1:]:
        if '=' in line:
            k, v = line.split('=', 1)
            tags[k.strip()] = v.strip()
    return tags


def last_block(text, header):
    i = text.rfind(header)
    return None if i < 0 else text[i:]


def eigenvalues(text, nspin):
    """Final eigenvalue/occupation table for every spin channel (single k-point)."""
    out = []
    for s in range(1, nspin + 1):
        hits = [m.end() for m in re.finditer(r'spin component %d\s*\n\s*\n\s*k-point\s+1 :[^\n]*\n\s*band No\.[^\n]*\n' % s, text)]
        rows = []
        for line in text[hits[-1]:].splitlines():
            parts = line.split()
            if len(parts) != 3 or not parts[0].isdigit():
                break
            rows.append([float(parts[1]), float(parts[2])])
        out.append(rows)
    return out


def core_potentials(text, nions):
    blk = last_block(text, 'average (electrostatic) potential at core')
    nums = re.findall(r'(\d+)\s*(-?\d+\.\d+)', blk.split('\n', 3)[3])
    pots = {}
    for idx, val in nums:
        idx = int(idx)
        if idx in pots or idx > nions:
            break
        pots[idx] = float(val)
    assert len(pots) == nions, (len(pots), nions)
    return [pots[i] for i in range(1, nions + 1)]


def site_moments(text, nions):
    """Final per-atom spin moments (LORBIT=10 'magnetization (x)' table, column tot), or None."""
    blk = last_block(text, 'magnetization (x)')
    if blk is None:
        return None
    rows = {}
    for line in blk.splitlines()[4:]:
        parts = line.split()
        if not parts or not parts[0].isdigit():
            break
        rows[int(parts[0])] = float(parts[-1])
    return [rows[i] for i in range(1, nions + 1)] if len(rows) == nions else None


def poscar(path):
    lines = path.read_text().splitlines()
    scale = float(lines[1].split()[0])
    cell = [[float(x) * scale for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    k = 7
    if lines[k].strip()[0] in 'sS':
        k += 1
    direct = lines[k].strip()[0] in 'dD'
    nat = sum(counts)
    pos = [[float(x) for x in lines[k + 1 + i].split()[:3]] for i in range(nat)]
    symbols = [s for s, c in zip(species, counts) for _ in range(c)]
    return {'cell_A': cell, 'symbols': symbols, 'frac': pos if direct else None}


def parse(folder):
    text = (folder/'OUTCAR').read_text(errors='ignore')
    tags = incar_tags(text)
    potcars = re.findall(r'POTCAR:\s+(PAW_PBE \S+ \S+)', text)
    nions = int(re.search(r'NIONS\s*=\s*(\d+)', text).group(1))
    nelect = float(re.search(r'NELECT\s*=\s*([\d.]+)\s+total number of electrons', text).group(1))
    nspin = int(re.search(r'ISPIN\s*=\s*(\d)', text).group(1))
    kline = text[text.index('k-points in reciprocal lattice and weights'):].splitlines()[1].split()
    energies = re.findall(r'free  energy   TOTEN  =\s+(-?[\d.]+)', text)
    e0 = re.findall(r'energy\(sigma->0\) =\s+(-?[\d.]+)', text)
    efermi = re.findall(r'E-fermi :\s+(-?[\d.]+)', text)
    mags = re.findall(r'number of electron\s+[\d.]+\s+magnetization\s+(-?[\d.]+)', text)
    # electronic steps per ionic step (from OSZICAR) to detect SCF that hit NELM
    osz = (folder/'OSZICAR').read_text().splitlines()
    nelm = int(tags.get('NELM', '60'))
    scf_counts, count = [], 0
    for line in osz:
        if re.match(r'\s*(DAV|RMM|CGA|CG):', line):
            count += 1
        elif re.match(r'\s*\d+ F=', line):
            scf_counts.append(count)
            count = 0
    contcar = poscar(folder/'CONTCAR')
    start = poscar(folder/'POSCAR')
    return {
        'incar': tags, 'potcar': sorted(set(potcars)), 'nions': nions, 'nelect': nelect, 'ispin': nspin,
        'kpoint': [float(x) for x in kline[:3]],
        'nupdown': float(re.search(r'NUPDOWN=\s*(-?[\d.]+)', text).group(1)),
        'ionic_steps': len(energies), 'reached_accuracy': 'reached required accuracy' in text,
        'scf_steps_last': scf_counts[-1] if scf_counts else None, 'nelm': nelm,
        'scf_hit_nelm': any(c >= nelm for c in scf_counts[-1:]),
        'toten_eV': float(energies[-1]), 'energy_sigma0_eV': float(e0[-1]), 'efermi_eV': float(efermi[-1]),
        'magnetization_muB': float(mags[-1]) if mags else 0.0,
        'cell_A': contcar['cell_A'], 'symbols': contcar['symbols'], 'frac_final': contcar['frac'],
        'frac_start': start['frac'], 'start_symbols': start['symbols'],
        'eigenvalues': eigenvalues(text, nspin), 'core_potentials_V': core_potentials(text, nions),
        'site_moments_muB': site_moments(text, nions),
    }


def main():
    raw = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT
    pool = {}
    for outcar in sorted(raw.rglob('OUTCAR')):
        rec = parse(outcar.parent)
        rec['name'] = str(outcar.parent.relative_to(raw))
        pool[rec['name']] = rec
    OUT.write_text(json.dumps(pool))
    print(len(pool), 'calculations ->', OUT, f'({OUT.stat().st_size / 1e6:.1f} MB)')


if __name__ == '__main__':
    main()
