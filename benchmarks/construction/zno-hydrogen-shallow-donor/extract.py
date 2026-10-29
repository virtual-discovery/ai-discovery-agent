"""Parse the released VASP runs of ZnO_ubiquitous into pool.json (energies, settings, eigenvalues,
core potentials, magnetization, structures).

    python3 -B extract.py [RAW]   (default ~/.cache/ai-theorist/zno-hydrogen-shallow-donor/raw)
"""
import json, os, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = Path(sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/.cache/ai-theorist/zno-hydrogen-shallow-donor/raw'))/'ZnO_ubiquitous'


def read_poscar(path):
    lines = Path(path).read_text().splitlines()
    scale = float(lines[1].split()[0])
    lat = [[float(x) * scale for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    i = 7
    if lines[i].strip()[0] in 'sS':
        i += 1
    direct = lines[i].strip()[0] in 'dD'
    pos = [[float(x) for x in lines[i + 1 + k].split()[:3]] for k in range(sum(counts))]
    if not direct:
        import numpy as np
        pos = (np.array(pos) @ np.linalg.inv(np.array(lat))).tolist()
    sym = [s for s, c in zip(species, counts) for _ in range(c)]
    return {'lattice_A': lat, 'species': sym, 'frac': [[round(v, 6) for v in p] for p in pos]}


def incar(path):
    out = {}
    for line in Path(path).read_text().splitlines():
        line = line.split('#')[0].split('!')[0].strip()
        if '=' in line:
            k, v = line.split('=', 1)
            out[k.strip()] = v.strip()
    return out


def parse_outcar(path):
    text = Path(path).read_text(errors='ignore')
    lines = text.splitlines()
    r = {}
    m = re.findall(r'free  energy   TOTEN  =\s+(-?[\d.]+)', text)
    r['toten_eV'] = float(m[-1])
    m = re.findall(r'energy\(sigma->0\) =\s+(-?[\d.]+)', text)
    r['energy_sigma0_eV'] = float(m[-1])
    r['ionic_steps'] = len(m)
    r['nelect'] = float(re.search(r'NELECT =\s+([\d.]+)', text).group(1))
    r['ispin'] = int(re.search(r'ISPIN  =\s+(\d)', text).group(1))
    r['encut_eV'] = float(re.search(r'ENCUT  =\s+([\d.]+)', text).group(1))
    aexx = re.search(r'AEXX\s*=\s*([\d.]+)', text)
    lhf = re.search(r'LHFCALC\s*=\s*([TF])', text)
    r['lhfcalc'] = bool(lhf and lhf.group(1) == 'T')
    r['aexx'] = float(aexx.group(1)) if aexx else None
    ldau = re.search(r'LDAU\s*=\s*([TF])', text)
    r['ldau'] = bool(ldau and ldau.group(1) == 'T')
    mags = re.findall(r'number of electron\s+[\d.]+\s+magnetization\s+(-?[\d.]+)', text)
    r['magnetization_muB'] = float(mags[-1]) if mags else None
    efs = re.findall(r'E-fermi :\s+(-?[\d.]+)', text)
    r['efermi_eV'] = float(efs[-1])
    # k-points and weights
    i = text.index('k-points in reciprocal lattice and weights')
    kp = []
    for line in text[i:].splitlines()[1:]:
        p = line.split()
        if len(p) != 4:
            break
        kp.append([float(x) for x in p])
    r['kpoints'] = [k[:3] for k in kp]
    r['kweights'] = [k[3] for k in kp]
    # last eigenvalue block
    j = max(n for n, l in enumerate(lines) if l.strip().startswith('E-fermi :'))
    eig = {}
    spin, k = 1, None
    for line in lines[j + 1:]:
        s = line.strip()
        if s.startswith('spin component'):
            spin = int(s.split()[2]); k = None; continue
        if s.startswith('k-point'):
            k = int(s.split()[1]) - 1; eig.setdefault(spin, {})[k] = []; continue
        if s.startswith('band No.'):
            continue
        p = s.split()
        if k is not None and len(p) == 3:
            try:
                eig[spin][k].append([float(p[1]), float(p[2])])
                continue
            except ValueError:
                pass
        if s.startswith('---') or s.startswith('soft charge') or 'total charge' in s or 'magnetization (x)' in s:
            break
        if not s and k is not None and len(eig[spin][k]) > 0 and spin == r['ispin'] and k == len(kp) - 1:
            break
    r['eigenvalues'] = {str(sp): [eig[sp][kk] for kk in sorted(eig[sp])] for sp in eig}
    # last core-potential block
    j = max(n for n, l in enumerate(lines) if 'average (electrostatic) potential at core' in l)
    vals = []
    for line in lines[j + 3:]:
        p = line.split()
        if not p:
            if vals:
                break
            continue
        try:
            vals += [float(p[q]) for q in range(1, len(p), 2)]
        except ValueError:
            break
    r['core_potentials_eV'] = vals
    # final magnetization per atom (tot column)
    js = [n for n, l in enumerate(lines) if l.strip() == 'magnetization (x)']
    if js:
        mom = []
        for line in lines[js[-1] + 4:]:
            p = line.split()
            if not p or p[0].startswith('---'):
                if mom: break
                continue
            if p[0] == 'tot':
                break
            mom.append(float(p[-1]))
        r['site_moments_muB'] = mom
    r['converged_ionic'] = 'reached required accuracy' in text
    return r


def parse_log(path):
    text = Path(path).read_text(errors='ignore')
    m = re.findall(r'F= ([-.E+\d]+) E0= ([-.E+\d]+)', text)
    mags = re.findall(r'mag=\s+(-?[\d.]+)', text)
    return {'toten_eV': float(m[-1][0]), 'energy_sigma0_eV': float(m[-1][1]), 'ionic_steps': len(m),
            'magnetization_muB': float(mags[-1]) if mags else None,
            'converged_ionic': 'reached required accuracy' in text}


def run_dir(d):
    d = Path(d)
    rec = {'path': str(d.relative_to(RAW))}
    if (d/'OUTCAR').exists():
        rec.update(parse_outcar(d/'OUTCAR'))
        rec['source_file'] = 'OUTCAR'
    else:
        rec.update(parse_log(d/'log'))
        rec['source_file'] = 'log (OUTCAR > 12 MiB not extracted)'
        inc = incar(d/'INCAR')
        rec['nelect'] = float(inc['NELECT']) if 'NELECT' in inc else None
    rec['incar'] = incar(d/'INCAR')
    rec['kpoints_file'] = (d/'KPOINTS').read_text().strip().splitlines()
    rec['potcar'] = [l.split('=')[1].strip() for l in (d/'POTCAR.spec').read_text().splitlines() if 'TITEL' in l]
    rec['initial'] = read_poscar(d/'POSCAR')
    rec['final'] = read_poscar(d/'CONTCAR')
    return rec


def main():
    pool = {}
    for d in sorted(RAW.glob('chem_pot/*')) + [RAW/'supercell_3_3_2'] + sorted(RAW.glob('supercell_3_3_2/*/q_*')):
        if not d.is_dir():
            continue
        for sub, level in ((d, 'PBE'), (d/'HSE06', 'HSE06')):
            if (sub/'INCAR').exists():
                key = str(d.relative_to(RAW)) + '/' + level
                pool[key] = run_dir(sub)
                print(key, pool[key]['energy_sigma0_eV'], pool[key].get('nelect'), pool[key]['ionic_steps'], pool[key]['source_file'])
    (HERE/'pool.json').write_text(json.dumps(pool))


if __name__ == '__main__':
    main()
