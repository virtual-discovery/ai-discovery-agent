"""Parse the unpacked archive into pool.json.gz (compact structures + energies).

    python3 -B extract.py [SRC]   (default ~/.cache/ai-theorist/sngeo2-coherent-miscibility/x)

Asserts the archive relations used by the task:
  * the strain-scan minima of each end member equal the SQS end-member energies (same settings);
  * the SQS file holds total energies per 6-atom cell for x = 0, 1 and, for 0 < x < 1, the formation
    energy per rutile primitive cell (2 cations): the subregular fit reproduces the paper's Table 1;
  * the HSE-setting PBE SQS totals (Alloy_PBE_Values) give per-cation formation energies equal to half
    of the released SQS values within 2 meV/cation at x = 0.5.
"""
import gzip, json, os, re, sys
from pathlib import Path
import numpy as np

SRC = Path(sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/.cache/ai-theorist/sngeo2-coherent-miscibility/x'))
HERE = Path(__file__).resolve().parent


def parse_poscar(text):
    L = text.strip('\n').split('\n')
    scale = float(L[1])
    A = np.array([[float(v) for v in L[i].split()] for i in (2, 3, 4)]) * scale
    el = L[5].split(); n = [int(v) for v in L[6].split()]
    el = [e.split('/')[0].split('_')[0] for e in el]
    N = sum(n)
    start = 8 if L[7].strip()[0] in 'DdCc' else 9
    fr = np.array([[float(v) for v in L[start + i].split()[:3]] for i in range(N)])
    sp = [e for e, k in zip(el, n) for _ in range(k)]
    return A, sp, fr


def compact(A, sp, fr):
    return {'lattice_A': np.round(A, 8).tolist(), 'species': sp, 'frac': np.round(fr, 7).tolist()}


def main():
    pool = {'orderings': {}, 'sqs': {}, 'strain': {}, 'hse_setting_pbe': {}}
    for fn in ('PBE', 'r2SCAN'):
        pool['orderings'][fn] = []
        for r in json.load(open(SRC/f'data_CE_{fn}.json')):
            A, sp, fr = parse_poscar(r['structure'])
            pool['orderings'][fn].append({**compact(A, sp, fr), 'energy_eV': r['energy']})
        key = {'PBE': 'PBE', 'r2SCAN': 'r2scan'}[fn]
        sq = json.load(open(SRC/f'data_SQS_{fn}.json'))[key]
        pool['sqs'][fn] = []
        for comp, v in sq.items():
            for i, r in enumerate(v if isinstance(v, list) else [v]):
                A, sp, fr = parse_poscar(r['structure'])
                pool['sqs'][fn].append({'x_Ge': float(comp), 'index': i, **compact(A, sp, fr), 'energy_as_released': r['energy']})
        for plane in ('001', '100', '101'):
            d = json.load(open(SRC/f'data_{plane}_{fn}.json'))
            for m in ('GeO2', 'SnO2'):
                pool['strain'][f'{fn}|{plane}|{m}'] = [{**compact(*parse_poscar(r['structure'])), 'energy_eV': r['energy']} for r in d[m]]
    # ordering order matches across functionals (same decorations)
    for a, b in zip(pool['orderings']['PBE'], pool['orderings']['r2SCAN']):
        assert a['species'] == b['species']
    # HSE-compatible PBE runs
    def outcar(folder):
        txt = (folder/'OUTCAR').read_text(errors='ignore')
        e = float(re.findall(r'energy\(sigma->0\) =\s+(-?\d+\.\d+)', txt)[-1])
        conv = 'reached required accuracy' in txt
        inc = {}
        for line in (folder/'INCAR').read_text().splitlines():
            line = line.split('#')[0]
            if '=' in line:
                k, v = line.split('=', 1)
                if k.strip() in ('ENCUT', 'ISMEAR', 'SIGMA', 'ISIF', 'EDIFF', 'EDIFFG', 'PREC'):
                    inc[k.strip()] = v.strip()
        kp = [l.strip() for l in (folder/'KPOINTS').read_text().splitlines() if l.strip()]
        pot = sorted(set(re.findall(r'POTCAR:\s+(PAW_PBE \S+ \S+)', txt)))
        A, sp, fr = parse_poscar((folder/'CONTCAR').read_text())
        return {**compact(A, sp, fr), 'energy_sigma0_eV': e, 'relaxation_converged': conv, 'incar': inc,
                'kpoints': ' '.join(kp[1:]), 'potcar_labels': pot}
    for folder in sorted((SRC/'Alloy_PBE_Values').iterdir()):
        pool['hse_setting_pbe'][f'Alloy_PBE_Values/{folder.name}'] = outcar(folder)
    for m in ('GeO2', 'SnO2'):
        pool['hse_setting_pbe'][f'{m}_PBE_values'] = outcar(SRC/f'{m}_PBE_values')
    # --- assertions ---
    for fn in ('PBE', 'r2SCAN'):
        ends = {r['x_Ge']: r['energy_as_released'] for r in pool['sqs'][fn] if r['x_Ge'] in (0.0, 1.0)}
        for m, x in (('GeO2', 1.0), ('SnO2', 0.0)):
            for plane in ('001', '100', '101'):
                emin = min(r['energy_eV'] for r in pool['strain'][f'{fn}|{plane}|{m}'])
                assert abs(emin - ends[x]) < 1e-5, (fn, plane, m, emin, ends[x])
        rows = [(r['x_Ge'], r['energy_as_released']) for r in pool['sqs'][fn] if 0 < r['x_Ge'] < 1]
        x = np.array([r[0] for r in rows]); y = np.array([r[1] for r in rows])
        a, b = np.linalg.lstsq(np.c_[x*(1-x), x*(1-x)*(1-2*x)], y, rcond=None)[0]
        paper = {'PBE': (0.915, -0.117), 'r2SCAN': (1.071, -0.219)}[fn]
        assert abs(a - paper[0]) < 1e-3 and abs(b - paper[1]) < 1e-3, (fn, a, b)
    H = pool['hse_setting_pbe']
    eG = H['GeO2_PBE_values']['energy_sigma0_eV'] / 2; eS = H['SnO2_PBE_values']['energy_sigma0_eV'] / 2
    rel = {r['index']: r['energy_as_released'] for r in pool['sqs']['PBE'] if r['x_Ge'] == 0.5}
    for i in range(3):
        h = H[f'Alloy_PBE_Values/Sn_50_{i+1}']
        nG, nS = h['species'].count('Ge'), h['species'].count('Sn')
        ef = (h['energy_sigma0_eV'] - nG*eG - nS*eS) / (nG + nS)
        assert abs(ef - rel[i]/2) < 0.002, (i, ef, rel[i])
    with gzip.open(HERE/'pool.json.gz', 'wt') as f:
        json.dump(pool, f, separators=(',', ':'))
    print('pool.json.gz written', (HERE/'pool.json.gz').stat().st_size, 'bytes')


if __name__ == '__main__':
    main()
