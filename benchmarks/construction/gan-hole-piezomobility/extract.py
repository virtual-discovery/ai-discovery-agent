"""Extract the released GaN strain archive into pool.json (mobilities, structures, metadata) and the staged band file.

    bash fetch.sh                      # downloads and unpacks into ~/.cache/ai-theorist/gan-hole-piezomobility/x
    python3 -B extract.py [ARCHIVE_ROOT]

For each of the 37 states (unstrained + compress/expand 1-3 % along e1..e6) it reads
  scf.in   : CELL_PARAMETERS (bohr) and ATOMIC_POSITIONS (crystal)
  epw2.in  : fermi_energy (the GW valence-band maximum used by the authors)
  gan.eig  : GW quasiparticle energies on the 10x10x10 coarse grid (bands 47-56 are staged)
  epw2.out : the zero-field SERTA table and the final table of the zero-field iterative BTE (IBTE) solve,
             with the iteration count, the last 'Max error' and the 'did not converge' flag.
The magnetic-field solve that follows in each epw2.out (Hall factor, B = 1e-10 T) is not used.
"""
import json
import re
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/gan-hole-piezomobility'
ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/gan-hole-piezomobility/x/cloud'
BANDS = (47, 56)          # 1-based band indices staged (top six valence bands, four lowest conduction bands)
NUM = r'[-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?'
ROW = re.compile(r'^\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(%s)\s+(%s)\s+(%s)\s+(%s)\s+(%s)\s*$' % ((NUM,) * 5))


def label_info(label):
    if label == 'unstrained':
        return {'component': None, 'sign': 0, 'archive_percent': 0}
    m = re.fullmatch(r'(compress|expand)_(\d)%_e(\d)', label)
    return {'component': int(m.group(3)), 'sign': -1 if m.group(1) == 'compress' else 1, 'archive_percent': int(m.group(2))}


def structure(path):
    lines = path.read_text().split('\n')
    i = next(k for k, l in enumerate(lines) if l.strip().startswith('CELL_PARAMETERS'))
    assert 'bohr' in lines[i]
    cell = [[float(x) for x in lines[i + j].split()] for j in (1, 2, 3)]
    j = next(k for k, l in enumerate(lines) if l.strip().startswith('ATOMIC_POSITIONS'))
    assert 'crystal' in lines[j]
    atoms = [[l.split()[0]] + [float(x) for x in l.split()[1:4]] for l in lines[j + 1:j + 5]]
    return cell, atoms


def fermi(path):
    return float(re.search(r'fermi_energy\s*=\s*(%s)' % NUM, path.read_text()).group(1))


def eig(path):
    lines = path.read_text().split('\n')
    nb, nk = map(int, re.search(r'bnbd=\s*40\s*\+\s*(\d+),\s*nks=\s*(\d+)', lines[0]).groups())
    vals = np.array([float(t) for t in ' '.join(lines[1:]).split()])
    arr = vals[:(43 + nb) * nk].reshape(nk, 43 + nb)
    k, e = arr[:, :3], arr[:, 3:]
    assert np.all(e[:, :40] == -100.0)
    return k, e[:, BANDS[0] - 1:BANDS[1]]


def table(lines, i):
    """Parse a 9-temperature mobility table starting at line i; return ({T: row}, next line index)."""
    out = {}
    while len(out) < 9:
        m = ROW.match(lines[i])
        if m:
            T = float(m.group(1))
            r1 = [float(m.group(k)) for k in (5, 6, 7)]
            r2 = [float(x) for x in lines[i + 1].split()[-3:]]
            r3 = [float(x) for x in lines[i + 2].split()[-3:]]
            out[T] = {'fermi_eV': float(m.group(2)), 'hole_density_cm3': float(m.group(3)), 'mu': [r1, r2, r3]}
            i += 3
        else:
            i += 1
    return out, i


def epw(path):
    lines = path.read_text().split('\n')
    s = next(k for k, l in enumerate(lines) if 'BTE in the self-energy relaxation time approximation (SERTA)' in l)
    serta, _ = table(lines, s)
    start = next(k for k, l in enumerate(lines) if 'Start solving iterative Boltzmann Transport Equation' in l)
    stop = next(k for k, l in enumerate(lines) if 'BTE in the SERTA with B-field' in l)
    assert start > s and stop > start
    iters = [k for k in range(start, stop) if re.match(r'\s+Iteration number:', lines[k])]
    blocks = [k for k in iters if any(ROW.match(lines[j]) for j in range(k, min(k + 12, stop)))]
    ibte, _ = table(lines, blocks[-1])
    errors = [float(lines[k].split()[0]) for k in range(start, stop) if lines[k].strip().endswith('Max error')]
    n_iter = int(lines[blocks[-1]].split()[-1])
    converged = not any('did not converge' in lines[k] for k in range(start, stop))
    sym = int(re.search(r'Symmetries of crystal:\s+(\d+)', '\n'.join(lines)).group(1))
    return serta, ibte, {'iterations': n_iter, 'last_max_change_cm2_per_Vs': errors[-1], 'converged_flag': converged,
                         'crystal_symmetries': sym}


def strain_from_cell(cell, cell0):
    """Rows are lattice vectors: A = A0 (I + eps) for a symmetric small strain eps."""
    F = np.linalg.solve(np.array(cell0), np.array(cell))
    eps = 0.5 * (F + F.T) - np.eye(3)
    return eps


def main():
    states = sorted(p.name for p in ROOT.iterdir() if (p/'epw2.out').exists())
    assert len(states) == 37, states
    cell0, _ = structure(ROOT/'unstrained/scf.in')
    pool, bands, kgrid = {}, {}, None
    for label in states:
        d = ROOT/label
        cell, atoms = structure(d/'scf.in')
        k, e = eig(d/'gan.eig')
        if kgrid is None:
            kgrid = k
        assert np.allclose(k, kgrid)
        serta, ibte, meta = epw(d/'epw2.out')
        eps = strain_from_cell(cell, cell0)
        voigt = [eps[0, 0], eps[1, 1], eps[2, 2], 2 * eps[1, 2], 2 * eps[0, 2], 2 * eps[0, 1]]
        info = label_info(label)
        pool[label] = {**info, 'cell_bohr': cell, 'atoms_crystal': atoms, 'vbm_eV': fermi(d/'epw2.in'),
                       'strain_voigt_from_cell': [round(v, 6) for v in voigt],
                       'serta': {str(int(T)): v for T, v in serta.items()},
                       'ibte': {str(int(T)): v for T, v in ibte.items()}, 'ibte_solver': meta}
        g = int(np.where(np.all(np.abs(k) < 1e-9, axis=1))[0][0])
        assert abs(e[g, :6].max() - pool[label]['vbm_eV']) < 2e-6, (label, e[g, :6].max(), pool[label]['vbm_eV'])
        bands[label] = np.round(e, 6).tolist()
    (HERE/'pool.json').write_text(json.dumps({'kgrid_crystal': kgrid.round(6).tolist(), 'band_indices': list(BANDS),
                                              'states': pool}, indent=1) + '\n')
    out = TASK/'environment/lab/bands.json'
    out.write_text(json.dumps({'band_indices': list(BANDS), 'energies_eV': bands}, separators=(',', ':')) + '\n')
    print(f'wrote {HERE/"pool.json"} and {out} ({out.stat().st_size / 1e6:.1f} MB)')


if __name__ == '__main__':
    main()
