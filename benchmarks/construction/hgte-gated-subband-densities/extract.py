"""Extract the archive members used by hgte-gated-subband-densities into one compact pool.

    bash run_in_image.sh extract.py            (archive cache mounted at /data, output /scratch/pool.npz)

Transport: the raw lock-in records (VG, B, VI, Vxx, Vxy) of the 26, 70 and 110(=107) nm Hall bars,
restricted to B <= 2.0 T, as float64 grids [gate, field, (VI, Vxx, Vxy)].
Calculations: every kdotpy run of FigureS3/calculations (main self-consistent runs and the two
hole-pocket wedge runs per density): the special energies of the .kpout file, the potential profile,
and the band-resolved integrated DOS (bands -3..6 and the total) on the archived 0.1 meV grid.
The 45 nm sample is not used (CompMat-Bench overlap).
"""
import glob
import json
import re
from pathlib import Path

import numpy as np

CALC = Path('/data/calc/data_repository/FigureS3/calculations')
MEAS = Path('/data/meas/data_repository/measurements')
OUT = Path('/scratch/pool.npz')
BANDS = [-3, -2, -1, 1, 2, 3, 4, 5, 6]
BMAX = 2.0
FOLDERS = {'26': 'HgTe26nm_FFT_04', '70': 'HgTe70nm_FFT_10', '107': 'HgTe107nm_FFT_03'}
MEASFILE = {'26': '26nm.csv', '70': '70nm.csv', '107': '110nm.csv'}


def kpout(path):
    out = {}
    for line in open(path):
        m = re.search(r'\[(\w+)\]:\s*(-?[\d.]+)', line)
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


def idos(path):
    with open(path) as fh:
        head = fh.readline().strip().split(',')
        units = fh.readline()
        assert units.startswith('meV')
        a = np.loadtxt(fh, delimiter=',')
    cols = ['IDOS(%d)' % b for b in BANDS] + ['IDOS(total)']
    data = np.zeros((len(a), len(cols)))
    for j, c in enumerate(cols):
        if c in head:
            data[:, j] = a[:, head.index(c)]
    # the per-band columns must add up to the archived total (all bands are listed in the file)
    allb = [i for i, c in enumerate(head) if c.startswith('IDOS(') and c != 'IDOS(total)']
    assert np.allclose(a[:, allb].sum(1), a[:, head.index('IDOS(total)')], atol=1e-6)
    return a[:, 0], data


def main():
    arrays, meta = {}, {'bands': BANDS, 'calcs': {}, 'transport': {}}
    for s, fname in MEASFILE.items():
        a = np.loadtxt(MEAS/fname, delimiter=',', skiprows=1)
        vg = np.unique(a[:, 0])
        b_all = np.unique(np.round(a[:, 1], 4))
        bsel = b_all[b_all <= BMAX + 1e-9]
        grid = np.full((len(vg), len(bsel), 3), np.nan)
        iv = np.searchsorted(vg, a[:, 0])
        ib = np.searchsorted(bsel, np.round(a[:, 1], 4))
        keep = np.round(a[:, 1], 4) <= BMAX + 1e-9
        grid[iv[keep], ib[keep]] = a[keep][:, 2:5]
        assert not np.isnan(grid).any()
        arrays[f'tr_{s}_vg'], arrays[f'tr_{s}_B'], arrays[f'tr_{s}_data'] = vg, bsel, grid
        meta['transport'][s] = {'file': fname, 'n_gate': len(vg), 'n_field': len(bsel)}
    for s, folder in FOLDERS.items():
        for kind, pattern in [('main', f'{folder}/EP*/n_*'), ('w105', f'{folder}_phi/EP*/varphi_105/n_*'),
                              ('w165', f'{folder}_phi/EP*/varphi_165/n_*')]:
            for d in sorted(glob.glob(str(CALC/pattern))):
                n = float(d.rsplit('n_', 1)[1])
                key = f'{s}_{kind}_{n:+.4f}'
                k = kpout(glob.glob(d + '/*.kpout')[0])
                E, data = idos(glob.glob(d + '/dos-integrated-byband*.csv')[0])
                pot = np.loadtxt(glob.glob(d + '/potential*.csv')[0], delimiter=',', skiprows=2)
                arrays[f'E_{key}'], arrays[f'I_{key}'], arrays[f'P_{key}'] = E, data, pot
                meta['calcs'][key] = {'sample': s, 'kind': kind, 'n_nm2': n, 'kpout': k,
                                      'kpout_text': open(glob.glob(d + '/*.kpout')[0]).read(),
                                      'folder': str(Path(d).relative_to(CALC))}
    np.savez_compressed(OUT, **arrays)
    (OUT.parent/'pool_meta.json').write_text(json.dumps(meta, indent=1))
    print(len(meta['calcs']), 'calculations;', {s: v['n_gate'] for s, v in meta['transport'].items()}, 'gate values')


if __name__ == '__main__':
    main()
