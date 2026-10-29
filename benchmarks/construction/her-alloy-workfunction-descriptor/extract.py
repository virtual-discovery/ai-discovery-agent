"""Extract the per-spot pool from the unpacked Zenodo 22031607 tarballs.

    python3 -B extract.py [RAW_DIR]     (default ~/.cache/ai-theorist/her-alloy-workfunction-descriptor/x)

Writes pool.json.gz next to this script: for each library and measurement area (MDI) the
EDX composition and stage position, the SECCM current densities at -100..-400 mV, the
DimeNet++ dilute (1/9 ML) adsorption-energy histogram (9,000 sites, 0.005 eV bins) and the
mean/std of the 1,000 DimeNet++ work functions; the same for the pure-metal linear-combination
surrogate. Also the DFT composition records, coverage paths and repulsion summary.
"""
import csv
import gzip
import json
import os
import sys
from pathlib import Path
import numpy as np

RAW = Path(sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/.cache/ai-theorist/her-alloy-workfunction-descriptor/x'))
HERE = Path(__file__).resolve().parent
LIBS = {'Lib-1': '0010402-au-ag-pt-pd', 'Lib-2': '0010728-au-ag-pd-pt', 'Lib-3': '0011298-pt-au-ag-pd'}
EDGES = np.round(np.arange(-1.2, 1.0 + 1e-9, 0.005), 4)     # bin edges (eV)
HADS = {'dimenet': 'hads-dimenet-r8-c128l3h1-agaupdpt-v1p1-14k-variant', 'lincomb': 'hads-lincomb-agaupdpt-v1p1-variant'}
WF = {'dimenet': 'wf-dimenet-r8-c128l3h1-agaupdpt-v1p1-1k-variant', 'lincomb': 'wf-lincomb-agaupdpt-v1p1-variant'}


def hist(e):
    e = e[np.isfinite(e)]
    assert e.min() > EDGES[0] and e.max() < EDGES[-1], (e.min(), e.max())
    c, _ = np.histogram(e, EDGES)
    nz = np.nonzero(c)[0]
    return {'first_bin': int(nz[0]), 'counts': c[nz[0]:nz[-1] + 1].tolist(), 'n': int(c.sum()),
            'mean': float(e.mean()), 'std': float(e.std())}


def main():
    pool = {'bin_edges_start_eV': float(EDGES[0]), 'bin_width_eV': 0.005, 'libraries': {}}
    for lab, d in LIBS.items():
        exp = RAW/'04-activity/experimental-data'/d
        comp = {int(r['MDI']): r for r in csv.DictReader(open(exp/'composition.csv'))}
        cur = {mv: {int(r['MDI']): float(r['I(A.cm-2)']) for r in csv.DictReader(open(exp/f'activity-{mv}mv.csv'))}
               for mv in (100, 200, 300, 400)}
        spots = []
        for mdi, r in sorted(comp.items()):
            try:
                c = {el: float(r[el]) / 100 for el in ('Ag', 'Au', 'Pd', 'Pt')}
            except ValueError:
                print(lab, 'MDI', mdi, 'has no EDX composition; skipped')
                continue
            s = {'mdi': mdi, 'x_mm': float(r['x']), 'y_mm': float(r['y']), 'composition': c,
                 'j_A_cm2': {str(mv): cur[mv][mdi] for mv in cur}}
            for sur in HADS:
                h = np.load(RAW/'03-inference/output'/d/HADS[sur]/f'{mdi:05d}.npz')
                w = np.load(RAW/'03-inference/output'/d/WF[sur]/f'{mdi:05d}.npz')
                assert np.allclose(h['fractions'], w['fractions'])
                assert np.allclose(h['fractions'], [c[k] for k in ('Ag', 'Au', 'Pd', 'Pt')], atol=0.002), (mdi, h['fractions'], c)
                wf = w['energies'][np.isfinite(w['energies'])]
                s[sur] = {'hads': hist(h['energies']), 'wf_mean_eV': float(wf.mean()), 'wf_std_eV': float(wf.std()),
                          'n_wf': int(len(wf))}
            spots.append(s)
        pool['libraries'][lab] = {'archive_dir': d, 'spots': spots}
        print(lab, len(spots), 'spots')
    ds = RAW/'01-dataset/data'
    pool['dft_adsorption'] = json.load(open(ds/'adsorption-energies.json'))
    for r in pool['dft_adsorption']:
        r.pop('h_slab_traj', None)
    pool['dft_work_functions'] = json.load(open(ds/'work-functions.json'))
    pool['dft_compositions'] = list(csv.DictReader(open(ds/'compositions.csv')))
    pool['coverage_isotherms'] = json.load(open(ds/'coverage-isotherms.json'))
    pool['repulsion'] = json.load(open(ds/'repulsion.json'))
    with gzip.open(HERE/'pool.json.gz', 'wt') as f:
        json.dump(pool, f)


if __name__ == '__main__':
    main()
