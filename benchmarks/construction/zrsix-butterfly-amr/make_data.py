"""Assemble the staged archive files and the lab data of zrsix-butterfly-amr from /data (fetch.py, make_model.py,
make_counts.py, make_cache.py outputs).

    bash run_in_image.sh make_data.py
"""
import json, shutil, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, '.')
import tb
from common import COMPOUNDS, GROUP, WINDOW, XGRID, SUBMESH, COUNT_SERVE as COUNT_RANGE, COUNT_MESHES, directions, bvec, dkey

SRC = Path('/data/ZrSiX data')
TASK = Path('/repo/benchmarks/tasks/zrsix-butterfly-amr')
LAB, AGENT = TASK/'environment/lab', TASK/'environment/agent'
HRFILE = {'ZrSiS': 'wannier90_hr.dat_nsymm16', 'ZrSiSe': 'wannier90_hr.dat_nsymm16', 'ZrSiTe': 'wannier90_hr.dat'}
TRANSPORT_COST = {40: 1, 80: 1, 160: 3}
MU_MARGIN = 0.020


def stage():
    arch = AGENT/'archive'
    if arch.exists():
        shutil.rmtree(arch)
    for p in sorted(SRC.rglob('*')):
        if p.is_dir() or p.name == 'POTCAR' or p.name.startswith('wannier90_hr'):
            continue
        dest = arch/p.relative_to(SRC)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, dest)
    for c in COMPOUNDS:
        f = SRC/'input_MR'/c/HRFILE[c]
        with open(f) as fh:
            head = [fh.readline() for _ in range(3)]
            nr, n = int(head[2]), 0
            while n < nr:
                line = fh.readline(); head.append(line); n += len(line.split())
            head.append(fh.readline())       # first matrix-element line shows the column layout
        (arch/'input_MR'/c/(HRFILE[c] + '.header')).write_text(''.join(head))
    pass   # the zip already contains the record README (README.txt)


def lab_data():
    d = LAB/'data'; d.mkdir(parents=True, exist_ok=True)
    spec = {'compounds': COMPOUNDS, 'files': {c: f'input_MR/{c}/{HRFILE[c]}' for c in COMPOUNDS},
            'transport_cost': {str(k): v for k, v in TRANSPORT_COST.items()}, 'count_meshes': COUNT_MESHES,
            'count_range': {c: list(COUNT_RANGE[c]) for c in COMPOUNDS},
            'transport_mu_range': {c: [round(WINDOW[c][0] + MU_MARGIN, 3), round(WINDOW[c][1] - MU_MARGIN, 3)] for c in COMPOUNDS},
            'transport_bands': {c: [g + 1 for g in GROUP[c]] for c in COMPOUNDS}, 'slice_bands': [10, 15],
            'field_menu': {}, 'dir_index': {}}
    for c in COMPOUNDS:
        m = tb.Model(f'/data/model/{c}.npz', c)
        np.savez_compressed(d/f'model_{c}.npz', R=m.R.astype(np.int16), H_re=m.H.real.astype(np.float32),
                            H_im=m.H.imag.astype(np.float32), nw=m.nw, A=m.A)
        shutil.copy2(f'/data/counts/{c}.npz', d/f'counts_{c}.npz')
        dirs = directions(c)
        Ks, Os, B = [], [], []
        for p, a in dirs:
            z = np.load(f'/data/cache/{c}/{dkey(p, a)}.npz')
            Ks.append(z['K']); Os.append(z['O']); B.append(z['bvec']); W = z['W']
            assert np.allclose(z['bvec'], bvec(p, a)) and np.allclose(z['X'], XGRID)
        np.savez_compressed(d/f'mr_{c}.npz', K=np.stack(Ks), O=np.stack(Os), W=W, X=XGRID, mesh=np.array(SUBMESH),
                            lo=WINDOW[c][0], ebin=float(z['ebin']), bands=np.array(GROUP[c]) + 1, bvec=np.array(B),
                            volume=abs(np.linalg.det(m.A)))
        menu = {}
        for p, a in dirs:
            menu.setdefault(p, []).append(a)
        spec['field_menu'][c] = menu
        spec['dir_index'][c] = {f'{p}:{a:.1f}': i for i, (p, a) in enumerate(dirs)}
    (d/'spec.json').write_text(json.dumps(spec, indent=1) + '\n')
    return spec


if __name__ == '__main__':
    stage()
    lab_data()
    print('staged and assembled')
