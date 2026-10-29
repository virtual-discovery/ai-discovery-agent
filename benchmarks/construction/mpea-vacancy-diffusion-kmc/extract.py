"""Map every archived hop onto the ideal 4x4x4 bcc supercell and write pool.json.

    ./run_in_image.sh python3 -B extract.py    (reads all_data.json from the download cache, mounted at /data)

Each archive entry holds pymatgen Structures (127 atoms + one vacancy) and DFT total energies (eV) of the two end
states and the transition state of one nearest-neighbour vacancy hop. Atoms are assigned to the nearest ideal site of
the 128-site conventional 4x4x4 bcc cell (site index ((i*4+j)*4+k)*2+b at (i+b/2, j+b/2, k+b/2) in units of a);
the assignment is asserted to be one-to-one with displacements below 0.25 a, the two end states must differ by one
vacancy-atom exchange between nearest-neighbour sites, and the transition-state atom closest to the hop midpoint
must be the hopping atom.
"""
import json
import os
from pathlib import Path
import numpy as np

SRC = Path(os.environ.get('MPEA_ALL_DATA', '/data/all_data.json' if Path('/data/all_data.json').exists()
                          else Path.home()/'.cache/ai-theorist/mpea-vacancy-diffusion-kmc/all_data.json'))
IDEAL = np.array([(i + b * .5, j + b * .5, k + b * .5) for i in range(4) for j in range(4) for k in range(4)
                  for b in (0, 1)]) / 4.0


def assign(structure):
    f = np.array([s['abc'] for s in structure['sites']])
    el = [s['species'][0]['element'] for s in structure['sites']]
    diff = f[:, None, :] - IDEAL[None]
    diff -= np.round(diff)
    dist = np.linalg.norm(diff, axis=2) * 4
    idx = dist.argmin(1)
    assert len(set(idx)) == len(idx) == 127 and dist.min(1).max() < 0.25
    occ = ['Va'] * 128
    for a, i in enumerate(idx):
        occ[i] = el[a]
    return occ, f, el


def main():
    data = json.loads(SRC.read_text())
    pool = []
    for n, h in enumerate(data):
        o1, _, _ = assign(h['end_state_1']['structure'])
        o2, _, _ = assign(h['end_state_2']['structure'])
        v1, v2 = o1.index('Va'), o2.index('Va')
        changed = [i for i in range(128) if o1[i] != o2[i]]
        assert sorted(changed) == sorted([v1, v2]) and o2[v1] == o1[v2]
        dv = IDEAL[v2] - IDEAL[v1]
        dv -= np.round(dv)
        assert abs(np.linalg.norm(dv) * 4 - np.sqrt(3) / 2) < 1e-6
        ts = h['transition_state']['structure']
        f = np.array([s['abc'] for s in ts['sites']])
        el = [s['species'][0]['element'] for s in ts['sites']]
        dd = f - (IDEAL[v1] + dv / 2)
        dd -= np.round(dd)
        assert el[int(np.linalg.norm(dd, axis=1).argmin())] == o1[v2]
        pool.append({'archive_index': n, 'occupation_initial': o1, 'vacancy_site_initial': v1, 'vacancy_site_final': v2,
                     'hopping_species': o1[v2],
                     'energy_eV': {'end_state_1': h['end_state_1']['energy'], 'end_state_2': h['end_state_2']['energy'],
                                   'transition_state': h['transition_state']['energy']},
                     'archived_cell_edge_A': {'end_states': h['end_state_1']['structure']['lattice']['a'],
                                              'transition_state': ts['lattice']['a']}})
    Path(__file__).with_name('pool.json').write_text(json.dumps(pool, separators=(',', ':')))
    print(len(pool), 'hops written to pool.json')


if __name__ == '__main__':
    main()
