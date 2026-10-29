"""Reference classical Monte Carlo for the graded CrI3 ordering temperature (and wrong-path variants).

    python3 -B mc_reference.py VARIANT OUT.json [sizes] [temps] [sweeps]

VARIANT: 'full' (the teacher model), 'nn_only' (J2 and J3 dropped). Runs every (L, T) job as an
independent chain started from the ferromagnet along z (seed 1000 L + 10 T, so a job is
reproducible on its own), on 4 worker processes, and writes the Binder cumulant, susceptibility
and specific heat with binning errors. OUT.json is updated after every job; jobs already in it
are skipped, so runs with different sizes or temperatures accumulate in one file.
"""
import os
import json
import sys
from multiprocessing import Pool
import numpy as np
import model as M
import mc


def table_for(variant, comp):
    tab = comp.bond_table()
    if variant == 'nn_only':
        tab = [x for x in tab if abs(np.linalg.norm(comp.basis[x[1]] + x[2] * comp.cell[0] + x[3] * comp.cell[1]
                                                    - comp.basis[x[0]]) / comp.L[0] - 1) < 0.01]
    return tab


def job(args):
    variant, name, L, T, sweeps, seed = args
    comp = M.Compound(name, M.load_pool()[name])
    lat = mc.Lattice(comp, L, L, table=table_for(variant, comp))
    rng = np.random.default_rng(seed)
    _, r = mc.run(lat, T, sweeps, rng, therm=sweeps // 5)
    r.update(L=L, variant=variant)
    print(json.dumps(r), flush=True)
    return r


def main():
    variant, out = sys.argv[1], sys.argv[2]
    sizes = [int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else '12,18,24,30').split(',')]
    temps = [float(x) for x in (sys.argv[4] if len(sys.argv) > 4 else '34,35,36,37,38,39,40,41').split(',')]
    sweeps = int(sys.argv[5]) if len(sys.argv) > 5 else 30000
    data = json.load(open(out)) if os.path.exists(out) else {'variant': variant, 'results': []}
    assert data['variant'] == variant
    done = {(r['L'], r['T_K']) for r in data['results']}
    jobs = [(variant, 'CrI3', L, T, sweeps, 1000 * L + int(10 * T)) for L in sorted(sizes, reverse=True) for T in temps
            if (L, T) not in done]
    with Pool(4) as p:
        for r in p.imap_unordered(job, jobs, chunksize=1):
            r['sweeps'] = sweeps
            data['results'].append(r)
            data['results'].sort(key=lambda x: (x['L'], x['T_K']))
            with open(out + '.tmp', 'w') as f:
                json.dump(data, f, indent=1)
            os.replace(out + '.tmp', out)


if __name__ == '__main__':
    main()
