"""Reference (hidden) diffusion quantities from the barrier model trained on every clean archived hop.

    ./run_in_image.sh python3 -B reference.py      -> reference.json  (~15 min on 4 cores)
"""
import json
import numpy as np
import kmclib as K
import common as C

STATS = dict(L=10, nsteps=50000, ntraj=600, nmc=150)


def main():
    ref = C.model(~K.PATHO)
    runs = {}
    for c in 'AB':
        runs[f'{c}@2000'] = K.simulate(C.xvec(c), C.T_ETA, ref, seed=11, **STATS)
    for T in C.T_ARRHENIUS:
        runs[f'C@{int(T)}'] = K.simulate(C.xvec('C'), T, ref, seed=13, **STATS)
    runs['B@2000:thermo_ideal'] = K.simulate(C.xvec('B'), C.T_ETA, C.thermo_ideal(ref), seed=17, **STATS)
    runs['B@2000:kra_simple'] = K.simulate(C.xvec('B'), C.T_ETA, C.kra_simple(ref), seed=19, **STATS)
    out = {k: {'log10_eta': float(r['log10_eta']), 'log10_Dtracer': [None if not np.isfinite(v) else float(np.log10(v)) for v in r['Dt']],
               'correlation_factor': [None if not np.isfinite(v) else float(v) for v in r['f']], 'f_vacancy': float(r['fv']),
               'Dv': float(r['Dv']), 'mean_accepted_barrier_eV': float(r['acc']), 'mean_encountered_barrier_eV': float(r['enc'])}
           for k, r in runs.items()}
    B = runs['B@2000']['Dt']
    fin = [i for i in range(6) if np.isfinite(B[i])]
    fast, slow = max(fin, key=lambda i: B[i]), min(fin, key=lambda i: B[i])
    Ccomp = runs['C@2000']['Dt']
    fastC = max([i for i in range(6) if np.isfinite(Ccomp[i])], key=lambda i: Ccomp[i])
    targets = {'log10_eta_vacancy': {c: out[f'{c}@2000']['log10_eta'] for c in 'ABC'},
               'fastest_tracer_B': K.EL[fast], 'slowest_tracer_B_info': K.EL[slow],
               'log10_tracer_ratio_B': float(np.log10(B[fast] / B[slow])),
               'fastest_tracer_C': K.EL[fastC],
               'activation_energy_C_eV': float(C.arrhenius(C.T_ARRHENIUS, [runs[f'C@{int(T)}']['Dt'][fastC] for T in C.T_ARRHENIUS])),
               'log10_eta_thermo_ideal_B': out['B@2000:thermo_ideal']['log10_eta'],
               'log10_eta_kra_simple_B': out['B@2000:kra_simple']['log10_eta']}
    json.dump({'stats': STATS, 'runs': out, 'targets': targets,
               'model': {'K0': ref['K0'].tolist(), 'W': ref['W'].tolist(), 'Q': ref['Q'].tolist(), 'V': ref['V'].tolist()}},
              open('reference.json', 'w'), indent=1)
    print(json.dumps(targets, indent=1))


if __name__ == '__main__':
    main()
