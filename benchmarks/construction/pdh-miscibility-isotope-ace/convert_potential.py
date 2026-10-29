"""Convert the released ACE B-basis potential (PdH_PBE_v03a.yaml) into a JSON C-tilde basis for the
numpy evaluator (ace_numpy.py). Uses python-ace (pyace 0.4.0rc1) on the build host only:

    $CACHE/venv/bin/python -B convert_potential.py $CACHE/potential_files/ace/PdH_PBE_v03a.yaml out.json

The conversion B-basis -> C-tilde basis is the library's own exact linear transformation (generalised
Clebsch-Gordan contraction); no parameter is changed.
"""
import json
import sys
import numpy as np
from pyace import ACEBBasisSet

b = ACEBBasisSet(sys.argv[1])
c = b.to_ACECTildeBasisSet()
rf = c.radial_functions
out = {'elements': list(c.elements_name), 'source': 'PdH_PBE_v03a.yaml (Mrovec et al. 2026, Zenodo 20591788, CC BY 4.0)',
       'inner_cutoff_type': rf.inner_cutoff_type, 'nradbase': c.nradbase, 'nradmax': c.nradmax, 'lmax': c.lmax,
       'cut': rf.cut, 'dcut': rf.dcut, 'cut_in': rf.cut_in, 'dcut_in': rf.dcut_in, 'prehc': rf.prehc,
       'lambhc': rf.lambhc, 'radbase': rf.radbasenameij,
       'crad': np.array(rf.crad).tolist(),
       'embeddings': {}, 'rank1': {}, 'functions': {}}
for mu, e in c.map_embedding_specifications.items():
    out['embeddings'][str(mu)] = {'npoti': e.npoti, 'ndensity': e.ndensity, 'FS_parameters': list(e.FS_parameters),
                                  'rho_core_cutoff': e.rho_core_cutoff, 'drho_core_cutoff': e.drho_core_cutoff}
for mu in range(len(c.elements_name)):
    out['rank1'][str(mu)] = [{'mus': list(f.mus), 'ns': list(f.ns), 'ctildes': np.array(f.ctildes).tolist()}
                             for f in c.basis_rank1[mu]]
    out['functions'][str(mu)] = [{'rank': f.rank, 'mus': list(f.mus), 'ns': list(f.ns), 'ls': list(f.ls),
                                  'ms_combs': np.array(f.ms_combs).tolist(), 'ctildes': np.array(f.ctildes).tolist()}
                                 for f in c.basis[mu]]
assert all(abs(x) < 1e-14 for x in getattr(c, 'E0vals', [0, 0])) if hasattr(c, 'E0vals') else True
json.dump(out, open(sys.argv[2], 'w'))
print('functions', {k: len(v) for k, v in out['functions'].items()}, 'rank1', {k: len(v) for k, v in out['rank1'].items()})
