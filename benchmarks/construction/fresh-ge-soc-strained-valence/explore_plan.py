"""Exploration batches (construction record): python3 explore_plan.py BATCHNAME > batch.json"""
import json
import sys


def req(a, c, soc=True, ecut=60, k=8, shift=0, occ='smearing', dg=0.005, bk=(), labels=False):
    return {'kind': 'pw', 'a_par_A': a, 'c_A': c, 'spin_orbit': soc, 'ecutwfc_Ry': ecut, 'kmesh': k, 'kshift': shift,
            'occupations': occ, 'degauss_Ry': dg if occ == 'smearing' else None,
            'band_kpoints_inverse_A': [list(x) for x in bk], 'symmetry_labels': labels}


G = (0.0, 0.0, 0.0)
B = {}
B['e1'] = {
    'e1_k8_bands': req(5.77, 5.77, bk=[G, (0, 0, 0.02), (0, 0, 0.04)], labels=True),
    'e1_k6': req(5.77, 5.77, k=6), 'e1_k10': req(5.77, 5.77, k=10), 'e1_k12': req(5.77, 5.77, k=12),
    'e1_ec40': req(5.77, 5.77, ecut=40), 'e1_ec50': req(5.77, 5.77, ecut=50), 'e1_ec80': req(5.77, 5.77, ecut=80),
    'e1_scalar_k8': req(5.77, 5.77, soc=False, bk=[G], labels=True),
    'e1_shift1': req(5.77, 5.77, shift=1), 'e1_fixed0': req(5.77, 5.77, occ='fixed'),
    'e1_fixed1': req(5.77, 5.77, occ='fixed', shift=1), 'e1_mv02': req(5.77, 5.77, dg=0.02),
    'e1_tet_k8': req(5.712, 5.813, bk=[G], labels=True),
}
if __name__ == '__main__':
    print(json.dumps([{'tag': t, 'request': r} for t, r in B[sys.argv[1]].items()], indent=1))
