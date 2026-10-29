"""Reference calculations for fresh-ge-soc-strained-valence (run through task_lab.run by run_reference.py).

Stage A fixes the protocol lattice parameter a0 (P = 0, spin-orbit, 80 Ry, Gamma-centred 10^3 mesh); stage B works
at strains defined relative to A0 (reference_lattice.json, written by make_task.py --stage-a).
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
G = [0.0, 0.0, 0.0]


def req(a, c, soc=True, ecut=60, k=8, shift=0, occ='smearing', dg=0.005, bk=(), labels=False):
    return {'kind': 'pw', 'a_par_A': round(a, 5), 'c_A': round(c, 5), 'spin_orbit': soc, 'ecutwfc_Ry': ecut, 'kmesh': k,
            'kshift': shift, 'occupations': occ, 'degauss_Ry': dg if occ == 'smearing' else None,
            'band_kpoints_inverse_A': [list(x) for x in bk], 'symmetry_labels': labels}


STAGE_A = {}
for a in (5.75, 5.77, 5.79):
    STAGE_A[f'A_soc_k10_ec80_{a}'] = req(a, a, ecut=80, k=10)
    STAGE_A[f'A_scalar_k10_ec80_{a}'] = req(a, a, soc=False, ecut=80, k=10)
for a in (5.75, 5.79):
    STAGE_A[f'A_soc_k8_ec60_{a}'] = req(a, a, ecut=60, k=8)
    STAGE_A[f'A_soc_k8_ec40_{a}'] = req(a, a, ecut=40, k=8)
STAGE_A['A_soc_k12_ec80_5.77'] = req(5.77, 5.77, ecut=80, k=12)

PLAN = dict(STAGE_A)
LAT = HERE/'reference_lattice.json'
if LAT.exists():
    L = json.loads(LAT.read_text())
    A0 = L['A0']
    STAGE_B = {}
    STAGE_B['B_a0_bands_k8'] = req(A0, A0, bk=[G], labels=True)
    STAGE_B['B_a0_bands_k12'] = req(A0, A0, k=12, bk=[G], labels=True)
    STAGE_B['B_a0_scalar_bands'] = req(A0, A0, soc=False, bk=[G], labels=True)
    STAGE_B['B_a0_shift1'] = req(A0, A0, shift=1)
    STAGE_B['B_a0_gamma_fixed'] = req(A0, A0, occ='fixed')
    STAGE_B['B_a0_k6_ec50'] = req(A0, A0, k=6, ecut=50, dg=0.01, bk=[G], labels=True)
    # c scans (sigma_zz = 0) at in-plane strains: spin-orbit 60 Ry Gamma-k8 with labelled Gamma levels; the scalar
    # runs give the cheap structural route and the spin-orbit-off ablation
    for e in L.get('scan_strains', []):
        for cz in L['scan_c'][str(e)]:
            STAGE_B[f'B_scan_soc_e{e}_c{cz}'] = req(A0 * (1 + e), A0 * (1 + cz), bk=[G], labels=True)
            STAGE_B[f'B_scan_scalar_e{e}_c{cz}'] = req(A0 * (1 + e), A0 * (1 + cz), soc=False, bk=[G], labels=True)
    e = -0.01
    STAGE_B['B_m1_c_a0'] = req(A0 * (1 + e), A0, bk=[G], labels=True)             # out-of-plane not relaxed
    STAGE_B['B_m1_cubic'] = req(A0 * (1 + e), A0 * (1 + e), bk=[G], labels=True)  # cubic (hydrostatic) compression
    STAGE_B['B_m1_k10_ec80'] = req(A0 * (1 + e), A0 * (1 + 0.76 * 0.01), k=10, ecut=80, bk=[G], labels=True)
    STAGE_B['B_m1_k8_ec80'] = req(A0 * (1 + e), A0 * (1 + 0.76 * 0.01), k=8, ecut=80, bk=[G], labels=True)
    PLAN.update(STAGE_B)
