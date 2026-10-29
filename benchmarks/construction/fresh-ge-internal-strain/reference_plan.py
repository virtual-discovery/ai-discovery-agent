"""Reference calculations for fresh-ge-internal-strain (requests only; run by run_reference.py through task_lab.run).

Stage A (equations of state) fixes the zero-pressure lattice parameter A0 and the 6 GPa lattice parameter A6;
stage B (everything at A0 and A6) uses the values derived from stage A, rounded to 1e-4 A and written in
reference_lattice.json by make_task.py --stage-a. Every request is an ordinary lab request.
Tags: <kind/geometry>_<a>_<mesh>_<cutoff/pseudo>[_variant].
"""
import json
from pathlib import Path

GBRV, PSL = 'ge_pbe_v1.4.uspp.F.UPF', 'Ge.pbe-n-rrkjus_psl.1.0.0.UPF'
ZERO = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]


def req(a, kind='scf', strain=ZERO, disp=(0, 0, 0), ecut=45, k=12, shift=0, occ='fixed', degauss=None, thr=None, pseudo=GBRV):
    return {'kind': kind, 'pseudopotential': pseudo, 'lattice_A': a, 'strain': strain, 'displacement_A': list(disp),
            'ecutwfc_Ry': ecut, 'kmesh': k, 'kshift': shift, 'occupations': occ, 'degauss_Ry': degauss,
            'forc_conv_thr_Ry_per_bohr': thr}


def shear(e4):          # engineering shear e4 = 2 eps_yz
    return [[0, 0, 0], [0, 0, e4 / 2], [0, e4 / 2, 0]]


def normal(exx):
    return [[exx, 0, 0], [0, 0, 0], [0, 0, 0]]


STAGE_A = {}
for a in (5.73, 5.75, 5.76, 5.77, 5.79):
    STAGE_A[f'cub_{a}_k12s0_e60'] = req(a, ecut=60)
for a in (5.73, 5.76, 5.79):
    STAGE_A[f'cub_{a}_k12s0_e45'] = req(a)
    STAGE_A[f'cub_{a}_k12s0_e35'] = req(a, ecut=35)
    STAGE_A[f'cub_{a}_k12s0_psl30'] = req(a, ecut=30, pseudo=PSL)
for a in (5.75, 5.77):
    STAGE_A[f'cub_{a}_k12s1_e60'] = req(a, ecut=60, shift=1)
for a in (5.58, 5.59, 5.60, 5.61, 5.62):
    STAGE_A[f'cub_{a}_k12s0_e60'] = req(a, ecut=60)
for a in (5.58, 5.60, 5.62):
    STAGE_A[f'cub_{a}_k12s0_e45'] = req(a)
    STAGE_A[f'cub_{a}_k12s0_e35'] = req(a, ecut=35)
    STAGE_A[f'cub_{a}_k12s0_psl30'] = req(a, ecut=30, pseudo=PSL)


def stage_b(a0, a6):
    B = {}
    # zero pressure: band ordering at Gamma (Gamma-containing meshes), and what a shifted mesh shows
    B['cub_a0_k8s0_e60'] = req(a0, ecut=60, k=8)
    B['cub_a0_k12s0_e60_mv002'] = req(a0, ecut=60, occ='smearing', degauss=0.002)
    B['cub_a0_k8s1_e45'] = req(a0, k=8, shift=1)
    # 6 GPa: clamped shear (internal-strain force and stress), k ladders for both mesh types
    for k in (8, 12, 16):
        B[f'e4_a6_k{k}s0'] = req(a6, strain=shear(0.01), k=k)
    for k in (8, 12):
        B[f'e4_a6_k{k}s1'] = req(a6, strain=shear(0.01), k=k, shift=1)
    B['e4_a6_k12s0_e60'] = req(a6, strain=shear(0.01), ecut=60)
    B['e4_a6_k12s0_mv002'] = req(a6, strain=shear(0.01), occ='smearing', degauss=0.002)
    B['e4_a6_k8s1_mv02'] = req(a6, strain=shear(0.01), k=8, shift=1, occ='smearing', degauss=0.02)
    B['e4m_a6_k12s0'] = req(a6, strain=shear(-0.01))
    B['e4s_a6_k12s0'] = req(a6, strain=shear(0.002))
    B['e4_a6_k12s0_psl30'] = req(a6, strain=shear(0.01), ecut=30, pseudo=PSL)
    B['cub_a6_k12s0'] = req(a6)
    B['cub_a6_k12s0_e60'] = req(a6, ecut=60)
    B['cub_a6_k8s1'] = req(a6, k=8, shift=1)
    B['cub_a6_k12s0_psl30'] = req(a6, ecut=30, pseudo=PSL)
    # zone-centre optical mode: displacement along a cube axis (no cubic term) and along [111]
    for k in (8, 12, 16):
        B[f'ph100_a6_k{k}s0'] = req(a6, disp=(0.01, 0, 0), k=k)
    for k in (8, 12):
        B[f'ph100_a6_k{k}s1'] = req(a6, disp=(0.01, 0, 0), k=k, shift=1)
    B['ph100_a6_k12s0_e60'] = req(a6, disp=(0.01, 0, 0), ecut=60)
    B['ph100_a6_k12s0_mv002'] = req(a6, disp=(0.01, 0, 0), occ='smearing', degauss=0.002)
    B['ph100_a6_k12s0_mv01'] = req(a6, disp=(0.01, 0, 0), occ='smearing', degauss=0.01)
    B['ph100_a6_k8s1_mv02'] = req(a6, disp=(0.01, 0, 0), k=8, shift=1, occ='smearing', degauss=0.02)
    B['ph100_a6_k12s0_psl30'] = req(a6, disp=(0.01, 0, 0), ecut=30, pseudo=PSL)
    B['ph100b_a6_k12s0'] = req(a6, disp=(0.03, 0, 0))
    for d in (0.02, -0.02, 0.05, -0.05):
        B[f'ph111_{d}_a6_k12s0'] = req(a6, disp=(d, d, d))
    # relaxed shear: displacement scan and relaxations
    B['e4u_a6_k12s0'] = req(a6, strain=shear(0.01), disp=(-0.012, 0, 0))
    B['e4u_a6_k8s1'] = req(a6, strain=shear(0.01), disp=(-0.012, 0, 0), k=8, shift=1)
    B['relax_e4_a6_k12s0_1e-5'] = req(a6, kind='relax', strain=shear(0.01), thr=1e-5)
    B['relax_e4s_a6_k12s0_1e-3'] = req(a6, kind='relax', strain=shear(0.002), thr=1e-3)
    B['relax_e4_a6_k8s1_1e-4'] = req(a6, kind='relax', strain=shear(0.01), k=8, shift=1, thr=1e-4)
    # C11, C12 (normal strain; no internal relaxation by symmetry)
    for e in (0.01, -0.01, 0.005, -0.005, 0.02):
        B[f'exx_{e}_a6_k12s0'] = req(a6, strain=normal(e))
    B['exx_0.01_a6_k8s1'] = req(a6, strain=normal(0.01), k=8, shift=1)
    B['exx_-0.01_a6_k8s1'] = req(a6, strain=normal(-0.01), k=8, shift=1)
    # Grueneisen: optical mode at +-1 % in a around A6
    for tag, f in (('m', 0.99), ('p', 1.01), ('mh', 0.995), ('ph', 1.005)):
        a = round(a6 * f, 4)
        B[f'ph100_a6{tag}_k12s0'] = req(a, disp=(0.01, 0, 0))
        B[f'cub_a6{tag}_k12s0'] = req(a)
    for tag, f in (('m', 0.99), ('p', 1.01)):
        a = round(a6 * f, 4)
        B[f'ph100_a6{tag}_k8s0'] = req(a, disp=(0.01, 0, 0), k=8)
        B[f'ph100_a6{tag}_k8s1_mv02'] = req(a, disp=(0.01, 0, 0), k=8, shift=1, occ='smearing', degauss=0.02)
    return B


def plan():
    lat = Path(__file__).with_name('reference_lattice.json')
    out = dict(STAGE_A)
    if lat.exists():
        L = json.loads(lat.read_text())
        out.update(stage_b(L['A0'], L['A6']))
    return out


PLAN = plan()
