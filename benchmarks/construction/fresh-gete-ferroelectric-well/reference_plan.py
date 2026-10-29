"""Reference requests for fresh-gete-ferroelectric-well (all through environment/lab/task_lab.run).

Stage A: rocksalt equation of state, k-point and cutoff ladders, and the energy/gradient grid around the R3m
minimum. Stage B (needs reference_stage_a.json written by make_task.py --stage-a): the R3m minimum itself on dense
meshes, the clamped-cubic polar scan at the rocksalt a0, the soft-mode series under compression, and the
spin-orbit comparison. Tags are unique; run_reference.py skips tags already present in reference_runs.jsonl.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
S2 = math.sqrt(2)


def req(a_rh, alpha=60.0, x=0.5, kmesh=16, kshift=0, ecut=50, occ='fixed', degauss=None, pp='sssp'):
    return {'kind': 'scf', 'pseudopotential_set': pp, 'a_rh_A': round(a_rh, 6), 'alpha_deg': round(alpha, 5),
            'x_Te': round(x, 6), 'ecutwfc_Ry': ecut, 'kmesh': kmesh, 'kshift': kshift, 'occupations': occ,
            'degauss_Ry': degauss}


def cub(a_cub, **kw):
    return req(a_cub / S2, **kw)


STAGE_A = {}
for a in (5.96, 5.98, 6.00, 6.02, 6.04, 6.06):
    for k in (12, 16, 20):
        STAGE_A[f'A_cub{a:.2f}_k{k}'] = cub(a, kmesh=k)
for a in (6.00, 6.02):
    STAGE_A[f'A_cub{a:.2f}_k24'] = cub(a, kmesh=24)
    STAGE_A[f'A_cub{a:.2f}_k12s1'] = cub(a, kmesh=12, kshift=1)
    STAGE_A[f'A_cub{a:.2f}_k8'] = cub(a, kmesh=8)
    for e in (40, 60, 75):
        STAGE_A[f'A_cub{a:.2f}_k12_e{e}'] = cub(a, kmesh=12, ecut=e)
    STAGE_A[f'A_cub{a:.2f}_k16_mv01'] = cub(a, kmesh=16, occ='smearing', degauss=0.01)
# R3m grid: centre and +- steps in a, alpha, x plus mixed points (energies + full gradients)
C = (4.36, 58.3, 0.470)
D = (0.03, 0.4, 0.005)
GRID = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1), (1, 1, 0), (1, 0, 1), (0, 1, 1)]
for i, (u, v, w) in enumerate(GRID):
    geo = dict(a_rh=C[0] + u * D[0], alpha=C[1] + v * D[1], x=C[2] + w * D[2])
    for k in (8, 12, 16):
        STAGE_A[f'A_r3m_g{i}_k{k}'] = req(**geo, kmesh=k)
STAGE_A['A_r3m_g0_k20'] = req(*C, kmesh=20)
STAGE_A['A_r3m_g0_k12_e60'] = req(*C, kmesh=12, ecut=60)
# a frequently used starting guess: relaxing only x in the experimental-like cell
for x in (0.470, 0.475, 0.480):
    STAGE_A[f'A_exp_cell_x{x:.3f}_k16'] = req(4.28, 58.36, x, kmesh=16)


def stage_b():
    ref = json.loads((HERE/'reference_stage_a.json').read_text())
    a0, (ra, ral, rx) = ref['cubic_a0_A'], ref['r3m_min']
    plan = {}
    for k in (16, 20, 24):
        plan[f'B_cub_a0_k{k}'] = cub(a0, kmesh=k)
        plan[f'B_r3m_min_k{k}'] = req(ra, ral, rx, kmesh=k)
    plan['B_r3m_min_k12'] = req(ra, ral, rx, kmesh=12)
    plan['B_cub_a0_k12'] = cub(a0, kmesh=12)
    # clamped cubic cell at a0: polar scan (gradient-bearing)
    for x in (0.475, 0.480, 0.485, 0.490):
        for k in (12, 16, 20):
            plan[f'B_clamp_x{x:.3f}_k{k}'] = cub(a0, x=x, kmesh=k)
    # same-volume comparison (cubic at the R3m volume) and R3m strain relaxed at the cubic volume
    vr = ref['r3m_min_volume_A3']
    a_same = (4 * vr) ** (1 / 3)
    plan['B_cub_at_r3m_volume_k20'] = cub(a_same, kmesh=20)
    # soft mode under compression: displacement 0.005 (and 0.0025, 0.02 checks) at x = 0.5 - dx
    for a in (5.80, 5.75, 5.70, 5.65):
        for dx in (0.005,):
            for k, occ, dg in ((20, 'smearing', 0.005), (24, 'smearing', 0.002), (20, 'smearing', 0.002),
                               (16, 'fixed', None), (12, 'smearing', 0.01)):
                plan[f'B_soft{a:.2f}_dx{dx}_k{k}_{occ}{dg or ""}'] = cub(a, x=0.5 - dx, kmesh=k, occ=occ, degauss=dg)
    for a in (5.80, 5.75, 5.70):
        plan[f'B_soft{a:.2f}_dx0.02_k20_smearing0.005'] = cub(a, x=0.48, kmesh=20, occ='smearing', degauss=0.005)
    # spin-orbit comparison at the protocol structures
    for k in (12, 16):
        for pp in ('pslibrary', 'pslibrary_fullrel'):
            plan[f'B_{pp}_cub_a0_k{k}'] = cub(a0, kmesh=k, pp=pp)
            plan[f'B_{pp}_r3m_min_k{k}'] = req(ra, ral, rx, kmesh=k, pp=pp)
    return plan
