"""Reference calculations for fresh-srtio3-soft-modes (requests only; run by run_reference.py through task_lab.run,
the code path the agent uses). Every request is built with the oracle's own request builders (solution/solve.py).

Stage A: needs nothing: frozen-phonon columns at a0 = A0 on 4/6/8 meshes (and 0.04 A amplitude, and at the
         experimental 3.905 A), dielectric ladder, tilt scan in the doubled cell (4x4x4) and its 6x6x6 check, cubic
         references, the strained c scan and the c = a0 (no Poisson relaxation) strained columns.
Stage B: needs stage A: equation of state (60 and 75 Ry), polar scans along the soft eigendisplacements of the
         4/6/8 force constants, the strained columns at the relaxed c (6x6x6 and 8x8x8), extra tilt angles.
A0 = 3.936 A is the oracle's rounded a0 (cubic fit of the 60 Ry, 6x6x6 EOS; make_task.py asserts it).
"""
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/fresh-srtio3-soft-modes'
spec = importlib.util.spec_from_file_location('solve', TASK/'solution/solve.py')
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)

A0 = 3.936
A_EXP = 3.905


def rows():
    path = HERE/'reference_runs.jsonl'
    if not path.exists():
        return {}
    return {r['tag']: r for r in (json.loads(l) for l in path.read_text().splitlines() if l.strip()) if 'result' in r}


def fp_tag(a, k, u, name, axis='z'):
    return f'fp_a{a}_k{k}_u{u}_{name}{axis}'


def stage_a():
    P = {}
    for k in (4, 6, 8):
        for name, r in S.fp_z_requests(A0, k).items():
            P[fp_tag(A0, k, 0.01, name)] = r
    P[fp_tag(A0, 6, 0.01, 'Sr')] = S.request('scf', 'perovskite_5', A0, A0, S.displaced([0], 2, 0.01), 6)
    for name, r in S.fp_z_requests(A0, 6, 0.04).items():
        P[fp_tag(A0, 6, 0.04, name)] = r
    for name, r in S.fp_z_requests(A_EXP, 6).items():
        P[fp_tag(A_EXP, 6, 0.01, name)] = r
    for k in (4, 6, 8):
        P[f'diel_a{A0}_k{k}'] = S.dielectric_request(A0, k)
    for th, r in S.tilt_requests(A0, 4, [0, 2, 3, 4, 5, 6]).items():
        P[f'tilt_a{A0}_k4_t{th}'] = r
    for th, r in S.tilt_requests(A0, 6, [0, 4]).items():
        P[f'tilt_a{A0}_k6_t{th}'] = r
    for k in (4, 6):
        P[f'cub_a{A0}_k{k}'] = S.cubic_request(A0, k)
    for c, r in zip(S.C_SCAN, S.c_scan_requests(A0)):
        P[f'strain_c{c}_k6'] = r
    for name, r in S.fp_x_requests(A0, A0).items():
        P[f'strainfp_c{A0}_k6_{name}x'] = r
    return P


def stage_b(R):
    P = {}
    for a, r in zip(S.EOS_A, S.eos_requests()):
        P[f'eos_a{a}_e60_k6'] = r
        P[f'eos_a{a}_e75_k6'] = dict(r, ecutwfc_Ry=75)
    for k in (4, 6, 8):
        fp = {n: R[fp_tag(A0, k, 0.01, n)]['result'] for n in S.COLS_Z}
        pattern = S.soft_pattern(S.z_block(fp))
        for kk in sorted({6, k, 8 if k == 8 else 6}):
            for t, r in zip(S.POLAR_TI_AMPS, S.polar_requests(A0, pattern, kk)):
                P[f'polar_pat{k}_k{kk}_t{t}'] = r
    P['cub_a3.936_k8'] = S.cubic_request(A0, 8)
    c_rel = round(S.c_relaxed([R[f'strain_c{c}_k6']['result'] for c in S.C_SCAN]), 3)
    for k in (6, 8):
        for name, r in S.fp_x_requests(A0, c_rel, k).items():
            P[f'strainfp_c{c_rel}_k{k}_{name}x'] = r
    for c, r in zip(S.C_SCAN, S.c_scan_requests(A0)):
        P[f'strain_c{c}_k6_e75'] = dict(r, ecutwfc_Ry=75)
    for th, r in S.tilt_requests(A0, 4, [1, 3.5, 4.5]).items():
        P[f'tilt_a{A0}_k4_t{th}'] = r
    # polar scans at the old (coarse, large) amplitudes are kept in reference_runs.jsonl as the coarse-scan variant
    for t in (0.08, 0.12):
        fp = {n: R[fp_tag(A0, 8, 0.01, n)]['result'] for n in S.COLS_Z}
        pattern = S.soft_pattern(S.z_block(fp))
        P[f'polar_pat8_k6_t{t}'] = S.request('scf', 'perovskite_5', A0, A0, (S.np.array(pattern) * t).tolist(), 6)
    return P


def plan(stage):
    if stage == 'A':
        return stage_a()
    return stage_b(rows())
