"""Reference calculations for fresh-mos2-strain-gap-crossover (requests only; run by run_reference.py through
task_lab.run, the code path the agent uses).

Stage A: scalar-relativistic relaxations at several lattice parameters and cutoffs (zero-stress lattice parameter).
Stage B: needs reference_lattice.json (make_task.py --stage-a): the oracle's spin-orbit ladder, a dense
         relaxed-ion ladder, the uniaxial-stress relaxations and clamped-ion cells.
Stage C: needs stage B results: spin-orbit runs at the relaxed geometries of the dense ladder, clamped-ion
         spin-orbit runs, and the spin-orbit runs of the uniaxial states (oracle, reference and wrong paths).
Every request is an ordinary lab request built with the oracle's own request builder (solution/solve.py).
"""
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/fresh-mos2-strain-gap-crossover'
spec = importlib.util.spec_from_file_location('solve', TASK/'solution/solve.py')
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
req, biaxial, uniaxial = S.request, S.biaxial, S.uniaxial

DENSE = (-0.03, -0.02, -0.015, -0.01, -0.005, 0.0, 0.005, 0.01, 0.015, 0.02, 0.03)
SOC_DENSE = (-0.02, -0.015, -0.01, -0.005, 0.0, 0.005, 0.01, 0.015)
CLAMPED = (-0.02, -0.01, 0.01, 0.02)
UNI_TRUTH = (0.004, 0.008, 0.012)
UNI_CLAMPED = (0.008, 0.012)
# dense sampling for the reference uniaxial spin-orbit run: K-valley grid and both families of Gamma-K lines
Q_FINE = (0.42, 0.46, 0.48, 0.50, 0.52, 0.54, 0.58)
UNI_FINE_POINTS = ([[0.0, 0.0]] + S.near_k(0.01, 2) + S.line(S.K_DIRS['K60'], Q_FINE) + S.line(S.K_DIRS['K180'], Q_FINE))

STAGE_A = {}
for a in (3.16, 3.17, 3.18, 3.19, 3.20):
    STAGE_A[f'relax_{a}_e45'] = req('relax', 'scalar', a, biaxial(0), S.H_START, bands=S.BIAXIAL_POINTS)
for a in (3.17, 3.18, 3.19):
    STAGE_A[f'relax_{a}_e55'] = req('relax', 'scalar', a, biaxial(0), S.H_START, ecut=55)
for a in (3.17, 3.18, 3.19):
    STAGE_A[f'relax_{a}_e70'] = req('relax', 'scalar', a, biaxial(0), S.H_START, ecut=70)
for a in (3.13, 3.15, 3.17):
    STAGE_A[f'relax_{a}_e35'] = req('relax', 'scalar', a, biaxial(0), S.H_START, ecut=35)
STAGE_A['vacuum_c15'] = req('scf', 'full', 3.18, biaxial(0), 1.5659, bands=S.BIAXIAL_POINTS, c=15)
STAGE_A['soc_k9_3.18'] = req('scf', 'full', 3.18, biaxial(0), 1.5659, bands=S.BIAXIAL_POINTS, kmesh=9)


def lattice():
    path = HERE/'reference_lattice.json'
    return json.loads(path.read_text()) if path.exists() else None


def rows():
    path = HERE/'reference_runs.jsonl'
    if not path.exists():
        return {}
    return {r['tag']: r for r in (json.loads(l) for l in path.read_text().splitlines() if l.strip()) if 'result' in r}


def stage_b(L):
    a0r, H = L['A0r'], L['hline']
    hline = S.linear_fit([tuple(p) for p in H])
    h0 = round(hline(a0r), 4)
    B = {}
    for tag, r in zip(('oracle_soc_-0.02', 'oracle_soc_-0.01', 'oracle_soc_+0.00', 'oracle_soc_+0.01'),
                      S.ladder_requests(a0r, hline)):
        B[tag] = r
    for e in DENSE:
        B[f'relax_bi_{e:+.3f}'] = req('relax', 'scalar', a0r, biaxial(e), S.H_START, bands=S.BIAXIAL_POINTS)
    for y, r in zip(S.UNI_EYY, S.uniaxial_relax_requests(a0r)):
        B[f'oracle_uni_relax_{y:.3f}'] = r
    for y in UNI_TRUTH:
        B[f'uni_relax_{y:.3f}'] = req('relax', 'scalar', a0r, uniaxial(y), S.H_START, bands=S.UNIAXIAL_POINTS)
    B['uni_relax_transverse_clamped'] = req('relax', 'scalar', a0r, uniaxial(0.0), S.H_START, bands=S.UNIAXIAL_POINTS)
    for y in UNI_CLAMPED:
        B[f'uni_clamped_ion_{y:.3f}'] = req('scf', 'scalar', a0r, uniaxial(y), h0)
    for e in CLAMPED:
        B[f'clamped_bi_{e:+.3f}'] = req('scf', 'scalar', a0r, biaxial(e), h0, bands=S.BIAXIAL_POINTS)
    return B


def stage_c(L, R):
    a0r, H = L['A0r'], L['hline']
    h0 = round(S.linear_fit([tuple(p) for p in H])(a0r), 4)
    C = {}
    for e in SOC_DENSE:
        rel = R[f'relax_bi_{e:+.3f}']['result']
        C[f'soc_bi_{e:+.3f}'] = req('scf', 'full', a0r, biaxial(e), rel['final_sulfur_height_A'], bands=S.BIAXIAL_POINTS)
    for e in CLAMPED:
        C[f'soc_clamped_bi_{e:+.3f}'] = req('scf', 'full', a0r, biaxial(e), h0, bands=S.BIAXIAL_POINTS)
    oracle_state = S.uniaxial_state([R[f'oracle_uni_relax_{y:.3f}']['result'] for y in S.UNI_EYY], L['stress_offset'])
    C['oracle_uni_soc'] = S.uniaxial_soc_request(a0r, oracle_state)
    C['uni_soc_fine'] = S.uniaxial_soc_request(a0r, truth_uniaxial_state(L, R), bands=UNI_FINE_POINTS)
    tc = R['uni_relax_transverse_clamped']['result']
    C['uni_soc_transverse_clamped'] = S.uniaxial_soc_request(
        a0r, {'eyy': 0.0, 'h': tc['final_sulfur_height_A'], 'shift': tc['final_sulfur_shift_A']}, bands=UNI_FINE_POINTS)
    clamped = [R[f'uni_clamped_ion_{y:.3f}']['result'] for y in UNI_CLAMPED]
    eyy_c = S.linear_root([(r['settings']['strain'][1][1], S.stress_2d(r)[1][1] - L['stress_offset']) for r in clamped])
    C['uni_soc_clamped_ion'] = S.uniaxial_soc_request(a0r, {'eyy': round(eyy_c, 6), 'h': h0, 'shift': [0.0, 0.0]},
                                                      bands=UNI_FINE_POINTS)
    raw = S.uniaxial_state([R[f'oracle_uni_relax_{y:.3f}']['result'] for y in S.UNI_EYY], 0.0)
    C['uni_soc_uncorrected_stress'] = S.uniaxial_soc_request(a0r, raw, bands=UNI_FINE_POINTS)
    return C


def uniaxial_relaxes(R):
    tags = [f'oracle_uni_relax_{y:.3f}' for y in S.UNI_EYY] + [f'uni_relax_{y:.3f}' for y in UNI_TRUTH]
    return sorted((R[t]['result'] for t in tags), key=lambda r: r['settings']['strain'][1][1])


def truth_uniaxial_state(L, R):
    """Reference uniaxial state: transverse strain at the minimum of the relaxed-ion energy E(eyy) (zero true
    transverse stress), S height and shift interpolated there."""
    rel = uniaxial_relaxes(R)
    eyy, _ = S.parabola_vertex([(r['settings']['strain'][1][1], r['total_energy_Ry']) for r in rel])
    near = sorted(rel, key=lambda r: abs(r['settings']['strain'][1][1] - eyy))[:2]
    h = S.linear_in([(r['settings']['strain'][1][1], r['final_sulfur_height_A']) for r in near], eyy)
    shift = [S.linear_in([(r['settings']['strain'][1][1], r['final_sulfur_shift_A'][i]) for r in near], eyy) for i in range(2)]
    return {'eyy': round(eyy, 6), 'h': h, 'shift': shift}


def plan(stage=None):
    L = lattice()
    out = dict(STAGE_A)
    if L and stage in (None, 'B', 'C'):
        out.update(stage_b(L))
    if L and stage in (None, 'C'):
        R = rows()
        try:
            out.update(stage_c(L, R))
        except KeyError:
            if stage == 'C':
                raise
    return out
