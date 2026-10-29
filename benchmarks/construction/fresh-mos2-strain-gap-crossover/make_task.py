"""Derive the hidden references, ablation values and calibration for fresh-mos2-strain-gap-crossover.

    python3 -B benchmarks/construction/fresh-mos2-strain-gap-crossover/make_task.py --stage-a   # reference_lattice.json
    python3 -B benchmarks/construction/fresh-mos2-strain-gap-crossover/make_task.py             # everything else

Reads reference_runs.jsonl (lab results of reference_plan made through task_lab.run by run_reference.py) and applies
the estimators of the task's own solution/solve.py. Writes environment/lab/private/hidden.json,
private/hidden_values.json, solution/ablation_values.json and the "calibration" block of design.json, and prints the
calibration table.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/fresh-mos2-strain-gap-crossover'
LAB = TASK/'environment/lab'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S = load('solve', TASK/'solution/solve.py')
sys.path.insert(0, str(HERE))
import reference_plan  # noqa: E402
T = load('task_lab_construction', LAB/'task_lab.py')
ROWS = {r['tag']: r for r in (json.loads(l) for l in (HERE/'reference_runs.jsonl').read_text().splitlines() if l.strip())}


def res(tag):
    row = ROWS[tag]
    if 'result' not in row:
        raise SystemExit(f'reference calculation {tag} failed: {row.get("error")}')
    return row['result']


def stage_a():
    rel = [res(f'relax_{a}_e45') for a in S.A_POINTS]
    a0, hline, offset = S.a0_and_height(rel)
    out = {'A0': a0, 'A0r': round(a0, 4), 'stress_offset': offset,
           'hline': [[r['settings']['lattice_A'], r['final_sulfur_height_A']] for r in rel]}
    (HERE/'reference_lattice.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))
    for cut in ('e45', 'e55', 'e70', 'e35'):
        pts = [(r['settings']['lattice_A'], S.stress_2d(r)[0][0], r['total_energy_Ry'], r['final_sulfur_height_A'])
               for t, r in ((t, ROWS[t].get('result')) for t in ROWS if t.startswith('relax_3.') and t.endswith(cut)) if r]
        if len(pts) < 2:
            continue
        pts.sort()
        a, s, e, h = map(np.array, zip(*pts))
        fit = np.polyfit(a, s, min(2, len(a) - 1))
        roots = [x.real for x in np.roots(fit) if abs(x.imag) < 1e-9 and 3.0 < x.real < 3.4]
        emin = S.parabola_vertex(list(zip(a, e)))[0] if len(a) >= 3 else None
        print(cut, 'stress-zero', roots, 'E-min', emin, 'points', [(x, round(y, 4)) for x, y in zip(a, s)])


TOL = {'a0_A': 0.003, 'crossover_lattice_A': 0.0035, 'k_inverse_A': 0.1,
       'deformation_potential_relative': {'valence_K': 0.08, 'conduction_K': 0.08, 'gap_K': 0.06},
       'transverse_strain_percent': 0.1, 'gap_eV': 0.015}
A_EXP = 3.16


def frame(a, strain, h, shift, c=20):
    """k-space symmetry frame of a cell: reciprocal vectors (1/A, with 2 pi) and the crystal point group on k."""
    cell, _, frac = T.geometry(a, strain, h, shift, c)
    _, ops = T.crystal_ops(cell, frac, ('Mo', 'S', 'S'))
    kops = []
    for M in ops:
        K = np.round(np.linalg.inv(M).T).astype(int)[:2, :2].tolist()
        if K not in kops:
            kops.append(K)
    b = 2 * np.pi * np.linalg.inv(cell).T
    return {'b_inverse_A': b[:2, :2].round(8).tolist(), 'k_ops': kops}


def ladder(tags, a0, **kw):
    return S.ladder_analysis([res(t) for t in tags], a0, **kw) if kw else S.ladder_analysis([res(t) for t in tags], a0)


def crossings_from_edges(rows, a0r, a0, fv, fc):
    """Range (percent, relative to a0) from per-strain edge dictionaries with functions for the two differences."""
    pts = [((1 + e) * a0r / a0 - 1, E) for e, E in rows]
    up = S.crossing([(s, fv(E)) for s, E in pts])
    lo = S.crossing([(s, fc(E)) for s, E in pts])
    return [100 * lo, 100 * up]


def quad_crossing(points):
    """Zero of the quadratic through the three points nearest the sign change."""
    pts = sorted(points)
    for i in range(len(pts) - 1):
        if pts[i][1] * pts[i + 1][1] <= 0:
            j = max(0, min(i - 1, len(pts) - 3))
            x, y = zip(*pts[j:j + 3])
            c = np.polyfit(x, y, 2)
            roots = [r.real for r in np.roots(c) if abs(r.imag) < 1e-12 and pts[i][0] - 1e-9 <= r.real <= pts[i + 1][0] + 1e-9]
            return roots[0]
    raise ValueError('no sign change')


def main():
    L = json.loads((HERE/'reference_lattice.json').read_text())
    a0r = L['A0r']
    h0 = round(S.linear_fit([tuple(p) for p in L['hline']])(a0r), 4)
    off = L['stress_offset']
    R = {t: r for t, r in ROWS.items() if 'result' in r}
    # ---------------- zero-stress lattice parameter: quadratic fit of the 2D stress of the 45 Ry relaxations
    def stress_zero(cut, lattices):
        pts = [(res(f'relax_{a}_{cut}')['settings']['lattice_A'], S.stress_2d(res(f'relax_{a}_{cut}'))[0][0]) for a in lattices]
        a, s = map(np.array, zip(*pts))
        c = np.polyfit(a, s, min(2, len(a) - 1))
        return float(min((r.real for r in np.roots(c) if abs(r.imag) < 1e-12), key=lambda r: abs(r - 3.18)))
    def energy_min(cut, lattices):
        return float(S.parabola_vertex([(res(f'relax_{a}_{cut}')['settings']['lattice_A'],
                                         res(f'relax_{a}_{cut}')['total_energy_Ry']) for a in lattices])[0])
    # ---------------- equilibrium lattice parameter: minimum of E(a) of the 45 Ry relaxations (five points); the
    # analytic stress at 45-55 Ry carries a Pulay offset (see README)
    a0 = energy_min('e45', (3.16, 3.17, 3.18, 3.19, 3.2))
    ref = {'a0_A': a0}
    # ---------------- biaxial: dense spin-orbit ladder at the relaxed heights of each strain (12x12)
    soc_tags = [f'soc_bi_{e:+.3f}' for e in reference_plan.SOC_DENSE]
    LD = ladder(soc_tags, a0)
    ref['direct_gap_strain_range_percent'] = LD['direct_gap_strain_range_percent']
    ref['edge_leaving_K_at_lower_end'] = LD['edge_leaving_K_at_lower_end']
    ref['edge_leaving_K_at_upper_end'] = LD['edge_leaving_K_at_upper_end']
    ref['deformation_potentials_eV'] = LD['deformation_potentials_eV']
    # ---------------- uniaxial zigzag compression
    rel = reference_plan.uniaxial_relaxes(R)
    y = np.array([r['settings']['strain'][1][1] for r in rel])
    syy = np.array([S.stress_2d(r)[1][1] - off for r in rel])
    cq = np.polyfit(y, syy, 2)
    eyy_quad = float(min((r.real for r in np.roots(cq) if abs(r.imag) < 1e-12), key=lambda r: abs(r - 0.008)))
    truth_state = reference_plan.truth_uniaxial_state(L, R)
    stress_state = S.uniaxial_state(sorted(rel, key=lambda r: abs(S.stress_2d(r)[1][1] - off))[:2], off)
    fine = res('uni_soc_fine')
    EU = S.band_edges(fine, lines=('K60', 'K180'))
    GU = S.gap_summary(EU)
    ref['uniaxial_zigzag_compression'] = {'transverse_strain_percent': 100 * truth_state['eyy'], 'gap_eV': GU['gap_eV'],
                                          'vbm_k_crystal': GU['vbm_k_crystal'], 'cbm_k_crystal': GU['cbm_k_crystal']}
    frames = {'biaxial': frame(a0r, S.biaxial(0.0), h0, (0.0, 0.0)),
              'uniaxial': frame(a0r, S.uniaxial(truth_state['eyy']), truth_state['h'], truth_state['shift'])}

    # ---------------- robustness: defensible alternatives that must pass
    Rb = {}
    Rb['a0_energy_minimum_55Ry'] = {'a0_A': energy_min('e55', (3.17, 3.18, 3.19))}
    Rb['a0_energy_minimum_70Ry'] = {'a0_A': energy_min('e70', (3.17, 3.18, 3.19))}
    Rb['a0_stress_zero_70Ry'] = {'a0_A': stress_zero('e70', (3.17, 3.18, 3.19))}
    Rb['a0_energy_minimum_three_points_45Ry'] = {'a0_A': S.a0_and_height([res(f'relax_{a}_e45') for a in S.A_POINTS])[0]}
    # the same band crossings reported against a different admissible a0 (strains rescaled consistently)
    for name, a0x in (('stress_70Ry_a0_with_consistent_strains', Rb['a0_stress_zero_70Ry']['a0_A']),):
        rng = [100 * ((1 + v / 100) * a0 / a0x - 1) for v in ref['direct_gap_strain_range_percent']]
        Rb[name] = {'a0_A': a0x, 'direct_gap_strain_range_percent': rng}
    Rb['quadratic_interpolation_of_crossings'] = {'direct_gap_strain_range_percent': [
        100 * quad_crossing([(s, E['CB_Q'] - E['CB_K']) for s, E in LD['rows']]),
        100 * quad_crossing([(s, E['VB_Gamma'] - E['VB_K']) for s, E in LD['rows']])]}
    coarse = [t for t in soc_tags if t.endswith(('-0.020', '-0.010', '+0.000', '+0.010'))]
    LC = ladder(coarse, a0)
    Rb['one_percent_soc_ladder'] = {'direct_gap_strain_range_percent': LC['direct_gap_strain_range_percent'],
                                    'deformation_potentials_eV': LC['deformation_potentials_eV']}
    # scalar-relativistic relaxed ladder + spin-orbit shifts of all four competing states measured at zero strain
    sr = {e: S.band_edges(res(f'relax_bi_{e:+.3f}')) for e in reference_plan.DENSE}
    so0, sr0 = S.band_edges(res('soc_bi_+0.000')), sr[0.0]
    d = {k: so0[k] - sr0[k] for k in ('VB_Gamma', 'VB_K', 'CB_K', 'CB_Q')}
    corr = [(e, {k: E[k] + d[k] for k in d}) for e, E in sorted(sr.items())]
    Rb['scalar_ladder_plus_soc_shifts_of_all_edges'] = {'direct_gap_strain_range_percent': crossings_from_edges(
        corr, a0r, a0, lambda E: E['VB_Gamma'] - E['VB_K'], lambda E: E['CB_Q'] - E['CB_K'])}
    half = [(s, E) for s, E in LD['rows'] if abs(abs(s) - 0.005) < 0.002]
    m, p = sorted(half, key=lambda x: x[0])
    Rb['deformation_potentials_half_percent'] = {'deformation_potentials_eV': {
        'valence_K': (p[1]['VB_K'] - m[1]['VB_K']) / (p[0] - m[0]), 'conduction_K': (p[1]['CB_K'] - m[1]['CB_K']) / (p[0] - m[0]),
        'gap_K': ((p[1]['CB_K'] - p[1]['VB_K']) - (m[1]['CB_K'] - m[1]['VB_K'])) / (p[0] - m[0])}}
    srm, srp = sr[-0.01], sr[0.01]
    Rb['deformation_potentials_scalar_relativistic'] = {'deformation_potentials_eV': {
        'valence_K': (srp['VB_K'] - srm['VB_K']) / 0.02, 'conduction_K': (srp['CB_K'] - srm['CB_K']) / 0.02,
        'gap_K': ((srp['CB_K'] - srp['VB_K']) - (srm['CB_K'] - srm['VB_K'])) / 0.02}}
    Rb['uniaxial_quadratic_fit_of_corrected_stress'] = {'uniaxial_zigzag_compression': {'transverse_strain_percent': 100 * eyy_quad}}
    Rb['uniaxial_offset_corrected_stress'] = {'uniaxial_zigzag_compression': {'transverse_strain_percent': 100 * stress_state['eyy']}}
    OU = S.band_edges(res('oracle_uni_soc'), lines=('K60', 'K180'))
    GO = S.gap_summary(OU)
    Rb['uniaxial_coarse_sampling'] = {'uniaxial_zigzag_compression': {'gap_eV': GO['gap_eV'], 'vbm_k_crystal': GO['vbm_k_crystal'],
                                                                       'cbm_k_crystal': GO['cbm_k_crystal']}}
    Rb['uniaxial_grid_point_edges'] = {'uniaxial_zigzag_compression': {
        'gap_eV': grid_gap(fine), 'cbm_k_crystal': [round(-1 / 3, 6), round(1 / 6, 6)]}}

    # ---------------- ablations: realistic wrong paths computed from the reference runs
    Ab = {}
    sr_rows = [(e, E) for e, E in sorted(sr.items())]
    Ab['scalar_relativistic_band_energies'] = ('dp1', {'direct_gap_strain_range_percent': crossings_from_edges(
        sr_rows, a0r, a0, lambda E: E['VB_Gamma'] - E['VB_K'], lambda E: E['CB_Q'] - E['CB_K'])})
    dK = so0['VB_K'] - sr0['VB_K']
    konly = [(e, dict(E, VB_K=E['VB_K'] + dK)) for e, E in sr_rows]
    Ab['spin_orbit_shift_of_K_valence_only'] = ('dp2', {'direct_gap_strain_range_percent': crossings_from_edges(
        konly, a0r, a0, lambda E: E['VB_Gamma'] - E['VB_K'], lambda E: E['CB_Q'] - E['CB_K'])})
    cl = [res(f'soc_clamped_bi_{e:+.3f}') for e in reference_plan.CLAMPED] + [res('soc_bi_+0.000')]
    orig = S.crossing
    S.crossing = mesh_crossing          # clamped ions: the compressive crossing lies beyond -2 %; extrapolate linearly
    try:
        LCl = S.ladder_analysis(cl, a0)
    finally:
        S.crossing = orig
    Ab['clamped_sulfur_height'] = ('dp3', {'direct_gap_strain_range_percent': LCl['direct_gap_strain_range_percent'],
                                           'deformation_potentials_eV': LCl['deformation_potentials_eV']})
    a35 = stress_zero('e35', (3.13, 3.15, 3.17))
    Ab['analytic_stress_zero_35Ry'] = ('dp4', {'a0_A': a35})
    Ab['analytic_stress_zero_45Ry'] = ('dp4', {'a0_A': stress_zero('e45', (3.16, 3.17, 3.18, 3.19, 3.2))})
    Ab['experimental_lattice_as_strain_reference'] = ('dp5', {'direct_gap_strain_range_percent': [
        100 * ((1 + v / 100) * a0 / A_EXP - 1) for v in ref['direct_gap_strain_range_percent']]})
    raw = [(s, S.band_edges(r, align=False)) for s, r in
           [((1 + r['settings']['strain'][0][0]) * a0r / a0 - 1, r) for r in (res(t) for t in soc_tags)]]
    rm = min(raw, key=lambda x: abs(x[0] + 0.01)); rp = min(raw, key=lambda x: abs(x[0] - 0.01))
    Ab['raw_eigenvalues_without_vacuum_alignment'] = ('dp6', {'deformation_potentials_eV': {
        'valence_K': (rp[1]['VB_K'] - rm[1]['VB_K']) / (rp[0] - rm[0]), 'conduction_K': (rp[1]['CB_K'] - rm[1]['CB_K']) / (rp[0] - rm[0])}})
    tc = S.gap_summary(S.band_edges(res('uni_soc_transverse_clamped'), lines=('K60', 'K180')))
    Ab['uniaxial_strain_instead_of_stress'] = ('dp7', {'uniaxial_zigzag_compression': {
        'transverse_strain_percent': 0.0, 'gap_eV': tc['gap_eV'], 'vbm_k_crystal': tc['vbm_k_crystal'], 'cbm_k_crystal': tc['cbm_k_crystal']}})
    one = S.gap_summary(S.band_edges(fine, lines=('K60',)))
    Ab['q_valleys_on_one_gamma_k_line'] = ('dp8', {'uniaxial_zigzag_compression': {
        'gap_eV': one['gap_eV'], 'vbm_k_crystal': one['vbm_k_crystal'], 'cbm_k_crystal': one['cbm_k_crystal']}})
    clamped_u = [res(f'uni_clamped_ion_{y:.3f}') for y in reference_plan.UNI_CLAMPED]
    eyy_c = S.linear_root([(r['settings']['strain'][1][1], S.stress_2d(r)[1][1] - off) for r in clamped_u])
    cu = S.gap_summary(S.band_edges(res('uni_soc_clamped_ion'), lines=('K60', 'K180')))
    raw_state = S.uniaxial_state([res(f'oracle_uni_relax_{y:.3f}') for y in S.UNI_EYY], 0.0)
    ru = S.gap_summary(S.band_edges(res('uni_soc_uncorrected_stress'), lines=('K60', 'K180')))
    Ab['uniaxial_zero_of_uncorrected_analytic_stress'] = ('dp9', {'uniaxial_zigzag_compression': {
        'transverse_strain_percent': 100 * raw_state['eyy'], 'gap_eV': ru['gap_eV'], 'vbm_k_crystal': ru['vbm_k_crystal'],
        'cbm_k_crystal': ru['cbm_k_crystal']}})
    Ab['uniaxial_clamped_ions'] = ('dp3', {'uniaxial_zigzag_compression': {
        'transverse_strain_percent': 100 * eyy_c, 'gap_eV': cu['gap_eV'], 'vbm_k_crystal': cu['vbm_k_crystal'],
        'cbm_k_crystal': cu['cbm_k_crystal']}})
    mesh_rows = []
    for s, r in [((1 + r['settings']['strain'][0][0]) * a0r / a0 - 1, r) for r in (res(t) for t in soc_tags)]:
        E = S.band_edges(r)
        nocc, vac = r['occupied_bands'], E['vacuum_eV']
        mesh_cb = min(p['energies_eV'][nocc] - vac for p in r['band_kpoints_results']
                      if any(abs(p['k_crystal'][0] - q) < 1e-5 and abs(p['k_crystal'][1] - q) < 1e-5 for q in (1 / 9, 2 / 9)))
        mesh_rows.append((s, mesh_cb - E['CB_K']))
    lo_mesh = mesh_crossing(mesh_rows)
    Ab['band_edges_from_9x9_mesh_points'] = ('dp10', {'direct_gap_strain_range_percent': [100 * lo_mesh, ref['direct_gap_strain_range_percent'][1]],
                                                      'edge_leaving_K_at_lower_end': {'band': 'conduction', 'k_crystal': [round(1 / 9, 6), round(1 / 9, 6)]}})
    return a0, a0r, ref, frames, Rb, Ab, LD, LC, LCl, sr_rows, konly, d, GU, GO, one, tc, cu, eyy_c, a35


def grid_gap(r):
    """Gap from the sampled grid points only (no parabolic refinement of the Q valleys)."""
    nocc = r['occupied_bands']
    vb = max(p['energies_eV'][nocc - 1] for p in r['band_kpoints_results'])
    cb = min(p['energies_eV'][nocc] for p in r['band_kpoints_results'])
    return cb - vb


def mesh_crossing(rows):
    """Strain where the lowest 9x9-mesh conduction state off K drops below CB(K); linear extrapolation beyond the
    sampled ladder if needed."""
    pts = sorted(rows)
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if y0 * y1 <= 0:
            return x0 - y0 * (x1 - x0) / (y1 - y0)
    (x0, y0), (x1, y1) = pts[0], pts[1]
    return x0 - y0 * (x1 - x0) / (y1 - y0)


def reference_payload(ref):
    p = {k: json.loads(json.dumps(v)) for k, v in ref.items()}
    p['evidence_ids'] = [1]
    return p


def calibration(ref, Rb, ablations):
    """Grade every robustness variation and ablation against the reference payload; tabulate errors per check."""
    import copy
    trace = [{'id': 1, 'status': 'complete'}]
    abl_mod = load('ablations_construction', TASK/'solution/ablations.py')

    def graded(changes):
        payload = reference_payload(ref)
        abl_mod.merge(payload, copy.deepcopy(changes))
        ok, rep = T.grade(payload, trace, [])
        return ok, rep
    rows = {}
    ok, base = graded({})
    assert ok, base
    for kind, items in (('robustness', Rb.items()), ('ablation', ((n, a['changes']) for n, a in ablations.items()))):
        for name, changes in items:
            ok, rep = graded(changes)
            for key, err in rep['errors'].items():
                if abs(err) > 1e-12:
                    rows.setdefault(key, {'robustness': {}, 'ablations': {}})[
                        'robustness' if kind == 'robustness' else 'ablations'][name] = round(err, 5)
            failing = [k for k, v in rep['checks'].items() if not v]
            rows.setdefault('_summary', {})[name] = {'kind': kind, 'passes': ok, 'failing_checks': failing}
    return rows


def build():
    a0, a0r, ref, frames, Rb, Ab, LD, LC, LCl, sr_rows, konly, d, GU, GO, one, tc, cu, eyy_c, a35 = main()
    hidden = {'reference': ref, 'tolerances': TOL, 'frames': frames,
              'provenance': {'runs': 'benchmarks/construction/fresh-mos2-strain-gap-crossover/reference_runs.jsonl '
                                     '(task_lab.run, QE 6.7, 2 MPI ranks)', 'A0r_of_stage_b': a0r}}
    (LAB/'private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    u = ref['uniaxial_zigzag_compression']
    values = [ref['a0_A'], *ref['direct_gap_strain_range_percent'], *ref['deformation_potentials_eV'].values(),
              u['transverse_strain_percent'], u['gap_eV'],
              a0 * (1 + ref['direct_gap_strain_range_percent'][0] / 100), a0 * (1 + ref['direct_gap_strain_range_percent'][1] / 100)]
    (LAB/'private/hidden_values.json').write_text(json.dumps(values) + '\n')
    ablations = {name: {'decision_point': dp, 'changes': ch} for name, (dp, ch) in Ab.items()}
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': ablations, 'robustness': Rb}, indent=1) + '\n')
    cal = calibration(ref, Rb, ablations)
    design = json.loads((TASK/'design.json').read_text())
    design['calibration'] = {'reference': ref, 'tolerances': TOL, 'errors_by_check': {k: v for k, v in cal.items() if k != '_summary'},
                             'variation_outcomes': cal['_summary']}
    (TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
    print(json.dumps({'reference': ref, 'calibration': cal}, indent=1, default=float))
    return ref, Rb, ablations


if __name__ == '__main__':
    if '--stage-a' in sys.argv:
        stage_a()
    else:
        build()
