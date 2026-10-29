"""Reference solution for fresh-mos2-strain-gap-crossover; runs as the Harbor oracle and in check_task.

It requests new pw.x calculations from the lab (11 of them). The functions at module level (request builders and
estimators) are also used by the construction scripts (benchmarks/construction/fresh-mos2-strain-gap-crossover/)
to plan the reference calculations and to apply the same estimators to them.
"""
import json
import math
import os
import subprocess
import time
from pathlib import Path

APP = Path(os.environ.get('APP_DIR', '/app'))
ECUT, KMESH, CELL_HEIGHT = 45, 12, 20
H_START = 1.565                         # starting S height (A) of the relaxations
A_POINTS = (3.17, 3.18, 3.19)           # lattice parameters around the energy minimum (A)
LADDER = (-0.02, -0.01, 0.0, 0.01)      # biaxial strains of the spin-orbit ladder (relative to a0)
UNI_EXX = -0.04                         # zigzag (x) strain of the uniaxial-stress state
UNI_EYY = (0.006, 0.010)                # trial transverse (armchair, y) strains bracketing zero transverse stress
SOC_KMESH = 9                           # spin-orbit scf mesh (agrees with 12x12 to 0.1 meV for the band edges)
K_DIRS = {'K60': (1 / 3, 1 / 3), 'K180': (-2 / 3, 1 / 3)}      # zone corners at 60 and 180 degrees from x


def r6(x):
    return round(x, 6)


def line(direction, ts):
    """Points t * K on the Gamma-K line towards the zone corner `direction` (crystal coordinates)."""
    return [[r6(t * direction[0]), r6(t * direction[1])] for t in ts]


def near_k(step, n):
    return [[r6(1 / 3 + i * step), r6(1 / 3 + j * step)] for i in range(-n, n + 1) for j in range(-n, n + 1)]


Q_TS = (0.40, 0.45, 0.50, 0.55, 0.60)
BIAXIAL_POINTS = [[0.0, 0.0], [r6(1 / 3), r6(1 / 3)], [r6(1 / 9), r6(1 / 9)], [r6(2 / 9), r6(2 / 9)]] + line(K_DIRS['K60'], Q_TS)
NEAR_K7 = [[r6(1 / 3 + d1), r6(1 / 3 + d2)] for d1, d2 in
           ((0.0, 0.0), (0.01, 0.0), (-0.01, 0.0), (0.0, 0.01), (0.0, -0.01), (0.01, -0.01), (-0.01, 0.01))]
UNIAXIAL_POINTS = [[0.0, 0.0]] + NEAR_K7 + line(K_DIRS['K60'], Q_TS) + line(K_DIRS['K180'], Q_TS)


# ---------------------------------------------------------------- requests
def request(kind, rel, a, strain, h, shift=(0.0, 0.0), bands=(), kmesh=KMESH, ecut=ECUT, c=CELL_HEIGHT):
    return {'kind': kind, 'relativistic': rel, 'lattice_A': round(a, 4),
            'strain': [[round(v, 6) for v in row] for row in strain], 'sulfur_height_A': round(h, 4),
            'sulfur_shift_A': [round(s, 4) for s in shift], 'cell_height_A': c, 'ecutwfc_Ry': ecut, 'kmesh': kmesh,
            'band_kpoints': [list(k) for k in bands]}


def biaxial(e):
    return [[e, 0.0], [0.0, e]]


def uniaxial(eyy):
    return [[UNI_EXX, 0.0], [0.0, eyy]]


def a0_requests():
    return [request('relax', 'scalar', a, biaxial(0.0), H_START, bands=BIAXIAL_POINTS) for a in A_POINTS]


def ladder_requests(a0r, hline):
    return [request('scf', 'full', a0r, biaxial(e), hline(a0r * (1 + e)), bands=BIAXIAL_POINTS, kmesh=SOC_KMESH)
            for e in LADDER]


def uniaxial_relax_requests(a0r):
    return [request('relax', 'scalar', a0r, uniaxial(y), H_START) for y in UNI_EYY]


def uniaxial_soc_request(a0r, state, bands=UNIAXIAL_POINTS, kmesh=SOC_KMESH):
    return request('scf', 'full', a0r, uniaxial(state['eyy']), state['h'], state['shift'], bands=bands, kmesh=kmesh)


# ---------------------------------------------------------------- estimators (pure functions of lab results)
def stress_2d(res):
    """In-plane Cauchy stress (N/m, tensile positive) from the lab's pw.x-convention stress (GPa) x cell height."""
    c = res['cell_vectors_A'][2][2]
    return [[-res['stress_GPa'][i][j] * c * 0.1 for j in range(2)] for i in range(2)]


def linear_root(points):
    """x where y = 0 on the straight line through the two points nearest zero (y)."""
    (x0, y0), (x1, y1) = sorted(points, key=lambda p: abs(p[1]))[:2]
    return x0 - y0 * (x1 - x0) / (y1 - y0)


def linear_in(points, x):
    (x0, y0), (x1, y1) = sorted(sorted(points, key=lambda p: abs(p[0] - x))[:2])
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def parabola_vertex(points):
    """Vertex (x, y) of the least-squares parabola through (x, y) points."""
    n = len(points)
    xm = sum(x for x, _ in points) / n
    X = [x - xm for x, _ in points]
    Y = [y for _, y in points]
    s = [sum(x ** k for x in X) for k in range(5)]
    t = [sum(y * x ** k for x, y in zip(X, Y)) for k in range(3)]
    # normal equations for y = c0 + c1 x + c2 x^2 (Cramer's rule)
    M = [[s[0], s[1], s[2]], [s[1], s[2], s[3]], [s[2], s[3], s[4]]]
    det = lambda A: (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                     + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))
    D = det(M)
    c = []
    for j in range(3):
        Mj = [row[:] for row in M]
        for i in range(3):
            Mj[i][j] = t[i]
        c.append(det(Mj) / D)
    xv = -c[1] / (2 * c[2])
    return xm + xv, c[0] + c[1] * xv + c[2] * xv ** 2


def linear_fit(points):
    n = len(points)
    xm, ym = sum(x for x, _ in points) / n, sum(y for _, y in points) / n
    slope = sum((x - xm) * (y - ym) for x, y in points) / sum((x - xm) ** 2 for x, _ in points)
    return lambda x: ym + slope * (x - xm)


def a0_and_height(relaxes):
    """Equilibrium lattice parameter (energy minimum of the relaxed-ion cells), the relaxed S height as a linear
    function of a, and the cell-averaged in-plane stress (N/m) the lab reports at a0.

    The energy route is used because the analytic stress at 45-55 Ry still carries a Pulay stress of about
    +0.25 N/m (it agrees between 45 and 55 Ry but moves at 70 Ry), while the energy minimum is the same at 45, 55
    and 70 Ry."""
    a0, _ = parabola_vertex([(r['settings']['lattice_A'], r['total_energy_Ry']) for r in relaxes])
    hline = linear_fit([(r['settings']['lattice_A'], r['final_sulfur_height_A']) for r in relaxes])
    stress_at_a0 = linear_fit([(r['settings']['lattice_A'], stress_2d(r)[0][0]) for r in relaxes])(a0)
    return a0, hline, stress_at_a0


def vacuum_level(res, width_A=1.5):
    """Plateau of the planar-averaged electrostatic potential far from the layer (the layer sits at c/2)."""
    pa = res['electrostatic_potential_planar_average']
    c = res['cell_vectors_A'][2][2]
    pts = [v for z, v in zip(pa['z_A'], pa['V_eV']) if min(z, c - z) <= width_A]
    return sum(pts) / len(pts)


def _close(k, ref, tol=1e-5):
    return abs(k[0] - ref[0]) < tol and abs(k[1] - ref[1]) < tol


def parabola_min(ts, es):
    """Minimum of the parabola through the lowest sampled point and its two neighbours."""
    i = min(range(len(es)), key=lambda j: es[j])
    if i == 0 or i == len(es) - 1:
        return es[i], ts[i]
    (x0, y0), (x1, y1), (x2, y2) = (ts[i - 1], es[i - 1]), (ts[i], es[i]), (ts[i + 1], es[i + 1])
    d = (x0 - x1) * (x0 - x2) * (x1 - x2)
    a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / d
    b = (x2 ** 2 * (y0 - y1) + x1 ** 2 * (y2 - y0) + x0 ** 2 * (y1 - y2)) / d
    c = y0 - a * x0 ** 2 - b * x0
    if a <= 0:
        return es[i], ts[i]
    t = -b / (2 * a)
    return c - b * b / (4 * a), t


def band_edges(res, lines=('K60',), align=True):
    """Band edges (eV; relative to the vacuum level unless align=False) and where they sit (crystal k).

    Gamma: the sampled Gamma point; K valley: extremum over the sampled points within 0.025 of K = (1/3, 1/3);
    Q valleys: parabolic minimum of the lowest conduction band along each sampled Gamma-K line."""
    vac = vacuum_level(res) if align else 0.0
    nocc = res['occupied_bands']
    pts = res['band_kpoints_results']
    E = {'vacuum_eV': vacuum_level(res)}
    for p in pts:
        if _close(p['k_crystal'], (0, 0)):
            E['VB_Gamma'], E['CB_Gamma'] = p['energies_eV'][nocc - 1] - vac, p['energies_eV'][nocc] - vac
    near = [p for p in pts if abs(p['k_crystal'][0] - 1 / 3) < 0.025 and abs(p['k_crystal'][1] - 1 / 3) < 0.025]
    top = max(near, key=lambda p: p['energies_eV'][nocc - 1])
    bottom = min(near, key=lambda p: p['energies_eV'][nocc])
    E['VB_K'], E['k_VB_K'] = top['energies_eV'][nocc - 1] - vac, top['k_crystal']
    E['VB_K_split'] = top['energies_eV'][nocc - 1] - top['energies_eV'][nocc - 2]
    E['CB_K'], E['k_CB_K'] = bottom['energies_eV'][nocc] - vac, bottom['k_crystal']
    for name in lines:
        d = K_DIRS[name]
        pl = sorted((p['k_crystal'][0] / d[0], p['energies_eV'][nocc] - vac) for p in pts
                    if abs(p['k_crystal'][0] * d[1] - p['k_crystal'][1] * d[0]) < 1e-6 and 0.2 < p['k_crystal'][0] / d[0] < 0.8)
        e, t = parabola_min([x for x, _ in pl], [y for _, y in pl])
        E['CB_Q_' + name], E['k_CB_Q_' + name] = e, [r6(t * d[0]), r6(t * d[1])]
    best = min(lines, key=lambda n: E['CB_Q_' + n])
    E['CB_Q'], E['k_CB_Q'] = E['CB_Q_' + best], E['k_CB_Q_' + best]
    return E


def gap_summary(E):
    vbm = max((E['VB_Gamma'], [0.0, 0.0]), (E['VB_K'], E['k_VB_K']), key=lambda x: x[0])
    cbm = min((E['CB_K'], E['k_CB_K']), (E['CB_Q'], E['k_CB_Q']), (E['CB_Gamma'], [0.0, 0.0]), key=lambda x: x[0])
    return {'gap_eV': cbm[0] - vbm[0], 'vbm_k_crystal': vbm[1], 'cbm_k_crystal': cbm[1]}


def crossing(points):
    """Strain where f changes sign, by linear interpolation between the bracketing ladder points."""
    pts = sorted(points)
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if y0 * y1 <= 0 and y0 != y1:
            return x0 - y0 * (x1 - x0) / (y1 - y0)
    raise ValueError('no sign change in the sampled range: %s' % pts)


def ladder_analysis(runs, a0):
    """Crossovers (percent, relative to a0), the edges that leave K, and deformation potentials (eV per unit strain)
    from spin-orbit runs at several biaxial strains of one reference lattice parameter."""
    rows = []
    for r in runs:
        e = r['settings']['strain'][0][0]
        strain = (1 + e) * r['settings']['lattice_A'] / a0 - 1          # strain relative to a0
        rows.append((strain, band_edges(r)))
    rows.sort(key=lambda x: x[0])
    upper = crossing([(s, E['VB_Gamma'] - E['VB_K']) for s, E in rows])
    lower = crossing([(s, E['CB_Q'] - E['CB_K']) for s, E in rows])
    q_at = min(rows, key=lambda x: abs(x[0] - lower))[1]['k_CB_Q']
    dps = {}
    for key, f in (('valence_K', lambda E: E['VB_K']), ('conduction_K', lambda E: E['CB_K']),
                   ('gap_K', lambda E: E['CB_K'] - E['VB_K'])):
        m = min(rows, key=lambda x: abs(x[0] + 0.01))
        p = min(rows, key=lambda x: abs(x[0] - 0.01))
        dps[key] = (f(p[1]) - f(m[1])) / (p[0] - m[0])
    return {'direct_gap_strain_range_percent': [100 * lower, 100 * upper],
            'edge_leaving_K_at_lower_end': {'band': 'conduction', 'k_crystal': q_at},
            'edge_leaving_K_at_upper_end': {'band': 'valence', 'k_crystal': [0.0, 0.0]},
            'deformation_potentials_eV': dps, 'rows': rows}


def uniaxial_state(relaxes, stress_offset):
    """Transverse strain at which the transverse stress equals the stress the lab reports for the unstrained cell at
    a0 (the Pulay offset of the analytic stress, isotropic in the plane), with the relaxed S height and shift there
    (linear in eyy)."""
    pts = [(r['settings']['strain'][1][1], stress_2d(r)[1][1] - stress_offset, r['final_sulfur_height_A'],
            r['final_sulfur_shift_A']) for r in relaxes]
    eyy = linear_root([(p[0], p[1]) for p in pts])
    h = linear_in([(p[0], p[2]) for p in pts], eyy)
    shift = [linear_in([(p[0], p[3][i]) for p in pts], eyy) for i in range(2)]
    return {'eyy': round(eyy, 6), 'h': h, 'shift': shift}


# ---------------------------------------------------------------- lab access
def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


EVIDENCE = []


def calc(req):
    job = lab('evaluate', req)
    while job['status'] == 'running':
        time.sleep(5)
        job = lab('job/%d' % job['id'])
    if job['status'] != 'complete':
        raise RuntimeError('calculation failed: %s' % job.get('error'))
    EVIDENCE.append(job['id'])
    return job['result']


def main():
    # Lattice parameter: zero in-plane stress with the S height relaxed (45 Ry: the stress agrees with 55 Ry, while
    # 35 Ry carries a large Pulay stress). The relaxed S height is linear in a over the strains used here.
    relax0 = [calc(r) for r in a0_requests()]
    a0, hline, offset = a0_and_height(relax0)
    a0r = round(a0, 4)
    # Spin-orbit band energies at relaxed-ion geometries (S height from the scalar-relativistic relaxations)
    ladder = [calc(r) for r in ladder_requests(a0r, hline)]
    L = ladder_analysis(ladder, a0)
    # Uniaxial stress along zigzag: transverse strain from zero transverse stress (relaxed ions), then spin-orbit
    # bands sampling both inequivalent families of Gamma-K lines (Q valleys along x and at +-60 degrees).
    urelax = [calc(r) for r in uniaxial_relax_requests(a0r)]
    state = uniaxial_state(urelax, offset)
    usoc = calc(uniaxial_soc_request(a0r, state))
    EU = band_edges(usoc, lines=('K60', 'K180'))
    G = gap_summary(EU)
    submission = {
        'a0_A': a0,
        'direct_gap_strain_range_percent': L['direct_gap_strain_range_percent'],
        'edge_leaving_K_at_lower_end': L['edge_leaving_K_at_lower_end'],
        'edge_leaving_K_at_upper_end': L['edge_leaving_K_at_upper_end'],
        'deformation_potentials_eV': L['deformation_potentials_eV'],
        'uniaxial_zigzag_compression': {'transverse_strain_percent': 100 * state['eyy'], 'gap_eV': G['gap_eV'],
                                        'vbm_k_crystal': G['vbm_k_crystal'], 'cbm_k_crystal': G['cbm_k_crystal']},
        'evidence_ids': list(EVIDENCE)}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    rows = '\n'.join(f"| {100 * s:+.2f} | {E['VB_Gamma'] - E['VB_K']:+.3f} | {E['CB_Q'] - E['CB_K']:+.3f} | "
                     f"{E['CB_K'] - E['VB_K']:.3f} | {E['VB_K']:.3f} | {E['CB_K']:.3f} |" for s, E in L['rows'])
    soc0 = min(L['rows'], key=lambda x: abs(x[0]))[1]
    lo, hi = L['direct_gap_strain_range_percent']
    dp = L['deformation_potentials_eV']
    (APP/'report.md').write_text(f"""# Monolayer MoS2 under strain (reference solution)

**Protocol.** PBE, PSLibrary 1.0.0 ultrasoft pseudopotentials, 45 Ry (density 360 Ry), 12x12 Gamma-centred
mesh, 20 A cell height. The lattice parameter, S height, in-plane S shift and transverse strain come from
scalar-relativistic calculations with relaxed internal coordinates; band energies come from fully relativistic
(spin-orbit) scf + band calculations at those geometries. Band energies are referred to the vacuum level, the
plateau of the planar-averaged electrostatic potential (the raw pw.x eigenvalues carry a reference that moves with
the cell area). Cutoff: band-energy differences agree to < 1 meV between 45 and 55 Ry, and the energy minimum
E(a) is the same at 45, 55 and 70 Ry; the analytic in-plane stress, however, still carries a Pulay stress of about
+0.25 N/m at 45-55 Ry (it agrees between 45 and 55 Ry but moves at 70 Ry; at 35 Ry it is off by ~2 N/m). The cell
height (15-30 A) changes vacuum-referenced edges by < 1 meV, and a 9x9 scf mesh reproduces 12x12 band energies to
0.1 meV.

**Lattice parameter.** a0 = {a0:.4f} A from the minimum of E(a) of relaxed-ion cells at 3.17, 3.18 and 3.19 A
(the analytic stress at 45 Ry would put the zero ~0.005 A lower). The relaxed S height falls linearly with a
(out-of-plane Poisson response), so each strained cell was given its relaxed height. The stress the lab reports
for the unstrained cell at a0, {offset:+.3f} N/m, is the Pulay offset; it is isotropic in the plane and was subtracted
in the zero-transverse-stress condition below.

**Biaxial strain (spin-orbit, relaxed ions, eV vs vacuum):**

| strain % | VB(Gamma)-VB(K) | CB(Q)-CB(K) | Eg(K-K) | VB(K) | CB(K) |
|---|---|---|---|---|---|
{rows}

At zero strain the gap is direct at K ({soc0['CB_K'] - soc0['VB_K']:.3f} eV). The spin-orbit splitting of the K
valence band ({soc0['VB_K_split']:.3f} eV) raises its upper branch by about half the splitting while the Gamma state is
hardly split, which keeps the VBM at K up to {hi:+.2f} %. Under tension the Gamma valence state (Mo dz2 + S pz,
out-of-plane, very sensitive to the S height) rises and the VBM moves K -> Gamma at {hi:+.2f} %. Under compression the
Q-valley conduction minimum (on Gamma-K, t ~ 0.51) falls and the CBM moves K -> Q at {lo:+.2f} %. The direct range is
therefore [{lo:+.2f}, {hi:+.2f}] %. Spin-orbit coupling also splits the Q conduction states (the lower branch drops by
~35 meV), so a correction applied only to the K valence band misplaces the compressive end; in scalar-relativistic
PBE the tensile crossover sits almost at zero strain.

Deformation potentials at zero strain (per unit biaxial strain, vacuum-referenced): VB(K) {dp['valence_K']:.2f} eV,
CB(K) {dp['conduction_K']:.2f} eV, Eg(K-K) {dp['gap_K']:.2f} eV. Clamping the S height changes them by 20-40 % and
moves both crossovers; the raw eigenvalues (no vacuum alignment) give a positive VB(K) deformation potential.

**Uniaxial stress along zigzag (eps_xx = -4 %).** Zero transverse stress at eps_yy = {100 * state['eyy']:+.3f} %
(relaxed ions, S shift along armchair {state['shift'][1]:+.4f} A). Uniaxial strain splits the six Q valleys into the pair
along k_x and the four at +-60 degrees; the pair along k_x is lower here, so the CBM is at {G['cbm_k_crystal']} (crystal) and
the VBM at {G['vbm_k_crystal']}; the gap is {G['gap_eV']:.3f} eV (indirect). A band path along the 60-degree Gamma-K line only would
miss the lowest valley.

**Limits.** PBE underestimates the quasiparticle gap by ~0.8-1 eV and has no excitonic effects (optical gap
~1.9 eV in experiment); crossover strains in experiments also depend on substrate, temperature and doping.
The strains are relative to this protocol's a0 (experiment: 3.16 A), static lattice, no phonon renormalisation.
""")
    lab('submit', submission)
    for _ in range(120):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
