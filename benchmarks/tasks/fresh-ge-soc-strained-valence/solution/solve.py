"""Reference solution for fresh-ge-soc-strained-valence; runs as the Harbor oracle and in check_task.

Route: a0 from spin-orbit pressures; labelled Gamma levels at a0; the out-of-plane relaxation at -1 % from two cheap
scalar-relativistic stress runs (sigma_zz = 0 by interpolation; the ratio D = -eps_perp/eps_par is then used along the
biaxial path); labelled spin-orbit Gamma levels at -1 %, at +-0.5 % (b from the symmetric difference, which cancels
the light-hole/split-off nonlinearity) and at -2 % / -2.5 % (s/split-off crossing, by parity labels).
"""
import json
import os
import subprocess
import time
from pathlib import Path

APP = Path(os.environ.get('APP_DIR', '/app'))
GAMMA = [0.0, 0.0, 0.0]


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def calc(a, c, soc=True, bands=False, ecut=60, k=8):
    req = {'kind': 'pw', 'a_par_A': round(a, 5), 'c_A': round(c, 5), 'spin_orbit': soc, 'ecutwfc_Ry': ecut, 'kmesh': k,
           'kshift': 0, 'occupations': 'smearing', 'degauss_Ry': 0.005,
           'band_kpoints_inverse_A': [GAMMA] if bands else [], 'symmetry_labels': bands}
    r = lab('evaluate', req)
    while r['status'] == 'running':
        time.sleep(10)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError('job %d failed: %s' % (r['id'], r.get('error')))
    return r['id'], r['result']


# ---------------------------------------------------------------- estimators (also applied by make_task.py)
def gamma_levels(result):
    """Identify the Gamma levels of a labelled spin-orbit run by character.

    Among spinor bands 23-30 (the s-like antibonding level and the p-like manifold, in a strain-dependent order) the
    odd-parity Kramers pair is the s-like level. Of the three even p-like levels, in the cubic crystal Gamma8+ is the
    fourfold one and Gamma7+ the split-off; in the tetragonal crystal the light hole and the split-off level carry the
    same label (they mix), so the heavy hole is the level whose label differs from the other two, the split-off is
    the lower of the two same-label levels and the light hole the upper."""
    band = next(b for b in result['band_kpoints'] if max(abs(v) for v in b['k_cartesian_inverse_A']) == 0)
    s, p = None, []
    for lev in band['symmetry']['levels']:
        i, j = lev['bands']
        if 23 <= i <= 30:
            name = lev['label'].split()[0]
            if name.endswith('-'):
                s = lev['energy_eV']
            else:
                p.append((lev['energy_eV'], name, (j - i + 1) // 2))
    out = {'s': s}
    if len(p) == 2:                                 # cubic: Gamma8+ (two Kramers pairs) and Gamma7+
        g8 = next(e for e, n, m in p if m == 2)
        g7 = next(e for e, n, m in p if m == 1)
        out.update(hh=g8, lh=g8, so=g7)
    else:
        names = [n for e, n, m in p]
        hh = next(e for e, n, m in p if names.count(n) == 1)
        same = sorted([e for e, n, m in p if names.count(n) == 2], reverse=True)
        out.update(hh=hh, lh=same[0], so=same[1])
    return out


def interp_root(x, y):
    """x where the straight line through two points (x0, y0), (x1, y1) crosses zero."""
    (x0, x1), (y0, y1) = x, y
    return x0 - y0 * (x1 - x0) / (y1 - y0)


def b_symmetric(dE_plus, X_plus, dE_minus, X_minus):
    """b from E(HH)-E(LH) = b X + O(X^2) at two strains of opposite sign (the quadratic term cancels)."""
    return (dE_plus - dE_minus) / (X_plus - X_minus)


def main():
    evidence = []
    # 1. zero-pressure lattice parameter: spin-orbit pressures bracketing a0 (cubic, Gamma-centred k8, 60 Ry)
    pts = []
    for a in (5.75, 5.79):
        i, r = calc(a, a)
        evidence.append(i)
        pts.append((a, r['pressure_GPa']))
    a0 = interp_root([p[0] for p in pts], [p[1] for p in pts])
    # 2. labelled levels at a0
    i, r = calc(a0, a0, bands=True)
    evidence.append(i)
    L0 = gamma_levels(r)
    delta0 = L0['hh'] - L0['so']
    s_level = L0['s'] - L0['hh']
    # 3. out-of-plane relaxation at eps_par = -1 %: sigma_zz = 0 from two scalar-relativistic stress runs
    e1 = -0.01
    szz = []
    for D in (0.70, 0.82):
        i, r = calc(a0 * (1 + e1), a0 * (1 - D * e1), soc=False)
        evidence.append(i)
        szz.append((-D * e1, r['stress_GPa'][2][2]))
    ez1 = interp_root([s[0] for s in szz], [s[1] for s in szz])
    D = -ez1 / e1

    def strained(e):
        ez = -D * e
        i, r = calc(a0 * (1 + e), a0 * (1 + ez), bands=True)
        evidence.append(i)
        L = gamma_levels(r)
        return ez, L, r
    ez, L1, r1 = strained(e1)
    c_over_a = (1 + ez) / (1 + e1)
    hh_lh = L1['hh'] - L1['lh']
    so_hh = L1['so'] - L1['hh']
    # 4. b in the small-strain limit from +-0.5 %
    split = {}
    for e in (0.005, -0.005):
        ezz, L, _ = strained(e)
        split[e] = (L['hh'] - L['lh'], 2 * e - 2 * ezz)
    b = b_symmetric(*split[0.005], *split[-0.005])
    # 5. s/split-off crossing on the relaxed path, tracked by parity labels
    track = {0.0: L0['s'] - L0['so'], e1: L1['s'] - L1['so']}
    for e in (-0.02, -0.025, -0.03):
        _, L, _ = strained(e)
        track[e] = L['s'] - L['so']
        if track[e] > 0:
            break
    xs = sorted(track)
    lo = max(x for x in xs if track[x] > 0) if any(track[x] > 0 for x in xs) else xs[0]
    hi = min(x for x in xs if x > lo)
    crossing = 100 * interp_root([lo, hi], [track[lo], track[hi]])
    submission = {'a0_A': a0, 'delta0_at_a0_eV': delta0, 's_level_at_a0_eV': s_level,
                  'c_over_a_at_minus1pct': c_over_a, 'hh_minus_lh_at_minus1pct_eV': hh_lh,
                  'so_minus_hh_at_minus1pct_eV': so_hh, 'b_shear_deformation_potential_eV': b,
                  's_so_crossing_strain_percent': crossing, 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    lin_b = hh_lh / (2 * e1 - 2 * ez)
    (APP/'report.md').write_text(f"""# Ge with spin-orbit coupling: strained valence levels (reference solution)

Protocol: PBE, SG15 fully relativistic ONCV file, noncollinear spin-orbit, 60 Ry, Gamma-centred 8x8x8 mesh,
Marzari-Vanderbilt 0.005 Ry (the protocol is a semimetal, so smearing is required; Gamma levels change by < 2 meV
between k6 and k12 and between 40 and 80 Ry). Band energies are read at Gamma itself (a shifted mesh never contains
Gamma; the inverted s band disperses by ~0.1 eV within 0.02 1/A of Gamma).

a0 = {a0:.4f} A from P = 0 (pressures {pts}); the stress needs >= 60 Ry (40 Ry gives a0 0.03 A too small). All
strains are measured from this a0, not from the experimental 5.658 A (which is about -2 % in this protocol).

Identification by bands.x double-group labels: at a0 the levels between the 3d bands and the conduction states are
Gamma8+ (fourfold, p3/2), Gamma7+ (split-off, p1/2) and the odd-parity Gamma7- s-like antibonding level, which lies
{-s_level:.3f} eV BELOW Gamma8+ (and below the split-off level): PBE inverts the gap, E_F cuts the Gamma8+ manifold
and the "first state above the valence top" (Gamma6-, p-like antibonding, ~2.3 eV) is not the s state. Delta0 =
{delta0:.4f} eV, s level = {s_level:.4f} eV.

Under (001) strain the light hole and the split-off level carry the same D4h label (they mix); the heavy hole is the
level with the other label (bands.x calls it G_6+, which is NOT the Koster Gamma7+ name of m_j = 3/2 states, so names
alone mislead). Out-of-plane relaxation: sigma_zz = 0 interpolated between two scalar-relativistic stress runs,
D = -eps_perp/eps_par = 2 C12/C11 = {D:.4f}; c/a at -1 % = {c_over_a:.5f}. At -1 %: E(HH)-E(LH) = {hh_lh:.4f} eV
(HH on top), E(SO)-E(HH) = {so_hh:.4f} eV.

b: E(HH)-E(LH) = b (exx + eyy - 2 ezz) at first order. The splitting at -1 % is not small compared with Delta0, and
the light hole is pushed by the split-off level (6x6 Bir-Pikus), so the linear formula at -1 % gives {lin_b:.3f} eV;
the symmetric difference of +-0.5 % (where HH and LH change order) cancels the quadratic term: b = {b:.3f} eV.

The s level rises fast under compression (large hydrostatic gap deformation potential) and, having odd parity, crosses
the even split-off level without repulsion. Tracking E(s) - E(SO) by parity: {track}; crossing at eps_par =
{crossing:.3f} %.

Experiment: Delta0 = 0.30 eV, E0 = +0.80-0.89 eV (direct gap, s level above Gamma8+), b = -2.2 to -2.9 eV. The PBE
s level is inverted (no hybrid or GW, unlike the paired HSE study); static lattice, bulk (no confinement); masses near
Gamma are unreliable in an inverted protocol and are not claimed.
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
