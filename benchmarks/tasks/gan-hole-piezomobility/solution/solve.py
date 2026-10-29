"""Reference solution for gan-hole-piezomobility; runs as the Harbor oracle and in check_task."""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
CELLS = {s['state']: np.array(s['cell_bohr']) for s in CATALOG['states']}
K = np.array(CATALOG['kgrid_crystal'])
KEYS = {'xx': (0, 0), 'yy': (1, 1), 'zz': (2, 2), 'yz': (1, 2), 'xz': (0, 2), 'xy': (0, 1)}
EVIDENCE = []


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(1)
        r = lab(f"job/{r['id']}")
    if r['status'] != 'complete':
        raise RuntimeError(r)
    EVIDENCE.append(r['id'])
    return r['result']


def tensor(result, T='300'):
    m = np.array(result['ibte_drift_mobility'][T])
    return 0.5 * (m + m.T)


def voigt(state):
    """Symmetric small strain of a state from its cell (rows = lattice vectors, A = A0 (I + eps))."""
    F = np.linalg.solve(CELLS['unstrained'], CELLS[state])
    e = 0.5 * (F + F.T) - np.eye(3)
    return [round(float(v), 4) + 0.0 for v in (e[0, 0], e[1, 1], e[2, 2], 2 * e[1, 2], 2 * e[0, 2], 2 * e[0, 1])]


def as_dict(m):
    return {k: float(m[i, j]) for k, (i, j) in KEYS.items()}


def geometric_next(m0, m1, m2):
    """Third point of a saturating series: increments shrink geometrically while they keep their sign and shrink;
    otherwise the last increment is continued."""
    d1, d2 = m1 - m0, m2 - m1
    out = m2 + d2
    same = (np.sign(d1) == np.sign(d2)) & (np.abs(d2) < np.abs(d1)) & (np.abs(d1) > 1e-9)
    out[same] = (m2 + d2 ** 2 / np.where(np.abs(d1) > 1e-9, d1, 1.0))[same]
    return out


def mirror(m, comp):
    out = m.copy()
    i, j = {4: (1, 2), 5: (0, 2), 6: (0, 1)}[comp]
    out[i, j] = out[j, i] = -m[i, j]
    return out


def gamma_levels(bands):
    e = np.array(bands['energies_eV'])[:, :6]
    g = int(np.where(np.all(np.abs(K) < 1e-9, axis=1))[0][0])
    z = int(np.where(np.all(np.abs(K - [0, 0, 0.1]) < 1e-9, axis=1))[0][0])
    vbm = e[g].max()
    return (np.sort(e[g])[::-1][0::2] - vbm) * 1000, (np.sort(e[z])[::-1][0::2] - vbm) * 1000


def main():
    # GW bands along the c-axis series: band order at Gamma and the light-z (crystal-field) level.
    levels = {s: gamma_levels(buy({'kind': 'bands', 'state': s}))
              for s in ('unstrained', 'compress_1%_e3', 'compress_2%_e3', 'compress_3%_e3', 'expand_3%_e3')}
    # Zero-field IBTE mobilities (the iterative solution, not SERTA).
    ibte = ['unstrained', 'compress_1%_e3', 'compress_2%_e3', 'expand_1%_e3', 'expand_2%_e3', 'expand_1%_e2',
            'expand_2%_e2', 'compress_3%_e4', 'expand_3%_e6', 'compress_2%_e4', 'compress_2%_e5', 'compress_2%_e6',
            'compress_2%_e1']
    res = {s: buy({'kind': 'mobility', 'state': s, 'solver': 'ibte'}) for s in ibte}
    mu = {s: tensor(r) for s, r in res.items()}
    flagged = sorted(s for s, r in res.items() if r['ibte_solver_log']['iterations'] >= 1000)
    worst = max(r['ibte_solver_log']['last_max_change_cm2_per_Vs'] for r in res.values())
    u = mu['unstrained']
    # Withheld states.
    y1, y2 = mu['compress_1%_e3'][2, 2], mu['compress_2%_e3'][2, 2]
    zz_c3 = float(np.sqrt(y2 * (2 * y2 - y1)))   # past the crossing: between saturation (y2) and no saturation
    c3 = geometric_next(u, mu['compress_1%_e3'], mu['compress_2%_e3'])
    c3[2, 2] = zz_c3
    preds = {'compress_3%_e3': as_dict(c3),
             'expand_3%_e3': as_dict(geometric_next(u, mu['expand_1%_e3'], mu['expand_2%_e3'])),
             'expand_3%_e2': as_dict(geometric_next(u, mu['expand_1%_e2'], mu['expand_2%_e2'])),
             'expand_3%_e4': as_dict(mirror(mu['compress_3%_e4'], 4)),
             'compress_3%_e6': as_dict(mirror(mu['expand_3%_e6'], 6))}
    assert set(preds) == set(CATALOG['to_predict'])
    # Design answers among the 2 % states (shear partners are mirror images, so one sign suffices).
    two = {s: m for s, m in mu.items() if '_2%_' in s}
    best_zz = max(two, key=lambda s: two[s][2, 2])
    inplane = {s: float(np.linalg.eigvalsh(m[:2, :2]).max()) for s, m in two.items()}
    best_in = max(inplane, key=inplane.get)
    w, v = np.linalg.eigh(two[best_in][:2, :2])
    angle = float(np.degrees(np.arctan2(v[1, -1], v[0, -1])) % 180)
    enh = {T: tensor(res[best_zz], T)[2, 2] / tensor(res['unstrained'], T)[2, 2] for T in ('100', '300', '500')}
    zz = {T: tensor(res[best_zz], T)[2, 2] for T in ('300', '500')}
    n_best = float(np.log(zz['300'] / zz['500']) / np.log(5 / 3))
    n_u = float(np.log(tensor(res['unstrained'], '300')[2, 2] / tensor(res['unstrained'], '500')[2, 2]) / np.log(5 / 3))
    submission = {'mobility_predictions_cm2_per_Vs': preds,
                  'best_mu_zz_2pct': {'strain_voigt': voigt(best_zz), 'mu_zz_cm2_per_Vs': float(two[best_zz][2, 2])},
                  'best_inplane_2pct': {'strain_voigt': voigt(best_in), 'mu_inplane_max_cm2_per_Vs': inplane[best_in]},
                  'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    lv = {s: ', '.join(f'{x:.0f}' for x in g) for s, (g, z) in levels.items()}
    lz = {s: ', '.join(f'{x:.0f}' for x in z) for s, (g, z) in levels.items()}
    (APP/'report.md').write_text(f"""# Strain control of the hole mobility tensor of wurtzite GaN (reference solution)

Conventions: zero-field IBTE drift mobility at 300 K, p = 1e18 cm^-3, symmetrised tensor in the cell frame (x || a1,
z || c). Strains are Voigt vectors computed from each cell; the archive's 'X% e4/e5/e6' states carry eps_ij = X%,
i.e. engineering shear 2X% (compress_2%_e6 has e6 = {voigt('compress_2%_e6')[5]}).

Bands (GW, top three Kramers levels at Gamma in meV relative to the VBM; same at kz = 0.1): {lv}; kz: {lz}. Unstrained
GaN has the two heavy-along-c levels (HH, LH) 7 meV apart, below kT, with the light-along-c crystal-field (CH) level
56 meV lower. Compressing c raises the CH level: it reaches the top between 1% and 2% (the deepest kz = 0.1 level stays
~270 meV below the VBM while the heavy pair falls 70 and 130 meV below the top at 2% and 3%). mu_zz jumps from
{y1:.0f} (1%) to {y2:.0f} (2%); tension along c leaves the order and mu unchanged. In-plane uniaxial and shear strains
split HH/LH by 60-170 meV and put a band that is light along the compressed principal direction on top; the response
is non-analytic at zero strain because HH/LH are quasi-degenerate, so a linear piezomobility tensor fitted to +-1%
fails already at 2-3% (one-sided for e3, saturating in-plane).

Predictions. compress_3%_e3: past the crossing the 2->3% increment of mu_zz lies between zero and the 1->2% increment;
I report the geometric mean of these bounds ({y2:.0f}-{2 * y2 - y1:.0f}), {zz_c3:.0f} cm2/Vs; populations of the heavy
levels fall from ~10% to ~1% and interband scattering into them weakens, so it keeps rising. expand_3%_e3 and
expand_3%_e2 come from saturating (geometric) extrapolation of the 0/1/2% series of the same component. The mirror
planes x->-x and y->-y of 6mm map +eps_yz, +eps_xz, +eps_xy onto their negatives, so the diagonal of mu is even and the
strained off-diagonal element odd in shear: expand_3%_e4 and compress_3%_e6 are the mirror images of compress_3%_e4
and expand_3%_e6 (sign-flipped off-diagonal). Mirror pairs in the archive agree to 1-3% (up to 6% at 1% e4), which is
the numerical noise of these runs. No exact symmetry relates e1 to e2 or e4 to e5 (90 degrees is not a 6mm rotation);
the x<->y swap is only a linear-order relation, so e2 targets were taken from the e2 series itself.

Design. Among the 2% states mu_zz is largest for c-axis compression ({best_zz}, {two[best_zz][2, 2]:.0f} cm2/Vs,
x{two[best_zz][2, 2] / u[2, 2]:.1f}); shear e4/e5 follow (~145). The largest in-plane directional mobility is for
in-plane shear ({best_in}, {inplane[best_in]:.0f} cm2/Vs along {angle:.0f} deg from x, near the compressive principal
axis), larger than uniaxial e1 compression ({inplane.get('compress_2%_e1', float('nan')):.0f}), because shear
compresses one diagonal and stretches the other. Between 300 and 500 K mu_zz falls as T^-{n_u:.1f} unstrained and
T^-{n_best:.1f} under 2% c-axis compression, so the enhancement shrinks with temperature: x{enh['100']:.1f} at 100 K,
x{enh['300']:.1f} at 300 K, x{enh['500']:.1f} at 500 K. SERTA underestimates mu by factors 2-3.5 that depend on the band
character and mis-ranks e4 above e3 for mu_zz, so only IBTE values were used. IBTE runs flagged as unconverged after
1000 iterations ({len(flagged)} bought) have last per-iteration changes <= {worst:.1e} cm2/Vs and are usable.

Not supported: a linear piezomobility tensor beyond |eps| << 1%; the exact strain of the band crossing (between 1 and
2%); combined strains (only single components were computed); quantitative mu_zz at 3% c-axis compression beyond
the bracket above; DFT-versus-GW ordering effects (only GW energies are available).
""")
    lab('submit', submission)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(2)


if __name__ == '__main__':
    main()
