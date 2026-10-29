"""Reference solution for tan-pressure-kappa-channels; runs as the Harbor oracle and in check_task.

Route: read the archived protocol (CONTROL: 16^3 grid, scalebroad 0.1); establish the row order of the archived
phonon-electron tables (ShengBTE irreducible order, third grid index fastest) from their frequency column and
validate the mapping against the lab's own phonon-electron channel at 60 GPa; interpolate the four-phonon and
phonon-electron rates mode by mode in pressure for the withheld pressures; transfer the four-phonon rates to
600 K with the mode-wise 600 K / 300 K ratio measured at 20 and 100 GPa; solve the iterative BTE with the
three-phonon scattering of the target pressure and the transferred rates as extra diagonal scattering.
"""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CAT = json.loads((APP/'catalog.json').read_text())
G16 = CAT['grid_16']
EVIDENCE = []


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def evaluate(req):
    r = lab('evaluate', req)
    while r['status'] == 'running':
        time.sleep(2)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(r)
    EVIDENCE.append(r['id'])
    return r['result']


def kappa(p, T, channels, extra=None, solver='iterative'):
    req = {'kind': 'kappa', 'pressure_GPa': p, 'mesh': 16, 'T_K': T, 'scalebroad': 0.1, 'channels': channels,
           'solver': solver}
    if extra is not None:
        req['extra_rates_ps'] = [[float('%.6g' % v) for v in row] for row in np.asarray(extra)]
    res = evaluate(req)
    return np.array(res['kappa_W_mK']), res


def fourph(p, T):
    r = evaluate({'kind': 'rates', 'pressure_GPa': p, 'mesh': 16, 'T_K': T, 'scalebroad': 0.1})
    return np.array(r['rate_4ph_pp_ps']) + np.array(r['rate_4ph_pm_ps']) + np.array(r['rate_4ph_mm_ps'])


def ratio(a, b):
    return np.where(b > 0, a / np.where(b > 0, b, 1), 1.0)


def orders():
    """Irreducible order produced by scanning the grid with a given index running fastest."""
    g2i = np.array(G16['grid_to_irreducible'])
    out = {}
    for name, fastest in (('i_fastest', 'i'), ('k_fastest', 'k')):
        seen, perm = set(), []
        for a in range(16):
            for b in range(16):
                for c in range(16):
                    i, j, k = (c, b, a) if fastest == 'i' else (a, b, c)
                    s = int(g2i[i + 16 * j + 256 * k])
                    if s not in seen:
                        seen.add(s)
                        perm.append(s)
        out[name] = np.array(perm)
    return out


def phel_table(text, perm):
    a = np.array([[float(x) for x in line.split()] for line in text.strip().split('\n')])
    n = len(perm)
    w, r = np.zeros((n, 6)), np.zeros((n, 6))
    w[perm] = a[:, 0].reshape(6, n).T
    r[perm] = a[:, 1].reshape(6, n).T
    return w, r


def main():
    control = evaluate({'kind': 'archived_file', 'pressure_GPa': 0, 'file': 'CONTROL'})['text']
    assert 'scalebroad=0.1' in control.replace(' ', '') and 'ngrid(:)=161616' in control.replace(' ', '')
    texts = {p: evaluate({'kind': 'archived_file', 'pressure_GPa': p, 'file': 'BTE.sr_phel_metal'})['text']
             for p in (40, 60, 100)}
    harm = evaluate({'kind': 'harmonic', 'pressure_GPa': 60, 'mesh': 16})
    w60 = np.array(harm['omega_rad_ps'])
    # Which irreducible order do the archived rows follow? Compare the archived frequency column (rad/ps).
    mism = {}
    for name, perm in orders().items():
        w, _ = phel_table(texts[60], perm)
        ok = w60 > 1
        mism[name] = float(np.median(np.abs(w[ok] - w60[ok]) / w60[ok]))
    best = min(mism, key=mism.get)
    perm = orders()[best]
    phel = {p: phel_table(t, perm)[1] for p, t in texts.items()}
    F = {(p, T): fourph(p, T) for (p, T) in ((20, 300), (60, 300), (100, 300), (20, 600), (100, 600))}
    # Validate the mapping and units against the lab's own channel at 60 GPa.
    k_lab, _ = kappa(60, 300, ['3ph', 'isotope', '4ph', 'phel'])
    k_mine, _ = kappa(60, 300, ['3ph', 'isotope', '4ph'], extra=phel[60])
    assert abs(k_mine[0, 0] / k_lab[0, 0] - 1) < 1e-3, (k_mine, k_lab)
    # Withheld pressures: mode-wise linear interpolation of the rates in pressure.
    f80 = 0.5 * (F[(60, 300)] + F[(100, 300)])
    e80 = 0.5 * (phel[60] + phel[100])
    f40 = 0.5 * (F[(20, 300)] + F[(60, 300)])
    k80, r80 = kappa(80, 300, ['3ph', 'isotope'], extra=f80 + e80)
    k80n, _ = kappa(80, 300, ['3ph', 'isotope'], extra=f80)
    k40, _ = kappa(40, 300, ['3ph', 'isotope', 'phel'], extra=f40)
    # 600 K: four-phonon rates follow the mode-wise 600/300 K ratio measured at 20 and 100 GPa.
    rT = 0.5 * (ratio(F[(20, 600)], F[(20, 300)]) + ratio(F[(100, 600)], F[(100, 300)]))
    k60h, _ = kappa(60, 600, ['3ph', 'isotope'], extra=F[(60, 300)] * rT + phel[60])
    k80rta, _ = kappa(80, 300, ['3ph', 'isotope'], extra=f80 + e80, solver='RTA')
    frac = 1 - k80[0, 0] / k80n[0, 0]
    sub = {'kappa_xx_80GPa_300K_W_mK': float(k80[0, 0]), 'kappa_zz_80GPa_300K_W_mK': float(k80[2, 2]),
           'kappa_xx_40GPa_300K_W_mK': float(k40[0, 0]), 'kappa_xx_60GPa_600K_W_mK': float(k60h[0, 0]),
           'phel_resistivity_fraction_80GPa_300K': float(frac), 'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(sub, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# theta-TaN under pressure: scattering channels (reference solution)

Protocol: the archived CONTROL files fix a 16x16x16 grid and scalebroad = 0.1; three-phonon scattering is solved
iteratively, four-phonon, isotope and phonon-electron scattering enter the diagonal (FourPhonon convention).

Archived phonon-electron tables: rows are branch-major over the ShengBTE irreducible list, which scans the grid with
the third index fastest (median frequency mismatch {mism['k_fastest']:.3f} against {mism['i_fastest']:.3f} for the
other order). The frequency column (rad/ps) comes from a different phonon calculation (EPW) and agrees only to
~1%, so rows are mapped by order, not by frequency. The mapping was validated against the lab's own phonon-electron
channel at 60 GPa (kappa_xx {k_mine[0, 0]:.1f} vs {k_lab[0, 0]:.1f} W/mK); the second column is 1/tau in ps^-1.

At 40 and 80 GPa the fourth-order constants (and at 80 GPa the phonon-electron table) are withheld. Rates were
interpolated mode by mode in pressure (60/100 GPa for 80 GPa; 20/60 GPa for the 40 GPa four-phonon rates) and
combined with the target pressure's own three-phonon scattering in one iterative solution; conductivities are not
combined by Matthiessen's rule. 600 K four-phonon rates at 60 GPa follow the 600/300 K ratio measured at 20 and
100 GPa (close to, but below, the classical T^2 factor).

Results (W/mK): kappa_xx(80 GPa) = {k80[0, 0]:.0f}, kappa_zz(80 GPa) = {k80[2, 2]:.0f} (RTA would give
{k80rta[0, 0]:.0f}); kappa_xx(40 GPa) = {k40[0, 0]:.0f}; kappa_xx(60 GPa, 600 K) = {k60h[0, 0]:.0f}; phonon-electron
share of the 80 GPa thermal resistivity = {frac:.3f}; released-pressure full-physics kappa_xx at 60 GPa =
{k_lab[0, 0]:.0f}. Four-phonon scattering (aaaa processes, made important because the large acoustic-optical gap
and acoustic bunching freeze out three-phonon channels) roughly halves kappa; with three-phonon scattering alone and
scalebroad 0.1 many modes have no allowed process and the iterative solution does not converge. Phonon-electron
scattering grows strongly above 20 GPa and produces the non-monotonic pressure dependence. The interpolated rates
carry a few-percent uncertainty; the pressure of the maximum is not resolved better than the 20 GPa spacing.
""")
    lab('submit', sub)
    for _ in range(120):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
