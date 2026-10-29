"""Reference solution for dense-hydrogen-llt-crossover; runs as the Harbor oracle and in check_task."""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
KB = 8.617333262e-5
EVID = []


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(req):
    r = lab('evaluate', req)
    while r['status'] == 'running':
        time.sleep(1)
        r = lab(f"job/{r['id']}")
    if r['status'] != 'complete':
        raise RuntimeError(str(r))
    EVID.append(r['id'])
    return r['result']


THERMO = {}


def thermo(n, t, rs, t0=10.0):
    """Mean pressure (GPa) and c_v/k_B after discarding the first t0 ps (start-up transient)."""
    r = buy({'kind': 'thermo', 'N': n, 'T_K': t, 'rs': rs})
    a = np.array(r['rows'])
    time_ps = (r['first_step'] + r['step_interval'] * np.arange(len(a))) * r['timestep_ps']
    sel = time_ps >= t0
    out = {'P': a[sel, 3].mean() / 1e4, 'cv': np.var(a[sel, 2]) / (n * KB * KB * t * t),
           'cv_all': np.var(a[:, 2]) / (n * KB * KB * t * t)}
    THERMO[(n, t, rs)] = out
    return out


def peak(xs, ys, npts=5):
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    i = int(np.argmax(ys))
    lo = max(0, min(i - npts // 2, len(xs) - npts))
    hi = min(len(xs), lo + npts)
    c = np.polyfit(xs[lo:hi], ys[lo:hi], 2)
    if c[0] < 0 and xs[lo] <= -c[1] / (2 * c[0]) <= xs[hi - 1]:
        return float(-c[1] / (2 * c[0]))
    return float(xs[i])


def flattest(x, p):
    x, p = np.asarray(x, float), np.asarray(p, float)
    s = np.diff(p) / np.diff(x)
    xm = 0.5 * (x[1:] + x[:-1])
    xs = peak(xm, s)
    return xs, float(np.interp(xs, x, p)), float(s.max())


def isotherm(n, t, rr):
    rows = [thermo(n, t, r) for r in rr]
    return [float(r) for r in rr], [x['P'] for x in rows], [x['cv'] for x in rows]


def main():
    # 1. Widom line: c_v maximum of the 2048-atom isotherms (transient removed).
    widom_rs = {'1100': ['1.44', '1.4425', '1.445', '1.4475', '1.45', '1.4525', '1.455', '1.4575'],
                '1200': ['1.45', '1.455', '1.4575', '1.46', '1.4625', '1.465', '1.4675', '1.47'],
                '1400': ['1.49', '1.495', '1.50', '1.505', '1.51', '1.515', '1.52']}
    widom, loops, notes = {}, {}, []
    for t, rr in widom_rs.items():
        x, p, c = isotherm(2048, int(t), rr)
        widom[t] = peak(x, c)
        loops[f'2048/{t}'] = flattest(x, p)[2]
    # 2. First order vs crossover: mechanical instability (dP/drs > 0) of the isotherm in the largest cell.
    sets = {900: ['1.41', '1.4125', '1.415', '1.4175', '1.42', '1.4225'],
            950: ['1.425', '1.4275', '1.43', '1.435', '1.44'],
            1000: ['1.425', '1.4275', '1.43', '1.4325', '1.435', '1.4375', '1.44']}
    cross = {}
    for t, rr in sets.items():
        x, p, c = isotherm(2048, t, rr)
        cross[t] = flattest(x, p)
        loops[f'2048/{t}'] = cross[t][2]
    for n in (128, 256):
        x, p, c = isotherm(n, 1000, ['1.43', '1.4325', '1.435', '1.4375'])
        loops[f'{n}/1000'] = flattest(x, p)[2]
    classification = {str(t): ('first_order' if loops[f'2048/{t}'] > 0 else 'crossover') for t in (900, 950, 1000, 1100)}
    # 3. Structure: long-range order of the molecular side at 900 K, compared with liquids of two sizes.
    sk = {}
    for n, t, rs in [(2048, 900, '1.425'), (2048, 900, '1.41'), (2048, 1100, '1.46'), (128, 1100, '1.46')]:
        r = buy({'kind': 'structure', 'N': n, 'T_K': t, 'rs': rs})
        late = [f[1] for f in r['frames'] if f[0] > 30]
        sk[(n, t, rs)] = (float(np.median(late)), r['msd'][-1] if r['msd'] else None)
    liquid_ref = sk[(2048, 1100, '1.46')][0]
    phase = 'crystal' if sk[(2048, 900, '1.425')][0] > 5 * liquid_ref else 'liquid'
    # 4. Nuclear quantum effects at 1000 K: crossover (maximal compressibility) of PIMD vs classical N = 2048.
    eos = []
    catalog = json.loads((APP/'catalog.json').read_text())
    pimd_rs = [d['rs'] for d in catalog['pimd_available']['1000K'] if 1.4599 <= float(d['rs']) <= 1.5101]
    for rs in pimd_rs:
        r = buy({'kind': 'pimd', 'T_label': '1000K', 'rs': rs})
        ps = [np.array(v)[len(v) // 5:, 5].mean() for v in r['replicas'].values()]
        eos.append((float(rs), float(np.mean(ps))))
    rs_q, p_q, _ = flattest(*zip(*eos))
    shift_1000 = p_q - cross[1000][1]
    # The 900K-labelled PIMD files: check what they contain before using them.
    temps, same = [], []
    for rs in ['1.48', '1.49']:
        r9 = buy({'kind': 'pimd', 'T_label': '900K', 'rs': rs})
        for v in r9['replicas'].values():
            temps.append(np.array(v)[:, 3].mean())
    shift_900 = None if abs(np.mean(temps) - 900) > 50 else 'unexpected'
    assert shift_900 is None
    submission = {'classification': classification, 'low_density_phase_900K': phase, 'widom_rs_N2048': widom,
                  'pimd_crossover_rs_1000K': rs_q, 'nqe_crossover_pressure_shift_GPa': {'1000': shift_1000, '900': shift_900},
                  'evidence_ids': EVID}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# Dense hydrogen: crossover vs first-order LLT (reference solution)

c_v/k_B = var(E_tot)/(N k_B^2 T^2) from the stitched NVT series after discarding the first 10 ps (the start-up
transient from the initial configuration inflates the whole-run variance and moves the peak; with whole runs the
c_v maximum shifts to lower rs). Widom line (N = 2048 c_v maximum): {', '.join(f'{t} K: rs = {v:.4f}' for t, v in widom.items())}.

First order vs crossover: maximal dP/drs along each isotherm (positive = mechanically unstable, van der Waals loop):
{json.dumps({k: round(v) for k, v in loops.items()})}. N = 128 and 256 show loops at 1000 K, which disappear at N = 2048,
so the small-cell loops are finite-size artefacts; only 900 K keeps a loop at N = 2048 -> first order at 900 K,
crossover at 950, 1000 and 1100 K (critical end point between 900 and 950 K, at the melting line).

900 K low-density side: max_k S(k)/N = {sk[(2048, 900, '1.425')][0]:.3f} at N = 2048 after 30 ps, vs {liquid_ref:.3f} for the
1100 K molecular liquid (N = 2048) and {sk[(128, 1100, '1.46')][0]:.3f} for N = 128 at 1100 K (small cells give large S/N even for
liquids, so only the N = 2048 comparison discriminates). The molecular side is a defective molecular crystal ({phase}) even
though atoms still diffuse (MSD at the longest lag {sk[(2048, 900, '1.425')][1]}): the 900 K transition is freezing/melting, not an LLT.

NQE at 1000 K: PIMD (N = 512) crossover at rs = {rs_q:.4f}, P = {p_q:.1f} GPa vs classical N = 2048 at rs = {cross[1000][0]:.4f},
P = {cross[1000][1]:.1f} GPa: shift {shift_1000:.1f} GPa (to lower pressure, lower density) and the PIMD EOS shows no plateau.
The files labelled 900K are byte-identical to 1000K files (Temp column {np.mean(temps):.0f} K), so no 900 K PIMD data exist and the
900 K quantum shift cannot be determined. AIMD was not used; the c_v peak height and the Widom-line pressures are
model (MACE-PBE) and protocol specific; the thermodynamic-limit critical point is only bracketed (900-950 K).
""")
    lab('submit', submission)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
