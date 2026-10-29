"""Reference solution for cnt-kohn-anomaly-nonadiabatic; runs as the Harbor oracle and in check_task."""
import json
import math
import os
import re
import subprocess
import tempfile
import time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
KB = 8.617333262145179e-05
RY = 13.605693122994
CM1 = 0.12398419843320026
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
        time.sleep(2)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(r.get('error'))
    EVIDENCE.append(r['id'])
    return r['result']


def table(mode, t):
    """Parse an archived self-energy table; identify the columns from the data, not from the header."""
    lines = buy({'kind': 'self_energy_table', 'mode': mode, 'temperature_K': t})['lines']
    d = np.array([l.split() for l in lines if not l.strip().startswith('#')], float)
    w, wq = d[:, 4] * 1000, d[0, 3] * 1000
    const = [c for c in (5, 6, 7) if np.ptp(d[:, c]) == 0]
    assert const == [6], 'expected exactly one frequency-independent column (the static part)'
    re_, im_ = d[:, 5], d[:, 7]
    assert abs(np.interp(0.0, w, im_)) < 0.05 * np.abs(im_).max()      # Im Pi vanishes at w = 0
    return {'w': w, 'wq': wq, 're': re_, 'im': im_, 'static': d[0, 6]}


def pole(w, re_, im_, wq, wb2):
    def rhs(o):
        g = -wq * np.interp(o, w, im_) / o
        return wb2 + g * g + 2 * wq * np.interp(o, w, re_), g
    grid = np.linspace(max(w[0], 1.0), w[-1], 2000)
    vals = [o * o - rhs(o)[0] for o in grid]
    roots = [brentq(lambda o: o * o - rhs(o)[0], grid[i], grid[i + 1]) for i in range(len(grid) - 1) if vals[i] * vals[i + 1] < 0]
    return roots[-1], rhs(roots[-1])[1]


def main():
    notes = {}
    # DFPT smearing (Ry in Quantum ESPRESSO) as a Fermi-Dirac temperature
    scf = buy({'kind': 'archive_file', 'file': 'dft_inputs/scf.in'})['text']
    degauss = float(re.search(r'degauss\s*=\s*([\d.dDeE+-]+)', scf).group(1).replace('d', 'e').replace('D', 'e'))
    smearing = re.search(r"smearing\s*=\s*'(\w+)'", scf).group(1)
    t_sigma = degauss * RY / KB
    released = CATALOG['requests']['self_energy_table']['released_temperatures_K']
    t_ref = min(released, key=lambda t: abs(t - t_sigma))
    assert abs(t_ref - t_sigma) < 0.1
    # broadening and sampling of the archived self-energy runs
    script = buy({'kind': 'archive_file', 'file': 'dft_inputs/specfun/input.py'})['text']
    eta_div = float(re.search(r'eta = temp \* KB_EV / (\d+)', script).group(1))
    nk_arch = int(float(re.search(r'NKF3 = int\(([\d.e]+)\)', script).group(1)))
    scan = buy({'kind': 'dfpt_scan'})
    # tables: the reference smearing for both modes, and the low-temperature table for validation
    ref20, ref1 = table(20, t_ref), table(1, t_ref)
    low20 = table(20, 50.0)
    pi_sigma = {20: ref20['static'], 1: ref1['static']}
    wq = {20: ref20['wq'], 1: ref1['wq']}
    wb2 = {nu: wq[nu] ** 2 - 2 * wq[nu] * pi_sigma[nu] for nu in (1, 20)}
    base = {'processes': ['k+q'], 'band_pairs': [[1, 1], [1, 2], [2, 1], [2, 2]]}
    eta = lambda t: KB * t / eta_div * 1000

    def evaluate(mode, t, omegas, nk):
        r = buy({'kind': 'self_energy', 'mode': mode, 'temperature_K': t, 'omega_meV': [float(x) for x in omegas],
                 'nk': nk, 'eta_meV': eta(t), **base})
        return np.array(r['re_pi_meV']), np.array(r['im_pi_meV'])
    # validation at 50 K against the archived table, and the effect of a coarse mesh
    idx = [int(np.argmin(np.abs(low20['w'] - x))) for x in (0.0, 115.0, 120.0, 125.0, 130.0)]
    wv = low20['w'][idx]
    r6, i6 = evaluate(20, 50.0, wv, nk_arch)
    r5, i5 = evaluate(20, 50.0, wv, 100000)
    dev6 = float(max(np.abs(r6 - low20['re'][idx]).max(), np.abs(i6 - low20['im'][idx]).max()))
    dev5 = float(max(np.abs(r5 - low20['re'][idx]).max(), np.abs(i5 - low20['im'][idx]).max()))
    assert dev6 < 0.3, dev6
    notes['validation_50K'] = {'nk_1e6_max_dev_meV': dev6, 'nk_1e5_max_dev_meV': dev5}
    # mode 20 at the target temperatures
    grid = np.r_[0.0, np.linspace(110.0, 135.0, 40)]
    omega_na, gamma, static = {}, {}, {}
    for t in (100.0, 300.0):
        r, i = evaluate(20, t, grid, nk_arch)
        static[t] = r[0]
        omega_na['%d' % t], gamma['%d' % t] = pole(grid[1:], r[1:], i[1:], wq[20], wb2[20])
    omega_ad = math.sqrt(wb2[20] + 2 * wq[20] * static[100.0])
    # RBM: check the evaluator against the archived static part at the reference smearing, then locate T*
    rs, _ = evaluate(1, t_ref, [0.0], 100000)
    assert abs(rs[0] - pi_sigma[1]) < 0.3
    ts, fs = [500.0, 520.0, 540.0, 560.0], []
    for t in ts:
        r, i = evaluate(1, t, [0.0, 0.5], 100000)
        gamma0 = -wq[1] * i[1] / 0.5                        # linewidth as Omega -> 0
        fs.append(wb2[1] + gamma0 ** 2 + 2 * wq[1] * r[0])   # Omega^2 in the limit Omega -> 0
    k = next(j for j in range(len(ts) - 1) if fs[j] < 0 <= fs[j + 1])
    tstar = ts[k] + (ts[k + 1] - ts[k]) * (-fs[k]) / (fs[k + 1] - fs[k])
    # threshold model: the archived program with main() parametrised; its parameter a is 1/m*
    src = buy({'kind': 'archive_file', 'file': 'Fig1/coupling_calculator.c'})['text']
    assert 'a = 1/m*' in src
    main_c = r'''int main(int argc, char **argv) {
    double a = atof(argv[1]), freq = atof(argv[2]); double omega_scan[OMEGA_POINTS]; double gamma = 1.0, g = 100.0;
    for (int i = 0; i < OMEGA_POINTS; i++) omega_scan[i] = freq * i / (OMEGA_POINTS - 1.0);
    for (int i = OMEGA_POINTS - 1; i >= 1; i--) {
        gamma = compute_linewidth(gamma, g, a, freq, omega_scan[i]);
        g = compute_coupling(gamma, g, a, N_KPOINTS, freq, omega_scan[i], TEMPERATURE, SIGMA);
        if (i == 1) printf("%.6f\n", g);
        if (isnan(gamma)) gamma = 1.0;
        if (isnan(g)) g = 100.0;
    }
    return 0;
}
'''
    g0 = {}
    with tempfile.TemporaryDirectory() as d:
        (Path(d)/'m.c').write_text(src[:src.index('int main(void)')] + main_c)
        subprocess.run(['gcc', '-O2', str(Path(d)/'m.c'), '-o', str(Path(d)/'m'), '-lm'], check=True)
        for nu in (1, 20):
            out = subprocess.run([str(Path(d)/'m'), '0.5', repr(wq[nu])], capture_output=True, text=True, check=True)
            g0[str(nu)] = float(out.stdout.split()[-1])
    # couplings at the Fermi points, for the comparison with g0
    gf = {}
    for nu in (1, 20):
        c = buy({'kind': 'coupling', 'mode': nu})
        g = np.array(c['k+q']['g'])
        kf = int(np.argmin(np.abs(np.array(c['k_crystal']) - 0.2875)))
        gf[nu] = (float(g[0, 1, kf]), float(g[1, 0, kf]))
    # DFPT smearing scan (Fermi-Dirac) used as a temperature proxy, for the report
    sig_t = [s * 1e-3 * RY / KB for s in scan['sigma_mRy']]
    om1 = [v * CM1 for v in scan['frequency_cm1']['1']]
    proxy = next(sig_t[j + 1] + (sig_t[j] - sig_t[j + 1]) * (-math.copysign(om1[j + 1] ** 2, om1[j + 1]))
                 / (om1[j] ** 2 - math.copysign(om1[j + 1] ** 2, om1[j + 1]))
                 for j in range(len(om1) - 1) if om1[j] > 0 > om1[j + 1])
    submission = {'omega_nonadiabatic_meV': omega_na, 'linewidth_meV': gamma, 'omega_adiabatic_100K_meV': omega_ad,
                  'rbm_instability_temperature_K': tstar, 'threshold_coupling_meV': g0, 'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# Kohn anomaly of the (3,3) nanotube (reference solution)

**Bare frequency.** The DFPT run behind omega_Q used degauss = {degauss} Ry ({smearing} smearing), i.e. a
Fermi-Dirac temperature of {t_sigma:.2f} K, which is the highest released table. DFPT frequencies already contain
the static adiabatic self-energy at that smearing (the archived Fermi-Dirac smearing scan softens the RBM from
{om1[0]:.1f} to {om1[-1]:.1f} meV between 10 and 0.5 mRy), so omega_b^2 = omega_Q^2 - 2 omega_Q Pi(0, T_sigma), with
Pi(0, T_sigma) = {pi_sigma[20]:.3f} meV (nu = 20) and {pi_sigma[1]:.3f} meV (RBM). In the archived tables the
w-independent column is the static part, although the header calls it Im Sigma; column 8 is Im Pi.

**Self-energies.** The evaluator was run with a single process (k -> k+Q covers the zone; adding k -> k-Q would
double count), all four energy-ordered band pairs (the Dirac crossing at E_F makes 2k_F backscattering the 1-2
pair), eta = k_B T/{eta_div:g} as in the archived runs, and nk = {nk_arch}. At 50 K it reproduces the archived table to
{dev6:.3f} meV; with nk = 1e5 the deviation is {dev5:.2f} meV (the Lorentzian of width k_B T/100 is not resolved),
so a high-temperature validation does not transfer to low T.

**Mode 20.** Nonadiabatic Omega = {omega_na['100']:.2f} / {omega_na['300']:.2f} meV at 100 / 300 K (gamma =
{gamma['100']:.2f} / {gamma['300']:.2f} meV), against the adiabatic {omega_ad:.2f} meV at 100 K and omega_Q =
{wq[20]:.2f} meV. The static self-energy diverges logarithmically (1D nesting), but for omega_Q >> k_B T the
dynamic self-energy at w ~ Omega saturates: the nonadiabatic softening is a few meV and non-monotonic in T
(minimum near omega_Q/4.6 k_B), while the adiabatic softening is ~24 meV at 100 K.

**RBM.** Omega^2 -> omega_b^2 + gamma_0^2 + 2 omega_Q Pi(0, T) as Omega -> 0, so T* = {tstar:.1f} K. Because the dynamic
self-energy reduces to the static one at Omega -> 0, nonadiabaticity only adds the linewidth term and shifts T* by
~10 K from the static estimate. T* is a small difference of large numbers: anchoring on the Fermi-Dirac DFPT
frequency at 10 mRy ({om1[0]:.1f} meV instead of the cold-smearing omega_Q = {wq[1]:.2f} meV) gives ~330 K, and the
DFPT smearing-as-temperature proxy gives {proxy:.0f} K. The scheme mismatch (cold-smearing omega_Q, Fermi-Dirac
subtraction) is an uncertainty larger than the nonadiabatic shift.

**Threshold model.** The archived program's parameter is a = 1/m*, so m* = 2 means a = 0.5: g0 = {g0['1']:.1f} meV
(RBM) and {g0['20']:.1f} meV (nu = 20). At the Fermi points |g| = {gf[1][0]:.0f}/{gf[1][1]:.0f} meV for the RBM (well
above g0: unstable) and {gf[20][0]:.0f}/{gf[20][1]:.0f} meV for nu = 20 (borderline; the full calculation keeps it stable).

**Not supported.** The adiabatic instability of nu = 20 lies far below 50 K and cannot be located from these data;
T* carries the reference-frequency uncertainty above; a Peierls distortion near room temperature for the RBM is
suggested, but its detection would have to contend with large linewidths (2 gamma ~ {2 * gamma['100']:.0f} meV for nu = 20).
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
