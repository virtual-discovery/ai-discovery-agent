"""Reference solution for hgte-gated-subband-densities; runs as the Harbor oracle and in check_task.

No literature value is used: the gate calibration comes from the purchased 70 nm Hall scan, subband
assignments from the purchased kdotpy runs, and densities from an SdH Fourier analysis of the purchased
raw transport blocks.
"""
import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import periodogram

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
E, H = 1.602176634e-19, 6.62607015e-34
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
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(r)
    EVIDENCE.append(r['id'])
    return r['result']


# ------------------------------------------------------------------ SdH analysis
def spectrum(B, y, bmax=1.0, nmax=8e11, nres=0.025e11, pad=4):
    dbi = 1 / (nmax * 1e4) * (E / H) / 2
    bmin = 1 / (1 / bmax + (nmax / nres) * dbi)
    ib = 1 / np.arange(1 / bmax, 1 / bmin + dbi, dbi)
    x = interp1d(B, y, bounds_error=False, fill_value=0)(ib)
    f, p = periodogram(x, fs=1 / dbi, window='hamming', detrend='linear', scaling='density', nfft=len(x) * pad)
    return E / H * f / 1e4 / 1e11, p          # n = e f / h per spin-resolved pocket, 1e11 cm^-2


class Block:
    def __init__(self, res, rpre):
        self.vg = np.array(res['VG_V'])
        B = np.array(res['B_T'])
        vi, vxx, vxy = (np.array(res[k]) for k in ('VI_V', 'Vxx_V', 'Vxy_V'))
        rho, rxy = vxx / vi * rpre / 3.0, vxy / vi * rpre          # L/W = 3
        sxx = rho / (rho ** 2 + rxy ** 2)
        spec = [spectrum(B, s) for s in sxx]
        self.n = spec[0][0]
        self.logp = np.log10(np.array([p for _, p in spec]))

    def averaged(self, center, hw=0.02):
        return self.logp[np.abs(self.vg - center) <= hw + 1e-9].mean(0)


def maxima(n, L, lo=0.15, hi=8.0):
    out = []
    for i in range(1, len(n) - 1):
        if lo <= n[i] <= hi and L[i] > L[i - 1] and L[i] >= L[i + 1]:
            a, b, c = L[i - 1], L[i], L[i + 1]
            d = 0.5 * (a - c) / (a - 2 * b + c) if a - 2 * b + c else 0.0
            out.append((n[i] + d * (n[1] - n[0]), b))
    return out


def contrast(L, n, half=0.25):
    """Log power relative to a running mean over +-half (1e11 cm^-2): peak prominence."""
    k = int(round(half / (n[1] - n[0])))
    padded = np.pad(L, (k, k), mode='edge')
    return L - np.convolve(padded, np.ones(2 * k + 1) / (2 * k + 1), mode='valid')


def strongest_near(n, L, guess, frac=0.25, width=None):
    w = frac * guess if width is None else width
    c = [p for p in maxima(n, L) if abs(p[0] - guess) <= w]
    return max(c, key=lambda p: p[1])[0]


# ------------------------------------------------------------------ model helpers
def band_at_fermi(calc, band, fermi=None):
    d = calc['integrated_dos_per_nm2']
    ef = fermi if fermi is not None else kp(calc)['ef']
    i = int(np.argmin(np.abs(np.array(d['E_meV']) - ef)))
    return d['IDOS(%d)' % band][i] * 1e3          # 1e11 cm^-2


def kp(calc):
    out = {}
    for line in calc['kpout'].splitlines():
        if '[' in line and ']:' in line:
            tag = line.split('[')[1].split(']')[0]
            out[tag] = float(line.split(']:')[1].split()[0])
    return out


def interp(x, pts):
    pts = sorted(pts)
    return float(np.interp(x, [p[0] for p in pts], [p[1] for p in pts]))


def main():
    rpre = CATALOG['samples']['70nm']['pre_resistor_ohm']
    # 1. Gate calibration of the 70 nm bar: single-carrier low-field Hall density on the electron side.
    hs = np.array(buy({'kind': 'hall_scan', 'sample': '70nm'})['rows'])
    nH = 1 / (hs[:, 1] * rpre * E) * 1e-15
    m = (hs[:, 0] > 0.05 + 1e-9) & (hs[:, 0] <= 0.29 + 1e-9)
    a, b = np.polyfit(hs[m, 0], nH[m], 1)
    ntot = lambda v: a * v + b                                   # 1e11 cm^-2, electrons positive

    # 2. Released self-consistent runs (main) and hole-pocket wedge runs for 70 nm.
    main_runs = {n: buy({'kind': 'calculation', 'sample': '70nm', 'n_nm2': n})
                 for n in (0.003, 0.004, 0.005, 0.006, 0.007, -0.001, -0.002, -0.003, -0.004)}
    dens = {n: {bd: band_at_fermi(c, bd) for bd in (-1, 1, 2, 3)} for n, c in main_runs.items()}
    model = lambda bd, v: interp(ntot(v), [(n * 1e3, dens[n][bd]) for n in dens])
    # Surface identification: the gate-side (top) surface is the band that follows the gate.
    def steeper(n1, n2):
        s1, s2 = (abs(dens[n1][bd] - dens[n2][bd]) for bd in (1, 2))
        return (1, 2) if s1 > s2 else (2, 1)
    top_e, bot_e = steeper(0.005, 0.006)
    top_h, bot_h = steeper(-0.001, -0.002)
    wedges = {}
    for n, w in ((-0.002, [165, 285]), (-0.003, [165, 285]), (-0.003, [105, 165])):
        c = buy({'kind': 'pocket_calculation', 'sample': '70nm', 'n_nm2': n, 'kphi_deg': w})
        width = w[1] - w[0]
        # kdotpy extrapolates a wedge to the full circle; the wedge's own E_F is not the system's E_F.
        wedges[(n, w[0])] = -band_at_fermi(c, -1, fermi=kp(main_runs[n])['ef']) / (360 / width)
    full = -dens[-0.003][-1]
    pair = 2 * (wedges[(-0.003, 165)] + wedges[(-0.003, 105)])
    model_hole = lambda v: interp(ntot(v), [(n * 1e3, wedges[(n, 165)]) for n in (-0.002, -0.003)])

    # 3. 70 nm SdH spectra at the graded gate voltages and along the hole ridge.
    blocks = {v: Block(buy({'kind': 'transport', 'sample': '70nm', 'vg_center': v}), rpre)
              for v in (0.30, 0.20, 0.05, -0.45, -0.55, -0.60, -0.65, -0.70, -0.75)}
    L = lambda v: blocks[v].averaged(v)
    n70 = blocks[0.30].n
    exp = {'top_surface_at_+0.20V': strongest_near(n70, L(0.20), model(top_e, 0.20)),
           'bottom_surface_at_+0.05V': strongest_near(n70, L(0.05), model(bot_e, 0.05)),
           'first_bulk_subband_at_+0.30V': strongest_near(n70, L(0.30), model(3, 0.30)),
           'top_surface_at_-0.45V': strongest_near(n70, L(-0.45), model(top_h, -0.45)),
           'bottom_surface_at_-0.45V': strongest_near(n70, L(-0.45), model(bot_h, -0.45))}
    # Hole ridge: the pocket density grows linearly as the gate depletes the film; harmonics and sum
    # frequencies of the surface states are stronger at single gate voltages, so track the ridge
    # (line search on peak prominence, since the absolute power falls steeply with frequency).
    hole_blocks = [-0.55, -0.60, -0.65, -0.70, -0.75]
    vg_all = np.concatenate([blocks[v].vg for v in hole_blocks])
    lp_all = np.array([contrast(l, blocks[-0.70].n) for v in hole_blocks for l in blocks[v].logp])
    best = None
    for sl in np.arange(2.0, 10.01, 0.05):
        for v0 in np.arange(-0.45, -0.0999, 0.005):
            nn = sl * (v0 - vg_all)
            ok = (nn > 0.9) & (nn < 4.0)
            if ok.sum() < 20:
                continue
            val = np.mean([np.interp(x, n70, lp_all[i]) for x, i in zip(nn[ok], np.where(ok)[0])])
            if best is None or val > best[0]:
                best = (val, sl, v0)
    _, sl, v0 = best
    ridge = [(v, strongest_near(n70, L(v), sl * (v0 - v), width=0.2)) for v in hole_blocks]
    exp['hole_pocket_at_-0.70V'] = dict(ridge)[-0.70]
    s2, i2 = np.polyfit([r[0] for r in ridge], [r[1] for r in ridge], 1)
    onset = -i2 / s2
    # Multiplicity from charge balance: holes = surface electrons - net density (gate line).
    nt = strongest_near(n70, L(-0.70), model(top_h, -0.70))
    nb = strongest_near(n70, L(-0.70), model(bot_h, -0.70))
    p_total = nt + nb - ntot(-0.70)
    multiplicity = int(round(p_total / exp['hole_pocket_at_-0.70V']))

    submission = {
        'experiment_1e11_cm2': exp,
        'model_1e11_cm2': {'bottom_surface_at_+0.20V': model(bot_e, 0.20),
                           'top_surface_at_-0.45V': model(top_h, -0.45),
                           'hole_pocket_at_-0.60V': model_hole(-0.60)},
        'hole_onset_gate_V': float(onset),
        'majority_pocket_multiplicity': multiplicity,
        'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# 70 nm HgTe film: gate-resolved subband densities (reference solution)

**Gate calibration.** Single-carrier low-field Hall density (|B| <= 0.2 T, L/W = 3, R_pre = {rpre} Ohm) fitted on the
electron side (0.05 < Vg <= 0.29 V): n_tot = {a:.3f} Vg + {b:.3f} (1e11 cm^-2). Its zero, {-b / a:.3f} V, is the charge
neutrality point of the net density; it is NOT where holes appear.

**SdH analysis.** sigma_xx = rho_xx/(rho_xx^2 + rho_xy^2), interpolated on a uniform 1/B grid for B <= 1 T, linear
detrend, Hamming window, periodogram; each peak is one spin-resolved Fermi pocket, n = e f / h (gated films break
inversion symmetry, so surface states, bulk subbands and hole pockets are all non-degenerate). Spectra are averaged
over +-0.02 V. Peaks were assigned with the released self-consistent kdotpy runs, mapped to Vg with the gate line.

**Surfaces.** The top (gate-side) surface is the band that follows the gate: band 1 on the electron side, band 2 on the
hole side (kdotpy orders bands by energy, and the two surface Dirac cones exchange order near n ~ 0.001 nm^-2). The
bottom surface is screened by the top surface and the bulk: its SdH density stays near 0.6-0.85e11 cm^-2 across the
whole gate range, while the top surface goes from {exp['top_surface_at_-0.45V']:.2f} (-0.45 V) to {exp['top_surface_at_+0.20V']:.2f}
(+0.20 V). A rigid-band filling of one run would put about half of the gate-induced charge on the bottom surface.
At +0.30 V the first bulk subband is at {exp['first_bulk_subband_at_+0.30V']:.2f}; harmonics of the strong bottom-surface
peak lie near 1.6.

**Hole pockets.** The released wedge runs (kphi 105-165 and 165-285 deg) show two inequivalent side-maximum pockets
along [1-10]-type and [110]-type directions (bulk inversion asymmetry lifts the four-fold degeneracy): main-run band -1
at n = -0.003 nm^-2 holds {full:.2f}, and 2 x (majority + minority) = {pair:.2f}. Wedge IDOS must be read at the main run's
E_F (the wedge run's own E_F assumes an isotropic extrapolation) and divided by 360/width. The SdH ridge of the majority
pocket was tracked across -0.75..-0.55 V (slope {sl:.2f} per V); at -0.70 V it is at {exp['hole_pocket_at_-0.70V']:.2f}, while
single spectra are dominated by harmonics and sum frequencies of the surface states. Charge balance gives
p_total = {p_total:.2f}, i.e. {p_total / exp['hole_pocket_at_-0.70V']:.2f} majority pockets -> {multiplicity}. The ridge extrapolates to
zero at {onset:.3f} V, well on the electron side of neutrality: bulk holes coexist with surface electrons (semimetal).

**Model vs experiment.** The self-consistent full-band Hartree runs reproduce the top surface and the first bulk
subband within a few per cent, but put the screened bottom surface 0.1-0.2e11 cm^-2 too high and the majority hole
pocket too low (-0.60 V: model {model_hole(-0.60):.2f}, SdH {dict(ridge)[-0.60]:.2f}).

**Not supported by these data.** The interface/trap charge and the dielectric constant are not separately determined
(only one archived parameter set: E_P = 21 eV, VBO = -620 meV, eps = 20.8, sample-specific boundary potential); a
lever-arm change on the hole side cannot be excluded; the minority pocket is not resolved by SdH.
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
