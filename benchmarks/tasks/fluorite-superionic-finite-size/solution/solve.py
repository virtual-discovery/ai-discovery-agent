"""Reference solution for fluorite-superionic-finite-size; runs as the Harbor oracle and in check_task."""
import json
import math
import os
import re
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
KB = 8.617333262e-5
MASS = {'CaF2': (40.078, 18.998403), 'PbF2': (207.2, 18.998403), 'UO2': (238.028909, 15.9994), 'AgI': (126.90447, 107.8682)}
RATIO = {'CaF2': 2.0, 'PbF2': 2.0, 'UO2': 2.0, 'AgI': 1.0}      # N(mobile type 2) / N(type 1)
FU_ATOMS = {'CaF2': 3, 'PbF2': 3, 'UO2': 3, 'AgI': 2}


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


EVIDENCE, RUNS = [], {}


def buy(material, n, T):
    r = lab('evaluate', {'kind': 'md_run', 'material': material, 'supercell': n, 'temperature_K': T})
    while r['status'] == 'running':
        time.sleep(1)
        r = lab(f"job/{r['id']}")
    EVIDENCE.append(r['id'])
    RUNS[(material, n, T)] = r['result']
    return r['result']


def frame_ps(deck):
    """Time between saved frames: timestep x dump interval, read from each run's own deck."""
    dt = float(re.search(r'variable\s+dt\s+equal\s+([\d.]+)', deck).group(1))
    every = int(re.search(r'variable\s+traj_print_interval\s+equal\s+(\d+)', deck).group(1))
    return dt * every


def msd_lab(run):
    rows = np.array(run['msd_table']['rows'])
    t = rows[:, 0] * frame_ps(run['input_deck'])
    m1, m2 = MASS[run['material']]
    a = RATIO[run['material']] * m2 / m1              # mobile sublattice relative to the immobile one
    return t, rows[:, 2] + a * (2 + a) * rows[:, 5]


def diffusion(run, window=(20.0, 100.0)):
    t, msd = msd_lab(run)
    sel = (t >= window[0]) & (t <= window[1])        # past the ballistic rise and vibrational plateau
    return float(np.polyfit(t[sel], msd[sel], 1)[0] / 6)


def plateau(run):
    t, msd = msd_lab(run)
    return float(msd[(t >= 1) & (t <= 3)].mean())


def diffusive(run):
    """Accept a run for the Arrhenius fit only if the MSD gains at least twice its vibrational plateau in 100 ps."""
    D = diffusion(run)
    return D > 0 and 600 * D >= 2 * plateau(run)


def heat_capacity_peak(runs):
    runs = sorted(runs, key=lambda r: r['temperature_K'])
    T = np.array([r['thermo_summary']['mean_temperature_K'] for r in runs])
    nfu = runs[0]['atoms'] / FU_ATOMS[runs[0]['material']]
    E = np.array([r['thermo_summary']['mean_total_energy_eV'] / nfu for r in runs])
    cv = np.gradient(E, T) / KB
    k = int(np.argmax(cv))
    if 0 < k < len(T) - 1:
        c = np.polyfit(T[k - 1:k + 2], cv[k - 1:k + 2], 2)
        return float(np.clip(-c[1] / (2 * c[0]), T[k - 1], T[k + 1])), dict(zip(T.round(0).tolist(), cv.round(2).tolist()))
    return float(T[k]), dict(zip(T.round(0).tolist(), cv.round(2).tolist()))


def arrhenius(Ts, Ds):
    slope, intercept = np.polyfit(1 / (KB * np.array(Ts, float)), np.log(Ds), 1)
    return float(-slope), float(math.exp(intercept))


def interp(T, D):
    lo, hi = max(t for t in D if t < T), min(t for t in D if t > T)
    f = (1 / T - 1 / lo) / (1 / hi - 1 / lo)
    return float(math.exp(math.log(D[lo]) + f * (math.log(D[hi]) - math.log(D[lo]))))


PLAN = {('CaF2', 4): [1100, 1200, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1750],
        ('PbF2', 4): [500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1100, 1250],
        ('PbF2', 5): [600, 800],
        ('UO2', 4): [2400, 2500, 2650, 2700, 2750, 2800, 2850, 2900, 3100, 3200],
        ('AgI', 2): [600, 900], ('AgI', 5): [600, 900],
        ('PbF2', 2): [600], ('PbF2', 3): [600], ('CaF2', 2): [1300], ('UO2', 3): [2750]}


def main():
    for (m, n), Ts in PLAN.items():
        for T in Ts:
            buy(m, n, T)
    big = {m: [RUNS[(m, 4, T)] for T in PLAN[(m, 4)]] for m in ('CaF2', 'PbF2', 'UO2')}
    D4 = {m: {r['temperature_K']: diffusion(r) for r in big[m]} for m in big}
    tc, cv = {}, {}
    for m in big:
        tc[m], cv[m] = heat_capacity_peak(big[m])
    ea, windows = {}, {}
    for m in ('CaF2', 'PbF2'):
        lo = [r['temperature_K'] for r in big[m] if r['temperature_K'] <= tc[m] - 100 and diffusive(r)]
        hi = [r['temperature_K'] for r in big[m] if r['temperature_K'] >= tc[m] + 100]
        ea[m] = {'ordered': arrhenius(lo, [D4[m][T] for T in lo])[0], 'superionic': arrhenius(hi, [D4[m][T] for T in hi])[0]}
        windows[m] = (lo, hi)
    D5 = {T: diffusion(RUNS[('PbF2', 5, T)]) for T in (600, 800)}
    pred = {'CaF2-4x4x4-1250K': interp(1250, D4['CaF2']),
            'UO2-4x4x4-2450K': interp(2450, D4['UO2']), 'UO2-4x4x4-3150K': interp(3150, D4['UO2']),
            # the 768-atom PbF2 cell is converged at 800 K but not at 600 K; the 1500-atom cell is used for N = 4116
            'PbF2-7x7x7-600K': D5[600], 'PbF2-7x7x7-800K': D5[800],
            # alpha-AgI: no measurable size dependence between 32 and 500 atoms
            'AgI-6x6x6-600K': diffusion(RUNS[('AgI', 5, 600)]), 'AgI-6x6x6-900K': diffusion(RUNS[('AgI', 5, 900)])}
    size600 = {n: diffusion(RUNS[('PbF2', n, 600)]) for n in (2, 3, 4, 5)}
    agi = {(n, T): diffusion(RUNS[('AgI', n, T)]) for n in (2, 5) for T in (600, 900)}
    uo2_3 = diffusion(RUNS[('UO2', 3, 2750)])
    caf2_2 = diffusion(RUNS[('CaF2', 2, 1300)])
    submission = {'transition_temperature_K': tc, 'activation_energy_eV': ea,
                  'diffusion_predictions_A2_per_ps': pred, 'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    frames = {m: frame_ps(RUNS[(m, 4, PLAN[(m, 4)][0])]['input_deck']) for m in big}
    (APP/'report.md').write_text(f"""# Fluorite superionics and cell size (reference solution)

**Diffusion coefficients.** D is the slope/6 of a straight-line fit to the mobile-ion MSD over 20-100 ps
lag, after the ballistic rise and the vibrational plateau, relative to the immobile sublattice
(MSD + a(2+a) MSD_cm with a = N_mob m_mob / (N_imm m_imm)). The lag column counts saved frames; the frame
interval is timestep x dump interval from each deck: {frames} ps (UO2 uses a 2 fs step). In the
zero-momentum frame the cation sublattice drifts (its centre of mass moves opposite to the anions), which is
not cation diffusion. Runs whose MSD gains less than twice its plateau in 100 ps are not diffusive and are
not used in Arrhenius fits.

**T_c.** c_V = d<E>/dT per formula unit from central differences of the mean total energy of the 768-atom
runs, peak by a parabola through the maximum: CaF2 {tc['CaF2']:.0f} K, PbF2 {tc['PbF2']:.0f} K, UO2
{tc['UO2']:.0f} K. The peaks of CaF2 and PbF2 are broad (c_V within noise over ~100 K: {cv['PbF2']}),
so T_c is uncertain by about +-50 K; the UO2 peak is sharp. Small cells shift the peak up: UO2 3x3x3 is
still ~200 K higher than 4x4x4, so the 768-atom UO2 value is not converged either.

**Activation energies (768 atoms).** CaF2: ordered {ea['CaF2']['ordered']:.2f} eV (T = {windows['CaF2'][0]}),
superionic {ea['CaF2']['superionic']:.2f} eV ({windows['CaF2'][1]}); PbF2: ordered {ea['PbF2']['ordered']:.2f} eV
({windows['PbF2'][0]}), superionic {ea['PbF2']['superionic']:.2f} eV ({windows['PbF2'][1]}). The 96-atom cells
show no such change of slope (one intermediate activation energy), so the order-disorder signature needs
hundreds of atoms. UO2 ln D(1/T) is curved on both sides of T_c; its activation energies depend on the window
and are not reported.

**Size effects.** At 600 K the PbF2 D is {size600} for n = 2, 3, 4, 5: the 96-atom cell overestimates D in
the ordered regime, while near and above T_c small cells underestimate it (CaF2 1300 K: {caf2_2:.3g} for 96
atoms vs {D4['CaF2'][1300]:.3g} for 768). The sign change and the saturation of D for n >= 5 rule out a
hydrodynamic Yeh-Hummer 1/L correction (the solid matrix does not flow) and pure statistical noise (block
errors are a few per cent at these temperatures). The size effect follows the collective Frenkel disorder
of the order-disorder transition, which a small periodic cell suppresses and smears. The 768-atom PbF2
cell agrees with 1500 atoms at 800 K ({D4['PbF2'][800]:.3g} vs {D5[800]:.3g}) but not at 600 K
({D4['PbF2'][600]:.3g} vs {D5[600]:.3g}), so the 4116-atom predictions use the 1500-atom cell. UO2 3x3x3 at
2750 K gives {uo2_3:.3g} vs {D4['UO2'][2750]:.3g} for 768 atoms.

**alpha-AgI.** D is the same within noise from 32 to 500 atoms ({agi}), with a single small activation
energy and no c_V peak: this type-I conductor already has a molten Ag sublattice, so there is no collective
disordering transition for the cell to suppress, and no size correction is applied.

**Not supported.** Thermodynamic-limit T_c values (the peak still moves between 324 and 768 atoms,
especially for UO2), activation energies of UO2, D in the ordered regime of cells larger than 1500 atoms
beyond the released points, and any comparison with experiment for these empirical force fields.
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
