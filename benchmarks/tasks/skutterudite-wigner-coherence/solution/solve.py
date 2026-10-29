"""Reference solution for skutterudite-wigner-coherence; runs as the Harbor oracle and in check_task."""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
RUNS = {r['run']: r for r in CATALOG['runs']}
THZ_TO_EV = 4.135667696923859e-3     # h * 1 THz
HBAR_EV_PS = 6.582119569e-4
EV = 1.602176634e-19
EVIDENCE = []


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(1)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(r.get('error'))
    EVIDENCE.append(r['id'])
    return r['result']


def volume(run):
    return abs(np.linalg.det(np.array(RUNS[run]['primitive_lattice_A'])))


def kappa(m, G, vol, gamma_iso=None):
    """SMA Wigner conductivity: populations (degenerate blocks) and coherences (non-degenerate pairs)."""
    f = np.array(m['frequency_THz']); w = np.array(m['weights'], float)
    g = np.array(m['gamma_THz']) + (0 if gamma_iso is None else gamma_iso)
    cv = np.array(m['heat_capacity_eV_per_K'])
    ok = f > 0.01                                     # phono3py cutoff: acoustic modes at Gamma excluded
    with np.errstate(divide='ignore', invalid='ignore'):
        tau = np.where(ok, 1.0 / (4 * np.pi * g), 0.0)   # ps; gamma = Im(self-energy) -> Gamma = 2 * 2 pi gamma
    conv = EV * 1e4 * 1e-12 / (vol * 1e-30)
    deg = np.abs(f[:, :, None] - f[:, None, :]) < 1e-4
    pair = ok[:, :, None] & ok[:, None, :]
    kp = ((np.where(deg & pair, G, 0).sum(2) * cv * tau).sum(1) * w).sum() * conv / w.sum()
    hw, hg = f * THZ_TO_EV, 2 * g * THZ_TO_EV        # full widths
    with np.errstate(divide='ignore', invalid='ignore'):
        s = hg[:, :, None] + hg[:, None, :]
        lor = 0.5 * s / ((hw[:, :, None] - hw[:, None, :]) ** 2 + 0.25 * s ** 2)
        x = cv / hw
        pref = 0.25 * (hw[:, :, None] + hw[:, None, :]) * (x[:, :, None] + x[:, None, :])
        term = np.where(pair & ~deg, G * pref * lor, 0.0)
    kc = (term.sum((1, 2)) * w).sum() * conv * HBAR_EV_PS / w.sum()
    return float(kp), float(kc)


def isotope(m):
    """Archived isotope rates, with numerically failed entries removed (rates far above the mode frequency)."""
    gi = m['gamma_isotope_THz']
    if gi is None:
        return None, []
    gi = np.array(gi)
    f = np.array(m['frequency_THz'])
    bad = (gi > 1.0) | ((gi > f) & (f > 0.01))
    return np.where(bad, 0.0, gi), [(int(q), round(float(f[q, s]), 3), float(gi[q, s])) for q, s in zip(*np.where(bad))]


def main():
    notes = {}
    # 1. Conventions: validate the evaluator on a released filled run against its phono3py log.
    run = 'YbCo4Sb12/hiphive'
    log = buy({'kind': 'archived_file', 'run': run, 'file': 'tc_YbCo4Sb12.out'})['thermal_conductivity_xx_W_mK']
    m = buy({'kind': 'mode_data', 'run': run, 'temperature_K': 300})
    G = np.array(buy({'kind': 'velocity_operator', 'run': run})['G_A2THz2'])
    raw = kappa(m, G, volume(run), np.array(m['gamma_isotope_THz']))
    gi, bad = isotope(m)
    fixed = kappa(m, G, volume(run), gi)
    notes['YbCo4Sb12'] = {'log_300K': (log['K_P']['300'], log['K_C']['300']), 'archived_isotope': raw,
                          'corrupted_isotope_removed': fixed, 'corrupted_entries': bad[:6], 'n_corrupted': len(bad)}
    # 2. Conflicting CoSb3 Wigner tables in the same folder: the phono3py log versus tc_data.dat.
    run = 'CoSb3/v_operator_mean_free_path_isotope'
    log = buy({'kind': 'archived_file', 'run': run, 'file': 'nohup.out'})['thermal_conductivity_xx_W_mK']
    tc = buy({'kind': 'archived_file', 'run': run, 'file': 'tc_data.dat'})['text']
    row = [l.split() for l in tc.splitlines() if l.strip().startswith('300')][0]
    notes['CoSb3'] = {'log_KP_300': log['K_P']['300'], 'tc_data_kP_300': float(row[1])}
    # 3. Isotope scattering matters for unfilled compounds at low T (FeSb3 log vs ph-ph only).
    run = 'FeSb3/hiphive'
    log = buy({'kind': 'archived_file', 'run': run, 'file': 'tc_FeSb3.out'})['thermal_conductivity_xx_W_mK']
    m = buy({'kind': 'mode_data', 'run': run, 'temperature_K': 100})
    Gd = np.einsum('qs,st->qst', np.array(m['v2_A2THz2']), np.eye(len(m['v2_A2THz2'][0])))
    gi, bad = isotope(m)
    notes['FeSb3'] = {'log_KP_100': log['K_P']['100'], 'with_isotope': kappa(m, Gd, volume(run), gi)[0],
                      'phph_only': kappa(m, Gd, volume(run))[0], 'n_corrupted': len(bad)}
    # 4. YbIr4Sb12: archived inputs, corrupted isotope rates removed.
    run = 'YbIr4Sb12/hiphive'
    G = np.array(buy({'kind': 'velocity_operator', 'run': run})['G_A2THz2'])
    res, corrupted = {}, None
    for T in (100, 300, 600, 650, 700):
        m = buy({'kind': 'mode_data', 'run': run, 'temperature_K': T})
        gi, bad = isotope(m)
        corrupted = bad
        res[T] = kappa(m, G, volume(run), gi)
        res[T] = res[T] + (kappa(m, G, volume(run), np.array(m['gamma_isotope_THz']))[0],)
    m_ybir = m
    temps = sorted(res)
    ratio = {T: res[T][1] / res[T][0] for T in temps}
    i = next(k for k, T in enumerate(temps) if ratio[T] >= 0.5)
    T0, T1 = temps[i - 1], temps[i]
    assert T1 - T0 == 50, 'crossover bracket must be adjacent grid temperatures'
    t_half = T0 + (0.5 - ratio[T0]) * (T1 - T0) / (ratio[T1] - ratio[T0])
    # 5. IrSb3: no corrupted rates; isotope scattering kept.
    run = 'IrSb3/hiphive'
    G = np.array(buy({'kind': 'velocity_operator', 'run': run})['G_A2THz2'])
    m50 = buy({'kind': 'mode_data', 'run': run, 'temperature_K': 50})
    m300 = buy({'kind': 'mode_data', 'run': run, 'temperature_K': 300})
    gi50, bad50 = isotope(m50)
    gi300, _ = isotope(m300)
    ir50 = kappa(m50, G, volume(run), gi50)
    ir300 = kappa(m300, G, volume(run), gi300)
    ir50_noiso = kappa(m50, G, volume(run))
    harm = {c: buy({'kind': 'harmonic', 'compound': c}) for c in ('YbIr4Sb12', 'IrSb3')}
    msd = {c: harm[c]['msd_A2']['300'] for c in harm}
    # Rattler frequency: the PDOS run used phonopy's QE calculator on eV/A^2 force constants, so its axis is
    # inflated by the ratio of the unit-conversion factors (check: the DOS extends far beyond the mode frequencies).
    h = harm['YbIr4Sb12']
    axis, yb = np.array(h['pdos_frequency_axis']), np.array(h['pdos_by_element']['Yb'])
    total = sum(np.array(v) for v in h['pdos_by_element'].values())
    top_axis = axis[total > 1e-6 * total.max()].max()
    top_modes = float(np.max(m_ybir['frequency_THz']))
    factor = h['pdos_run']['phonopy_yaml']['frequency_unit_conversion_factor'] / RUNS['YbIr4Sb12/hiphive']['frequency_unit_conversion_factor']
    rattler = float(axis[np.argmax(yb)] / factor)
    kB, amu = 1.380649e-23, 1.66053907e-27
    m_yb = [a['mass_amu'] for a in RUNS['YbIr4Sb12/hiphive']['atoms'] if a['symbol'] == 'Yb'][0]
    einstein = np.sqrt(3 * kB * 300 / (m_yb * amu * msd['YbIr4Sb12']['Yb'] * 1e-20)) / (2 * np.pi) / 1e12
    submission = {'YbIr4Sb12': {'kappa_P_W_mK': {'300': res[300][0], '700': res[700][0]},
                                'kappa_C_W_mK': {'100': res[100][1], '300': res[300][1], '700': res[700][1]},
                                'T_half_K': float(t_half), 'rattler_frequency_THz': rattler},
                  'IrSb3': {'kappa_P_50K_W_mK': ir50[0], 'kappa_TOT_300K_W_mK': ir300[0] + ir300[1]},
                  'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    y = notes['YbCo4Sb12']
    (APP/'report.md').write_text(f"""# Boltzmann-to-Wigner crossover in Yb-filled skutterudites (reference solution)

**Framework.** kappa_P is the single-mode-relaxation-time populations conductivity (Peierls-Boltzmann with
degenerate subspaces treated as blocks); kappa_C is the Wigner coherences conductivity between non-degenerate
bands, (w_s+w_s')/4 (C_s/w_s + C_s'/w_s') |V_ss'|^2 times a Lorentzian of full width (Gamma_s+Gamma_s')/2, with
Gamma = 1/tau = 4 pi gamma (phono3py gamma is the imaginary part of the self-energy). Coherences matter when
linewidths become comparable to interband spacings (tau ~ 1/Delta w_avg).

**Validation.** With the archived inputs my evaluator reproduces the YbCo4Sb12 phono3py log at 300 K
(log K_P, K_C = {y['log_300K']}; mine {y['archived_isotope'][0]:.4f}, {y['archived_isotope'][1]:.4f}). The CoSb3 folder holds two
inconsistent Wigner tables: the phono3py log (K_P = {notes['CoSb3']['log_KP_300']} at 300 K) and tc_data.dat
(kP = {notes['CoSb3']['tc_data_kP_300']:.2f}), written by a post-processing script whose populations prefactor is 0.5
instead of 0.25 (and which uses omega_1 for omega_2); I trust the phono3py log. For FeSb3 at 100 K, ph-ph linewidths
alone give {notes['FeSb3']['phph_only']:.3f} W/mK, adding the archived isotope rates gives {notes['FeSb3']['with_isotope']:.3f}
(log {notes['FeSb3']['log_KP_100']}): the archived gamma excludes isotope scattering, which must be added.

**Corrupted isotope rates.** In the Yb-filled runs the archived gamma_isotope has {y['n_corrupted']} (YbCo4Sb12) and
{len(corrupted)} (YbIr4Sb12) entries of 10^3-10^7 THz, all at the q-point nearest Gamma (examples q, f, gamma_iso:
{corrupted[:3]}), including the acoustic branches, i.e. linewidths 10^3-10^7 times the mode frequency. Every other
isotope rate is below 0.01 THz (mostly 1e-5 to 1e-3 THz) and mass-variance scattering of long-wavelength acoustic
phonons is weaker still, so this is a numerical failure of the isotope integration, not physics. Removing these entries raises kappa_P of YbIr4Sb12 at 300 K from {res[300][2]:.3f} to {res[300][0]:.3f} W/mK
(YbCo4Sb12: {y['archived_isotope'][0]:.3f} to {y['corrupted_isotope_removed'][0]:.3f}); the archived YbCo4Sb12 table
reproduces only with the corrupted rates, so it is affected too. Discarding the q-point instead would bias the BZ sum.

**Results (YbIr4Sb12).** kappa_P = {res[300][0]:.3f} (300 K), {res[700][0]:.3f} (700 K) W/mK; kappa_C = {res[100][1]:.4f}
(100 K), {res[300][1]:.4f} (300 K), {res[700][1]:.4f} (700 K). kappa_C/kappa_P = {', '.join('%d K: %.3f' % (T, ratio[T]) for T in temps)};
it reaches 0.5 at T_half = {t_half:.0f} K (linear interpolation between {T0} and {T1} K). kappa_P falls roughly as 1/T while
kappa_C rises and saturates, so wave-like tunnelling becomes comparable to particle-like propagation at the
temperatures where skutterudite thermoelectrics operate.

**IrSb3 (unfilled).** kappa_P(50 K) = {ir50[0]:.1f} W/mK with isotope scattering ({ir50_noiso[0]:.1f} without);
kappa_TOT(300 K) = {ir300[0] + ir300[1]:.2f} W/mK with kappa_C = {ir300[1]:.3f}, i.e. Boltzmann regime (no corrupted rates).

**Rattling.** Yb MSD at 300 K is {msd['YbIr4Sb12']['Yb']:.3f} A^2 against {msd['YbIr4Sb12']['Sb']:.3f} (Sb) and
{msd['YbIr4Sb12']['Ir']:.3f} (Ir). The archived projected DOS extends to {top_axis:.1f} on its own axis while the largest
mode frequency is {top_modes:.2f} THz: the PDOS run used phonopy's QE calculator (conversion factor
{h['pdos_run']['phonopy_yaml']['frequency_unit_conversion_factor']}) on force constants in eV/A^2 (phono3py factor
{RUNS['YbIr4Sb12/hiphive']['frequency_unit_conversion_factor']}), so its axis is inflated by {factor:.4f}. Corrected, the
Yb-projected DOS peaks at {rattler:.3f} THz, consistent with the Einstein frequency from the Yb MSD ({einstein:.3f} THz)
and the lowest optical Gamma mode; the Yb band (0.65-1.25 THz) lies among the acoustic branches and hybridises with them. Filling shortens lifetimes across the spectrum, pushing them toward the Wigner limit; it is
not only a resonant scattering channel, since kappa_C/kappa_P grows by an order of magnitude relative to IrSb3.

**Competing explanations and limits.** A glass-like (Cahill/Allen-Feldman) picture would make kappa_C dominant at all
T; here kappa_P still dominates below ~650 K. Magnetism: only harmonic AFM (DFT+U) data are archived (Gamma
frequencies similar to NM), so its effect on kappa cannot be assessed here. The fc3 source matters for YbFe4Sb12
(finite differences vs hiPhive change kappa_P by a factor ~2.4 but kappa_C by ~10%). Results are SMA on a 9x9x9 mesh
with three-phonon and isotope scattering only (no four-phonon, no full BTE solution, no mesh extrapolation).
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
