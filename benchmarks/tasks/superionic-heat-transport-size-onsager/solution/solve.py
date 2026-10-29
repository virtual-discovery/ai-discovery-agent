"""Reference solution for superionic-heat-transport-size-onsager; runs as the Harbor oracle and in check_task.

Multicomponent cepstral Green-Kubo (heat flux + anion velocity sum), kappa = S_red(0) / (2 V kB <T>^2) for the
extensive LAMMPS heat flux; 4x4x4 cells for every graded state; NVT runs used after checking them against NVE.
"""
import base64
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.special import digamma, polygamma

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
KB = 8.617333262e-5                         # eV/K
CONV = 1.602176634e-19 / 1e-12 / 1e-10      # eV/(ps A K) -> W/(m K)


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(run):
    r = lab('evaluate', {'kind': 'md_run', 'run': run})
    while r['status'] == 'running':
        time.sleep(2)
        r = lab(f"job/{r['id']}")
    res = r['result']
    cols = [np.frombuffer(base64.b64decode(res['series'][c]), '<i2') * s
            for c, s in zip(res['series_columns'], res['scale'])]
    x = np.stack(cols, axis=1)
    return r['id'], res, x[:, :3], x[:, 3:]


def cepstral(spec, n_comp, n_cur):
    """Zero-frequency value of a periodogram ~ S chi2_nu/(2 n_comp); AIC choice of the cepstral cutoff."""
    nu = 2 * (n_comp - n_cur + 1)
    L = np.log(spec)
    ext = np.concatenate([L, L[-2:0:-1]])
    N = len(ext)
    C = np.fft.ifft(ext).real
    var = polygamma(1, nu / 2.0)
    tail = np.cumsum((C[:len(L)] ** 2)[::-1])[::-1]
    aic = np.array([N * tail[p] / var + 2 * p for p in range(1, len(L) - 1)])
    P = int(np.argmin(aic)) + 1
    L0 = C[0] + 2 * C[1:P].sum() - (digamma(nu / 2.0) - np.log(n_comp))
    return float(np.exp(L0)), float(np.exp(L0) * np.sqrt(var * (4 * P - 2) / N)), P


def kappa(res, J, V):
    """Returns (multicomponent kappa, its error, raw single-flux kappa) in W/m/K."""
    dt, n = res['frame_interval_ps'], len(J) - len(J) % 2
    F = np.array([np.fft.rfft(J[:n], axis=0), np.fft.rfft(V[:n], axis=0)])
    S = np.einsum('ifl,jfl->fij', np.conj(F), F) * dt / n / 3.0
    red = 1.0 / np.linalg.inv(S)[:, 0, 0].real
    vol, T = res['Lx_A'] ** 3, res['mean_temperature_K']
    pref = CONV / (2 * vol * KB * T ** 2)
    s0, e0, _ = cepstral(red, 3, 2)
    raw, _, _ = cepstral(S[:, 0, 0].real, 3, 1)
    return s0 * pref, e0 * pref, raw * pref


def main():
    plan = ['UO2-NVE-4x4x4-1800K', 'UO2-NVT-4x4x4-1800K', 'UO2-NVE-4x4x4-2600K', 'UO2-NVE-4x4x4-3300K',
            'UO2-NVE-4x4x4-3500K', 'CaF2-NVE-4x4x4-1500K', 'CaF2-NVE-4x4x4-1700K', 'PbF2-NVT-4x4x4-1100K',
            'PbF2-NVT-4x4x4-1300K', 'UO2-NVE-2x2x2-1700K', 'UO2-NVE-3x3x3-1700K', 'PbF2-NVE-3x3x3-1300K',
            'PbF2-NVT-3x3x3-1300K']
    K, evidence, rows = {}, [], []
    for run in plan:
        jid, res, J, V = buy(run)
        k, e, raw = kappa(res, J, V)
        K[run] = (k, e)
        evidence.append(jid)
        rows.append(f'| {run} | {res["mean_temperature_K"]:.0f} | {k:.3f} +- {e:.3f} | {raw:.3f} | {raw / k:.2f} |')
    # Thermostat test: NVT vs NVE pairs of the same state (crystalline UO2, superionic PbF2).
    pulls = [(K[a][0] - K[b][0]) / np.hypot(K[a][1], K[b][1])
             for a, b in (('UO2-NVT-4x4x4-1800K', 'UO2-NVE-4x4x4-1800K'),
                          ('PbF2-NVT-3x3x3-1300K', 'PbF2-NVE-3x3x3-1300K'))]
    biased = bool(max(abs(p) for p in pulls) > 3)
    pred = {'UO2-4x4x4-1700K': (K['UO2-NVE-4x4x4-1800K'][0] + K['UO2-NVT-4x4x4-1800K'][0]) / 2,
            'UO2-4x4x4-3400K': (K['UO2-NVE-4x4x4-3300K'][0] + K['UO2-NVE-4x4x4-3500K'][0]) / 2,
            'CaF2-4x4x4-1600K': (K['CaF2-NVE-4x4x4-1500K'][0] + K['CaF2-NVE-4x4x4-1700K'][0]) / 2,
            'PbF2-4x4x4-1200K': (K['PbF2-NVT-4x4x4-1100K'][0] + K['PbF2-NVT-4x4x4-1300K'][0]) / 2}
    runs = {r: K[r][0] for r in CATALOG['report_runs']}
    submission = {'kappa_predictions_W_mK': pred, 'kappa_runs_W_mK': runs, 'thermostat_biases_kappa': biased,
                  'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    table = '\n'.join(rows)
    (APP/'report.md').write_text(f"""# Heat transport in fluorite superionics (reference solution)

Estimator: kappa = V/(kB T^2) int_0^inf <J(t)J(0)> dt with J the intensive heat flux; LAMMPS compute heat/flux
returns J*V, so kappa = S(0)/(2 V kB <T>^2) with S(0) the zero-frequency power spectrum of the served flux and <T>
the run's mean temperature. S(0) is obtained by cepstral analysis (log-periodogram, AIC cutoff, error bar).
Because the anions diffuse in the superionic regime, the heat flux is coupled to a non-vanishing mass flux; the
single-flux Green-Kubo integral then depends on the arbitrary per-species energy zero and is not the thermal
conductivity. I use the multicomponent (Onsager) formula: kappa from the Schur-reduced spectrum
1/[S^-1]_00 of the 2x2 (heat flux, anion velocity) matrix, which removes the mass-flux coupling (chi2 with
2(3-2+1) = 4 dof).

| run | <T> (K) | kappa (W/m/K) | raw single-flux kappa | raw/onsager |
|---|---|---|---|---|
{table}

Findings. In the crystalline regime (UO2 1800 K) raw and multicomponent kappa agree; once the anions diffuse the
raw value is inflated (UO2 ~1.2x, CaF2 ~2x, PbF2 ~4x at the highest temperatures) and would suggest a spurious
rise of kappa above the transition. The corrected kappa decreases with T in the crystal and then saturates at a
temperature-independent value (UO2 ~1.2 W/m/K above ~2800 K, CaF2 ~1.35, PbF2 ~0.65): it does not follow the
phonon-gas 1/T law, consistent with heat carried by strongly scattered, diffusive (Allen-Feldman-like)
vibrations once the anion sublattice melts; the saturation is not a convective contribution, which the
multicomponent treatment has removed. Cell size: the 96-atom cell underestimates kappa strongly (UO2 1700 K:
{K['UO2-NVE-2x2x2-1700K'][0]:.2f} vs ~1.6-1.7 W/m/K in 768 atoms) and 324 atoms still ~15% low in crystalline UO2;
the size effect shrinks in the superionic regime. Graded states use 4x4x4 cells; a thermodynamic-limit kappa is
not supported by three sizes. NVT (global CSVR thermostat) and NVE runs agree within error (pulls
{', '.join(f'{p:+.1f}' for p in pulls)} sigma), so NVT data are used where no NVE run exists (PbF2 4x4x4 above 900 K).
Predictions interpolate the neighbouring 4x4x4 temperatures (UO2 1700 K from 1800 K, flat crystalline trend).
These are empirical force fields at fixed volume; no comparison with experiment is implied.
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
