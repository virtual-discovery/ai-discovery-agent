"""Reference solution for metal-resistivity-te-vs-anharmonicity; runs as the Harbor oracle and in check_task.

Route: archived EPW iterative-BTE resistivities for the purchasable runs, plus a Ziman/Bloch-Gruneisen phonon
surrogate validated on released same-lattice pairs,
    K(ifc, T) = sum_{q nu} |q|^2 [x / sinh x]^2 / omega^2,  x = omega / 2 k_B T   (q folded into the first BZ),
whose ratio between two force-constant sets at one lattice estimates the resistivity ratio.
"""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
KB = 0.08617333262          # meV / K
GRID = 16
EVIDENCE = []


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def evaluate(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(1)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(r)
    EVIDENCE.append(r['id'])
    return r['result']


RHO, SERTA, KERN, META = {}, {}, {}, {}


def transport(run):
    res = evaluate({'kind': 'transport', 'run': run})
    META[run] = res
    for row in res['results']:
        RHO[(run, row['T_K'])] = 1e6 / (np.trace(np.array(row['ibte_final_sigma_per_ohm_cm'])) / 3)
        SERTA[(run, row['T_K'])] = 1e6 / (np.trace(np.array(row['serta_sigma_per_ohm_cm'])) / 3)
    return res


def phon(ifc):
    res = evaluate({'kind': 'phonons', 'ifc': ifc, 'grid': GRID})
    q = np.array(res['grid']['q_cart_2pi_over_alat'])
    w = np.array(res['grid']['frequencies_meV'])
    KERN[ifc] = (np.sum(q ** 2, axis=1)[:, None] * np.ones_like(w), w, res['moments'])
    return res


def K(ifc, T):
    q2, w, _ = KERN[ifc]
    m = w > 0.05
    x = w[m] / (2 * KB * T)
    return float(np.sum(q2[m] * (x / np.sinh(x)) ** 2 / w[m] ** 2))


def interp(run, T, order=2):
    temps = np.array(sorted(t for (r, t) in RHO if r == run))
    idx = np.sort(np.argsort(np.abs(temps - T))[:order + 1])
    return float(np.polyval(np.polyfit(temps[idx], [RHO[(run, t)] for t in temps[idx]], order), T))


def main():
    dfpt = lambda m, T: f'{m}:DFPT:a({T}K)'
    tdep = lambda m, T, lat=None: f'{m}:TDEP{T}K:a({T if lat is None else lat}K)'
    PB_HARM, AL_HARM, NB_HARM = 'Pb/4-epw/1-DFPT/epw', 'Al/4-epw/1-DFPT/ibte.out-mesh50', 'Nb/4-epw/1-DFPT/ibte.out-mesh50'
    pb_te = lambda T: f'Pb/4-epw/2-with_TE/{T}K'
    pb_ta = lambda T: f'Pb/4-epw/4-TE_and_anharm/{T}K'
    # --- purchases: transport (2 units each) ---
    al_te = transport('Al/4-epw/2-with_TE/481K')
    mesh = al_te['fine_k_mesh'][0]
    harm_al = transport(f'Al/4-epw/1-DFPT/ibte.out-mesh{mesh}')       # same fine meshes as the TE run
    assert harm_al['fine_k_mesh'] == al_te['fine_k_mesh'] and harm_al['fine_q_mesh'] == al_te['fine_q_mesh']
    transport(PB_HARM)
    for T in (210, 310, 410, 510):
        transport(pb_te(T))
    for T in (210, 310):
        transport(pb_ta(T))
    for run in ('Nb/4-epw/2-with_TE/1000K', 'Nb/4-epw/2-with_TE/2000K', 'Nb/4-epw/3-TE_and_anharm/1000K', NB_HARM):
        transport(run)
    # --- purchases: phonons (1 unit each) ---
    for T in (0, 210, 310, 410, 510, 610):
        phon(dfpt('Pb', T))
    for T in (210, 310, 510):
        phon(tdep('Pb', T))
    for T in (410, 610):
        phon(tdep('Pb', T, 0))
    for T in (1000, 2000):
        phon(dfpt('Nb', T))
        phon(tdep('Nb', T))
    S = lambda a, b, T: K(a, T) / K(b, T)
    # 1. Al: thermal-expansion change at 481 K against the harmonic run with identical fine meshes.
    al_change = 100 * (RHO[('Al/4-epw/2-with_TE/481K', 481.0)] / interp(AL_HARM, 481) - 1)
    # Surrogate validation on released same-lattice Pb pairs (TE+anh vs TE); 10 K is left out (grid-limited,
    # Bloch-Gruneisen regime: the EPW 10 K values jump by a factor ~2 between lattices 0.1% apart).
    harm = {T: RHO[(PB_HARM, float(T))] for T in (10, 110, 210, 310, 410, 510, 610)}
    c_pb = {T: S(tdep('Pb', T), dfpt('Pb', T), T) / (RHO[(pb_ta(T), float(T))] / RHO[(pb_te(T), float(T))])
            for T in (210, 310)}
    c = float(np.mean(list(c_pb.values())))
    inv2 = {T: KERN[tdep('Pb', T)][2]['inverse_square_mean_meV-2'] / KERN[dfpt('Pb', T)][2]['inverse_square_mean_meV-2']
            / (RHO[(pb_ta(T), float(T))] / RHO[(pb_te(T), float(T))]) for T in (210, 310)}
    # 2. Pb thermal expansion at 610 K: phonon softening (surrogate across lattices) times an electronic /
    #    deformation-potential residual e(T) = R_EPW / S calibrated at 410 and 510 K and extrapolated linearly.
    R = {T: RHO[(pb_te(T), float(T))] / harm[T] for T in (210, 310, 410, 510)}
    S_te = {T: S(dfpt('Pb', T), dfpt('Pb', 0), T) for T in (310, 410, 510, 610)}
    e = {T: R[T] / S_te[T] for T in (310, 410, 510)}
    te610 = S_te[610] * (2 * e[510] - e[410]) * harm[610]
    te610_linear = 2 * RHO[(pb_te(510), 510.0)] - RHO[(pb_te(410), 410.0)]
    # 3. Pure anharmonic change at the 0 K lattice (TDEP at the 0 K lattice vs harmonic DFPT at the 0 K lattice).
    anh = {str(T): 100 * (S(tdep('Pb', T, 0), dfpt('Pb', 0), T) / c - 1) for T in (410, 610)}
    anh_expanded = {T: 100 * (RHO[(pb_ta(T), float(T))] / RHO[(pb_te(T), float(T))] - 1) for T in (210, 310)}
    anh_expanded[510] = 100 * (S(tdep('Pb', 510), dfpt('Pb', 510), 510) / c - 1)
    # 4. Pb TE + anharmonicity at 510 K: TE run times the same-lattice surrogate ratio.
    ta510 = RHO[(pb_te(510), 510.0)] * S(tdep('Pb', 510), dfpt('Pb', 510), 510) / c
    sep510 = RHO[(pb_te(510), 510.0)] * (1 + 0.5 * (anh['410'] + anh['610']) / 100)
    # 5. Nb TE + anharmonicity at 2000 K, calibrated on Nb's own 1000 K pair (d-band coupling differs from Pb/Al).
    te1, te2 = RHO[('Nb/4-epw/2-with_TE/1000K', 1000.0)], RHO[('Nb/4-epw/2-with_TE/2000K', 2000.0)]
    ta1 = RHO[('Nb/4-epw/3-TE_and_anharm/1000K', 1000.0)]
    c_nb = S(tdep('Nb', 1000), dfpt('Nb', 1000), 1000) / (ta1 / te1)
    nb2000 = te2 * S(tdep('Nb', 2000), dfpt('Nb', 2000), 2000) / c_nb
    nb_harm = {T: RHO[(NB_HARM, float(T))] for T in (1000, 2000)}
    submission = {'al_te_change_percent_481K': al_change, 'pb_te_resistivity_610K_uohm_cm': te610,
                  'pb_anharmonic_change_percent_0K_lattice': anh, 'pb_te_anharm_resistivity_510K_uohm_cm': ta510,
                  'nb_te_anharm_resistivity_2000K_uohm_cm': nb2000, 'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    ser = SERTA[(pb_ta(310), 310.0)] / RHO[(pb_ta(310), 310.0)]
    (APP/'report.md').write_text(f"""# Thermal expansion vs anharmonicity in the phonon-limited resistivity (reference solution)

**Conventions.** rho = 1/(trace sigma / 3) of the converged iterative BTE (last IBTE iteration), micro-Ohm cm.
SERTA overestimates rho by {100 * (ser - 1):.1f}% for Pb at 310 K (TE+anh) and is not used. Effects are ratios of runs
at the same temperature and identical numerical settings.

**Al, thermal expansion at 481 K.** TE runs use {mesh}^3 fine meshes; the harmonic reference is the {mesh}^3 run
(quadratic interpolation 400-600 K): +{al_change:.2f}%. Mixing in the 70^3 reference would add about +2 points of
pure mesh error, comparable to the effect.

**Surrogate.** For two force-constant sets at one lattice (same electrons, same DFPT deformation potentials) the
resistivity ratio is estimated by K = sum |q|^2 (x/sinh x)^2 / omega^2 (Ziman transport factor (1 - cos theta) ~ q^2,
|g|^2 ~ 1/omega, Bose occupation). On released Pb TE+anh/TE pairs it overestimates by c = {c_pb[210]:.4f} (210 K) and
{c_pb[310]:.4f} (310 K); the lambda-like <omega^-2> moment printed by the lab is off by {inv2[210]:.3f} and {inv2[310]:.3f}
and, being dominated by near-Gamma modes, transfers badly between lattices. The 10 K EPW values are not used:
rho(10 K) is 1.57 (0 K lattice) vs 2.95 micro-Ohm cm (10 K lattice, 0.1% larger), far above Bloch-Gruneisen
expectations, i.e. grid-limited.

**Pb, pure anharmonicity at the 0 K lattice** (TDEP at the 0 K lattice vs harmonic DFPT): {anh['410']:.1f}% at 410 K and
{anh['610']:.1f}% at 610 K. At the expanded lattice the correction is {anh_expanded[210]:.1f}% (210 K) and
{anh_expanded[310]:.1f}% (310 K) in the released EPW pairs and {anh_expanded[510]:.1f}% at 510 K (surrogate): TDEP
stiffening matters more where expansion has softened the DFPT phonons, so the anharmonic correction is not a
property of T alone and cannot be transferred between lattices.

**Pb, thermal expansion at 610 K.** R = rho_TE/rho_harm = {R[310]:.3f}, {R[410]:.3f}, {R[510]:.3f} at 310-510 K. Phonon
softening alone (surrogate DFPT(a_T) vs DFPT(a_0)) accounts for {S_te[410]:.3f}, {S_te[510]:.3f}, {S_te[610]:.3f}; the
residual e = {e[410]:.3f}, {e[510]:.3f} (electronic structure and deformation potentials at the expanded lattice) grows with
T and was extrapolated linearly: rho_TE(610 K) = {te610:.1f} micro-Ohm cm (harmonic {harm[610]:.2f}). Linear
extrapolation of rho_TE ({te610_linear:.1f}) and phonon-only scaling ({S_te[610] * harm[610]:.1f}) miss the
acceleration as Pb approaches melting.

**Pb, TE + anharmonicity at 510 K:** {ta510:.2f} micro-Ohm cm vs harmonic {harm[510]:.2f} and TE-only
{RHO[(pb_te(510), 510.0)]:.2f}. Combining separately estimated effects would give {sep510:.1f}: the effects do not
multiply or add, because the anharmonic correction depends on the lattice. In Pb the two nearly cancel.

**Nb, TE + anharmonicity at 2000 K:** {nb2000:.1f} micro-Ohm cm, calibrated on Nb's own 1000 K pair
(c = {c_nb:.3f}, vs ~{c:.3f} for Pb: d-band coupling is not free-electron-like). Unlike Pb and Al the effects do not
cancel: TE+anh/harmonic = {ta1 / nb_harm[1000]:.3f} at 1000 K, and the harmonic value {nb_harm[2000]:.1f} would be ~12% high.

**Uncertainty and limits.** Surrogate transfer errors are 1-3% for Pb/Al pairs and ~5% for Nb before calibration; the
TE extrapolation carries ~3-5%. These data cannot separate electronic from deformation-potential changes under
expansion, say nothing below ~100 K (grid-limited), and the Nb 1500 K TE+anh output is missing from the archive.
""")
    lab('submit', submission)
    for _ in range(120):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(2)


if __name__ == '__main__':
    main()
