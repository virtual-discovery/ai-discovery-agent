"""Estimators for metal-resistivity-te-vs-anharmonicity, evaluated on pool.json (construction side).

Every graded quantity is an archived EPW iterative-BTE resistivity (or a ratio of two) that the lab withholds.
This module computes, from the released data only, the reference route the oracle follows, the defensible
alternatives (robustness) and the realistic wrong paths (ablations). The phonon surrogate is

    K_w(ifc, T) = sum_{q, nu} w(q, nu) [x / sinh x]^2 / omega_{q nu}^2,   x = hbar omega / 2 k_B T,

summed over a Gamma-centred n^3 grid (modes below 0.05 meV dropped). It is the Bloch-Gruneisen/Ziman
high-temperature kernel with |g|^2 ~ 1/omega at fixed deformation potentials; w = |q|^2 (q folded into the first
Brillouin zone) is the (1 - cos theta) transport factor of normal processes, w = 1 gives the lambda-like
<omega^-2> moment. For two force-constant sets at the same lattice (same electrons and deformation potentials)
rho_X / rho_Y ~ K_w(X) / K_w(Y).
"""
import json
from functools import lru_cache
from pathlib import Path
import numpy as np
import phonons as P

HERE = Path(__file__).resolve().parent
POOL = json.loads((HERE/'pool.json').read_text())
RUNS = {r['id']: r for r in POOL['runs']}
IFCS = {r['id']: r for r in POOL['ifcs']}
GRID = 16                      # grid the oracle requests from the lab


def rho(run_id, T, tensor='ibte'):
    """Resistivity in micro-Ohm cm at an archived temperature: 1e6 / (trace(sigma) / 3)."""
    for row in RUNS[run_id]['results']:
        if abs(row['T_K'] - T) < 1e-6:
            s = row['ibte_final_sigma_per_ohm_cm'] if tensor == 'ibte' else row['serta_sigma_per_ohm_cm']
            return 1e6 / (np.trace(np.array(s)) / 3)
    raise KeyError((run_id, T))


def rho_interp(run_id, T, tensor='ibte', order=2):
    """Polynomial interpolation of a multi-temperature run through the nearest archived temperatures."""
    temps = np.array([row['T_K'] for row in RUNS[run_id]['results']])
    vals = np.array([rho(run_id, t, tensor) for t in temps])
    idx = np.sort(np.argsort(np.abs(temps - T))[:order + 1])
    return float(np.polyval(np.polyfit(temps[idx], vals[idx], order), T))


@lru_cache(maxsize=None)
def grid(ifc_id, n=GRID, vectors=False):
    rec = IFCS[ifc_id]
    ph = P.from_record(rec)
    q = P.folded_grid(rec['at_alat'], n)
    w2, e = np.linalg.eigh(ph.dynmat(q))
    w = np.sign(w2) * np.sqrt(np.abs(w2)) * P.RY_MEV
    return q, w, e


def kernel(ifc_id, T, weight='q2', n=GRID):
    q, w, e = grid(ifc_id, n)
    q2 = np.sum(q ** 2, axis=1)[:, None] * np.ones_like(w)
    if weight == 'q2':
        wt = q2
    elif weight == 'uniform':
        wt = np.ones_like(w)
    elif weight == 'long_q2':
        qh = q / np.maximum(np.linalg.norm(q, axis=1, keepdims=True), 1e-12)
        wt = np.abs(np.einsum('qa,qan->qn', qh, e)) ** 2 * q2
    else:
        raise ValueError(weight)
    m = w > 0.05
    x = w[m] / (2 * P.KB_MEV * T)
    return float(np.sum(wt[m] * (x / np.sinh(x)) ** 2 / w[m] ** 2))


def lab_moment(ifc_id, key, n=GRID):
    return P.report(IFCS[ifc_id], n)['moments'][key]


def ratio(a, b, T, weight='q2'):
    if weight == 'inv2_moment':   # the lab's printed <omega^-2> (no occupation, no transport factor)
        return lab_moment(a, 'inverse_square_mean_meV-2') / lab_moment(b, 'inverse_square_mean_meV-2')
    if weight == 'omega_log':
        return (lab_moment(b, 'omega_log_meV') / lab_moment(a, 'omega_log_meV')) ** 2
    return kernel(a, T, weight) / kernel(b, T, weight)


# ----- Pb, Nb, Al run and force-constant names ------------------------------------------------------------------
PB_HARM = 'Pb/4-epw/1-DFPT/epw'
pb_te = lambda T: f'Pb/4-epw/2-with_TE/{T}K'
pb_ta = lambda T: f'Pb/4-epw/4-TE_and_anharm/{T}K'
pb_a0 = lambda T: f'Pb/4-epw/3-anharm/{T}K'
dfpt = lambda m, T: f'{m}:DFPT:a({T}K)'
tdep = lambda m, T, lat=None: f'{m}:TDEP{T}K:a({T if lat is None else lat}K)'


def pb_calibration(Ts=(210, 310), weight='q2', tensor='ibte'):
    """Surrogate over-estimate factor c = S / EPW on purchasable Pb TE+anh vs TE pairs (same lattice)."""
    cs = [ratio(tdep('Pb', T), dfpt('Pb', T), T, weight) / (rho(pb_ta(T), T, tensor) / rho(pb_te(T), T, tensor))
          for T in Ts]
    return float(np.mean(cs))


# ----- graded quantities -----------------------------------------------------------------------------------------
def al_te_change_481(ref='Al/4-epw/1-DFPT/ibte.out-mesh50', tensor='ibte', order=2):
    return 100 * (rho('Al/4-epw/2-with_TE/481K', 481, tensor) / rho_interp(ref, 481, tensor, order) - 1)


def pb_te_610(method='phonon_e_linear', weight='q2', tensor='ibte'):
    harm = {T: rho(PB_HARM, T, tensor) for T in (10, 110, 210, 310, 410, 510, 610)}
    te = {T: rho(pb_te(T), T, tensor) for T in (10, 110, 210, 310, 410, 510)}
    R = {T: te[T] / harm[T] for T in te}
    if method.startswith('phonon'):
        S = {T: ratio(dfpt('Pb', T), dfpt('Pb', 0), T, weight) for T in (310, 410, 510, 610)}
        e = {T: R[T] / S[T] for T in (310, 410, 510)}
        if method == 'phonon_e_linear':          # electronic/deformation-potential residual extrapolated linearly
            e610 = 2 * e[510] - e[410]
        elif method == 'phonon_e_const':
            e610 = e[510]
        elif method == 'phonon_e_quadratic':
            e610 = 3 * e[510] - 3 * e[410] + e[310]
        elif method == 'phonon_only':            # wrong: thermal expansion acts only through the phonons
            e610 = 1.0
        return S[610] * e610 * harm[610]
    if method == 'R_quadratic':
        return (3 * R[510] - 3 * R[410] + R[310]) * harm[610]
    if method == 'rho_cubic':
        Ts = np.array([210, 310, 410, 510])
        return float(np.polyval(np.polyfit(Ts, [te[t] for t in Ts], 3), 610))
    if method == 'rho_quadratic':
        Ts = np.array([310, 410, 510])
        return float(np.polyval(np.polyfit(Ts, [te[t] for t in Ts], 2), 610))
    if method == 'rho_linear':                   # wrong: linear extrapolation of the TE series
        return 2 * te[510] - te[410]
    raise ValueError(method)


def pb_anh0K(T, weight='q2', calibrate=True, calib_T=(210, 310), mode='surrogate', tensor='ibte'):
    """Pure anharmonic change (percent) at the 0 K lattice."""
    c = pb_calibration(calib_T, weight, tensor) if calibrate else 1.0
    if mode == 'surrogate':
        r = ratio(tdep('Pb', T, 0), dfpt('Pb', 0), T, weight) / c
    elif mode == 'T_lattice_transfer':            # wrong: the expanded-lattice correction used at the 0 K lattice
        r = ratio(tdep('Pb', T), dfpt('Pb', T), T, weight) / c
    else:
        raise ValueError(mode)
    return 100 * (r - 1)


def pb_te_anh_510(method='surrogate', weight='q2', calib_T=(210, 310), tensor='ibte'):
    T = 510
    if method == 'surrogate':
        c = pb_calibration(calib_T, weight, tensor)
        return rho(pb_te(T), T, tensor) * ratio(tdep('Pb', T), dfpt('Pb', T), T, weight) / c
    if method == 'R_harm_linear':                # TE+anh / harmonic ratio extrapolated from 210 and 310 K
        R = {t: rho(pb_ta(t), t, tensor) / rho(PB_HARM, t, tensor) for t in (210, 310)}
        return (R[310] + 2 * (R[310] - R[210])) * rho(PB_HARM, T, tensor)
    if method.startswith('linear_fit'):           # linear fit of the TE+anh series
        Ts = (10, 110, 210, 310) if method == 'linear_fit_with_10K' else (110, 210, 310)
        return float(np.polyval(np.polyfit(Ts, [rho(pb_ta(t), t, tensor) for t in Ts], 1), T))
    if method in ('separable_multiplicative', 'separable_additive'):
        # wrong: combine the separately estimated effects; the 0 K-lattice anharmonic change at 510 K is
        # interpolated between the surrogate estimates at 410 and 610 K
        a = 0.5 * (pb_anh0K(410, weight, tensor=tensor) + pb_anh0K(610, weight, tensor=tensor)) / 100
        harm, te = rho(PB_HARM, T, tensor), rho(pb_te(T), T, tensor)
        return te * (1 + a) if method == 'separable_multiplicative' else te + a * harm
    raise ValueError(method)


def nb_te_anh_2000(method='surrogate_nb', weight='q2', tensor='ibte'):
    te1, te2 = rho('Nb/4-epw/2-with_TE/1000K', 1000, tensor), rho('Nb/4-epw/2-with_TE/2000K', 2000, tensor)
    ta1 = rho('Nb/4-epw/3-TE_and_anharm/1000K', 1000, tensor)
    harm = lambda T: rho('Nb/4-epw/1-DFPT/ibte.out-mesh50', T, tensor)
    S2 = ratio(tdep('Nb', 2000), dfpt('Nb', 2000), 2000, weight)
    if method == 'surrogate_nb':                 # calibrated on Nb's own 1000 K pair
        c = ratio(tdep('Nb', 1000), dfpt('Nb', 1000), 1000, weight) / (ta1 / te1)
        return te2 * S2 / c
    if method == 'surrogate_uncalibrated':       # wrong: surrogate validated on Pb/Al used as is
        return te2 * S2
    if method == 'surrogate_pb_calibrated':      # wrong: calibration factor transferred from Pb
        return te2 * S2 / pb_calibration(weight=weight, tensor=tensor)
    if method == 'ratio_to_harmonic':            # TE+anh / harmonic ratio held at its 1000 K value
        return ta1 / harm(1000) * harm(2000)
    if method == 'proportional_T':
        return ta1 * 2000 / 1000
    if method == 'cancellation':                 # wrong: TE and anharmonicity cancel (as for Pb, Al)
        return harm(2000)
    if method == 'constant_anharmonic_fraction':  # wrong: TE+anh / TE held at its 1000 K value
        return te2 * ta1 / te1
    raise ValueError(method)


def truths():
    return {'al_te_change_percent_481K': al_te_change_481(),
            'pb_te_resistivity_610K': rho(pb_te(610), 610),
            'pb_anharmonic_change_percent_0K_lattice': {str(T): 100 * (rho(pb_a0(T), T) / rho(PB_HARM, T) - 1)
                                                        for T in (410, 610)},
            'pb_te_anharm_resistivity_510K': rho(pb_ta(510), 510),
            'nb_te_anharm_resistivity_2000K': rho('Nb/4-epw/3-TE_and_anharm/2000K', 2000)}


if __name__ == '__main__':
    print(json.dumps(truths(), indent=1))
    print('Al TE change 481:', {k: round(al_te_change_481(**kw), 3) for k, kw in
                                {'mesh50_quad': {}, 'mesh50_lin': dict(order=1), 'mesh70': dict(ref='Al/4-epw/1-DFPT/ibte.out-mesh70'),
                                 'mesh60': dict(ref='Al/4-epw/1-DFPT/ibte.out-mesh60'), 'serta': dict(tensor='serta')}.items()})
    print('Pb TE 610:', {m: round(pb_te_610(m), 3) for m in ('phonon_e_linear', 'phonon_e_const', 'phonon_e_quadratic', 'phonon_only',
                                                         'R_quadratic', 'rho_cubic', 'rho_quadratic', 'rho_linear')})
    print('Pb TE 610 uniform weight:', {m: round(pb_te_610(m, weight='uniform'), 3) for m in ('phonon_e_linear', 'phonon_only')})
    for T in (410, 610):
        print('Pb anh0K', T, {k: round(pb_anh0K(T, **kw), 3) for k, kw in {
            'q2_cal': {}, 'q2_uncal': dict(calibrate=False), 'q2_cal_110_310': dict(calib_T=(110, 210, 310)),
            'q2_cal_with_10K': dict(calib_T=(10, 110, 210, 310)), 'long_q2_cal': dict(weight='long_q2'),
            'uniform_cal': dict(weight='uniform'), 'inv2_moment_uncal': dict(weight='inv2_moment', calibrate=False),
            'inv2_moment_cal': dict(weight='inv2_moment'), 'omega_log_uncal': dict(weight='omega_log', calibrate=False),
            'omega_log_cal': dict(weight='omega_log'),
            'T_lattice_transfer': dict(mode='T_lattice_transfer'), 'serta': dict(tensor='serta')}.items()})
    print('Pb TE+anh 510:', {k: round(pb_te_anh_510(**kw), 3) for k, kw in {
        'q2': {}, 'long_q2': dict(weight='long_q2'), 'uniform': dict(weight='uniform'), 'inv2_moment': dict(weight='inv2_moment'),
        'omega_log': dict(weight='omega_log'), 'R_harm_linear': dict(method='R_harm_linear'),
        'linear_fit': dict(method='linear_fit'), 'linear_fit_with_10K': dict(method='linear_fit_with_10K'),
        'sep_mult': dict(method='separable_multiplicative'), 'sep_add': dict(method='separable_additive'),
        'serta': dict(tensor='serta'), 'q2_cal_with_10K': dict(calib_T=(10, 110, 210, 310))}.items()})
    print('Nb TE+anh 2000:', {m: round(nb_te_anh_2000(m), 3) for m in ('surrogate_nb', 'surrogate_uncalibrated', 'surrogate_pb_calibrated',
                                                                    'ratio_to_harmonic', 'proportional_T', 'cancellation', 'constant_anharmonic_fraction')})
    print('Nb TE+anh 2000 other weights:', {w: round(nb_te_anh_2000('surrogate_nb', weight=w), 3) for w in ('long_q2', 'uniform', 'inv2_moment', 'omega_log')},
          'serta', round(nb_te_anh_2000('surrogate_nb', tensor='serta'), 3), 'long_q2 uncal', round(nb_te_anh_2000('surrogate_uncalibrated', weight='long_q2'), 3))
