"""Compute reference values, robustness variants and ablation values from pool.json. Imported by make_task.py."""
import numpy as np
import model as M

TS_CV = [100, 200, 300, 500, 600]       # converged 10k runs bracketing 300 K (TaSe3 has no 400 K matrices)


def compute(pool):
    W = pool['wires']
    out = {'ref': {}, 'rob': {}, 'abl': {}, 'info': {}}
    # ---- stiffness (TaSe3)
    d1 = M.dyn_from_texts(W['TaSe3']['dfpt']['1d'])
    d0 = M.dyn_from_texts({k: v for k, v in W['TaSe3']['dfpt']['none'].items() if not k.endswith('dyn0')})
    k_ref = float(np.mean([M.axial_stiffness(d1, 6, f) for f in (0.005, 0.01, 0.02)]))
    out['ref']['k'] = k_ref
    out['rob']['k_q005'] = M.axial_stiffness(d1, 6, 0.05)
    out['rob']['k_T0_hessian'] = M.axial_stiffness(M.dyn_from_texts(W['TaSe3']['runs']['0K/10000']['hessian']), 6, 0.01)
    out['abl']['k_plain_asr'] = M.axial_stiffness(d0, 6, 0.01)
    # finite strain secant from the archived strained calculation
    import re
    st = pool['strained_TaSe3']
    e_s = float(re.findall(r'!\s+total energy\s+=\s+(\S+)', st['pw_output'])[-1])
    e_0 = float(re.findall(r'!\s+total energy\s+=\s+(\S+)', W['TaSe3']['dft_files']['scf.out'])[-1])
    c0 = d1['cell'][2, 2] * M.BOHR_A
    c_s = float(re.findall(r'CELL_PARAMETERS angstrom\s*\n.*\n.*\n\s*\S+\s+\S+\s+(\S+)', st['pw_input'])[0])
    eps = c_s / c0 - 1
    out['abl']['k_finite_strain'] = 2 * (e_s - e_0) * M.RY_EV / eps ** 2
    out['info']['strain'] = eps
    out['info']['k_T300_hessian'] = M.axial_stiffness(M.dyn_from_texts(W['TaSe3']['runs']['300K/10000']['hessian']), 6, 0.01)
    # ---- stress slope (TaSe3, 300 K)
    V = M.cell_volume_A3(d1)
    out['info']['V_cell_A3'] = V
    runs = W['TaSe3']['runs']
    T = np.array([0, 50, 100, 200, 300, 400, 500, 600])
    P = np.array([M.stress_zz(runs[f'{t}K/10000']['log']['final_excerpt'])[0] for t in T])
    A = np.array([M.stress_zz(runs[f'{t}K/10000']['log']['final_excerpt'], 'average')[0] for t in T])
    loc = np.isin(T, [200, 300, 400])
    dP = M.slope(T[loc], P[loc], 300)
    dA = M.slope(T[loc], A[loc], 300)
    out['info'].update(dPdT=dP, dAdT=dA)
    alpha = M.alpha_per_K(dP, V, k_ref) * 1e6
    out['ref']['alpha'] = alpha
    wide = np.isin(T, [100, 200, 300, 400, 500])
    out['rob']['alpha_wide_quadratic'] = M.alpha_per_K(M.slope(T[wide], P[wide], 300, 2), V, k_ref) * 1e6
    out['rob']['alpha_global_cubic_T0_hessian_k'] = M.alpha_per_K(M.slope(T, P, 300, 3), V, out['rob']['k_T0_hessian']) * 1e6
    out['rob']['alpha_k_q005'] = M.alpha_per_K(dP, V, out['rob']['k_q005']) * 1e6
    out['abl']['alpha_sign'] = -alpha
    out['abl']['alpha_average_stress'] = M.alpha_per_K(dA, V, k_ref) * 1e6
    out['abl']['alpha_plain_asr'] = M.alpha_per_K(dP, V, out['abl']['k_plain_asr']) * 1e6
    out['abl']['alpha_finite_strain'] = M.alpha_per_K(dP, V, out['abl']['k_finite_strain']) * 1e6
    out['abl']['alpha_half_k'] = 2 * alpha
    out['abl']['alpha_supercell_volume'] = 6 * alpha
    out['info']['alpha_T300_hessian_k'] = M.alpha_per_K(dP, V, out['info']['k_T300_hessian']) * 1e6
    # ---- heat capacity (CuC2 3 zero modes: straight chain; TaSe3 4: translations + axial rotation)
    for wire, nz in (('CuC2', 3), ('TaSe3', 4)):
        rw = W[wire]['runs']
        temps = [t for t in (0, 50, 100, 200, 300, 400, 500, 600, 800, 1000) if f'{t}K/10000' in rw and rw[f'{t}K/10000']['auxiliary']]
        F = {t: M.grid_frequencies(M.dyn_from_texts(rw[f'{t}K/10000']['auxiliary'])) for t in temps}
        H = {t: M.grid_frequencies(M.dyn_from_texts(rw[f'{t}K/10000']['hessian'])) for t in temps if rw[f'{t}K/10000']['hessian']}
        key = wire.lower()
        out['ref'][f'cv_{key}'] = M.cv_ratio(F, 300, TS_CV, 2, nz)
        out['rob'][f'cv_{key}_cubic'] = M.cv_ratio(F, 300, TS_CV, 3, nz)
        out['rob'][f'cv_{key}_linear_200_500'] = M.cv_ratio(F, 300, [200, 300, 500], 1, nz)
        conv = [t for t in temps if M.converged(rw[f'{t}K/10000']['log']['final_excerpt'])]
        out['rob'][f'cv_{key}_all_converged_cubic'] = M.cv_ratio(F, 300, conv, 3, nz)
        out['abl'][f'cv_{key}_paper_sign'] = M.cv_ratio(F, 300, TS_CV, 2, nz, sign=+1)
        out['abl'][f'cv_{key}_harmonic_sscha_freq'] = M.cv_ratio(F, 300, TS_CV, 2, nz, derivative=False)
        out['abl'][f'cv_{key}_per_atom_3N'] = M.cv_ratio(F, 300, TS_CV, 2, nz, per='atoms')
        hz = [t for t in TS_CV if t in H]
        Hpos = {t: (H[t][0], np.abs(H[t][1])) for t in hz}
        out['abl'][f'cv_{key}_hessian_freq'] = M.cv_ratio(Hpos, 300, hz, 2, nz)
        if wire == 'CuC2':
            out['abl']['cv_cuc2_with_unconverged_400K'] = M.cv_ratio(F, 300, [100, 200, 300, 400, 500, 600], 3, nz)
            out['info']['cv_cuc2_with_unconverged_400K_quadratic'] = M.cv_ratio(F, 300, [100, 200, 300, 400, 500, 600], 2, nz)
            out['info']['cuc2_converged_temps'] = conv
        # DFPT harmonic (T-independent) heat capacity: the paper's C_harm
        dd = M.dyn_from_texts(W[wire]['dfpt']['1d'])
        Fd = {300: M.grid_frequencies(dd)}
        out['abl'][f'cv_{key}_dfpt_harmonic'] = M.cv_ratio(Fd, 300, [300], 0, nz, derivative=False)
    # ---- TaSe3 instability (Hessian, 4 rigid-body modes at Gamma excluded)
    soft = {}
    for key, r in runs.items():
        if r['hessian']:
            qs, Wh = M.grid_frequencies(M.dyn_from_texts(r['hessian']))
            soft[key] = (M.softest_hessian(qs, Wh, 4), float(Wh.min()))
    out['info']['tase3_softest'] = soft
    t100 = M.crossing(650, soft['650K/100000'][0], 700, soft['700K/100000'][0])
    t10 = M.crossing(650, soft['650K/10000'][0], 700, soft['700K/10000'][0])
    out['info'].update(T_inst_100k=t100, T_inst_10k=t10)
    out['ref']['T_inst'] = 0.5 * (t100 + t10)
    out['rob']['T_inst_released_600_800_10k'] = M.crossing(600, soft['600K/10000'][0], 800, soft['800K/10000'][0])
    out['rob']['T_inst_released_600_800_100k'] = M.crossing(600, soft['600K/10000'][0], 800, soft['800K/100000'][0])
    # quadratic fit of signed omega^2 over released 500, 600, 800, 1000 K (10k)
    tt = np.array([500, 600, 800, 1000.])
    ss = np.array([np.sign(soft[f'{int(t)}K/10000'][0]) * soft[f'{int(t)}K/10000'][0] ** 2 for t in tt])
    roots = np.roots(np.polyfit(tt, ss, 2))
    out['rob']['T_inst_quadratic_fit'] = float(min(r.real for r in roots if abs(r.imag) < 1e-9 and 500 < r.real < 1000))
    out['abl']['T_inst_rotation_mode'] = 0.0
    out['abl']['T_inst_first_released_unstable'] = 800.0
    return out


if __name__ == '__main__':
    import json, sys
    from pathlib import Path
    pool = M.load_pool(Path.home()/'.cache/ai-theorist/wires-1d-anharmonic-thermal-contraction/pool.json')
    print(json.dumps(compute(pool), indent=1, default=float))
