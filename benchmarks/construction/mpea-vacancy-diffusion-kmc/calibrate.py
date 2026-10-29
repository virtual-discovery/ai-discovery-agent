"""Robustness variants and ablation (wrong-path) values for mpea-vacancy-diffusion-kmc.

    ./run_in_image.sh python3 -B calibrate.py     -> calibration.json and solution/ablation_values.json (~25 min)

Variants are defensible alternative analyses trained only on purchasable hops (the 100 hidden hops excluded).
Ablations are the realistic wrong paths of design.json, computed from the archived data with the same engine.
"""
import json
import math
import numpy as np
import kmclib as K
import common as C

STATS = dict(L=10, nsteps=25000, ntraj=120, nmc=100)
HID = C.hidden_indices()
AVAIL = np.ones(K.NH, bool)
AVAIL[HID] = False
CLEAN = AVAIL & ~K.PATHO
IDS = ['H%04d' % n for n in HID]
B1 = K.ET - K.E1


def subset(frac, seed):
    m = CLEAN & (np.random.default_rng(seed).random(K.NH) < frac)
    m[:35] = CLEAN[:35]
    return m


def random_plan(n, seed):
    """The pure/dilute hops plus a random sample of alloy hops, n in total (a budget-sized purchase)."""
    rng = np.random.default_rng(seed)
    m = np.zeros(K.NH, bool)
    m[:35] = True
    alloy = np.where(AVAIL & (np.arange(K.NH) >= 35))[0]
    m[rng.choice(alloy, n - 35, replace=False)] = True
    return m & ~K.PATHO


def barrier_preds(m, kra_only=False, species_mean=False):
    if species_mean:
        km = np.array([K.KRA[CLEAN & (K.SP == s)].mean() for s in range(6)])[K.SP[HID]]
    else:
        km = K.predict_kra(m['K0'], m['W'], m['Q'] if m['Q'] is not None else np.zeros((6, 6, 6)), HID)
    de = 0 if kra_only else K.predict_de(m['V'], HID)
    return {i: float(v) for i, v in zip(IDS, km + 0.5 * de)}


def fields(m, seed, stats=STATS, with_variants=True):
    """All kMC-derived graded fields for a model."""
    out, runs = {}, {}
    for c in 'ABC':
        runs[c] = K.simulate(C.xvec(c), C.T_ETA, m, seed=seed, **stats)
    out['log10_eta_vacancy'] = {c: float(runs[c]['log10_eta']) for c in 'ABC'}
    Bt = runs['B']['Dt']
    fin = [i for i in range(6) if np.isfinite(Bt[i])]
    fast, slow = max(fin, key=lambda i: Bt[i]), min(fin, key=lambda i: Bt[i])
    out['fastest_tracer_B'], out['slowest_tracer_B'] = K.EL[fast], K.EL[slow]
    out['log10_tracer_ratio_B'] = float(np.log10(Bt[fast] / Bt[slow]))
    Bj = runs['B']['Dj']
    out['_jump_count_ratio_B'] = float(np.log10(Bj[fast] / Bj[slow]))
    Ct = runs['C']['Dt']
    fc = max([i for i in range(6) if np.isfinite(Ct[i])], key=lambda i: Ct[i])
    Ds = []
    for T in C.T_ARRHENIUS:
        r = runs['C'] if T == C.T_ETA else K.simulate(C.xvec('C'), T, m, seed=seed + 1, **stats)
        Ds.append(r['Dt'][fc])
    out['activation_energy_C_eV'] = float(C.arrhenius(C.T_ARRHENIUS, Ds))
    out['_mean_accepted_barrier_C_eV'] = float(runs['C']['acc'])
    if with_variants and m['V'] is not None:
        out['log10_eta_thermo_ideal_B'] = float(K.simulate(C.xvec('B'), C.T_ETA, C.thermo_ideal(m), seed=seed + 2, **stats)['log10_eta'])
        out['log10_eta_kra_simple_B'] = float(K.simulate(C.xvec('B'), C.T_ETA, C.kra_simple(m), seed=seed + 3, **stats)['log10_eta'])
    return out


def main():
    cal = {'stats': STATS, 'variants': {}, 'ablations': {}}
    variants = {'all_purchasable_hops': C.model(CLEAN), 'half_of_hops_seed0': C.model(subset(0.5, 0)),
                'half_of_hops_seed1': C.model(subset(0.5, 1)), 'random_800_hops_seed0': C.model(random_plan(800, 0)),
                'random_800_hops_seed1': C.model(random_plan(800, 1)), 'random_800_hops_seed2': C.model(random_plan(800, 2)),
                'with_pathological_entries': C.model(AVAIL), 'linear_kra_no_quadratic': C.model(CLEAN, quad=False),
                'two_shell_pair_hamiltonian': C.model(CLEAN, nsh=2), 'four_orbit_kra': C.model(CLEAN, n_orb=4)}
    for k, (name, m) in enumerate(variants.items()):
        f = fields(m, seed=100 + 10 * k)
        f['barrier_predictions_eV'] = barrier_preds(m)
        cal['variants'][name] = f
        print(name, json.dumps({a: b for a, b in f.items() if a != 'barrier_predictions_eV'}), flush=True)
    ref = C.model(~K.PATHO)
    av = variants['all_purchasable_hops']
    A = cal['ablations']
    ti = fields(C.thermo_ideal(ref), seed=300, with_variants=False)
    A['thermo_ideal_as_real'] = {'log10_eta_vacancy': ti['log10_eta_vacancy'], 'log10_tracer_ratio_B': ti['log10_tracer_ratio_B'],
                                 'activation_energy_C_eV': ti['activation_energy_C_eV'],
                                 'barrier_predictions_eV': barrier_preds(av, kra_only=True)}
    sm_model = {'K0': np.array([K.KRA[CLEAN & (K.SP == s)].mean() for s in range(6)]), 'W': np.zeros((6, 0, 6)), 'Q': None,
                'V': av['V']}
    sm = fields(sm_model, seed=310, with_variants=False)
    A['species_mean_barriers'] = {'log10_eta_vacancy': sm['log10_eta_vacancy'], 'log10_tracer_ratio_B': sm['log10_tracer_ratio_B'],
                                  'fastest_tracer_B': sm['fastest_tracer_B'], 'activation_energy_C_eV': sm['activation_energy_C_eV'],
                                  'barrier_predictions_eV': barrier_preds(av, species_mean=True)}
    rf = fields(ref, seed=320)
    A['jump_count_tracers'] = {'log10_tracer_ratio_B': rf['_jump_count_ratio_B']}
    K0d, Wd, Qd = K.train_kra(CLEAN, target=B1)
    A['direct_forward_barrier_fit'] = {'barrier_predictions_eV': {i: float(v) for i, v in zip(IDS, K.predict_kra(K0d, Wd, Qd, HID))}}
    info = cal['information_only'] = {}
    info['random_configurations'] = {'log10_eta_vacancy': {c: float(K.simulate(C.xvec(c), C.T_ETA, ref, seed=330, **dict(STATS, nmc=0))['log10_eta'])
                                                        for c in 'ABC'}}
    info['small_cell_short_runs'] = {'log10_eta_vacancy': {c: float(K.simulate(C.xvec(c), C.T_ETA, ref, seed=340,
                                                                            **dict(STATS, L=4, nsteps=2000, ntraj=12))['log10_eta'])
                                                        for c in 'ABC'}}
    lin = variants['linear_kra_no_quadratic']
    pure_lin = lin['K0'] + np.array([lin['W'][s][:, s].sum() for s in range(6)])
    shift = {c: float(np.dot(np.array(C.xvec(c)) / sum(C.xvec(c)), pure_lin - K.PURE_KRA) / (K.KB * C.T_ETA * math.log(10)))
             for c in 'ABC'}
    A['rom_from_linear_model_pure_limit'] = {'log10_eta_vacancy': {c: rf['log10_eta_vacancy'][c] + shift[c] for c in 'ABC'}}
    xC = np.array(C.xvec('C')) / sum(C.xvec('C'))
    A['ea_from_rom_barrier'] = {'activation_energy_C_eV': float(np.dot(xC, K.PURE_KRA))}
    A['ea_from_mean_accepted_barrier'] = {'activation_energy_C_eV': rf['_mean_accepted_barrier_C_eV']}
    A['untested_attribution'] = {'log10_eta_thermo_ideal_B': rf['log10_eta_vacancy']['B'],
                                 'log10_eta_kra_simple_B': rf['log10_eta_vacancy']['B']}
    cal['reference_model_moderate_stats'] = rf
    cal['pure_limit_linear_model_eV'] = dict(zip(K.EL, pure_lin.tolist()))
    cal['pure_kra_dft_eV'] = dict(zip(K.EL, K.PURE_KRA.tolist()))
    json.dump(cal, open(K.HERE/'calibration.json', 'w'), indent=1)
    strip = lambda d: {k: v for k, v in d.items() if not k.startswith('_')}
    (K.TASK/'solution/ablation_values.json').write_text(json.dumps(
        {'ablations': A, 'robustness': {k: strip(v) for k, v in cal['variants'].items()}}, indent=1) + '\n')
    for k, v in A.items():
        print(k, json.dumps({a: b for a, b in v.items() if a != 'barrier_predictions_eV'}))


if __name__ == '__main__':
    main()
