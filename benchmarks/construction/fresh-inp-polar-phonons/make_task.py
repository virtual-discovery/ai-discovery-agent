"""Build the private reference, ablation and robustness values for fresh-inp-polar-phonons.

    python3 -B benchmarks/construction/fresh-inp-polar-phonons/make_task.py

Reads the reference calculations made by reference_runs.sh (each one a call of the task's own task_lab.run,
stored under ~/.cache/ai-theorist/fresh-inp-polar-phonons/runs/<name>/result.json), copies their results into
reference_results.json (kept in this folder), and writes
  environment/lab/private/hidden.json, environment/lab/private/hidden_values.json,
  solution/ablation_values.json
and prints the calibration table. The analysis functions are the oracle's own (solution/solve.py).
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/fresh-inp-polar-phonons'
RUNS = Path.home()/'.cache/ai-theorist/fresh-inp-polar-phonons/runs'
spec = importlib.util.spec_from_file_location('solve', TASK/'solution/solve.py')
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)

EOS = ['eos_a5.88_k6s1', 'eos_a5.92_k6s1', 'eos_a5.96_k6s1', 'eos_a6.00_k6s1', 'eos_a6.04_k6s1']
EOS65 = [n + '_e65' for n in EOS]
NEEDED = EOS + EOS65 + [
    'ref_a0_dfpt_k6s1', 'ref_a0_diel_k8s1', 'ref_a0_diel_k6s1', 'ref_a0_diel_k4s1', 'ref_a0_diel_k10s1',
    'ref_a0_diel_k8s0', 'ref_a0_dfpt_k4s1', 'ref_a0_dfpt_k8s0',
    'ref_am_dfpt_k6s1', 'ref_am_diel_k8s1', 'ref_am_diel_k4s1', 'ref_am_diel_k8s0', 'ref_am_diel_k10s1',
    'ref_ap_dfpt_k6s1', 'ref_ap_diel_k8s1', 'ref_ap_diel_k4s1', 'ref_ap_diel_k8s0', 'ref_ap_diel_k10s1',
    'ref_amh_dfpt_k6s1', 'ref_amh_diel_k8s1', 'ref_aph_dfpt_k6s1', 'ref_aph_diel_k8s1',
    'ref_strain_dfpt_k6s1', 'ref_strain_dfpt_k4s1', 'ref_aexp_dfpt_k6s1']
A0 = 5.9624
TOL = {'absolute': {'a0_A': 0.010, 'omega_TO_cm-1': 2.0, 'omega_LO_cm-1': 2.0, 'born_charge_In': 0.020,
                    'eps_inf': 0.20, 'eps_0': 0.35,
                    'strained_omega_LO_q_z_cm-1': 2.0, 'strained_omega_LO_q_x_cm-1': 2.0,
                    'strained_LO_anisotropy': 0.8},
       'relative': {'dlnZ_dlnV': 0.15, 'dln_eps_inf_dlnV': 0.20, 'gruneisen_TO': 0.10}}
CATEGORICAL = ['LOTO_splitting_under_compression']


def load_results():
    store = HERE/'reference_results.json'
    have = json.loads(store.read_text()) if store.exists() else {}
    for name in NEEDED:
        path = RUNS/name/'result.json'
        if path.exists():
            have[name] = json.loads(path.read_text())
    missing = [n for n in NEEDED if n not in have]
    if missing:
        sys.exit('missing reference runs: %s' % missing)
    store.write_text(json.dumps({n: have[n] for n in NEEDED}, indent=1) + '\n')
    return have


def tensors(r):
    """(Z* tensors with the sum rule imposed, eps_inf tensor)."""
    return S.asr_charges(r), np.array(r['dielectric_tensor_eps_inf'])


def zscalar(r, raw=False):
    return float(r['born_charges_raw']['In'][0][0]) if raw else float(S.asr_charges(r)[0][0][0])


def eps_scalar(r):
    return float(r['dielectric_tensor_eps_inf'][0][0])


def point(ph, diel=None, zfun=zscalar):
    """Cubic-cell values: TO from the dfpt run `ph`, eps and Z* from `diel` (default: ph itself)."""
    diel = ph if diel is None else diel
    z, e = zfun(diel), eps_scalar(diel)
    to, lo = S.to_lo(ph, [np.eye(3) * z, -np.eye(3) * z], np.eye(3) * e)
    return {'TO': to, 'LO': lo, 'Z': z, 'eps': e, 'V': ph['cell_volume_A3']}


def derivatives(lo_pt, hi_pt):
    dlnv = np.log(hi_pt['V'] / lo_pt['V'])
    d = lambda key, f=lambda p, k: p[k]: float(np.log(f(hi_pt, key) / f(lo_pt, key)) / dlnv)
    split = lambda p, k: p['LO'] - p['TO']
    osc = lambda p, k: p['LO'] ** 2 - p['TO'] ** 2
    out = {'dlnZ_dlnV': d('Z'), 'dln_eps_inf_dlnV': d('eps'), 'gruneisen_TO': -d('TO'),
           'dln_splitting_dlnV': d('S', split), 'dln_oscillator_dlnV': d('O', osc)}
    out['LOTO_splitting_under_compression'] = 'grows' if out['dln_splitting_dlnV'] < 0 else 'shrinks'
    return out


GRADED_DERIV = ['dlnZ_dlnV', 'dln_eps_inf_dlnV', 'gruneisen_TO', 'LOTO_splitting_under_compression']


def submission(a0, p0, deriv, strained):
    out = {'a0_A': a0, 'omega_TO_cm-1': p0['TO'], 'omega_LO_cm-1': p0['LO'], 'born_charge_In': p0['Z'],
           'eps_inf': p0['eps'], 'eps_0': p0['eps'] * (p0['LO'] / p0['TO']) ** 2}
    out.update({k: deriv[k] for k in GRADED_DERIV})
    out['strained_omega_LO_q_z_cm-1'], out['strained_omega_LO_q_x_cm-1'] = strained
    return out


def strained_scalar_tensors(r):
    """Wrong path: strained TO modes, but Z* and eps_inf replaced by their isotropic averages."""
    z, e = tensors(r)
    zs = np.array([np.eye(3) * np.trace(z[0]) / 3, np.eye(3) * np.trace(z[1]) / 3])
    es = np.eye(3) * np.trace(e) / 3
    return float(S.modes(r, (0, 0, 1), zs, es)[0][-1]), float(S.modes(r, (1, 0, 0), zs, es)[0][-1])


def strained_swapped(r):
    """Wrong path: scalar LO formula with q||z attached to the doublet TO and q||x to the singlet."""
    z, e = tensors(r)
    w, vec = S.modes(r)
    opt = [(w[i], vec[:, i]) for i in range(3, 6)]
    zpol = max(opt, key=lambda t: abs(t[1][2]) + abs(t[1][5]))[0]
    xy = [t[0] for t in opt if t[0] != zpol][0]
    omega = r['cell_volume_A3'] / S.BOHR_A ** 3
    mu = 1 / (1 / (r['atomic_masses_amu']['In'] * S.AMU_RY) + 1 / (r['atomic_masses_amu']['P'] * S.AMU_RY))
    nac = lambda zz, ee: 4 * np.pi * 2 / omega * zz ** 2 / ee / mu * S.RY_CM ** 2
    return (float(np.sqrt(xy ** 2 + nac(z[0][2][2], e[2][2]))), float(np.sqrt(zpol ** 2 + nac(z[0][0][0], e[0][0]))))


def main():
    R = load_results()
    a_grid = [R[n]['lattice_A'] for n in EOS]
    a0_fit, b0 = S.birch_murnaghan_a0(a_grid, [R[n]['total_energy_Ry'] for n in EOS])
    a0_65, _ = S.birch_murnaghan_a0(a_grid, [R[n]['total_energy_Ry'] for n in EOS65])
    p = [R[n]['pressure_GPa'] for n in EOS]
    a0_p0 = float(np.roots(np.polyfit(a_grid, p, 3))[-1].real)
    assert abs(a0_fit - A0) < 5e-5, (a0_fit, A0)

    def mean_diel(*names):
        rs = [R[n] for n in names]
        fake = json.loads(json.dumps(rs[0]))
        fake['dielectric_tensor_eps_inf'] = np.mean([r['dielectric_tensor_eps_inf'] for r in rs], axis=0).tolist()
        for el in ('In', 'P'):
            fake['born_charges_raw'][el] = np.mean([r['born_charges_raw'][el] for r in rs], axis=0).tolist()
        return fake

    # Reference: TO and force constants from 6x6x6-shifted DFPT; eps_inf and Z* from the shifted-mesh plateau
    # (mean of the 8x8x8 and 10x10x10 shifted electric-field responses); centred +-1% differences in a.
    plateau = {s: mean_diel(f'ref_{s}_diel_k8s1', f'ref_{s}_diel_k10s1') for s in ('a0', 'am', 'ap')}
    ref_pts = {s: point(R[f'ref_{s}_dfpt_k6s1'], plateau[s]) for s in ('a0', 'am', 'ap')}
    ref_deriv = derivatives(ref_pts['am'], ref_pts['ap'])
    st = R['ref_strain_dfpt_k6s1']
    ref = submission(A0, ref_pts['a0'], ref_deriv, S.strained_lo(st))
    ref['strained_LO_anisotropy_cm-1'] = ref['strained_omega_LO_q_x_cm-1'] - ref['strained_omega_LO_q_z_cm-1']

    def with_(**changes):
        out = {k: v for k, v in ref.items() if k != 'strained_LO_anisotropy_cm-1'}
        out.update(changes)
        return out

    # ---- robustness (defensible alternatives; must pass)
    rob = {}
    s8 = {s: point(R[f'ref_{s}_dfpt_k6s1'], R[f'ref_{s}_diel_k8s1']) for s in ('a0', 'am', 'ap')}
    rob['oracle_like_k8_tensors'] = submission(A0, s8['a0'], derivatives(s8['am'], s8['ap']), S.strained_lo(st))
    s6 = {s: point(R[f'ref_{s}_dfpt_k6s1']) for s in ('a0', 'am', 'ap')}
    rob['all_from_k6_dfpt'] = submission(A0, s6['a0'], derivatives(s6['am'], s6['ap']), S.strained_lo(st))
    s10 = {s: point(R[f'ref_{s}_dfpt_k6s1'], R[f'ref_{s}_diel_k10s1']) for s in ('a0', 'am', 'ap')}
    rob['k10_tensors'] = submission(A0, s10['a0'], derivatives(s10['am'], s10['ap']), S.strained_lo(st))
    half = {s: point(R[f'ref_{s}_dfpt_k6s1'], R[f'ref_{s}_diel_k8s1']) for s in ('amh', 'aph')}
    dh = derivatives(half['amh'], half['aph'])
    rob['half_percent_steps'] = with_(**{k: dh[k] for k in GRADED_DERIV})
    d1 = derivatives(s8['am'], s8['a0'])
    rob['one_sided_compression'] = with_(**{k: d1[k] for k in GRADED_DERIV})
    rob['a0_from_zero_pressure'] = with_(a0_A=a0_p0)
    rob['a0_from_65Ry_eos'] = with_(a0_A=a0_65)
    raw = {s: point(R[f'ref_{s}_dfpt_k6s1'], R[f'ref_{s}_diel_k8s1'], lambda r: zscalar(r, raw=True))
           for s in ('a0', 'am', 'ap')}
    rob['raw_born_charge_no_sum_rule'] = submission(A0, raw['a0'], derivatives(raw['am'], raw['ap']), S.strained_lo(st))
    k4ph = point(R['ref_a0_dfpt_k4s1'], plateau['a0'])
    rob['TO_from_k4_dfpt_and_k4_strained'] = with_(**{'omega_TO_cm-1': k4ph['TO'], 'omega_LO_cm-1': k4ph['LO']},
                                                   **dict(zip(['strained_omega_LO_q_z_cm-1', 'strained_omega_LO_q_x_cm-1'],
                                                              S.strained_lo(R['ref_strain_dfpt_k4s1']))))

    # ---- ablations (realistic wrong paths; must fail)
    abl = {}
    g = {'a0': point(R['ref_a0_dfpt_k8s0'])}
    g.update({s: point(R[f'ref_{s}_dfpt_k6s1'], R[f'ref_{s}_diel_k8s0']) for s in ('am', 'ap')})
    abl['gamma_centred_k8'] = submission(A0, g['a0'], derivatives(g['am'], g['ap']), S.strained_lo(st))
    k4 = point(R['ref_a0_dfpt_k6s1'], R['ref_a0_diel_k4s1'])
    abl['k4_shifted_tensors'] = with_(**{'omega_LO_cm-1': k4['LO'], 'born_charge_In': k4['Z'], 'eps_inf': k4['eps'],
                                        'eps_0': k4['eps'] * (k4['LO'] / k4['TO']) ** 2})
    mixed = {s: point(R[f'ref_{s}_dfpt_k6s1'], R[f'ref_{s}_diel_k4s1']) for s in ('am', 'ap')}
    dm = derivatives(mixed['am'], ref_pts['a0'])
    abl['side_volume_on_cheaper_mesh'] = with_(**{k: dm[k] for k in GRADED_DERIV})
    dm2 = derivatives(ref_pts['a0'], mixed['ap'])
    abl['side_volume_on_cheaper_mesh_expansion'] = with_(**{k: dm2[k] for k in GRADED_DERIV})
    # informational (not an ablation): 4x4x4 at both side volumes, centred; the k error largely cancels
    dk4 = derivatives(mixed['am'], mixed['ap'])
    ex = point(R['ref_aexp_dfpt_k6s1'])
    abl['experimental_lattice'] = with_(a0_A=5.8687, **{'omega_TO_cm-1': ex['TO'], 'omega_LO_cm-1': ex['LO'],
                                                        'born_charge_In': ex['Z'], 'eps_inf': ex['eps'],
                                                        'eps_0': ex['eps'] * (ex['LO'] / ex['TO']) ** 2})
    p0 = ref_pts['a0']
    abl['analytic_frequency_as_LO'] = with_(**{'omega_LO_cm-1': p0['TO'], 'eps_0': p0['eps']})
    eps0 = ref['eps_0']
    lo_eps0 = float(np.sqrt(p0['TO'] ** 2 + (p0['LO'] ** 2 - p0['TO'] ** 2) * p0['eps'] / eps0))
    abl['nonanalytic_term_with_eps0'] = with_(**{'omega_LO_cm-1': lo_eps0, 'eps_0': p0['eps'] * (lo_eps0 / p0['TO']) ** 2})
    abl['strained_scalar_tensors'] = with_(**dict(zip(['strained_omega_LO_q_z_cm-1', 'strained_omega_LO_q_x_cm-1'],
                                                      strained_scalar_tensors(st))))
    abl['strained_swapped_TO_pairing'] = with_(**dict(zip(['strained_omega_LO_q_z_cm-1', 'strained_omega_LO_q_x_cm-1'],
                                                          strained_swapped(st))))
    abl['strained_as_cubic'] = with_(**{'strained_omega_LO_q_z_cm-1': p0['LO'], 'strained_omega_LO_q_x_cm-1': p0['LO']})
    # Z*-only reasoning: d ln(w_LO^2 - w_TO^2)/d ln V = 2 d ln Z* - 1 with eps_inf held fixed. With
    # w_LO^2 = w_TO^2 + O:  dS/dlnV = O (dlnO/dlnV) / (2 w_LO) + gamma_TO w_TO (1 - w_TO/w_LO), S = w_LO - w_TO.
    z_only = 2 * ref_deriv['dlnZ_dlnV'] - 1
    w_lo, w_to = ref['omega_LO_cm-1'], ref['omega_TO_cm-1']
    ds = z_only * (w_lo ** 2 - w_to ** 2) / (2 * w_lo) + ref_deriv['gruneisen_TO'] * w_to * (1 - w_to / w_lo)
    abl['splitting_trend_from_born_charge_only'] = with_(LOTO_splitting_under_compression='grows' if ds < 0 else 'shrinks')
    print('Z*-only dS/dlnV (cm^-1):', ds)

    hidden = {'reference': ref, 'tolerances': TOL,
              'provenance': {'a0_birch_murnaghan_A': a0_fit, 'B0_GPa': b0, 'a0_zero_pressure_A': a0_p0,
                             'a0_65Ry_A': a0_65, 'dln_oscillator_dlnV': ref_deriv['dln_oscillator_dlnV'],
                             'dln_splitting_dlnV': ref_deriv['dln_splitting_dlnV'],
                             'runs': 'benchmarks/construction/fresh-inp-polar-phonons/reference_results.json'}}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = [v for v in ref.values() if not isinstance(v, str)] + [a0_fit, b0, ref_deriv['dln_oscillator_dlnV'],
                                                                   ref_deriv['dln_splitting_dlnV']]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1) + '\n')

    # ---- calibration table
    sys.path.insert(0, str(TASK/'environment/lab'))
    fields = [k for k in ref if k not in ('strained_LO_anisotropy_cm-1', *CATEGORICAL)] + ['strained_LO_anisotropy']

    def err(sub, field):
        if field == 'strained_LO_anisotropy':
            v = sub['strained_omega_LO_q_x_cm-1'] - sub['strained_omega_LO_q_z_cm-1']
            return v - ref['strained_LO_anisotropy_cm-1'], v
        e = sub[field] - ref[field]
        return (e / abs(ref[field]) if field in TOL['relative'] else e), sub[field]

    table = {}
    for field in fields:
        tol = TOL['relative'].get(field, TOL['absolute'].get(field))
        row = {'reference': ref['strained_LO_anisotropy_cm-1'] if field == 'strained_LO_anisotropy' else ref[field],
               'tolerance': ('relative ' if field in TOL['relative'] else 'absolute ') + str(tol),
               'robustness': {}, 'ablations_failing_this_field': {}}
        for name, sub in rob.items():
            e, v = err(sub, field)
            if abs(e) > 1e-12:
                row['robustness'][name] = {'value': v, 'error': e}
        for name, sub in abl.items():
            e, v = err(sub, field)
            if abs(e) > tol:
                row['ablations_failing_this_field'][name] = {'value': v, 'error': e}
        worst = max([abs(x['error']) for x in row['robustness'].values()], default=0.0)
        nearest = min([abs(x['error']) for x in row['ablations_failing_this_field'].values()], default=None)
        row['max_robust_error'], row['nearest_failing_ablation_error'] = worst, nearest
        table[field] = row
        flag = '' if worst <= tol else '  <-- ROBUSTNESS FAILS'
        print(f"{field:30s} ref={row['reference']:.5g} tol={tol} max_robust={worst:.4g} nearest_abl={nearest}{flag}")
    for c in CATEGORICAL:
        table[c] = {'reference': ref[c], 'tolerance': 'exact',
                    'robustness': {n: sub[c] for n, sub in rob.items()},
                    'ablations_failing_this_field': {n: sub[c] for n, sub in abl.items() if sub[c] != ref[c]}}
        print(c, table[c])
    for name, sub in abl.items():
        failing = [f for f in fields if abs(err(sub, f)[0]) > TOL['relative'].get(f, TOL['absolute'].get(f))]
        failing += [c for c in CATEGORICAL if sub[c] != ref[c]]
        print(f'ablation {name:40s} fails {failing}')
    print('informational: 4x4x4 meshes at both side volumes (centred):', {k: dk4[k] for k in GRADED_DERIV})
    print('dln_splitting_dlnV reference / robustness:', ref_deriv['dln_splitting_dlnV'],
          {n: derivatives(*pair)['dln_splitting_dlnV'] for n, pair in
           {'k8': (s8['am'], s8['ap']), 'k6': (s6['am'], s6['ap']), 'k10': (s10['am'], s10['ap']),
            'half': (half['amh'], half['aph']), 'one_sided': (s8['am'], s8['a0'])}.items()})
    (HERE/'calibration.json').write_text(json.dumps(table, indent=1) + '\n')
    print(json.dumps({'a0_fit': a0_fit, 'a0_p0': a0_p0, 'a0_65': a0_65, 'B0_GPa': b0,
                      'dln_oscillator_dlnV': ref_deriv['dln_oscillator_dlnV']}, indent=1))


if __name__ == '__main__':
    main()
