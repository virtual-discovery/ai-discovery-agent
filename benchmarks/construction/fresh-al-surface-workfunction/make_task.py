"""Generate the private reference data, ablation values and calibration for fresh-al-surface-workfunction.

    python3 -B benchmarks/construction/fresh-al-surface-workfunction/make_task.py

Inputs are construction runs: lab requests executed through the task's own QE code path (slab_qe.calculate,
called by qe_batch.py inside the runtime image on the builder's cpuset; reference_batch.json lists every request).
They are read from ~/.cache/ai-theorist/fresh-al-surface-workfunction/runs when present, otherwise from
reference_runs.json.gz in this folder (the parsed results with the vacuum levels pre-extracted). The script writes
environment/lab/private/hidden.json, private/hidden_values.json, solution/ablation_values.json and
calibration.json, refreshes reference_runs.json.gz, and prints the calibration table.
"""
import gzip
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/fresh-al-surface-workfunction'
RUNS = Path.home()/'.cache/ai-theorist/fresh-al-surface-workfunction/runs'
STORE = HERE/'reference_runs.json.gz'
J = 16.0217663                       # eV/A^2 -> J/m^2
FACES = ('111', '100')
REF_N = list(range(7, 15)) + [18]    # thick-slab average: all slabs from 7 layers up, two beyond the lab limit
ORACLE_N = list(range(7, 13))        # what the reference solution buys (lab limit at kpoints 18 is 12)
HIDDEN_N = 18
ANCHOR = -537.46                     # eV/atom; graded thick-slab energies are F18 - 18 * ANCHOR
RELAX_REF = [9, 10, 11, 12]
KZ_REF = {'111': 24, '100': 27}
PRIM_K = 40                          # densest primitive-cell mesh, the "converged" separate bulk reference
TOL = {'sigma_J_per_m2': 0.008, 'phi_eV': 0.03, 'd12_percent': 0.3, 'thick_energy_eV': {'111': 0.006, '100': 0.006},
       'thick_phi_eV': 0.04}
VKEYS = ('electrostatic_potential_eV', 'total_local_potential_eV')


def plateau(r, key, width):
    z, v = np.array(r['potential_z_A']), np.array(r[key])
    c, planes = r['cell_height_A'], np.array(r['atom_z_A'])
    dist = np.min(np.abs(((z[:, None] - planes[None, :]) + c / 2) % c - c / 2), axis=1)
    return float(v[dist >= dist.max() - width].mean())


def load_all():
    if RUNS.is_dir():
        data = {}
        for d in sorted(RUNS.iterdir()):
            if (d/'result.json').exists():
                r = json.loads((d/'result.json').read_text())
                if 'potential_z_A' in r:
                    r['vacuum_levels'] = {k: {str(w): plateau(r, k, w) for w in (0.5, 1.0)} for k in VKEYS}
                    for k in ('potential_z_A', 'electrostatic_potential_eV', 'total_local_potential_eV', 'pw_input'):
                        r.pop(k, None)
                data[d.name] = r
        STORE.write_bytes(gzip.compress(json.dumps(data, sort_keys=True).encode(), mtime=0))
        return data
    return json.loads(gzip.decompress(STORE.read_bytes()))


R = load_all()


def phi(r, key='electrostatic_potential_eV', width=0.5):
    return r['vacuum_levels'][key][str(width)] - r['fermi_energy_eV']


def slab(face, n, k=18, prefix='s'):
    return R[f'{prefix}{face}_N{n}_k{k}']


def eb_oriented(face, k=18, kz=None):
    kz = kz or (KZ_REF[face] if k == 18 else (12 if face == '111' else 18))
    return R[f'bulk_or{face}_k{k}_kz{kz}']['free_energy_eV_per_atom']


def area(face):
    return slab(face, 7)['surface_cell_area_A2']


def sigma_of(r, eb):
    return (r['free_energy_eV'] - r['layers'] * eb) / (2 * area(r['face'])) * J


def ladder(face, ns, eb, k=18, prefix='s', key='electrostatic_potential_eV', width=0.5):
    rows = [slab(face, n, k, prefix) for n in ns]
    return float(np.mean([sigma_of(r, eb) for r in rows])), float(np.mean([phi(r, key, width) for r in rows]))


def thick_pred(face, eb, sig):
    return HIDDEN_N * (eb - ANCHOR) + 2 * area(face) * sig / J


def fit(face, ns, k=18, prefix='s'):
    e = [slab(face, n, k, prefix)['free_energy_eV'] for n in ns]
    slope, icpt = np.polyfit(ns, e, 1)
    return float(slope), float(icpt / (2 * area(face)) * J)


def d12(n):
    r = R[f'r111_N{n}_k18_both']
    sp, d0 = r['interlayer_spacings_A'], r['bulk_interlayer_spacing_A']
    return 50.0 * ((sp[0] + sp[-1]) / d0 - 2.0)


def estimate(sig, ph, e18, ph18, d, order='phi100>phi111'):
    """A complete set of graded fields (evidence ids are added by the ablation mutator)."""
    return {'surface_energy_J_per_m2': {f: sig[f] for f in FACES}, 'work_function_eV': {f: ph[f] for f in FACES},
            'work_function_ordering': order, 'relaxation_d12_percent_111': d,
            'thick_slab_predictions': {f'{f}_N{HIDDEN_N}': {'free_energy_offset_eV': e18[f], 'work_function_eV': ph18[f]}
                                       for f in FACES}}


def standard(ns=ORACLE_N, kz=18, k=18, width=0.5, d=None, key='electrostatic_potential_eV'):
    sig, ph, e18 = {}, {}, {}
    for f in FACES:
        eb = eb_oriented(f, k, kz if k == 18 else None)
        sig[f], ph[f] = ladder(f, ns[f] if isinstance(ns, dict) else ns, eb, k, key=key, width=width)
        e18[f] = thick_pred(f, eb_oriented(f, 18, kz) if k == 18 else eb_oriented(f, 18, 18), sig[f])
    return estimate(sig, ph, e18, dict(ph), float(np.mean([d12(9), d12(10)])) if d is None else d)


def main():
    ref = {'sigma': {}, 'phi': {}, 'thick': {}}
    detail = {}
    for f in FACES:
        eb = eb_oriented(f)
        ref['sigma'][f], ref['phi'][f] = ladder(f, [n for n in REF_N if n != HIDDEN_N], eb)
        h = slab(f, HIDDEN_N, prefix='h')
        rows = [slab(f, n) for n in REF_N if n != HIDDEN_N] + [h]
        ref['sigma'][f] = float(np.mean([sigma_of(r, eb) for r in rows]))
        ref['phi'][f] = float(np.mean([phi(r) for r in rows]))
        ref['thick'][f'{f}_N{HIDDEN_N}'] = {'free_energy_offset_eV': h['free_energy_eV'] - HIDDEN_N * ANCHOR,
                                            'work_function_eV': phi(h)}
        detail[f] = {'eb_oriented': eb, 'sigma_N': {r['layers']: round(sigma_of(r, eb), 4) for r in rows},
                     'phi_N': {r['layers']: round(phi(r), 4) for r in rows},
                     'thick_residual_meV': 1000 * (h['free_energy_eV'] - HIDDEN_N * ANCHOR - thick_pred(f, eb, ref['sigma'][f]))}
    ref['d12_111'] = float(np.mean([d12(n) for n in RELAX_REF]))
    ref['ordering'] = 'phi100>phi111' if ref['phi']['100'] > ref['phi']['111'] else 'phi111>phi100'
    detail['d12_N'] = {n: round(d12(n), 3) for n in [7] + RELAX_REF}

    # Defensible variations (must pass) -------------------------------------------------------------------------
    robust = {
        'oracle_plan_N7_12': standard(),
        'window_N8_12': standard(ns=list(range(8, 13))),
        'window_N9_12_kz9_or_12': None,
        'k24_ladders_for_sigma_phi': None,
        'wide_vacuum_window_2A': standard(width=1.0),
        'd12_from_N10_only': standard(d=d12(10)),
    }
    low = standard(ns=list(range(9, 13)))           # coarser kz for the oriented bulk: 12 for (111), 9 for (100)
    for f, kz in (('111', 12), ('100', 9)):
        eb = eb_oriented(f, 18, kz)
        s9 = float(np.mean([sigma_of(slab(f, n), eb) for n in range(9, 13)]))
        low['surface_energy_J_per_m2'][f] = s9
        low['thick_slab_predictions'][f'{f}_N{HIDDEN_N}']['free_energy_offset_eV'] = thick_pred(f, eb, s9)
    robust['window_N9_12_kz9_or_12'] = low
    k24 = standard(ns={'111': list(range(7, 13)), '100': list(range(7, 12))}, k=24)
    k18 = robust['oracle_plan_N7_12']
    k24['thick_slab_predictions'] = k18['thick_slab_predictions']
    robust['k24_ladders_for_sigma_phi'] = k24

    # Wrong paths (must fail) ----------------------------------------------------------------------------------
    abl = {}
    d_ok = float(np.mean([d12(9), d12(10)]))
    # dp1: pw.x default band count for thick slabs
    sig, ph, e18 = {}, {}, {}
    for f in FACES:
        eb = eb_oriented(f, 18, 18)
        sig[f], ph[f] = ladder(f, ORACLE_N, eb, prefix='dflt_s' if f == '111' else 'dflt2_s')
        e18[f] = thick_pred(f, eb, sig[f])
    abl['default_bands_ladder'] = estimate(sig, ph, e18, dict(ph), d_ok)
    sig2, e182 = {}, {}
    for f in FACES:
        slope, sfit = fit(f, ORACLE_N, prefix='dflt_s' if f == '111' else 'dflt2_s')
        sig2[f], e182[f] = sfit, HIDDEN_N * (slope - ANCHOR) + 2 * area(f) * sfit / J
    abl['default_bands_slope_fit'] = estimate(sig2, ph, e182, dict(ph), d_ok)
    # dp2: separately converged primitive-cell bulk reference (Boettger divergence)
    ep = R[f'bulk_prim_k{PRIM_K}']['free_energy_eV_per_atom']
    sig, ph, e18 = {}, {}, {}
    for f in FACES:
        sig[f], ph[f] = ladder(f, ORACLE_N, ep)
        e18[f] = thick_pred(f, ep, sig[f])
    abl['primitive_bulk_reference'] = estimate(sig, ph, e18, dict(ph), d_ok)
    # dp3: straight-line fit of E_N over the purchasable, quantum-size-oscillating window
    sig, e18 = {}, {}
    for f in FACES:
        slope, sfit = fit(f, ORACLE_N)
        sig[f], e18[f] = sfit, HIDDEN_N * (slope - ANCHOR) + 2 * area(f) * sfit / J
    abl['slope_fit_window_N7_12'] = estimate(sig, k18['work_function_eV'], e18, dict(k18['work_function_eV']), d_ok)
    sig, e18 = {}, {}
    for f in FACES:
        slope, sfit = fit(f, list(range(8, 13)))
        sig[f], e18[f] = sfit, HIDDEN_N * (slope - ANCHOR) + 2 * area(f) * sfit / J
    abl['slope_fit_window_N8_12'] = estimate(sig, k18['work_function_eV'], e18, dict(k18['work_function_eV']), d_ok)
    # dp4: vacuum level from the total local potential (includes exchange-correlation)
    abl['total_local_potential_vacuum_level'] = standard(key='total_local_potential_eV')
    # dp5: quantum-size effects: thin slabs only, or convergence judged from two thin slabs
    abl['thin_slab_ladder_N3_6'] = standard(ns=[3, 4, 5, 6])
    abl['single_7_layer_slab'] = standard(ns=[7])
    # dp6: relaxation from a thin slab
    abl['relaxation_from_7_layers'] = standard(d=d12(7))
    # dp7: cheap 12x12 mesh taken as converged for sigma and phi
    sig, ph, e18 = {}, {}, {}
    for f in FACES:
        sig[f], ph[f] = ladder(f, ORACLE_N, eb_oriented(f, 12), k=12)
        e18[f] = thick_pred(f, eb_oriented(f, 18, 18), sig[f])
    abl['k12_ladder_for_sigma_phi'] = estimate(sig, ph, e18, dict(ph), d_ok)
    # dp8: one surface per slab in the normalisation
    one = standard()
    one['surface_energy_J_per_m2'] = {f: 2 * v for f, v in one['surface_energy_J_per_m2'].items()}
    abl['area_counted_once'] = one

    hidden = {'reference': ref, 'tolerances': TOL}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = [ref['sigma'][f] for f in FACES] + [ref['phi'][f] for f in FACES] + [ref['d12_111']]
    values += [ref['thick'][k][q] for k in ref['thick'] for q in ('free_energy_offset_eV', 'work_function_eV')]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': robust}, indent=1) + '\n')

    def errors(p):
        e = {f'sigma_{f}': p['surface_energy_J_per_m2'][f] - ref['sigma'][f] for f in FACES}
        e.update({f'phi_{f}': p['work_function_eV'][f] - ref['phi'][f] for f in FACES})
        e['d12'] = p['relaxation_d12_percent_111'] - ref['d12_111']
        for k, v in p['thick_slab_predictions'].items():
            e[f'E_{k}_meV'] = 1000 * (v['free_energy_offset_eV'] - ref['thick'][k]['free_energy_offset_eV'])
            e[f'phi_{k}'] = v['work_function_eV'] - ref['thick'][k]['work_function_eV']
        return {k: round(v, 4) for k, v in e.items()}

    calib = {'reference': ref, 'detail': detail, 'tolerances': TOL,
             'robustness': {k: errors(v) for k, v in robust.items()},
             'ablations': {k: errors(v) for k, v in abl.items()}}
    (HERE/'calibration.json').write_text(json.dumps(calib, indent=1, default=float) + '\n')
    design_path = TASK/'design.json'
    design = json.loads(design_path.read_text())
    design['calibration'] = calibration_table(ref, calib)
    design_path.write_text(json.dumps(design, indent=1) + '\n')
    print(json.dumps(calib, indent=1, default=float))


def calibration_table(ref, calib):
    """Per graded field: reference, robustness errors, nearest failing ablation errors, tolerance."""
    fields = {'sigma_111': ('surface_energy_J_per_m2.111', ref['sigma']['111'], TOL['sigma_J_per_m2']),
              'sigma_100': ('surface_energy_J_per_m2.100', ref['sigma']['100'], TOL['sigma_J_per_m2']),
              'phi_111': ('work_function_eV.111', ref['phi']['111'], TOL['phi_eV']),
              'phi_100': ('work_function_eV.100', ref['phi']['100'], TOL['phi_eV']),
              'd12': ('relaxation_d12_percent_111', ref['d12_111'], TOL['d12_percent']),
              'E_111_N18_meV': ('thick 111_N18 free_energy_offset_eV', ref['thick']['111_N18']['free_energy_offset_eV'],
                                1000 * TOL['thick_energy_eV']['111']),
              'E_100_N18_meV': ('thick 100_N18 free_energy_offset_eV', ref['thick']['100_N18']['free_energy_offset_eV'],
                                1000 * TOL['thick_energy_eV']['100']),
              'phi_111_N18': ('thick 111_N18 work_function_eV', ref['thick']['111_N18']['work_function_eV'], TOL['thick_phi_eV']),
              'phi_100_N18': ('thick 100_N18 work_function_eV', ref['thick']['100_N18']['work_function_eV'], TOL['thick_phi_eV'])}
    table = {}
    for key, (label, value, tol) in fields.items():
        rob = {k: v[key] for k, v in calib['robustness'].items()}
        failing = {k: v[key] for k, v in calib['ablations'].items() if abs(v[key]) > tol}
        nearest = dict(sorted(failing.items(), key=lambda kv: abs(kv[1]))[:3])
        table[label] = {'reference': round(value, 4), 'robustness_errors': rob, 'nearest_failing_ablations_errors': nearest,
                        'tolerance': f'+-{tol}'}
    table['work_function_ordering'] = {'reference': ref['ordering'], 'nearest_ablation': 'unresolved (fails)',
                                       'margin': f"phi(100) - phi(111) = {ref['phi']['100'] - ref['phi']['111']:.3f} eV"}
    return table


if __name__ == '__main__':
    main()
