"""Reference values, ablation and robustness values, private data and calibration for ge-pressure-polymorph-pathways.

    python3 -B make_task.py        (after build_records.py and reference.py)

Truth (released data + released potential, nothing from the literature):
  * near-equilibrium frames of every archived record of the five phases (common.clean_frames);
  * for ST12, BC8 and R8 the archived compressed/expanded structures keep their zero-pressure internal
    coordinates (ISIF = 6 relaxations and scaled single points); each such frame's DFT energy is corrected by the
    released potential's fixed-cell ionic relaxation energy (reference.py, lab code path);
  * third-order Birch-Murnaghan E(V) fit per phase; equal-enthalpy pressures; E0 differences; E0 ordering.
"""
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline

from common import POOL, PHASES, METASTABLE, EV_A3_GPA, all_frames, clean_frames, fit, crossing, enthalpy

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/ge-pressure-polymorph-pathways'
REF = json.loads((HERE/'ref_dp.json').read_text())
CORR = REF['relax_corrections']
TOL = {'transition_pressures_GPa': {'CD_to_bSn': 0.2, 'bSn_to_ST12': 0.35, 'bSn_to_BC8': 0.2},
       'zero_pressure_energy_differences_meV_per_atom': {'ST12_minus_BC8': 1.0, 'R8_minus_BC8': 1.0}}
ORACLE_RECORDS = {'CD': [35, 36, 37, 38, 98], 'bSn': [25, 26, 28, 30, 61], 'ST12': [1, 48, 49, 50, 51],
                  'BC8': [30, 31, 32, 33, 61], 'R8': [29, 30, 31, 32]}


def corr(o):
    c = CORR[o['tag']]
    return (c['dp_relaxed_eV'] - c['dp_unrelaxed_eV']) / c['natoms']


def payload_from_fits(F, form='bm3'):
    c = lambda a, b: crossing(F[a], F[b], form)
    return {'transition_pressures_GPa': {'CD_to_bSn': c('CD', 'bSn'), 'bSn_to_ST12': c('bSn', 'ST12'), 'bSn_to_BC8': c('bSn', 'BC8')},
            'zero_pressure_energy_differences_meV_per_atom': {'ST12_minus_BC8': (F['ST12'][0] - F['BC8'][0]) * 1e3,
                                                              'R8_minus_BC8': (F['R8'][0] - F['BC8'][0]) * 1e3},
            'metastable_order_at_0GPa': sorted(METASTABLE, key=lambda g: F[g][0]),
            'extra': {'bSn_to_R8': c('bSn', 'R8'), 'ST12_to_BC8': c('ST12', 'BC8'),
                      'E0_BC8_minus_CD_meV': (F['BC8'][0] - F['CD'][0]) * 1e3,
                      'fits': {g: [float(x) for x in F[g]] for g in F}}}


def analyse(sel, correct=METASTABLE, form='bm3'):
    F = {}
    for g in PHASES:
        fs = sel[g]
        F[g] = fit([o['V'] for o in fs], [o['E'] + (corr(o) if g in correct else 0.0) for o in fs], form)
    return payload_from_fits(F, form)


def from_H_functions(Hf):
    """Payload from per-phase enthalpy functions H(P) (eV/atom)."""
    def c(a, b):
        Ps = np.linspace(-1.9, 11.9, 691)
        d = np.array([Hf[a](p) - Hf[b](p) for p in Ps])
        i = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
        return float(Ps[i[0]] - d[i[0]] * (Ps[i[0] + 1] - Ps[i[0]]) / (d[i[0] + 1] - d[i[0]])) if len(i) else float('nan')
    return {'transition_pressures_GPa': {'CD_to_bSn': c('CD', 'bSn'), 'bSn_to_ST12': c('bSn', 'ST12'), 'bSn_to_BC8': c('bSn', 'BC8')},
            'zero_pressure_energy_differences_meV_per_atom': {'ST12_minus_BC8': float(Hf['ST12'](0) - Hf['BC8'](0)) * 1e3,
                                                              'R8_minus_BC8': float(Hf['R8'](0) - Hf['BC8'](0)) * 1e3},
            'metastable_order_at_0GPa': sorted(METASTABLE, key=lambda g: float(Hf[g](0))),
            'extra': {'bSn_to_R8': c('bSn', 'R8'), 'ST12_to_BC8': c('ST12', 'BC8')}}


def strip(p):
    return {k: v for k, v in p.items() if k != 'extra'}


def main():
    full = {g: clean_frames(g) for g in PHASES}
    oracle = {g: [o for o in full[g] if int(o['key'].split('_')[-1]) in ORACLE_RECORDS[g] and o['nf'] > 1] for g in PHASES}
    reference = analyse(full)
    # ---------------------------------------------------------------- robustness (defensible variants)
    relax_only = {g: [o for o in full[g] if o['nf'] > 1] for g in PHASES}
    robust = {
        'vinet_eos': analyse(full, form='vinet'),
        'oracle_records_only': analyse(oracle),
        'oracle_records_without_highest_pressure_run': analyse(
            {g: [o for o in oracle[g] if o['key'] != max(oracle[g], key=lambda o: -o['V'])['key']] for g in PHASES}),
        'correct_all_phases': analyse(full, correct=PHASES),
        'stricter_hydrostatic_filter': analyse({g: clean_frames(g, dev=0.2) for g in PHASES}),
        'relaxation_runs_only': analyse(relax_only),
    }
    # harmonic variant: correction from DFT forces with the DP Hessian is computed in reference_harmonic (README);
    # ---------------------------------------------------------------- ablations (wrong paths from the real data)
    abl = {}
    # dp1: rattled training configurations included (all frames of all records, no filter)
    allf = {g: all_frames(g) for g in PHASES}
    abl['rattled_frames_included'] = analyse(allf, correct=[])
    # dp2: enthalpy from the printed stress of each frame, H_i = E_i + P_i V_i, quadratic H(P) per phase
    Hf = {}
    for g in PHASES:
        fs = full[g]
        P = np.array([o['P'] for o in fs])
        H = np.array([o['E'] + (corr(o) if g in METASTABLE else 0) + o['P'] / EV_A3_GPA * o['V'] for o in fs])
        Hf[g] = np.poly1d(np.polyfit(P, H, 2))
    abl['enthalpy_from_printed_stress'] = from_H_functions(Hf)
    # dp3: linear estimate P = dE0/dV0 (no enthalpy / common tangent)
    Fr = {g: np.array(reference['extra']['fits'][g]) for g in PHASES}
    lin = lambda a, b: (Fr[b][0] - Fr[a][0]) / (Fr[a][1] - Fr[b][1]) * EV_A3_GPA
    p3 = json.loads(json.dumps(strip(reference)))
    p3['transition_pressures_GPa'] = {'CD_to_bSn': lin('CD', 'bSn'), 'bSn_to_ST12': lin('ST12', 'bSn'), 'bSn_to_BC8': lin('BC8', 'bSn')}
    abl['linear_dE_over_dV'] = p3
    # dp4: archived metastable structures used as they are (no internal relaxation)
    abl['archived_structures_uncorrected'] = analyse(full, correct=[])
    # dp5: surrogate in place of DFT
    curves = REF['dp_relaxed_curves']
    dpH = {g: CubicSpline([q['P_GPa'] for q in curves[g]], [q['H_eV_atom'] for q in curves[g]]) for g in PHASES}
    abl['surrogate_only'] = from_H_functions(dpH)
    dftH = {g: (lambda P, p=Fr[g]: enthalpy(p, P)[0]) for g in PHASES}
    abl['surrogate_for_metastable_phases'] = from_H_functions({**{g: dftH[g] for g in ['CD', 'bSn']},
                                                               **{g: dpH[g] for g in METASTABLE}})
    # dp6: residual-force filter applied to the metastable phases (compressed frames discarded, EOS extrapolated)
    sel6 = {g: ([o for o in full[g] if o['F'] < 0.05] if g in METASTABLE else full[g]) for g in PHASES}
    F6 = {g: fit([o['V'] for o in sel6[g]], [o['E'] for o in sel6[g]]) for g in PHASES}
    p6 = payload_from_fits(F6)
    for k in p6['transition_pressures_GPa']:
        if not np.isfinite(p6['transition_pressures_GPa'][k]):
            p6['transition_pressures_GPa'][k] = crossing(F6['bSn'], F6[k.split('_to_')[1]], lo=-10, hi=60)
    abl['force_filter_discards_compressed_frames'] = p6
    # dp7: PSTRESS taken as the pressure of the final relaxation step (H = E + PSTRESS V), linear interpolation
    def pst_points(g):
        pts = {}
        for o in all_frames(g):
            if o['nf'] > 1 and o['j'] == o['nf'] - 1:
                pts[o['pst']] = o['E'] + (corr(o) if g in METASTABLE and o['tag'] in CORR else 0) + o['pst'] / EV_A3_GPA * o['V']
        Ps = sorted(pts)
        return lambda P, Ps=Ps, Hs=[pts[p] for p in Ps]: float(np.interp(P, Ps, Hs))
    abl['pstress_as_pressure'] = from_H_functions({g: pst_points(g) for g in PHASES})
    # A wrong path that finds no crossing in range would still report something; fields without a finite
    # value (and the zero-pressure fields of the pressure-only ablation dp7) are set to the reference so that
    # each ablation is judged only on the quantities its error actually changes.
    for name, p in abl.items():
        for sec in ('transition_pressures_GPa', 'zero_pressure_energy_differences_meV_per_atom'):
            for k, v in p[sec].items():
                if not np.isfinite(v) or (name == 'pstress_as_pressure' and sec.startswith('zero')):
                    p[sec][k] = reference[sec][k]
    # ---------------------------------------------------------------- write
    ref = strip(reference)
    hidden = {'reference': ref, 'tolerances': TOL, 'extra_reference': reference['extra']}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    vals = [*ref['transition_pressures_GPa'].values(), *ref['zero_pressure_energy_differences_meV_per_atom'].values(),
            reference['extra']['bSn_to_R8'], reference['extra']['ST12_to_BC8']]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps([round(v, 6) for v in vals]) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps(
        {'ablations': {k: strip(v) for k, v in abl.items()}, 'robustness': {k: strip(v) for k, v in robust.items()}},
        indent=1) + '\n')
    # ---------------------------------------------------------------- calibration printout
    def line(name, p):
        t, d = p['transition_pressures_GPa'], p['zero_pressure_energy_differences_meV_per_atom']
        ok = all(abs(t[k] - ref['transition_pressures_GPa'][k]) <= TOL['transition_pressures_GPa'][k] for k in t) and \
            all(abs(d[k] - ref['zero_pressure_energy_differences_meV_per_atom'][k]) <= TOL['zero_pressure_energy_differences_meV_per_atom'][k] for k in d) and \
            p['metastable_order_at_0GPa'] == ref['metastable_order_at_0GPa']
        print('%-46s Pt %6.3f  ST12 %6.3f  BC8 %6.3f | dST %6.2f dR8 %6.2f | %-18s %s' % (
            name, t['CD_to_bSn'], t['bSn_to_ST12'], t['bSn_to_BC8'], d['ST12_minus_BC8'], d['R8_minus_BC8'],
            ','.join(p['metastable_order_at_0GPa']), 'PASS' if ok else 'fail'))
    line('REFERENCE', reference)
    print('   extra:', {k: v for k, v in reference['extra'].items() if k != 'fits'})
    for k, v in robust.items():
        line('robust: ' + k, v)
    for k, v in abl.items():
        line('ablation: ' + k, v)
    calib = {'reference': ref, 'tolerances': TOL,
             'robustness': {k: strip(v) for k, v in robust.items()}, 'ablations': {k: strip(v) for k, v in abl.items()}}
    (HERE/'calibration.json').write_text(json.dumps(calib, indent=1) + '\n')


if __name__ == '__main__':
    main()
