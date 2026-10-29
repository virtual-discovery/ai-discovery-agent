"""Generate the staged, lab, private and ablation data for co3sn2s2-dmi-diffuse.

Run inside the runtime image (numpy 1.21 / scipy 1.8, as the gate grades there), from this folder:

    docker run --rm --cpuset-cpus <set> --user $(id -u):$(id -g) -e HOME=/tmp -e PYTHONDONTWRITEBYTECODE=1 \
        -e CO3SN2S2_CACHE=/cache -v <cache>:/cache:ro -v <repo>:/repo -w /repo/benchmarks/construction/co3sn2s2-dmi-diffuse \
        ai-theorist-benchmark-runtime:20261007 python3 -B make_task.py

Writes environment/agent/catalog.json and archive/, environment/lab/data.json, private/hidden.json,
private/hidden_values.json, solution/ablation_values.json, both public.json files, and calibration.json here.
"""
import json
import shutil
from pathlib import Path
import numpy as np
import dms
import hamiltonian as H
import predict as P

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/co3sn2s2-dmi-diffuse'
BUDGET = 16
COSTS = {'narrow_scan': 2, 'wide_region': 1}
HALF = 0.48
REFLECTIONS = [(0, 3), (1, 1), (-1, 2), (1, -2), (1, 4)]          # (h, l) of allowed (h 0 l) inside the wide coverage
EXTRA = [(0, 1.5), (1, 3.5), (-1, 3.5), (0, 2)]                    # half-integer L (AFM stacking) and centring-forbidden
POSITIONS = REFLECTIONS + EXTRA
NARROW_HIDDEN = [170, 175, 200]
WIDE_HIDDEN = (180, (1, 1))
THRESHOLDS = {'narrow_profile_chi2r': 12.0, 'wide_101_profile_chi2r': 15.0,
              'canting_angle_tol_deg': 0.5, 'gap_rel_tol': 0.06}
ORACLE_NARROW = [140, 178, 185, 192, 220]
ORACLE_WIDE = [(165, (0, 3)), (180, (0, 3)), (165, (1, 1)), (300, (1, 1))]


def ints(a):
    return [[int(v) for v in row] for row in np.asarray(a)]


def floats(a, sig=8):
    return [[None if not np.isfinite(v) else float('%.*g' % (sig, v)) for v in row] for row in np.asarray(a)]


def served(scan):
    for k in ('counts', 'monitor', 'normcount'):
        assert np.allclose(scan[k], np.round(scan[k]))
    return {'qx_edges': [float('%.10g' % v) for v in scan['qx_edges']],
            'qy_edges': [float('%.10g' % v) for v in scan['qy_edges']],
            'counts': ints(scan['counts']), 'monitor': ints(scan['monitor']),
            'normalization': floats(scan['normalization']), 'normcount': ints(scan['normcount'])}


def as_scan(d):
    return {k: (np.array(v, dtype=float) if k != 'normalization' else
                np.array([[np.nan if x is None else x for x in row] for row in v], dtype=float))
            for k, v in d.items() if k in ('qx_edges', 'qy_edges', 'counts', 'monitor', 'normalization', 'normcount')}


def main():
    narrow = {T: dms.load('narrow', run) for T, run in dms.NARROW.items()}
    wide = {T: dms.load('wide', run) for T, run in dms.WIDE.items()}
    # ---------- lab data (what purchases return) ----------
    lab = {'costs': COSTS, 'narrow': {}, 'wide': {}}
    for T, scan in narrow.items():
        if T not in NARROW_HIDDEN:
            lab['narrow'][str(T)] = {'run': dms.NARROW[T], **served(scan)}
    for T, scan in wide.items():
        for h, l in POSITIONS:
            if (T, (h, l)) == WIDE_HIDDEN:
                continue
            c0, c1 = dms.nominal(h, l), dms.nominal(*WIDE_HIDDEN[1])
            if T == WIDE_HIDDEN[0]:          # no served 180 K pixel may lie inside the withheld (1 0 1) square
                assert abs(c0[0] - c1[0]) > 2*HALF + 0.03 or abs(c0[1] - c1[1]) > 2*HALF + 0.03, (h, l)
            region = dms.crop(scan, c0, HALF)
            lab['wide']['%g|%g,%g' % (T, h, l)] = {'run': dms.WIDE[T], 'hkl': [h, 0, l], **served(region)}
    # ---------- profiles (graded definition) ----------
    prof = {}
    for T, scan in narrow.items():
        v, e, n, pk = dms.profile(scan, dms.nominal(0, 3))
        prof[('n', T)] = (v, e)
    for T, scan in wide.items():
        for h, l in REFLECTIONS:
            v, e, n, pk = dms.profile(scan, dms.nominal(h, l))
            prof[('w', T, (h, l))] = (v, e)
            if (T, (h, l)) != WIDE_HIDDEN:     # served region must reproduce the full-map profile exactly
                rv, re_, rn, rpk = dms.profile(as_scan(lab['wide']['%g|%g,%g' % (T, h, l)]), dms.nominal(h, l))
                assert np.allclose(rv, v, rtol=1e-6) and np.allclose(re_, e, rtol=1e-6), (T, h, l)
    for T in lab['narrow']:
        rv, _, _, _ = dms.profile(as_scan(lab['narrow'][T]), dms.nominal(0, 3))
        assert np.allclose(rv, prof[('n', int(T))][0], rtol=1e-6)
    truth_n = {str(T): prof[('n', T)] for T in NARROW_HIDDEN}
    truth_w = prof[('w', 180, (1, 1))]
    # ---------- reference recipe and variants ----------
    def p003(narrow_T, wide_T):
        d = {T: prof[('n', T)] for T in narrow_T}
        d.update({T: prof[('w', T, (0, 3))] for T in wide_T})
        return d
    oracle = p003(ORACLE_NARROW, [165, 180])
    g = P.global_fit(oracle)
    p101 = {T: prof[('w', T, (1, 1))] for T in (165, 300, 1.5)}
    crit180 = P.critical(g, 180)
    ref_n = {str(T): P.predict(g, T) for T in NARROW_HIDDEN}
    ref_w = P.predict_101(p101, crit180)

    def scores(pn=None, pw=None):
        out = {}
        if pn is not None:
            out.update({'chi2r_%s' % T: P.chi2r(pn[T], *truth_n[T]) for T in truth_n})
        if pw is not None:
            out['chi2r_101'] = P.chi2r(pw, *truth_w)
        return out
    lst = lambda d: {T: [float(x) for x in v] for T, v in d.items()}
    abl, rob, calib = {}, {}, {}
    calib['reference_recipe'] = scores(ref_n, ref_w)
    # robustness
    allp = p003(ORACLE_NARROW + [250], [165, 180, 300])
    g_all = P.global_fit(allp)
    rob['narrow_more_data'] = lst({str(T): P.predict(g_all, T) for T in NARROW_HIDDEN})
    rob['narrow_bin_interpolation'] = lst({str(T): P.interp_bins(oracle, T) for T in NARROW_HIDDEN})
    par = P.independent_fits(oracle)
    rob['narrow_independent_fits'] = lst({str(T): P.predict_independent(par, T) for T in NARROW_HIDDEN})
    rob['wide_101_background_only'] = [float(x) for x in P.predict_101(p101, crit180, ratio=0.0)]
    rob['wide_101_xy_ratio'] = [float(x) for x in P.predict_101(p101, crit180, ratio=(1 + 1/9)/18)]
    rob['wide_101_raw_background'] = [float(x) for x in P.predict_101(p101, crit180, smooth=False)]
    calib['robustness'] = {'narrow_more_data': scores(rob['narrow_more_data']),
                           'narrow_bin_interpolation': scores(rob['narrow_bin_interpolation']),
                           'narrow_independent_fits': scores(rob['narrow_independent_fits']),
                           'wide_101_background_only': scores(pw=rob['wide_101_background_only']),
                           'wide_101_xy_ratio': scores(pw=rob['wide_101_xy_ratio']),
                           'wide_101_raw_background': scores(pw=rob['wide_101_raw_background'])}
    calib['robustness']['wide_101_background_165K_only'] = None   # filled below
    # ablations (diffuse scattering)
    g_narrow = P.global_fit(p003(ORACLE_NARROW, []))
    abl['narrow_scans_only'] = lst({str(T): P.predict(g_narrow, T) for T in NARROW_HIDDEN})
    abl['narrow_scans_only_bin_interpolation'] = lst({str(T): P.interp_bins(p003(ORACLE_NARROW, []), T) for T in NARROW_HIDDEN})
    abl['narrow_scans_only_independent_fits'] = lst({str(T): P.predict_independent(P.independent_fits(p003(ORACLE_NARROW, [])), T)
                                                     for T in NARROW_HIDDEN})
    abl['mirror_about_tc'] = lst({str(T): P.predict_mirror(g, T) for T in NARROW_HIDDEN})
    g_flat = P.global_fit(oracle, fix_b1=0.0)
    abl['temperature_independent_background'] = lst({str(T): P.predict(g_flat, T) for T in NARROW_HIDDEN})
    abl['wide_101_equal_to_003'] = [float(x) for x in P.predict_101(p101, crit180, ratio=1.0)]
    abl['wide_101_c_axis_polarisation_only'] = [float(x) for x in P.predict_101(p101, crit180, ratio=0.889)]
    abl['wide_101_in_plane_afm_correlations'] = [float(x) for x in P.predict_101(p101, crit180, ratio=4*(1 + 1/9)/2)]
    abl['wide_101_background_from_1p5K'] = [float(x) for x in P.linear_q(*p101[1.5]) + crit180/9]
    abl['wide_101_background_from_300K'] = [float(x) for x in P.linear_q(*p101[300]) + crit180/9]
    rob['wide_101_background_165K_only'] = [float(x) for x in P.linear_q(*p101[165]) + crit180/9]
    calib['robustness']['wide_101_background_165K_only'] = scores(pw=rob['wide_101_background_165K_only'])
    calib['ablations'] = {k: scores(abl[k]) for k in ('narrow_scans_only', 'narrow_scans_only_bin_interpolation',
                                                      'narrow_scans_only_independent_fits', 'mirror_about_tc',
                                                      'temperature_independent_background')}
    calib['ablations'].update({k: scores(pw=abl[k]) for k in ('wide_101_equal_to_003', 'wide_101_c_axis_polarisation_only',
                                                             'wide_101_in_plane_afm_correlations', 'wide_101_background_from_1p5K',
                                                             'wide_101_background_from_300K')})
    calib['global_fit'] = {'b0': g['b0'], 'b1': g['b1'], 's': g['s'], 'A': {str(k): v for k, v in g['A'].items()},
                           'kappa': {str(k): v for k, v in g['k'].items()}, 'chi2r': g['chi2r']}
    # ---------- Hamiltonian ----------
    ref = H.gap_and_tilt(H.terms(H.dm_authors()))
    ll = H.landau_lifshitz_k0(np.array(ref['spins']), H.terms(H.dm_authors()))
    assert abs(ll[0] - ref['gap_meV']) < 1e-3
    variants = {'reversed_generator_bond': H.gap_and_tilt(H.terms(H.dm_authors(), gen=H.reverse(H.GEN))),
                'copy_matrix_to_all_bonds': H.gap_and_tilt(H.terms(H.dm_authors(), copy=True)),
                'out_of_plane_dm_only': H.gap_and_tilt(H.terms(H.dm_authors() - np.array([[0, 0, 15.], [0, 0, 0], [-15., 0, 0]]))),
                'in_plane_dm_only': H.gap_and_tilt(H.terms(np.array([[0, 0, 15.], [0, 0, 0], [-15., 0, 0]]))),
                'no_interlayer': H.gap_and_tilt(H.terms(H.dm_authors(), jc=0.0))}
    S9 = np.array(ref['spins'])
    T0 = H.terms(H.dm_authors())
    gap_s015 = float(H.lswt(S9[[0, 3, 6]], [(H.PRIM[i], H.PRIM[j], d, M/3) for i, j, d, M in T0], [np.zeros(3)], spin=0.15)[0][0])
    abl['tilt_reversed_bond'] = variants['reversed_generator_bond']['tilt_deg']
    abl['gap_reversed_bond'] = variants['reversed_generator_bond']['gap_meV']
    abl['tilt_copy_matrix'] = variants['copy_matrix_to_all_bonds']['tilt_deg']
    abl['gap_copy_matrix'] = variants['copy_matrix_to_all_bonds']['gap_meV']
    abl['tilt_dz_only'] = variants['out_of_plane_dm_only']['tilt_deg']
    abl['gap_dz_only'] = variants['out_of_plane_dm_only']['gap_meV']
    abl['gap_in_plane_only'] = variants['in_plane_dm_only']['gap_meV']
    abl['tilt_in_plane_only'] = variants['in_plane_dm_only']['tilt_deg']
    abl['gap_spin_from_ordered_moment'] = gap_s015
    abl['gap_bonds_counted_twice'] = 2*ref['gap_meV']
    rob['tilt_rounded_ground_state'] = round(ref['tilt_deg'], 2)
    rob['gap_from_conventional_cell'] = float(np.min(H.lswt(S9, T0, [np.zeros(3)])[0]))
    calib['hamiltonian'] = {'reference': {k: v for k, v in ref.items() if k != 'spins'}, 'landau_lifshitz_k0_meV': ll.tolist(),
                            'variants': {k: {kk: vv for kk, vv in v.items() if kk != 'spins'} for k, v in variants.items()},
                            'gap_spin_0p15_meV': gap_s015}
    # ---------- write files ----------
    tc = {'narrow': {str(T): {'run': dms.NARROW[T]} for T in NARROW_HIDDEN}}
    catalog = {
        'material': 'Co3Sn2S2 single crystal (kagome Weyl ferromagnet), diffuse neutron scattering on DMC (PSI), '
                    'released by D. W. Tam et al., Materials Cloud Archive 2026 (CC BY 4.0)',
        'structure': {'space_group': 'R-3m (hexagonal setting)', 'a_A': dms.A_LAT, 'c_A': dms.C_LAT,
                      'Co_site': 'Wyckoff 9d, (1/2, 0, 1/2) and R-3m images (9 Co per hexagonal cell)',
                      'curie_temperature_K_literature': 177},
        'scattering_plane': '(H 0 L). Qx is along a* and Qy along c*, both in 1/A; the nominal position of (h 0 l) is '
                            '(h*4*pi/(sqrt(3)*a), l*2*pi/c) with the a and c above.',
        'pixels': 'Square bins of 0.015 1/A. qx_edges and qy_edges are the bin edges; arrays are indexed [ix][iy].',
        'channels': {'counts': 'detector counts', 'monitor': 'monitor counts', 'normalization': 'vanadium normalisation '
                     '(null where undefined)', 'normcount': 'number of contributing detector frames'},
        'normalized_intensity': 'counts * normcount / (monitor * normalization), as in the archived readcombined.py',
        'archive_files': ['archive/README.md', 'archive/readcombined.py', 'archive/SpinW_Co3Sn2S2_DMI.m'],
        'purchase_costs': COSTS, 'budget': BUDGET,
        'narrow_scans_available': [{'temperature_K': T, 'run': dms.NARROW[T]} for T in sorted(narrow) if T not in NARROW_HIDDEN],
        'wide_regions_available': [{'temperature_K': T, 'run': dms.WIDE[T], 'hkl': [h, 0, l],
                                    'nominal_center_A-1': [round(dms.nominal(h, l)[0], 4), round(dms.nominal(h, l)[1], 4)]}
                                   for T in sorted(wide) for h, l in POSITIONS if (T, (h, l)) != WIDE_HIDDEN],
        'wide_region_note': 'Regions are centred on the nominal position of the listed (h 0 l); positions with half-integer l '
                            'or with -h + l not a multiple of 3 are not nuclear Bragg reflections of R-3m.',
        'wide_region_half_width_A-1': HALF,
        'narrow_scans_coverage': 'the narrow scans cover an arc around (0 0 3); they are served whole',
        'withheld': {'narrow_scans': [{'temperature_K': T, 'run': dms.NARROW[T]} for T in NARROW_HIDDEN],
                     'wide_region': {'temperature_K': 180, 'run': dms.WIDE[180], 'hkl': [1, 0, 1]}},
        'profile_annuli_pixels': dms.EDGES_PIX,
    }
    agent = TASK/'environment/agent'
    (agent/'catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (agent/'archive').mkdir(exist_ok=True)
    for name in ('README.md', 'readcombined.py', 'SpinW_Co3Sn2S2_DMI.m'):
        shutil.copy(dms.CACHE/'raw'/name, agent/'archive'/name)
    (TASK/'environment/lab/data.json').write_text(json.dumps(lab, separators=(',', ':')) + '\n')
    hidden = {'narrow_profiles': {T: {'values': v.tolist(), 'errors': e.tolist()} for T, (v, e) in truth_n.items()},
              'wide_101_profile_180K': {'values': truth_w[0].tolist(), 'errors': truth_w[1].tolist()},
              'canting_angle_deg': ref['tilt_deg'], 'gap_meV': ref['gap_meV'], 'thresholds': THRESHOLDS}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = [x for T in truth_n for x in truth_n[T][0].tolist()] + truth_w[0].tolist() + [ref['tilt_deg'], ref['gap_meV']]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob}, indent=1) + '\n')
    for path, extra in ((TASK/'environment/lab/public.json', {'reserved_validation_charge': 0, 'requests': {
            'narrow_scan': {'kind': 'narrow_scan', 'temperature_K': '<from narrow_scans_available>', 'cost': COSTS['narrow_scan']},
            'wide_region': {'kind': 'wide_region', 'temperature_K': '<1.5|165|180|300>', 'hkl': '[h, 0, l] from wide_regions_available',
                            'cost': COSTS['wide_region']}}, 'catalog': '/app/catalog.json'}),
                        (TASK/'tests/public.json', {'min_report_chars': 400})):
        path.write_text(json.dumps({'task_id': 'co3sn2s2-dmi-diffuse', 'solver_budget': BUDGET, **extra}, indent=1) + '\n')
    calib['truth'] = {'narrow': {T: v[0].tolist() for T, v in truth_n.items()}, 'wide_101_180K': truth_w[0].tolist()}
    calib['visible_profiles'] = {('%s_%s' % (k[0], k[1]) if k[0] == 'n' else '%s_%s_%d%d' % (k[0], k[1], *k[2])): v[0].tolist()
                                 for k, v in prof.items()}
    (HERE/'calibration.json').write_text(json.dumps(calib, indent=1, default=float) + '\n')
    print(json.dumps({k: calib[k] for k in ('reference_recipe', 'robustness', 'ablations')}, indent=1, default=float))
    print(json.dumps({'tilt': ref['tilt_deg'], 'gap': ref['gap_meV'], 'LL': ll[:2].tolist(),
                      'variants': {k: (round(v['tilt_deg'], 3), round(v['gap_meV'], 3)) for k, v in variants.items()},
                      'gap_S0.15': gap_s015}, indent=1))


if __name__ == '__main__':
    main()
