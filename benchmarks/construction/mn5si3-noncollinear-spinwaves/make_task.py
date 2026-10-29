"""Write the catalog, lab data and private targets of mn5si3-noncollinear-spinwaves from CACHE/pool.json.

    python3 -B make_task.py CACHE        # then calibrate.py writes solution/ablation_values.json
"""
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/mn5si3-noncollinear-spinwaves'
CACHE = Path(sys.argv[1])
BUDGET = 25
THRESHOLDS = {'l_rms_rlu': 0.065, 'optical_window_meV': [8.0, 15.5], 'canting_window_deg': [9.0, 19.0]}
COSTS = {'exchange_set': 1, 'dispersion_table': 1, 'constant_energy_cut': 2, 'energy_scan': 1, 'gap_table': 1,
         'q_scan': 1}
EXPENSIVE_SCANS = {}
SOURCE = 'retrieved author-archived data (no new measurement or calculation)'


def r4(v):
    return [float('%.4g' % x) for x in v]


def main():
    P = json.loads((CACHE/'pool.json').read_text())
    labels = P['structure']['model_index_label']
    a = P['a_A']
    lat = [round(a, 4), round(a * np.sqrt(3), 4), round(a * P['c_over_a'], 4)]
    # ------------------------------------------------------------ exchange sets (labels instead of indices)
    notes = {'noncollinear_r4.5': 'KKR (JuKKR, Liechtenstein formula) couplings between Mn2a moments computed in the noncollinear '
                                  'reference state; all pairs up to 4.5 A.',
             'collinear_AFMz_r4.5': 'KKR couplings between Mn2a moments computed in a collinear reference state with the Mn2a '
                                    'moments along +-c; all pairs up to 4.5 A.',
             'noncollinear_r17': 'KKR couplings between Mn2a moments in the noncollinear reference state; all pairs up to 17 A.',
             'mn1_couplings': 'Mn1-Mn2a (2.96 A) and Mn1-Mn1 (2.41 A, along c) couplings that the authors used for a model in '
                              'which Mn1 carries a moment. They were not computed with KKR.'}
    exch = {}
    for name, bonds in P['exchange'].items():
        exch[name] = {'name': name, 'description': notes[name], 'units': 'meV', 'sign': 'J > 0 ferromagnetic',
                      'listing': 'every bond appears once from each end (i->j and j->i)',
                      'bonds': [{'i': labels[b['i']], 'j': labels[b['j']], 'vector_A': b['vec_A'], 'J_meV': b['J_meV']}
                                for b in bonds]}
    # ------------------------------------------------------------ dispersion tables (high-energy l rows withheld)
    disp, hidden_l, optical = {}, [], None
    qname = {'h': '(Qh,1,0)', 'k': '(0,Qk,1)', 'l': '(0,3,Ql)'}
    for key, rows in P['dispersion'].items():
        pts = []
        for r in rows:
            q, qe, E = r[:3]
            if key == 'h' and E > 9:          # the (2,1,0) zone-centre excitation from the withheld energy scan
                optical = {'E_meV': E, 'E_err_meV': r[3]}
                continue
            if key == 'l':
                hidden_l.append((E, q, qe))
                continue
            pts.append({'q_rlu': q, 'q_err_rlu': round(qe, 5), 'E_meV': E})
        if key != 'l':
            disp[key] = {'direction': key, 'Q': qname[key], 'T_K': 10, 'points': pts,
                         'method': 'Branch positions from fits of constant-energy scans with two Gaussians of common width '
                                   'on a constant background (q_err = 0 marks a zone-centre energy from an energy scan).'}
    # ------------------------------------------------------------ cuts (Fig. 4d = l at 11-20 meV withheld)
    cuts = {}
    for k, c in P['cuts'].items():
        if c['direction'] == 'l':
            continue
        cuts[k] = {'direction': c['direction'], 'Q': qname[c['direction']], 'energy_meV': c['energy_meV'], 'T_K': 10,
                   'x_rlu': c['x'], 'counts': c['counts'], 'error': c['error'],
                   'note': 'archived as plotted: a constant may have been added to the counts for display'}
    escans = {}
    for k, s in P['escans'].items():
        if k == 'optical_210_10K':
            continue
        escans[k] = {kk: vv for kk, vv in s.items() if kk != 'T_measured_K'}
    qscans = {}
    for k, s in P['qscans'].items():
        qscans[k] = {kk: vv for kk, vv in s.items() if kk != 'T_measured_K'}
    data = {'costs': COSTS, 'expensive_scans': EXPENSIVE_SCANS, 'exchange': exch, 'dispersion': disp, 'cuts': cuts,
            'escans': escans, 'gap_table': P['gap_table'], 'qscans': qscans}
    (TASK/'environment/lab/data.json').write_text(json.dumps(data))
    # ------------------------------------------------------------ hidden targets
    hidden_l.sort()
    targets = {}
    for E, q, qe in hidden_l:
        targets.setdefault('%g' % E, {})['low' if q < 1 else 'high'] = {'q_rlu': q, 'q_err_rlu': qe}
    assert all(len(v) == 2 for v in targets.values()) and optical is not None
    hidden = {'l_targets': targets, 'optical_210': optical, 'canting_reference_deg': P['archived_canting_deg'],
              'mn1_carries_moment': True, 'thresholds': THRESHOLDS,
              'withheld': {'l_cuts': {k: c for k, c in P['cuts'].items() if c['direction'] == 'l'},
                           'optical_scan': P['escans']['optical_210_10K']}}
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    values = [t['q_rlu'] for v in targets.values() for t in v.values()] + [optical['E_meV'], P['archived_canting_deg']]
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    HIDDEN_E = sorted(float(E) for E in targets)
    # ------------------------------------------------------------ catalog (agent)
    sites = [{'label': s['label'], 'element': s['element'], 'frac': s['frac']} for s in P['structure']['sites']]
    catalog = {
        'material': 'Mn5Si3, AFM1 phase (T = 10 K, zero field); hexagonal P6_3/mcm above T_N, described here in the '
                    'orthorhombic cell a x sqrt(3)a x c that contains the AFM1 magnetic cell',
        'lattice_A': {'a': lat[0], 'b': lat[1], 'c': lat[2], 'axes': 'x || a, y || b, z || c (orthogonal)'},
        'sites': sites,
        'site_notes': ['Mn1 and Mn2 are the two Mn Wyckoff sites of the hexagonal structure.',
                       'Mn2b sites carry no ordered moment in either AFM phase (neutron diffraction); the KKR exchange tables '
                       'therefore contain only Mn2a.',
                       'KKR local moments in the noncollinear reference state: about 2.44 muB on Mn2a and about 1.2 muB on Mn1.'],
        'momentum_convention': 'Q = (H, K, L) in reciprocal-lattice units of the orthorhombic cell',
        'phase_diagram_B_T': {'description': 'Archived H-T phase boundaries (field along c; temperatures in K, fields in T)',
                              **P['phase_diagram']},
        'budget': BUDGET,
        'costs': {'exchange_set': 1, 'dispersion_table': 1, 'constant_energy_cut': 2, 'energy_scan': 1,
                  'gap_table': 1, 'q_scan': 1},
        'exchange_sets': {k: v['description'] for k, v in exch.items()},
        'dispersion_tables': {k: {'Q': v['Q'], 'energies_meV': sorted({p['E_meV'] for p in v['points']})}
                              for k, v in disp.items()},
        'constant_energy_cuts': [{'direction': c['direction'], 'Q': c['Q'], 'energy_meV': c['energy_meV']}
                                 for c in cuts.values()],
        'energy_scans': {k: {kk: s[kk] for kk in ('Q_rlu', 'T_K', 'B_T')} for k, s in escans.items()},
        'gap_table': 'Fitted gap-mode energies at Q = (1 + q_h, 2, 0) for q_h <= 0.06 at 10 K (AFM1) and 80 K (AFM2, two modes)',
        'q_scans': {k: {kk: s[kk] for kk in ('Q', 'E_meV', 'T_K', 'B_T')} for k, s in qscans.items()},
        'to_predict': {'l_dispersion': {'Q': '(0,3,Ql)', 'T_K': 10, 'energies_meV': HIDDEN_E},
                       'zone_centre_excitation': {'Q_rlu': [2, 1, 0], 'T_K': 10, 'energy_window_meV': [5, 20]}},
    }
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    public = {'task_id': 'mn5si3-noncollinear-spinwaves', 'solver_budget': BUDGET, 'reserved_validation_charge': 0,
              'requests': {'exchange_set': {'kind': 'exchange_set', 'name': '<name from exchange_sets>', 'cost': 1},
                           'dispersion_table': {'kind': 'dispersion_table', 'direction': 'h | k', 'cost': 1},
                           'constant_energy_cut': {'kind': 'constant_energy_cut', 'direction': 'h | k',
                                                   'energy_meV': '<from constant_energy_cuts>', 'cost': 2},
                           'energy_scan': {'kind': 'energy_scan', 'name': '<name from energy_scans>', 'cost': 1},
                           'gap_table': {'kind': 'gap_table', 'cost': 1},
                           'q_scan': {'kind': 'q_scan', 'name': '<name from q_scans>', 'cost': 1}},
              'catalog': '/app/catalog.json'}
    (TASK/'environment/lab/public.json').write_text(json.dumps(public, indent=1) + '\n')
    (TASK/'tests/public.json').write_text(json.dumps({'task_id': 'mn5si3-noncollinear-spinwaves', 'solver_budget': BUDGET,
                                                      'min_report_chars': 400}, indent=1) + '\n')
    worth = len(exch) + len(disp) + 2 * len(cuts) + len(escans) + 2 + 1 + len(qscans)
    print('catalog worth', worth, 'units; budget', BUDGET, '; hidden', len(values), 'values; optical', optical)


if __name__ == '__main__':
    main()
