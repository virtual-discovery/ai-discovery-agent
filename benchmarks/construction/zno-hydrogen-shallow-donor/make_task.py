"""Generate the catalog, lab data, private targets and ablation values for zno-hydrogen-shallow-donor.

    python3 -B make_task.py          (reads pool.json written by extract.py)

Prints the calibration table (reference, robustness variants, ablations).
"""
import copy
import json
import math
import numpy as np
from pathlib import Path
import model as M

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/zno-hydrogen-shallow-donor'
P = json.loads((HERE/'pool.json').read_text())
T = 1000.0
LEVEL_DEFECTS = ('OH_O', 'H_O', 'V_O')
EIG_WINDOW = (-5.0, 4.0)   # eV around each run's Fermi energy

THRESHOLDS = {'band_gap_tol_eV': 0.03, 'delta_mu_tol_eV': 0.03, 'V_O_level_tol_eV': 0.05,
              'shallow_level_window_eV': [3.00, 3.30], 'fermi_level_tol_eV': 0.025,
              'log10_tol': {'V_O': 0.15, 'OH_O^+': 0.12, 'OH_O^0': 0.7},
              'ignore_levels_outside_eV': [0.0, 3.55], 'interior_margin_eV': 0.3}


def run_id(path):
    """Archive folder path -> run id (the PBE relaxation is the folder, the hybrid single point its HSE06 subfolder)."""
    return path


def pool_key_to_id(key):
    base, level = key.rsplit('/', 1)
    return base if level == 'PBE' else base + '/HSE06'


def lab_run(rec):
    out = {'archive_folder': rec['path'], 'source_file_parsed': rec['source_file'], 'INCAR': rec['incar'],
           'KPOINTS': rec['kpoints_file'], 'POTCAR_TITEL': rec['potcar'],
           'final_energy_TOTEN_eV': rec['toten_eV'], 'final_energy_sigma0_eV': rec['energy_sigma0_eV'],
           'total_magnetization_muB': rec['magnetization_muB'], 'initial_structure': rec['initial'],
           'final_structure': rec['final']}
    if rec['source_file'] == 'OUTCAR':
        out.update({'NELECT': rec['nelect'], 'ISPIN': rec['ispin'], 'ENCUT_eV': rec['encut_eV'],
                    'LHFCALC': rec['lhfcalc'], 'AEXX': rec['aexx'], 'E_fermi_eV': rec['efermi_eV'],
                    'kpoints_reciprocal': rec['kpoints'], 'kpoint_weights': rec['kweights'],
                    'core_potentials_eV': rec['core_potentials_eV'],
                    'site_moments_muB': rec.get('site_moments_muB'),
                    'relaxation_reached_required_accuracy': rec['converged_ionic'] if 'IBRION' in rec['incar'] else None})
    else:
        out.update({'NELECT': float(rec['incar']['NELECT']) if 'NELECT' in rec['incar'] else None,
                    'ionic_steps_in_log': rec['ionic_steps'], 'relaxation_reached_required_accuracy': rec['converged_ionic'],
                    'note': 'Only the run log of this calculation is served: energies, magnetization, structures.'})
    return out


def lab_eigen(rec):
    lo, hi = rec['efermi_eV'] + EIG_WINDOW[0], rec['efermi_eV'] + EIG_WINDOW[1]
    spins = {}
    for sp, ks in rec['eigenvalues'].items():
        spins['spin_' + sp] = [[[n + 1, e, o] for n, (e, o) in enumerate(b) if lo <= e <= hi] for b in ks]
    return {'E_fermi_eV': rec['efermi_eV'], 'kpoints_reciprocal': rec['kpoints'], 'kpoint_weights': rec['kweights'],
            'window_eV': [round(lo, 4), round(hi, 4)], 'columns': ['band', 'energy_eV', 'occupation'], 'eigenvalues': spins}


# ------------------------------------------------------------------ analysis pipeline (all variants)
def pipeline(level='HSE06', mu=None, energy='energy_sigma0_eV', gap_from='supercell', carriers='fd', shallow_at_cbm=False,
             **kw):
    pool = P if energy == 'energy_sigma0_eV' else {k: dict(v, energy_sigma0_eV=v['toten_eV']) for k, v in P.items()}
    mu = mu or M.chem_pots(pool, level)
    F, vbm, cbm = M.formation_energies(pool, level, mu=mu, **kw)
    gap = cbm - vbm
    if gap_from == 'bulk':
        bv, bc = M.band_edges(pool['chem_pot/ZnO/' + level])
        gap = bc - bv
    E0 = {k: v['E0'] for k, v in F.items()}
    if shallow_at_cbm:   # neutral shallow donors as D+ plus an electron at the CBM
        for d in ('OH_O', 'H_O'):
            E0[(d, 0)] = E0[(d, 1)] + gap
    V = abs(np.linalg.det(pool['supercell_3_3_2/' + level]['final']['lattice_A'])) * 1e-24
    if carriers == 'boltzmann':
        orig = M.carriers
        def boltz(ef, g, T_, me=M.M_E, mh=M.M_H):
            kT = M.KB * T_
            nc = lambda m: 2.0 * (m * M.M0C2 * kT / (2 * math.pi * M.HBARC ** 2)) ** 1.5
            return nc(me) * math.exp((ef - g) / kT), nc(mh) * math.exp(-ef / kT)
        M.carriers = boltz
        ef, conc, n, p = M.sc_fermi(E0, gap, T, 36 / V)
        M.carriers = orig
    else:
        ef, conc, n, p = M.sc_fermi(E0, gap, T, 36 / V)
    levels = {d: M.stable_levels({q: E0[(d, q)] for q in M.DEFECTS[d]}, gap)[0] for d in LEVEL_DEFECTS}
    vo_stable = M.stable_levels({q: E0[('V_O', q)] for q in (0, 1, 2)}, gap)[1]
    lg = lambda x: math.log10(x)
    return {
        'band_gap_eV': gap,
        'delta_mu_eV': {'O': mu['dmu_O'], 'H': mu['dmu_H']},
        'transition_levels': {d: [{'from': x['from'], 'to': x['to'], 'level_eV': x['level_eV']}
                                  for x in levels[d] if 0.0 <= x['level_eV'] <= gap + 1e-9]
                              for d in LEVEL_DEFECTS},
        'V_O_plus_ever_stable': 1 in vo_stable,
        'fermi_level_eV': ef,
        'log10_concentration_cm3': {'V_O': lg(sum(conc[('V_O', q)] for q in (0, 1, 2))),
                                    'OH_O^+': lg(conc[('OH_O', 1)]), 'OH_O^0': lg(conc[('OH_O', 0)])},
        '_E0': {'%s|%d' % k: v for k, v in E0.items()},
        '_mu': mu,
    }


def fields(res):
    return {k: v for k, v in res.items() if not k.startswith('_')}


def main():
    ref = pipeline()
    refE0 = ref['_E0']
    print('REFERENCE', json.dumps(fields(ref), indent=1))
    # ------------------------------------------------------------- catalog and lab data
    runs, eig = {}, {}
    catalog_runs = []
    for key, rec in P.items():
        rid = pool_key_to_id(key)
        runs[rid] = lab_run(rec)
        if 'eigenvalues' in rec:
            eig[rid] = lab_eigen(rec)
        parts = rid.split('/')
        entry = {'id': rid, 'eigenvalues_available': rid in eig}
        if parts[0] == 'chem_pot':
            entry.update(group='reference phase', phase=parts[1])
        elif len(parts) <= 2:
            entry.update(group='host supercell')
        else:
            entry.update(group='defect supercell', defect=parts[1], charge_folder=parts[2])
        catalog_runs.append(entry)
    catalog_runs.sort(key=lambda e: e['id'])
    catalog = {
        'material': 'wurtzite ZnO; VASP calculations released by Aicher, Berger, Aschauer and Diwald (2025), Materials Cloud Archive record c91vp-p9h03 (CC BY 4.0)',
        'archive_description': ('Authors\' description (paper): spin-polarized DFT with VASP; structures relaxed with PBE, '
                                'followed by HSE06 single-point energy calculations on the PBE structures; 3x3x2 supercell of '
                                'the hexagonal unit cell (72 atoms) with a 2x2x2 Gamma-centred k-mesh; reference phases Zn, ZnO, '
                                'Zn(OH)2 and the molecules O2, H2, H2O in their own cells.'),
        'folder_layout': ('Each folder holds a relaxation; its HSE06 subfolder holds the hybrid single point on the relaxed '
                          'structure. Defect folders are named <defect>/q_<charge>, where the charge is the authors\' label. '
                          'Defect names are the authors\'; the structures show what each one is.'),
        'experiment': ('EPR of ZnO particles from gas-phase synthesis (particle formation at about 1000 K), annealed in vacuum '
                       'at 873 K, shows a resonance at g = 1.96 attributed to a paramagnetic shallow donor.'),
        'purchases': {'run': 'settings (INCAR, KPOINTS, POTCAR titles, electron count), final energies, total magnetization and '
                             'site moments, Fermi energy, initial and final structures, k-points with weights, and the '
                             'per-atom average electrostatic potential at the cores, as printed by VASP (1 unit)',
                      'eigenvalues': 'final Kohn-Sham eigenvalues and occupations of both spin channels at every k-point, '
                                     'within -5 to +4 eV of the run\'s Fermi energy, with the k-point weights (1 unit)'},
        'budget': 35,
        'runs': catalog_runs,
    }
    lab_dir = TASK/'environment/lab'
    (lab_dir/'data.json').write_text(json.dumps({'runs': runs, 'eigenvalues': eig}))
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    # ------------------------------------------------------------- hidden targets
    reference = fields(ref)
    hidden = {'reference': reference, 'thresholds': THRESHOLDS,
              'categorical': {'epr_centre': {'defect': 'OH_O', 'charge': 0}, 'oh_o_site': 'BC_parallel',
                              'V_O_plus_ever_stable': False}}
    (lab_dir/'private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    hv = [reference['band_gap_eV'], reference['delta_mu_eV']['O'], reference['delta_mu_eV']['H'], reference['fermi_level_eV']]
    hv += [x['level_eV'] for d in LEVEL_DEFECTS for x in reference['transition_levels'][d]]
    hv += list(reference['log10_concentration_cm3'].values())
    (lab_dir/'private/hidden_values.json').write_text(json.dumps(hv) + '\n')
    # ------------------------------------------------------------- robustness and ablations
    rob = {
        'band_filling_unaligned_cbm': pipeline(bf_align=False),
        'shallow_donors_at_cbm': pipeline(shallow_at_cbm=True),
        'toten_energies': pipeline(energy='toten_eV'),
        'bulk_zno_gap': pipeline(gap_from='bulk'),
        'boltzmann_carriers': pipeline(carriers='boltzmann'),
        'alignment_radius_4A': pipeline(rexcl=4.0),
    }
    pbe = pipeline(level='PBE')
    abl = {
        'pbe_pipeline': pbe,
        'hydrogen_reservoir_H2': pipeline(mu=M.chem_pots(P, 'HSE06', h_limit=['H2', 'H2O'])),
        'no_band_filling': pipeline(bf=False),
        'band_filling_equal_k_weights': pipeline(bf_equal_weights=True),
    }
    hse = ref
    # consecutive-adjacent levels for V_O (every q/q-1 pair) instead of the stable ladder
    E = {q: refE0['V_O|%d' % q] for q in (0, 1, 2)}
    adjacent = [{'from': 2, 'to': 1, 'level_eV': E[1] - E[2]}, {'from': 1, 'to': 0, 'level_eV': E[0] - E[1]}]
    # HSE defects with PBE reference phases (mixed Hamiltonians)
    mixed_mu = M.chem_pots(P, 'PBE')
    out = {'robustness': {k: fields(v) for k, v in rob.items()},
           'ablations': {k: fields(v) for k, v in abl.items()},
           'V_O_adjacent_levels': adjacent}
    (TASK/'solution/ablation_values.json').write_text(json.dumps(out, indent=1) + '\n')
    # ------------------------------------------------------------- calibration printout
    def row(name, r):
        f = fields(r)
        lv = {d: [(x['from'], x['to'], round(x['level_eV'], 3)) for x in f['transition_levels'][d]] for d in LEVEL_DEFECTS}
        print(f"{name:30s} gap={f['band_gap_eV']:.4f} dmuO={f['delta_mu_eV']['O']:.3f} dmuH={f['delta_mu_eV']['H']:.3f} "
              f"EF={f['fermi_level_eV']:.4f} logVO={f['log10_concentration_cm3']['V_O']:.3f} "
              f"logOH+={f['log10_concentration_cm3']['OH_O^+']:.3f} logOH0={f['log10_concentration_cm3']['OH_O^0']:.3f} "
              f"VO+stable={f['V_O_plus_ever_stable']} {lv}")
    row('reference', ref)
    for k, v in rob.items():
        row('rob:' + k, v)
    for k, v in abl.items():
        row('abl:' + k, v)
    print('V_O adjacent levels', adjacent)
    print('mixed PBE references dmu', mixed_mu['dmu_O'], mixed_mu['dmu_H'])


if __name__ == '__main__':
    main()
