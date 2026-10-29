"""Write the staged, lab, private and ablation data for gan-irradiation-defect-levels.

    cd benchmarks/construction/gan-irradiation-defect-levels && python3 -B make_task.py

Inputs: pool.json (extract.py) and dielectric.json (make_dielectric.py). Run in the runtime image.
"""
import copy
import os
import json
from pathlib import Path
import numpy as np
import model as M
import reference as R

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/gan-irradiation-defect-levels'
BUDGET = 43
THRESHOLDS = {'gap_tol_eV': 0.03, 'level_tol_eV': 0.12, 'level_window_eV': 0.3, 'merge_width_eV': 0.1, 'binding_tol_eV': 0.15,
              'frenkel_tol_eV': {'N': 0.12, 'Ga': 0.25}}
ARCHIVE_DESCRIPTION = (
    'Includes VASP structure files (POSCAR & CONTCAR) and select output files (OUTCAR & OSZICAR) of hybrid '
    'functional (HSE06) density-functional theory calculations for gallium nitride (GaN). Files are included for '
    'the perfect structure (96-atom supercell, no defects), native point defects (Ga/N vacancies, interstitials, '
    'and antisites), a Ga/N divacancy (VGa-VN), and a carbon substitutional atom on a N site (CN). Various charge '
    'states are included for each defect type. Information about the INCAR and KPOINT settings are detailed in the '
    'OUTCAR files.')


def r6(x):
    return [round(float(v), 6) for v in x]


def lab_record(rec):
    nspin = len(rec['eigenvalues'])
    return {
        'run_settings': rec['incar'], 'potcar_titles': rec['potcar'],
        'kpoint_reciprocal_coordinates': rec['kpoint'], 'nelect': rec['nelect'], 'nions': rec['nions'],
        'nupdown': rec['nupdown'], 'ionic_steps': rec['ionic_steps'],
        'reached_required_accuracy': rec['reached_accuracy'],
        'energy_sigma_to_0_eV': rec['energy_sigma0_eV'], 'free_energy_TOTEN_eV': rec['toten_eV'],
        'fermi_energy_eV': rec['efermi_eV'], 'total_magnetization_muB': rec['magnetization_muB'],
        'cell_A': rec['cell_A'], 'symbols': rec['symbols'],
        'initial_fractional_positions': [r6(p) for p in rec['frac_start']],
        'final_fractional_positions': [r6(p) for p in rec['frac_final']],
        'eigenvalues_eV_and_occupations': {('spin_up', 'spin_down')[s]: rec['eigenvalues'][s] for s in range(nspin)},
        'core_potentials_eV': rec['core_potentials_V'],
        'core_potentials_note': "OUTCAR block 'average (electrostatic) potential at core' of the final step, per atom in the order of symbols",
        'site_spin_moments_muB': rec['site_moments_muB'],
    }


def payload_from(r):
    return {'band_gap_eV': r['gap_eV'],
            'transition_levels': {d: r['levels'][d] for d in R.GRADED_LEVEL_DEFECTS},
            'divacancy_binding_energy_eV': r['binding_energy_eV'],
            'frenkel_pair_energy_eV': dict(r['frenkel_pair_energy_eV'])}


def adjacent_levels(E, vbm, qs):
    """Naive ladder: (q/q-1) for every consecutive pair, ignoring which states are ever stable."""
    return [{'from': q, 'to': q - 1, 'level_eV': E[q - 1] - E[q] - vbm} for q in qs]


def main():
    pool = M.load_pool()
    diel = json.loads((HERE/'dielectric.json').read_text())
    assert not diel.get('placeholder') or os.environ.get('ALLOW_PLACEHOLDER'), 'run make_dielectric.py first'
    ids = sorted(pool)
    calculations = {n: lab_record(pool[n]) for n in ids}
    data = {'calculations': calculations,
            'dielectric': {'eps_inf': diel['eps_inf'], 'eps_0': diel['eps_0'],
                           'axes': 'Cartesian axes of the supercell (x, y, z); z is the wurtzite c axis',
                           'description': {'eps_inf': 'high-frequency (electronic) dielectric tensor',
                                           'eps_0': 'static dielectric tensor (electronic + ionic)'},
                           'method': diel['method']}}
    catalog = {
        'archive': 'A. S. Hauck, M. Jin, B. R. Tuttle, "Structures and energies of GaN native point defects, divacancy, '
                   'and carbon substitutional from hybrid functional calculations", Zenodo (2026), '
                   'doi:10.5281/zenodo.21893631, CC BY 4.0; used in J. Appl. Phys. 139, 245702 (2026)',
        'archive_description': ARCHIVE_DESCRIPTION,
        'calculations': ids,
        'naming': 'Folder names as released: perfect and its subfolders are defect-free supercells; defect folders '
                  'are <defect><site label>_<charge>, e.g. Va_Ga1_-2 (Ga vacancy, q = -2), Va_N1 (N vacancy), Ga_i1/'
                  'Ga_i2 and N_i1/N_i2 (two interstitial configurations each), Ga_N1 and N_Ga1 (antisites), VV_GaN1 '
                  '(V_Ga-V_N divacancy), C_N1 (C on a N site).',
        'requests': {'calculation': {'kind': 'calculation', 'id': '<id from calculations>', 'cost': 1},
                     'dielectric': {'kind': 'dielectric', 'cost': 1}},
        'budget': BUDGET,
    }
    # ---------------- reference and variants
    ref = R.analyse()
    reference = {'band_gap_eV': ref['gap_eV'], 'vbm_eV': ref['vbm'], 'cbm_eV': ref['cbm'],
                 'levels': {d: ref['levels'][d] for d in R.GRADED_LEVEL_DEFECTS},
                 'divacancy_binding_energy_eV': ref['binding_energy_eV'],
                 'frenkel_pair_energy_eV': ref['frenkel_pair_energy_eV'],
                 'binding_states': ref['binding_states'], 'frenkel_states': ref['frenkel_states'],
                 'other_levels': {d: ref['levels'][d] for d in ('VV', 'Ga_N', 'N_Ga')}}
    hidden = {'reference': reference, 'thresholds': THRESHOLDS}
    values = [ref['gap_eV'], ref['binding_energy_eV'], ref['frenkel_pair_energy_eV']['N'],
              ref['frenkel_pair_energy_eV']['Ga']]
    values += [l['level_eV'] for d in R.GRADED_LEVEL_DEFECTS for l in ref['levels'][d]]
    variants = {
        'edges_R_point': dict(edges='R_point'), 'edges_efermi': dict(edges='efermi'),
        'edges_Ga_d': dict(edges='gamma_Ga_d'), 'host_gamma_only': dict(host='perfect/kpt_gam'),
        'eps_inf': dict(eps='eps_inf'), 'no_correction': dict(scheme='none'), 'pc_only': dict(scheme='pc_only'), 'alignment_wrong_sign': dict(scheme='efnv_wrong_sign'),
        'single_site_i1': dict(sites={'Ga_i': 'Ga_i1', 'N_i': 'N_i1'}),
        'single_site_i2': dict(sites={'Ga_i': 'Ga_i2', 'N_i': 'N_i2'}),
        'neutral_states': dict(binding_charge_states='neutral', fp_charge_states='neutral'),
        'edges_aligned_N_sites': dict(edges='gamma_aligned_N_sites'),
        # defensible variations
        'edges_aligned_all_sites': dict(edges='gamma_aligned_all_sites'), 'eps_0_isotropic': dict(eps='eps_0_iso'),
        'region_beyond_5.5A': dict(region=5.5), 'region_beyond_4.0A': dict(region=4.0),
    }
    results = {k: R.analyse(**kw) for k, kw in variants.items()}
    abl = {k: payload_from(results[k]) for k in ['edges_R_point', 'edges_efermi', 'edges_Ga_d', 'host_gamma_only',
                                                  'eps_inf', 'no_correction', 'pc_only', 'single_site_i1',
                                                  'alignment_wrong_sign']}
    E = {d: {int(q): v for q, v in ref['energies'][d].items()} for d in ('V_N', 'Ga_i')}
    abl['adjacent_levels'] = {'V_N': adjacent_levels(E['V_N'], ref['vbm'], [3, 2, 1, 0]),
                              'Ga_i': adjacent_levels(E['Ga_i'], ref['vbm'], [3, 2, 1])}
    robust = {k: payload_from(results[k]) for k in ['edges_aligned_all_sites', 'eps_0_isotropic',
                                                     'region_beyond_5.5A', 'region_beyond_4.0A']}
    # ---------------- write
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps(data, separators=(',', ':')) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': robust,
                                                                   'oracle_expected': payload_from(ref)}, indent=1) + '\n')
    for p in (TASK/'environment/lab/public.json', TASK/'tests/public.json'):
        pub = json.loads(p.read_text())
        pub['solver_budget'] = BUDGET
        if p.parent.name == 'lab':
            pub['requests'] = catalog['requests']
            pub['catalog'] = '/app/catalog.json'
        p.write_text(json.dumps(pub, indent=1) + '\n')
    # ---------------- calibration printout
    def brief(r):
        return {'gap': round(r['gap_eV'], 4), 'Eb': round(r['binding_energy_eV'], 3),
                'FP_N': round(r['frenkel_pair_energy_eV']['N'], 3), 'FP_Ga': round(r['frenkel_pair_energy_eV']['Ga'], 3),
                'levels': {d: [(l['from'], l['to'], round(l['level_eV'], 3)) for l in r['levels'][d]]
                           for d in R.GRADED_LEVEL_DEFECTS}}
    out = {'reference': brief(ref), **{k: brief(v) for k, v in results.items()}}
    (HERE/'calibration.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out['reference'], indent=1))


if __name__ == '__main__':
    main()
