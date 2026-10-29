"""Write benchmarks/tasks/mnte-altermagnet-ahc-doping/design.json from calibration.json (make_task.py) and
oracle_run.json (record_oracle.py, optional).

    python3 -B benchmarks/construction/mnte-altermagnet-ahc-doping/write_design.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/mnte-altermagnet-ahc-doping'
CAL = json.loads((HERE/'calibration.json').read_text())
ORACLE = json.loads((HERE/'oracle_run.json').read_text()) if (HERE/'oracle_run.json').exists() else {}
A, R, TG = CAL['ablations'], CAL['robustness'], CAL['targets']
PM = CAL['per_mesh']
r = lambda x, n=1: round(x, n)
e = lambda x: float('%.3g' % x)

DPS = [
    {'id': 'dp1', 'fork': 'k-mesh convergence of sigma_xy: Berry-curvature hot spots where the SOC-split top valence bands nearly '
                          'touch make coarse meshes noisy (spikes, spurious sign changes); which mesh and smearing are converged',
     'wrong_choice': 'Stop at a 100, 150 or 200 mesh (1-3 units) with 1 meV smearing',
     'consequence': f"sign changes {A['mesh_100']['sigma_xy_sign_changes']} (100) or {A['mesh_150']['sigma_xy_sign_changes']} (150) "
                    f"instead of {TG['sigma_xy_sign_changes']}; sigma_xy(2e21) {r(A['mesh_100']['sigma_xy_S_per_cm']['2.0e21'])} (100) or "
                    f"{r(A['mesh_200']['sigma_xy_S_per_cm']['2.0e21'])} (200) instead of {r(TG['sigma_xy_S_per_cm']['2.0e21'])} S/cm",
     'ablation': 'mesh_200'},
    {'id': 'dp2', 'fork': 'Where the valence-band maximum is: off the high-symmetry lines (k_frac ~ (0.02, 0, 0.42)), above the '
                          'band-path maximum and the maximum on the DFT scf grid',
     'wrong_choice': 'VBM from the band path (as in bands_ori.dat) or from the DFT 6x6x8 eigenvalues',
     'consequence': f"hole density 0.30 eV below the VBM {e(A['vbm_from_band_path'])} or {e(A['vbm_from_dft_grid'])} instead of "
                    f"{e(TG['hole_density_0p30_below_vbm_cm3'])} cm^-3",
     'ablation': 'vbm_from_band_path'},
    {'id': 'dp3', 'fork': 'Energy reference: the VASP E-fermi (6.2098 eV, OUTCAR) lies in the gap 78 meV above the VBM; wt.in E_FERMI '
                          '(6.1301 eV) is the authors\' zero, close to the VBM',
     'wrong_choice': 'Take the OUTCAR Fermi energy as the band edge',
     'consequence': f"hole density 0.30 eV below the 'edge' {e(A['vbm_from_outcar_efermi'])} instead of "
                    f"{e(TG['hole_density_0p30_below_vbm_cm3'])} cm^-3",
     'ablation': 'vbm_from_outcar_efermi'},
    {'id': 'dp4', 'fork': 'Counting holes in spinor bands: each of the 32 Wannier bands holds one electron per k-point',
     'wrong_choice': 'Multiply by a spin degeneracy of 2',
     'consequence': 'densities doubled; sigma_xy read at half the requested density ('
                    + ', '.join(f"{k}: {r(v)}" for k, v in A['spin_degeneracy_twice']['sigma_xy_S_per_cm'].items())
                    + ' S/cm instead of ' + ', '.join(f"{k}: {r(v)}" for k, v in TG['sigma_xy_S_per_cm'].items()) + ')',
     'ablation': 'spin_degeneracy_twice'},
    {'id': 'dp5', 'fork': 'Whether the nonzero sigma_yz and sigma_zx of the archived (unsymmetrized) Hamiltonian are physical: the '
                          'authors symmetrized it with WannSymm before WannierTools; C2z forbids both',
     'wrong_choice': 'Report the service values',
     'consequence': f"max |sigma_yz|, |sigma_zx| {r(A['raw_forbidden_components']['yz'], 2)} and "
                    f"{r(A['raw_forbidden_components']['zx'], 2)} S/cm instead of 0",
     'ablation': 'raw_forbidden_components'},
    {'id': 'dp6', 'fork': 'Magnetic symmetry from the structure and moments versus the archived symmetries.dat, whose matrices '
                          '(sixfold, threefold, C2z without c/2) are not symmetries of the structure and contradict the labels '
                          'in wannsymm.out; and the Neel-along-x case (unitary C2x and C2y t(c/2))',
     'wrong_choice': 'Use symmetries.dat as written (C3 T kills sigma_xy), or treat in-plane Neel directions as equivalent',
     'consequence': 'allowed components [] for Neel along y, or ["xy"] for Neel along x',
     'ablation': 'symmetries_dat_literal'},
    {'id': 'dp7', 'fork': 'Whether the Wannier model describes the conduction bands: the frozen window ends at 6.7 eV, below the '
                          'DFT CBM; on the fitting grid the Wannier conduction bands are off by up to 1.6 eV',
     'wrong_choice': 'Compute electron-doped sigma_xy from the model',
     'consequence': f"sigma_xy {r(A['electron_doped_from_model'], 1)} S/cm reported instead of null",
     'ablation': 'electron_doped_from_model'},
    {'id': 'dp8', 'fork': 'Occupation smearing: the zero-temperature limit (the archived run used 1 meV, Eta_Arc in wt.in) versus a '
                          'large smearing that hides the mesh noise',
     'wrong_choice': 'Smooth the noisy curves with a room-temperature smearing (kT = 25 meV)',
     'consequence': f"sigma_xy(2e21) {r(A['room_temperature_smearing']['sigma_xy_S_per_cm']['2.0e21'])} instead of "
                    f"{r(TG['sigma_xy_S_per_cm']['2.0e21'])} S/cm",
     'ablation': 'room_temperature_smearing'},
]

design = {
    'task_id': 'mnte-altermagnet-ahc-doping', 'catalog_id': 'mnte-altermagnet-ahc-doping', 'start': 'released',
    'source_group': 'kluczyk_2024_mnte_ahe', 'title': 'MnTe altermagnet: anomalous Hall conductivity under hole doping',
    'decision_points': DPS,
    'graded_fields': [
        'sigma_xy_S_per_cm at 2.0e21 and 6.0e21 cm^-3 within 15% of the archived WannierTools output',
        f"sigma_xy_sign_changes == {TG['sigma_xy_sign_changes']} (5 meV grid, mu from VBM-1.00 to VBM-0.25 eV)",
        'sign_change_hole_density_cm3 within 8%',
        'hole_density_0p30_below_vbm_cm3 within 5%',
        'forbidden_components_max_abs_S_per_cm: |yz|, |zx| <= 0.5',
        'allowed_hall_components == {neel_y: [xy], neel_x: []}',
        'electron_doped_sigma_xy_S_per_cm is null'],
    'length_estimate': {
        'oracle_lab_evaluations': ORACLE.get('evaluations'), 'oracle_charge': ORACLE.get('charge_used'), 'budget': 32,
        'catalog_units': 'meshes 48-400 cost 1+1+1+2+3+5+8 = 21 units per Neel direction (42 for both), plus any number of '
                         'bands / berry_curvature / electron_count requests (1 unit each)',
        'tool_calls': '45-75', 'gpt6astra_cost_usd': '2.5-4.5',
        'why': 'Locating the off-path VBM, checking the Wannier model against the DFT eigenvalues, a mesh-convergence study of a '
               'noisy Berry-curvature integral under a budget that covers two dense meshes, the p -> mu -> sigma cascade, a '
               'magnetic-symmetry analysis that must not trust the archived symmetries.dat, and judging which model outputs '
               '(forbidden components, conduction bands) are artefacts.',
        'oracle_errors': ORACLE.get('errors')},
    'runtime_dependencies': ['numpy and scipy (runtime image); no extra layer'],
    'data_provenance': {
        'archive_url': 'https://zenodo.org/records/18863304', 'doi': '10.5281/zenodo.18863304',
        'paper_doi': '10.1103/PhysRevB.110.155201', 'files_used': ['dft_data.zip', 'fig6_script.zip'],
        'md5': {'dft_data.zip': 'cef0196e73cec5deefad7e707137b940', 'fig6_script.zip': '15e2f717c3d2c28f11c19025122769ce'},
        'members': 'benchmarks/construction/mnte-altermagnet-ahc-doping/members.md5', 'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': ('Adapted from the candidate design (see the construction README). Graded sigma_xy values are the archived '
              'WannierTools output (801 x 801 x 401 mesh, WannSymm-symmetrized Hamiltonian that is not archived) read at the '
              'chemical potentials of the requested hole densities; the lab evaluates the archived unsymmetrized Hamiltonian '
              'with the same protocol, and its converged 300/400-mesh values lie 3-7% from the archive at the graded points. '
              'Targets near the band edge (realistic dopings < 1e21 cm^-3; the design\'s E-E_F = -0.2 eV) are not graded: '
              'they are not converged on any affordable mesh and the archived curve itself has 1-3 meV spikes there. '
              'Lucky passes exist: the 150 mesh gives one sign change for some VBM offsets, and the 200 mesh with step '
              'occupation passes sigma_xy(2e21); the 100 mesh and the 200 mesh with 0.5-2 meV smearing fail.'),
    'calibration': {
        'sigma_xy_S_per_cm': {
            'reference_archive': TG['sigma_xy_S_per_cm'],
            'mu_below_vbm_eV': {k: r(v, 4) for k, v in CAL['mu_of_p_eV_below_vbm'].items()},
            'robustness': {'mesh_400_model': PM['400']['sigma_xy_at_p'], 'mesh_300_model': PM['300']['sigma_xy_at_p'],
                           'mesh_200_model': PM['200']['sigma_xy_at_p'],
                           'step_occupation_400': CAL['smearing_400']['0.0']['sigma_xy_at_p']},
            'nearest_ablation': {'spin_degeneracy_twice': A['spin_degeneracy_twice']['sigma_xy_S_per_cm'],
                                 'mesh_by_mesh': {n: PM[n]['sigma_xy_at_p'] for n in PM},
                                 'room_temperature_smearing_400': A['room_temperature_smearing']['sigma_xy_S_per_cm'],
                                 'not_ablation_wt_out_100K_smearing_400': A['wt_out_100K_smearing']['sigma_xy_S_per_cm']},
            'tolerance': '+-15%'},
        'sigma_xy_sign_changes': {
            'reference_archive': TG['sigma_xy_sign_changes'],
            'by_mesh': {n: PM[n]['sign_changes'] for n in PM},
            'crossings_below_vbm_by_mesh_eV': {n: PM[n]['crossings_eV_below_vbm'] for n in PM},
            'by_smearing_400': {k: v['sign_changes'] for k, v in CAL['smearing_400'].items()},
            'tolerance': 'exact'},
        'sign_change_hole_density_cm3': {
            'reference': TG['sign_change_hole_density_cm3'],
            'robustness': {n: PM[n]['sign_change_density'] for n in ('200', '300', '400')},
            'nearest_ablation': {'spin_degeneracy_twice': A['spin_degeneracy_twice']['sign_change_hole_density_cm3']},
            'tolerance': '+-8%'},
        'hole_density_0p30_below_vbm_cm3': {
            'reference': TG['hole_density_0p30_below_vbm_cm3'],
            'robustness': {'mesh_300': PM['300']['hole_density_0p30'], 'mesh_200': PM['200']['hole_density_0p30'],
                           'mesh_100': PM['100']['hole_density_0p30'], 'wt_in_efermi_as_vbm': R['wt_in_efermi_as_vbm']},
            'nearest_ablation': {'vbm_from_band_path': A['vbm_from_band_path'], 'vbm_from_dft_grid': A['vbm_from_dft_grid'],
                                 'vbm_from_outcar_efermi': A['vbm_from_outcar_efermi'],
                                 'spin_degeneracy_twice': A['spin_degeneracy_twice']['hole_density_0p30_below_vbm_cm3']},
            'tolerance': '+-5%'},
        'forbidden_components_max_abs_S_per_cm': {
            'reference': {'yz': 0, 'zx': 0}, 'archive_max_abs': '< 1e-7',
            'nearest_ablation': {'raw_400': A['raw_forbidden_components'],
                                 'raw_by_mesh': {n: [r(PM[n]['max_abs_yz'], 2), r(PM[n]['max_abs_zx'], 2)] for n in PM}},
            'tolerance': '<= 0.5 S/cm'},
        'allowed_hall_components': {'reference': TG['allowed_hall_components'],
                                    'nearest_ablation': {'symmetries_dat_literal': A['symmetries_dat_literal'],
                                                         'inplane_isotropic': A['inplane_isotropic']}, 'tolerance': 'exact'},
        'electron_doped_sigma_xy_S_per_cm': {'reference': None, 'nearest_ablation': {'from_model': A['electron_doped_from_model']},
                                             'evidence': CAL['wannier_vs_dft_on_grid_eV'], 'tolerance': 'must be null'},
        'energy_references_eV': {'vbm': CAL['vbm_eV'], 'vbm_k_frac': CAL['vbm_k_frac'], 'band_path_max': CAL['vbm_path_eV'],
                                 'dft_grid_max': CAL['vbm_dft_grid_eV'], 'wt_in_E_FERMI': CAL['wt_in_E_FERMI'],
                                 'outcar_E_fermi': CAL['outcar_E_fermi']}},
}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print('wrote design.json with', len(DPS), 'decision points')
