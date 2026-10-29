"""Write design.json (decision points, calibration table) from reference.json and calibration.json.

    python3 -B write_design.py
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/mpea-vacancy-diffusion-kmc'
cal = json.loads((HERE/'calibration.json').read_text())
ref = json.loads((HERE/'reference.json').read_text())['targets']
priv = json.loads((TASK/'environment/lab/private/hidden.json').read_text())
hb, tol = priv['hidden_barriers_eV'], priv['tolerances']
V, A = cal['variants'], cal['ablations']
r3 = lambda x: round(x, 3)


def rmse(p):
    return r3(math.sqrt(sum((p[k] - v) ** 2 for k, v in hb.items()) / len(hb)))


def row(field, getter, ablations, tolerance, reference):
    rob = {n: r3(getter(v)) for n, v in V.items()}
    return {'reference': r3(reference), 'robustness': rob,
            'robustness_max_deviation': r3(max(abs(x - reference) for x in rob.values())) if field != 'barrier_rmse_eV' else None,
            'nearest_ablations': ablations, 'tolerance': tolerance}


calibration = {
    'barrier_rmse_eV (100 hidden hops)': {
        'reference_full_data_model': rmse(V['all_purchasable_hops']['barrier_predictions_eV']),
        'robustness': {n: rmse(v['barrier_predictions_eV']) for n, v in V.items()},
        'nearest_ablations': {n: rmse(A[n]['barrier_predictions_eV']) for n in
                              ('species_mean_barriers', 'direct_forward_barrier_fit', 'thermo_ideal_as_real')},
        'tolerance': f"<= {tol['barrier_rmse_eV']}"}}
for c in 'ABC':
    calibration[f'log10_eta_vacancy:{c}'] = row('eta', lambda v, c=c: v['log10_eta_vacancy'][c],
        {n: r3(A[n]['log10_eta_vacancy'][c]) for n in ('thermo_ideal_as_real', 'species_mean_barriers', 'rom_from_linear_model_pure_limit')},
        f"+-{tol['log10_eta']}", ref['log10_eta_vacancy'][c])
calibration['log10_tracer_ratio_B'] = row('ratio', lambda v: v['log10_tracer_ratio_B'],
    {n: r3(A[n]['log10_tracer_ratio_B']) for n in ('jump_count_tracers', 'species_mean_barriers', 'thermo_ideal_as_real')},
    f"+-{tol['log10_tracer_ratio_B']}", ref['log10_tracer_ratio_B'])
calibration['activation_energy_C_eV'] = row('ea', lambda v: v['activation_energy_C_eV'],
    {n: r3(A[n]['activation_energy_C_eV']) for n in ('ea_from_rom_barrier', 'ea_from_mean_accepted_barrier', 'thermo_ideal_as_real', 'species_mean_barriers')},
    f"+-{tol['activation_energy_eV']} eV", ref['activation_energy_C_eV'])
for k in ('log10_eta_thermo_ideal_B', 'log10_eta_kra_simple_B'):
    calibration[k] = row(k, lambda v, k=k: v[k], {'untested_attribution': r3(A['untested_attribution'][k])}, f"+-{tol[k]}", ref[k])
calibration['fastest_tracer_B'] = {'reference': ref['fastest_tracer_B'], 'robustness': {n: v['fastest_tracer_B'] for n, v in V.items()},
                                   'tolerance': 'exact'}
calibration['information_only (not graded)'] = {
    'random_configurations_log10_eta': {c: r3(x) for c, x in cal['information_only']['random_configurations']['log10_eta_vacancy'].items()},
    'small_cell_L4_short_runs_log10_eta': {c: r3(x) for c, x in cal['information_only']['small_cell_short_runs']['log10_eta_vacancy'].items()},
    'pure_limit_of_linear_kra_model_eV': {e: r3(x) for e, x in cal['pure_limit_linear_model_eV'].items()},
    'pure_element_kra_dft_eV': {e: r3(x) for e, x in cal['pure_kra_dft_eV'].items()}}

design = {
    'task_id': 'mpea-vacancy-diffusion-kmc', 'catalog_id': 'mpea-vacancy-diffusion-kmc', 'start': 'released',
    'source_group': 'lee_2026_mpea_diffusion_elce',
    'title': 'Vacancy diffusion in refractory multi-principal-element alloys from released hop barriers',
    'decision_points': [
        {'id': 'dp1', 'kind': 'competing explanation / model class (detailed balance, discovered site-energy disorder)',
         'fork': 'Whether hop barriers carry the end-state energy difference (KRA + dE/2 with an end-state Hamiltonian) or are modelled '
                 'as direction-independent / directly fitted forward barriers',
         'wrong_choice': 'Use the KRA alone (equal forward and backward barriers: the thermo-ideal model) or fit E_TS - E_initial directly '
                         'on symmetric pair descriptors, assuming (as the paper reports) that thermodynamics is minor',
         'consequence': f"Archived |E2 - E1| has median 0.43 eV; hidden-barrier RMSE {rmse(A['thermo_ideal_as_real']['barrier_predictions_eV'])} eV "
                        f"(direct fit {rmse(A['direct_forward_barrier_fit']['barrier_predictions_eV'])} eV); log10 eta of A off by "
                        f"{r3(A['thermo_ideal_as_real']['log10_eta_vacancy']['A'] - ref['log10_eta_vacancy']['A'])}",
         'ablation': 'thermo_ideal_as_real'},
        {'id': 'dp2', 'kind': 'descriptor range / model class',
         'fork': 'Local-environment barrier model vs composition-averaged (species-mean) barriers',
         'wrong_choice': 'One barrier per hopping element (alloy mean of the purchased KRA), with the end-state Hamiltonian',
         'consequence': f"Hidden-barrier RMSE {rmse(A['species_mean_barriers']['barrier_predictions_eV'])} eV; tracer ratio of B "
                        f"{r3(A['species_mean_barriers']['log10_tracer_ratio_B'])} vs {r3(ref['log10_tracer_ratio_B'])}; "
                        f"log10 eta of A {r3(A['species_mean_barriers']['log10_eta_vacancy']['A'])} vs {r3(ref['log10_eta_vacancy']['A'])}",
         'ablation': 'species_mean_barriers'},
        {'id': 'dp3', 'kind': 'estimator (tracer vs jump frequency; correlation)',
         'fork': 'Tracer diffusivities from mean-squared displacements vs from hop counts',
         'wrong_choice': 'D*_i from the number of hops times the hop length squared (uncorrelated random walk)',
         'consequence': f"log10 tracer ratio of B {r3(A['jump_count_tracers']['log10_tracer_ratio_B'])} vs {r3(ref['log10_tracer_ratio_B'])}: "
                        'the fast species (Nb, f ~ 0.11) is far more correlated than the slow ones (f ~ 0.35-0.4)',
         'ablation': 'jump_count_tracers'},
        {'id': 'dp4', 'kind': 'evidence discovered by validation (model extrapolation vs released data); cascade into eta',
         'fork': 'Source of the pure-element migration barriers in the rule-of-mixtures benchmark',
         'wrong_choice': 'Take E_m from the pure-element limit of a fitted linear local-environment model instead of the released pure-element hops',
         'consequence': 'The linear model misses the pure limits by up to +0.33 eV (Nb 0.70 vs 0.37 eV); log10 eta shifts by '
                        f"+{r3(A['rom_from_linear_model_pure_limit']['log10_eta_vacancy']['A'] - ref['log10_eta_vacancy']['A'])} (A) and "
                        f"+{r3(A['rom_from_linear_model_pure_limit']['log10_eta_vacancy']['C'] - ref['log10_eta_vacancy']['C'])} (C)",
         'ablation': 'rom_from_linear_model_pure_limit'},
        {'id': 'dp5', 'kind': 'regime judgement (effective vs microscopic activation energy)',
         'fork': 'Activation energy from the temperature dependence of the simulated tracer diffusivity vs from a barrier average',
         'wrong_choice': 'Report the rule-of-mixtures barrier, or the mean barrier of accepted hops, as the activation energy',
         'consequence': f"{r3(A['ea_from_rom_barrier']['activation_energy_C_eV'])} or {r3(A['ea_from_mean_accepted_barrier']['activation_energy_C_eV'])} eV "
                        f"instead of {r3(ref['activation_energy_C_eV'])} eV: trapping and correlated slow paths add ~0.35 eV",
         'ablation': 'ea_from_rom_barrier'},
        {'id': 'dp6', 'kind': 'discriminating calculation (mechanism attribution)',
         'fork': 'Whether thermodynamic (site-energy) or kinetic (barrier-distribution) disorder controls sluggishness, tested by the variants',
         'wrong_choice': 'Attribute without running the variants (variants assumed to behave like the real alloy)',
         'consequence': f"Thermo-ideal B is {r3(ref['log10_eta_thermo_ideal_B'])} and KRA-simple B {r3(ref['log10_eta_kra_simple_B'])} against "
                        f"{r3(ref['log10_eta_vacancy']['B'])} for the real alloy: both matter, the barrier distribution most",
         'ablation': 'untested_attribution'}],
    'graded_fields': [
        f"barrier_predictions_eV: RMSE over 100 hidden archived hops <= {tol['barrier_rmse_eV']} eV",
        f"log10_eta_vacancy A, B, C (2000 K) within +-{tol['log10_eta']}",
        'fastest_tracer_B exact',
        f"log10_tracer_ratio_B within +-{tol['log10_tracer_ratio_B']}",
        f"activation_energy_C_eV (fastest tracer of C, 1700-2300 K) within +-{tol['activation_energy_eV']} eV",
        f"log10_eta_thermo_ideal_B and log10_eta_kra_simple_B within +-{tol['log10_eta_thermo_ideal_B']}"],
    'length_estimate': {},
    'runtime_dependencies': ['numpy (runtime image); the lab compiles its kMC kernel (kmc.c, ~230 lines of C) with the image gcc at build '
                             'time; no extra packages'],
    'data_provenance': {'archive_url': 'https://archive.materialscloud.org/records/bhv3r-a8q53', 'doi': '10.24435/materialscloud:5g-95',
                        'files_used': ['all_data.json', 'files_description.md'],
                        'md5': {'all_data.json': '15df9136ed668320057131c8c59cad2a', 'files_description.md': '8849580fa4e1bb9f72bca8d1e74e571f'},
                        'license': 'CC BY 4.0'},
    'exceptions': {},
    'notes': 'Diffusion targets are the builders\' reconstruction (reference barrier model trained on all 1,799 clean hops, lab kMC engine), '
             'not values from the paper; the authors\' separate DFT energy cluster expansion is not archived, so the end-state energetics '
             'come from the released hop energies and the thermodynamic contribution is larger than in the paper. The SRO-equilibration and '
             'finite-size forks of the design are not graded: alloys where they matter (Cr-Nb-Mo-Ta near its miscibility gap) were too '
             'sensitive to the purchased subset (log10 eta from -2.0 to -3.4 across 800-hop random subsets); they remain rubric items.',
    'calibration': calibration}
design['length_estimate'] = {
    'oracle_lab_evaluations': 47, 'oracle_charge': 54, 'budget': 60,
    'catalog_size': '1,845 purchasable hops (93 units at 20 per unit) plus unlimited kMC runs (1-8 units each)',
    'oracle_plan': '40 hop purchases (800 hops: 35 pure/dilute + 765 random alloy hops) + 7 kMC runs of 2 units',
    'oracle_container_seconds': 35, 'harbor_oracle_wall_minutes': 1.3,
    'tool_calls': '50-80', 'wall_time_h': '1.5-3', 'gpt6astra_cost_usd': '3-5', 'pilot': 'pending (coordinator)'}
(TASK/'design.json').write_text(json.dumps(design, indent=1) + '\n')
print(json.dumps(calibration, indent=1)[:3000])
