"""Write the task package files for superionic-heat-transport-size-onsager (run in the runtime image).

    docker run --rm --cpuset-cpus 6,7,22,23 --user $(id -u):$(id -g) -v $PWD:/c \
        -v ~/.cache/ai-theorist/superionic-heat-transport-size-onsager:/cache \
        -v ~/Desktop/workspace/AI-Theorist/benchmarks/tasks/superionic-heat-transport-size-onsager:/task \
        ai-theorist-benchmark-runtime:20261007 python3 -B /c/make_task.py

Inputs: /cache/pool.npz, /cache/pool_meta.json (extract.py), /cache/kappa_table.json (analysis.py).
Outputs: catalog.json, lab data (data.json + series.npz without the withheld runs), private targets,
hidden_values.json and solution/ablation_values.json; prints the calibration numbers.
"""
import json
import numpy as np

TASK = '/task'
BUDGET = 48
COST = {2: 1, 3: 2, 4: 4}
# Withheld states (4x4x4 cell, nominal temperature): every released run of the state is withheld.
PREDICT = {'UO2-4x4x4-1700K': ('UO2', 4, 1700), 'UO2-4x4x4-3400K': ('UO2', 4, 3400),
           'CaF2-4x4x4-1600K': ('CaF2', 4, 1600), 'PbF2-4x4x4-1200K': ('PbF2', 4, 1200)}
# Served runs whose kappa must be reported.
REPORT_RUNS = ['UO2-NVE-4x4x4-3500K', 'CaF2-NVE-4x4x4-1700K', 'PbF2-NVT-4x4x4-1300K']
TOL = {'prediction_rel': 0.12, 'run_rel': 0.12}
ANION = {'UO2': 'O', 'CaF2': 'F', 'PbF2': 'F'}
COLUMNS = {'UO2': ['c_myHeatFlux[1]', 'c_myHeatFlux[2]', 'c_myHeatFlux[3]', 'c_VOx', 'c_VOy', 'c_VOz'],
           'CaF2': ['c_myHeatFlux[1]', 'c_myHeatFlux[2]', 'c_myHeatFlux[3]', 'c_VFx', 'c_VFy', 'c_VFz'],
           'PbF2': ['c_myHeatFlux[1]', 'c_myHeatFlux[2]', 'c_myHeatFlux[3]', 'c_VFx', 'c_VFy', 'c_VFz']}


def r(x, n=6):
    return float('%.*g' % (n, x))


def main():
    pool = np.load('/cache/pool.npz')
    meta = json.load(open('/cache/pool_meta.json'))
    kt = json.load(open('/cache/kappa_table.json'))

    def kap(m, e, n, T, f='mc'):
        key = f'{m}-{e}-{n}x{n}x{n}-{T}K'
        return kt[key][f] if key in kt else None

    def state(m, n, T, f='mc', ens=('NVE', 'NVT')):
        v = [kap(m, e, n, T, f) for e in ens if kap(m, e, n, T, f) is not None]
        return float(np.mean(v)) if v else None

    withheld = {k for k, v in meta.items() for (m, n, T) in PREDICT.values()
                if v['material'] == m and v['supercell'] == n and v['temperature_K'] == T}
    served = [k for k in meta if k not in withheld]
    assert not set(REPORT_RUNS) & withheld

    # ---- lab data -------------------------------------------------------------------------------
    data = {}
    for k in served:
        v = meta[k]
        data[k] = {kk: v[kk] for kk in ('material', 'ensemble', 'supercell', 'temperature_K', 'atoms', 'Lx_A',
                                        'mean_temperature_K', 'std_temperature_K', 'mean_pressure_bar',
                                        'total_energy_drift_eV_per_ns', 'scale', 'deck')}
        data[k]['cost'] = COST[v['supercell']]
    with open(f'{TASK}/environment/lab/data.json', 'w') as fh:
        json.dump(data, fh, separators=(',', ':'))
    np.savez(f'{TASK}/environment/lab/series.npz', **{k: pool[k] for k in served})

    # ---- catalog --------------------------------------------------------------------------------
    runs = [{'run': k, 'material': meta[k]['material'], 'ensemble': meta[k]['ensemble'],
             'supercell': meta[k]['supercell'], 'temperature_K': meta[k]['temperature_K'],
             'atoms': meta[k]['atoms'], 'cost': COST[meta[k]['supercell']]}
            for k in sorted(served, key=lambda k: (meta[k]['material'], meta[k]['ensemble'], meta[k]['supercell'],
                                                    meta[k]['temperature_K']))]
    catalog = {
        'description': ('Released classical MD runs (LAMMPS, rigid-ion/EAM force fields, fixed volume) of three '
                        'fluorite-structure ionic conductors, recorded for heat-transport analysis. Each run is '
                        '0.8 ns of production after equilibration at the target temperature; NVE runs are '
                        'microcanonical, NVT runs keep a stochastic velocity-rescaling (CSVR) thermostat on during '
                        'production. The input deck of every run is returned with its data.'),
        'materials': {m: {'mobile_anion': ANION[m], 'atoms_per_conventional_cell': 12,
                          'series_columns': COLUMNS[m]} for m in ANION},
        'request': {'kind': 'md_run', 'run': '<run id from runs_available>'},
        'example': {'kind': 'md_run', 'run': 'PbF2-NVE-3x3x3-700K'},
        'costs': {'2x2x2 cell (96 atoms)': 1, '3x3x3 cell (324 atoms)': 2, '4x4x4 cell (768 atoms)': 4},
        'budget': BUDGET,
        'returned_data': {
            'series': ('Six time series sampled every frame_interval_ps = 0.08 ps (10000 frames, 800 ps), in the '
                       'order of series_columns: the three Cartesian components of the deck\'s compute '
                       'myHeatFlux and of the summed anion velocities, in LAMMPS metal units, exactly as '
                       'written to the thermo log except that the 8 fs thermo output has been band-limited to '
                       '6.25 THz (ideal low-pass) and resampled to 0.08 ps.'),
            'encoding': ('series[column] is base64 of little-endian int16 values; physical value = int16 * '
                         'scale[column]. In Python: np.frombuffer(base64.b64decode(s), "<i2") * scale.'),
            'metadata': ('cell edge Lx_A, time-averaged temperature, its standard deviation, mean pressure, '
                         'total-energy drift over the production run and the LAMMPS input deck')},
        'runs_available': runs,
        'predict': [{'id': k, 'material': m, 'supercell': n, 'temperature_K': T} for k, (m, n, T) in PREDICT.items()],
        'report_runs': REPORT_RUNS,
    }
    with open(f'{TASK}/environment/agent/catalog.json', 'w') as fh:
        json.dump(catalog, fh, indent=1)

    # ---- hidden targets -------------------------------------------------------------------------
    targets = {k: state(m, n, T) for k, (m, n, T) in PREDICT.items()}
    run_targets = {k: kt[k]['mc'] for k in REPORT_RUNS}
    pairs = [kt[k.replace('NVE', 'NVT')]['mc'] / kt[k]['mc'] for k in kt if '-NVE-' in k and k.replace('NVE', 'NVT') in kt]
    hidden = {'prediction_targets_W_mK': targets, 'run_targets_W_mK': run_targets,
              'thermostat_biases_kappa': False, 'thresholds': TOL,
              'reference': {'withheld_runs': {k: kt[k]['mc'] for k in sorted(withheld)},
                            'nvt_over_nve_pairs': {'n': len(pairs), 'mean': float(np.mean(pairs)),
                                                   'std': float(np.std(pairs))}}}
    with open(f'{TASK}/environment/lab/private/hidden.json', 'w') as fh:
        json.dump(hidden, fh, indent=1)
    hv = list(targets.values()) + list(run_targets.values()) + [kt[k]['mc'] for k in sorted(withheld)]
    with open(f'{TASK}/environment/lab/private/hidden_values.json', 'w') as fh:
        json.dump(hv, fh)

    # ---- ablations and robustness (from the served-format data) ---------------------------------
    def oracle_like(f='mc', ens=None):
        """The oracle's prediction scheme with another estimator / ensemble choice."""
        s = (lambda m, n, T: state(m, n, T, f, ens)) if ens else (lambda m, n, T: state(m, n, T, f))
        return {'UO2-4x4x4-1700K': state('UO2', 4, 1800, f, ens or ('NVE', 'NVT')),
                'UO2-4x4x4-3400K': (kap('UO2', 'NVE', 4, 3300, f) + kap('UO2', 'NVE', 4, 3500, f)) / 2,
                'CaF2-4x4x4-1600K': (kap('CaF2', 'NVE', 4, 1500, f) + kap('CaF2', 'NVE', 4, 1700, f)) / 2,
                'PbF2-4x4x4-1200K': (kap('PbF2', 'NVT', 4, 1100, f) + kap('PbF2', 'NVT', 4, 1300, f)) / 2}

    def at_cell(n, ens=('NVE', 'NVT')):
        return {k: state(m, n, T, 'mc', ens) for k, (m, n0, T) in PREDICT.items()}

    def inv_T(m, Ts, T):
        """Phonon-gas law kappa = A/T fitted to crystalline 4x4x4 runs, extrapolated to T."""
        A = np.mean([state(m, 4, t) * t for t in Ts])
        return A / T

    base = oracle_like()
    runs_plain = {k: kt[k]['plain'] for k in REPORT_RUNS}
    ab = {
        'plain_flux': {'pred': oracle_like('plain'), 'runs': runs_plain},
        'smallest_cells': {'pred': at_cell(2)},
        'factor_two': {'pred': {k: 2 * v for k, v in base.items()}, 'runs': {k: 2 * v for k, v in run_targets.items()}},
        'phonon_gas': {'pred': {**base, 'UO2-4x4x4-3400K': inv_T('UO2', (1800, 1900), 3400),
                                'CaF2-4x4x4-1600K': inv_T('CaF2', (1000, 1100), 1600),
                                'PbF2-4x4x4-1200K': inv_T('PbF2', (500, 600), 1200)}},
        'cell3_converged': {'pred': {**base, 'UO2-4x4x4-1700K': state('UO2', 3, 1700)}},
        'thermostat_bias': {'pred': {**oracle_like(ens=('NVE',)), 'PbF2-4x4x4-1200K': kap('PbF2', 'NVE', 3, 1200)},
                            'bias': True},
    }
    rob = {
        'lowpass_3THz': {'pred': oracle_like('mc3'), 'runs': {k: kt[k]['mc3'] for k in REPORT_RUNS}},
        'nve_only': {'pred': oracle_like(ens=('NVE',))},
        'single_neighbour': {'pred': {'UO2-4x4x4-1700K': kap('UO2', 'NVT', 4, 1800),
                                      'UO2-4x4x4-3400K': kap('UO2', 'NVE', 4, 3300),
                                      'CaF2-4x4x4-1600K': kap('CaF2', 'NVE', 4, 1500),
                                      'PbF2-4x4x4-1200K': kap('PbF2', 'NVT', 4, 1300)}},
        'linear_extrapolation_both_ensembles': {'pred': {'UO2-4x4x4-1700K': 2 * state('UO2', 4, 1800) - state('UO2', 4, 1900),
                                               'UO2-4x4x4-3400K': (state('UO2', 4, 3300) + state('UO2', 4, 3500)) / 2,
                                               'CaF2-4x4x4-1600K': (state('CaF2', 4, 1500) + state('CaF2', 4, 1700)) / 2,
                                               'PbF2-4x4x4-1200K': base['PbF2-4x4x4-1200K']}},
    }
    with open(f'{TASK}/solution/ablation_values.json', 'w') as fh:
        json.dump({'ablations': ab, 'robustness': rob}, fh, indent=1)

    # ---- calibration print ----------------------------------------------------------------------
    print('targets', {k: round(v, 4) for k, v in targets.items()})
    print('run targets', {k: round(v, 4) for k, v in run_targets.items()})
    print('NVT/NVE pairs', len(pairs), round(float(np.mean(pairs)), 4), round(float(np.std(pairs)), 4))

    def rel(d, ref):
        return {k: round(v / ref[k] - 1, 3) for k, v in d.items()}
    print('oracle-like', rel(base, targets))
    for name, v in {**ab, **rob}.items():
        print(name, 'pred', rel(v['pred'], targets), 'runs', rel(v['runs'], run_targets) if 'runs' in v else '')
    print('cell3 all', rel(at_cell(3), targets))
    print('served', len(served), 'withheld', sorted(withheld))
    print('catalog total cost', sum(x['cost'] for x in runs))


if __name__ == '__main__':
    main()
