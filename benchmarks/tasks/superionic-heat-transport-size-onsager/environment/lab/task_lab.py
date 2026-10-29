"""Lab logic for superionic-heat-transport-size-onsager: retrieval of author-released MD flux series and grading.

A purchase returns one released run: its LAMMPS deck, thermodynamic metadata and the heat-flux and anion-velocity
time series (band-limited to 6.25 THz, 0.08 ps frames). The runs of the states listed under `predict` are not in
this lab at all; their Green-Kubo conductivities, and those of the `report_runs`, grade the submission.
"""
import base64
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).parent
DATA = json.loads((HERE/'data.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
SERIES = np.load(HERE/'series.npz')
PREDICT = set(PRIVATE['prediction_targets_W_mK'])
USAGE = 'Use {"kind": "md_run", "run": "<run id from runs_available in /app/catalog.json>"}'
SOURCE = 'retrieved author-released MD run (F. Grasselli 2022, Materials Cloud yzhmw-vje16); no new simulation'


def validate(request):
    if request.get('kind') != 'md_run' or set(request) != {'kind', 'run'} or not isinstance(request['run'], str):
        raise ValueError(USAGE)
    run = request['run']
    if run not in DATA:
        parts = run.split('-')
        if len(parts) == 4 and f'{parts[0]}-{parts[2]}-{parts[3]}' in PREDICT:
            raise ValueError('This state is withheld: it is one of the states you must predict')
        raise ValueError('Unknown run. ' + USAGE)
    return DATA[run]['cost']


def run(request, workdir):
    d = DATA[request['run']]
    q = np.ascontiguousarray(SERIES[request['run']]).astype('<i2')
    columns = {'UO2': ['c_myHeatFlux[1]', 'c_myHeatFlux[2]', 'c_myHeatFlux[3]', 'c_VOx', 'c_VOy', 'c_VOz']}.get(
        d['material'], ['c_myHeatFlux[1]', 'c_myHeatFlux[2]', 'c_myHeatFlux[3]', 'c_VFx', 'c_VFy', 'c_VFz'])
    return {'run': request['run'], 'material': d['material'], 'ensemble': d['ensemble'], 'supercell': d['supercell'],
            'nominal_temperature_K': d['temperature_K'], 'atoms': d['atoms'], 'Lx_A': d['Lx_A'],
            'mean_temperature_K': round(d['mean_temperature_K'], 2), 'std_temperature_K': round(d['std_temperature_K'], 2),
            'mean_pressure_bar': round(d['mean_pressure_bar'], 1),
            'total_energy_drift_eV_per_ns': round(d['total_energy_drift_eV_per_ns'], 4),
            'frame_interval_ps': 0.08, 'n_frames': int(q.shape[0]), 'series_columns': columns,
            'scale': d['scale'], 'series': {c: base64.b64encode(q[:, i].tobytes()).decode() for i, c in enumerate(columns)},
            'input_deck': d['deck'], 'source': SOURCE}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def grade(payload, trace, validation):
    t = PRIVATE['thresholds']
    fields = ['kappa_predictions_W_mK', 'kappa_runs_W_mK', 'thermostat_biases_kappa', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    pt, rt = PRIVATE['prediction_targets_W_mK'], PRIVATE['run_targets_W_mK']
    pp, rp = payload['kappa_predictions_W_mK'], payload['kappa_runs_W_mK']
    if not isinstance(pp, dict) or set(pp) != set(pt) or not all(finite(v) and v > 0 for v in pp.values()):
        return False, {'error': 'kappa_predictions_W_mK must give a positive finite kappa for every id in predict'}
    if not isinstance(rp, dict) or set(rp) != set(rt) or not all(finite(v) and v > 0 for v in rp.values()):
        return False, {'error': 'kappa_runs_W_mK must give a positive finite kappa for every run in report_runs'}
    if not isinstance(payload['thermostat_biases_kappa'], bool):
        return False, {'error': 'thermostat_biases_kappa must be true or false'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    checks = {'evidence_executed': evidence}
    report = {'relative_error': {}}
    for k, ref in pt.items():
        e = pp[k] / ref - 1
        report['relative_error'][k] = e
        checks[f'predict_{k}'] = abs(e) <= t['prediction_rel']
    for k, ref in rt.items():
        e = rp[k] / ref - 1
        report['relative_error'][k] = e
        checks[f'run_{k}'] = abs(e) <= t['run_rel']
    checks['thermostat_verdict'] = payload['thermostat_biases_kappa'] == PRIVATE['thermostat_biases_kappa']
    report['checks'] = checks
    return all(checks.values()), report
