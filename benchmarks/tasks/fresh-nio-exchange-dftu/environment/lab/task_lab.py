"""Lab logic for fresh-nio-exchange-dftu: new Quantum ESPRESSO 6.7 PBE+U SCF calculations for rocksalt NiO.

validate(request) checks a structured SCF request against the schema in public.json (no defaults) and returns
its charge; run(request, workdir) writes the pw.x input, runs `mpirun -np CPUS pw.x` (threads pinned to 1) in the
job folder, parses pw.out and the XML and returns the physical results; pw.in and pw.out stay in the job folder
as evidence. Nothing hidden is ever served: the graded values come from construction-time reference
calculations made with this same run() and stored in private/.
"""
import json
import math
import os
import re
import shutil
import subprocess
import time
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).parent
PUBLIC = json.loads((HERE/'public.json').read_text())
SPEC = PUBLIC['scf_request']
FIXED = SPEC['fixed']
CPUS = int(os.environ.get('LAB_MPI_PROCS', '2'))
HA_EV = 27.211386245988
NI_UPF, O_UPF = FIXED['pseudopotentials']['Ni'], FIXED['pseudopotentials']['O']
NI_MASS, O_MASS = 58.6934, 15.999
FIELDS = ['kind', 'cell', 'starting_magnetization', 'hubbard_U_eV', 'kpoints', 'kshift', 'ecutwfc_Ry', 'smearing',
          'degauss_Ry']
CELLS = {name: (c['vectors_in_units_of_a'], c['Ni_crystal'], c['O_crystal']) for name, c in SPEC['cells'].items()}


def _number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def normalize(request):
    """Return the canonical request, or raise ValueError with a usable message. Every field is required."""
    if not isinstance(request, dict) or request.get('kind') != 'scf':
        raise ValueError('Use {"kind": "scf", ...}; GET /task (python3 /app/lab_client.py task) lists the schema')
    unknown = set(request) - set(FIELDS)
    if unknown:
        raise ValueError(f'Unknown fields {sorted(unknown)}; allowed: {FIELDS}')
    missing = [f for f in FIELDS if f not in request]
    if missing:
        raise ValueError(f'Missing fields {missing}; every field is required (no defaults)')
    cell = request['cell']
    if cell not in CELLS:
        raise ValueError(f'cell must be one of {sorted(CELLS)}')
    n_ni = len(CELLS[cell][1])
    mag = request['starting_magnetization']
    if mag is not None:
        if not isinstance(mag, list) or len(mag) != n_ni or not all(_number(m) and -1 <= m <= 1 for m in mag):
            raise ValueError(f'starting_magnetization must be null (non-spin-polarized) or a list of {n_ni} values in '
                             f'[-1, 1], one per Ni of cell {cell} in the order of its Ni_crystal list')
        mag = [round(float(m), 3) for m in mag]
        if all(m == 0 for m in mag):
            raise ValueError('starting_magnetization: give at least one nonzero value, or null for a '
                             'non-spin-polarized calculation')
    U = request['hubbard_U_eV']
    lo, hi = SPEC['hubbard_U_eV_range']
    if not _number(U) or not lo <= U <= hi or abs(round(U, 1) - U) > 1e-9:
        raise ValueError(f'hubbard_U_eV must be a number in [{lo}, {hi}] eV in steps of 0.1')
    k = request['kpoints']
    kmax = SPEC['max_kpoints_per_direction'][cell]
    if not isinstance(k, list) or len(k) != 3 or not all(type(n) is int and 1 <= n <= kmax for n in k):
        raise ValueError(f'kpoints must be [n1, n2, n3], integers from 1 to {kmax} for cell {cell} '
                         '(Monkhorst-Pack mesh along the three cell vectors)')
    for key in ('kshift', 'ecutwfc_Ry', 'smearing', 'degauss_Ry'):
        value = request[key]
        if isinstance(value, bool) or value not in SPEC['allowed'][key]:
            raise ValueError(f'{key} must be one of {SPEC["allowed"][key]}')
    return {'kind': 'scf', 'cell': cell, 'starting_magnetization': mag, 'hubbard_U_eV': round(float(U), 1),
            'kpoints': list(k), 'kshift': request['kshift'], 'ecutwfc_Ry': request['ecutwfc_Ry'],
            'smearing': request['smearing'], 'degauss_Ry': request['degauss_Ry']}


def charge_of(req):
    """Units: ceil(w_cell x n1 n2 n3 / 1000 x (ecutwfc / 50)^1.5), at least 1 (about one unit per minute)."""
    rule = SPEC['charge']
    nk = req['kpoints'][0] * req['kpoints'][1] * req['kpoints'][2]
    units = rule['cell_weight'][req['cell']] * nk / 1000 * rule['ecut_factor'][str(req['ecutwfc_Ry'])]
    return max(1, math.ceil(units - 1e-9))


def validate(request):
    return charge_of(normalize(request))


def write_input(req, folder):
    vectors, ni, ox = CELLS[req['cell']]
    a = FIXED['lattice_parameter_A']
    mag = req['starting_magnetization']
    n_fu = len(ni)
    if mag is None:
        values, labels = [0.0], ['Ni1'] * n_fu
    else:   # one Ni species per distinct starting magnetization
        values = sorted(set(mag), reverse=True)
        labels = ['Ni%d' % (values.index(m) + 1) for m in mag]
    species = ['Ni%d' % (i + 1) for i in range(len(values))]
    system = [f'ibrav=0, nat={2 * n_fu}, ntyp={len(species) + 1}',
              f"ecutwfc={req['ecutwfc_Ry']}, ecutrho={FIXED['ecutrho_over_ecutwfc'] * req['ecutwfc_Ry']}",
              f"nbnd={FIXED['bands_per_formula_unit'] * n_fu + 4}",
              f"occupations='smearing', smearing='{req['smearing']}', degauss={req['degauss_Ry']}"]
    if mag is None:
        system.append('nspin=1')
    else:
        system.append('nspin=2')
        system += [f'starting_magnetization({i + 1})={v}' for i, v in enumerate(values)]
    if req['hubbard_U_eV'] > 0:
        system += ['lda_plus_u=.true., lda_plus_u_kind=0', "U_projection_type='atomic'"]
        system += [f"Hubbard_U({i + 1})={req['hubbard_U_eV']}" for i in range(len(species))]
    lines = ['&CONTROL', "calculation='scf', prefix='nio', outdir='./work'", f"pseudo_dir='{HERE}'",
             "disk_io='low'", '/', '&SYSTEM', *system, '/',
             '&ELECTRONS', f"conv_thr={FIXED['conv_thr_Ry_per_formula_unit'] * n_fu:.1e}, electron_maxstep=200",
             f"mixing_beta={FIXED['mixing_beta']}", '/',
             'ATOMIC_SPECIES', *[f'{s} {NI_MASS} {NI_UPF}' for s in species], f'O {O_MASS} {O_UPF}',
             'CELL_PARAMETERS angstrom', *['%.10f %.10f %.10f' % tuple(a * x for x in v) for v in vectors],
             'ATOMIC_POSITIONS crystal', *['%s %.6f %.6f %.6f' % (labels[i], *p) for i, p in enumerate(ni)],
             *['O %.6f %.6f %.6f' % tuple(p) for p in ox],
             'K_POINTS automatic', '%d %d %d %d %d %d' % (*req['kpoints'], *([req['kshift']] * 3)), '']
    (folder/'pw.in').write_text('\n'.join(lines))


def band_gap(root):
    """Gap (eV) on the k mesh when every spin channel holds an integer number of electrons, else 0 (metal)."""
    bands = root.find('./output/band_structure')
    nelec = float(bands.findtext('nelec'))
    eig = [[float(v) * HA_EV for v in k.findtext('eigenvalues').split()] for k in bands.findall('ks_energies')]
    if bands.findtext('lsda').strip() == 'true':
        nup = int(bands.findtext('nbnd_up'))
        magnet = float(root.findtext('./output/magnetization/total'))
        channels = [([e[:nup] for e in eig], (nelec + magnet) / 2), ([e[nup:] for e in eig], (nelec - magnet) / 2)]
    else:
        channels = [(eig, nelec / 2)]
    top, bottom = [], []
    for levels, filling in channels:
        nocc = round(filling)
        if abs(filling - nocc) > 0.02:
            return 0.0
        if nocc > 0:
            top.append(max(e[nocc - 1] for e in levels))
        bottom.append(min(e[nocc] for e in levels))
    return round(max(0.0, min(bottom) - max(top)), 4)


def hubbard_occupations(text, n_ni):
    """Ni 3d occupation-matrix traces and eigenvalues per spin, from the last write_ns block."""
    block = text[text.rfind('--- enter write_ns ---'):text.rfind('--- exit write_ns ---')]
    out = []
    for part in re.split(r'\natom\s+\d+\s+Tr\[ns\(na\)\]', block)[1:n_ni + 1]:
        traces = [float(x) for x in re.match(r'\s*\(up, down, total\)\s*=\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)',
                                             part).groups()]
        eigen = [[float(x) for x in m.split()] for m in re.findall(r'eigenvalues:\s*\n\s*([-\d. ]+)\n', part)]
        out.append({'trace_up': traces[0], 'trace_down': traces[1], 'trace_total': traces[2],
                    'eigenvalues_up': eigen[0], 'eigenvalues_down': eigen[1]})
    return out


def parse(folder, req):
    text = (folder/'pw.out').read_text()
    vectors, ni, ox = CELLS[req['cell']]
    n_fu = len(ni)
    energies = re.findall(r'^!\s+total energy\s+=\s+([-\d.]+)\s+Ry', text, re.M)
    if not energies or 'JOB DONE.' not in text or 'convergence has been achieved' not in text:
        raise ValueError('SCF did not converge within 200 iterations; the job is charged')
    root = ET.parse(folder/'work/nio.save/data-file-schema.xml').getroot()
    free = float(energies[-1])
    if abs(2 * float(root.findtext('./output/total_energy/etot')) - free) > 1e-6:
        raise RuntimeError('XML and printed energies disagree')
    ts = float(re.findall(r'smearing contrib\. \(-TS\)\s+=\s+([-\d.]+)\s+Ry', text)[-1])
    result = {'cell': req['cell'], 'formula_units_per_cell': n_fu,
              'free_energy_Ry_per_cell': free,
              'smearing_contribution_minus_TS_Ry_per_cell': ts,
              'internal_energy_Ry_per_cell': round(free - ts, 8),
              'fermi_energy_eV': round(float(root.findtext('./output/band_structure/fermi_energy')) * HA_EV, 4),
              'band_gap_eV': band_gap(root),
              'scf_iterations': int(re.findall(r'convergence has been achieved in\s+(\d+)\s+iterations', text)[-1]),
              'irreducible_k_points': int(re.findall(r'number of k points=\s+(\d+)', text)[-1])}
    if req['starting_magnetization'] is not None:
        block = text[text.rfind('Magnetic moment per site'):]
        sites = [float(v) for v in re.findall(r'magn:\s+([-\d.]+)', block)][:2 * n_fu]
        result.update(Ni_moments_muB=sites[:n_fu], O_moments_muB=sites[n_fu:],
                      total_magnetization_muB_per_cell=round(float(root.findtext('./output/magnetization/total')), 5),
                      absolute_magnetization_muB_per_cell=round(float(root.findtext('./output/magnetization/absolute')), 5))
    if req['hubbard_U_eV'] > 0:
        result['Ni_3d_occupations'] = hubbard_occupations(text, n_fu)
    return result


def run(request, workdir):
    req = normalize(request)
    folder = Path(workdir)
    folder.mkdir(parents=True, exist_ok=True)
    write_input(req, folder)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               OMPI_MCA_btl_vader_single_copy_mechanism='none')
    start = time.monotonic()
    try:
        with (folder/'pw.out').open('w') as out, (folder/'pw.err').open('w') as err:
            try:
                proc = subprocess.run(['mpirun', '--allow-run-as-root', '--oversubscribe', '-np', str(CPUS), 'pw.x',
                                       '-npool', str(CPUS), '-in', 'pw.in'], cwd=folder, stdout=out, stderr=err,
                                      env=env, timeout=FIXED['timeout_s'])
            except subprocess.TimeoutExpired:
                raise ValueError('pw.x exceeded the time limit; the job is charged') from None
        if 'convergence NOT achieved' in (folder/'pw.out').read_text(errors='ignore'):
            raise ValueError('SCF did not converge within 200 iterations; the job is charged')
        if proc.returncode:
            raise ValueError('pw.x exited with an error; the job is charged')
        result = parse(folder, req)
    finally:   # keep pw.in/pw.out as evidence, always drop wavefunctions and scratch
        shutil.rmtree(folder/'work', ignore_errors=True)
    result['wall_seconds'] = round(time.monotonic() - start, 1)
    result['request'] = req
    return result


def finite(x):
    return _number(x)


KEYED = {'J1_meV': ('4.0', '8.0'), 'J2_meV': ('4.0', '8.0'), 'neel_temperature_mean_field_K': ('4.0', '8.0')}
SCALAR = ('pbe_fm_minus_afm2_meV_per_Ni', 'fm_gap_closing_U_eV')


def grade(payload, trace, validation):
    private = json.loads((HERE/'private/hidden.json').read_text())
    ref, tol = private['reference'], private['tolerances']
    missing = [f for f in [*KEYED, *SCALAR, 'evidence_ids'] if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    for field, keys in KEYED.items():
        value = payload[field]
        if not isinstance(value, dict) or set(value) != set(keys) or not all(finite(v) for v in value.values()):
            return False, {'error': f'{field} must map exactly the keys {list(keys)} to finite numbers'}
    bad = [f for f in SCALAR if not finite(payload[f])]
    if bad:
        return False, {'error': f'fields must be finite numbers: {bad}'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    items = [(f'{f}[{k}]', payload[f][k], ref[f][k], f) for f, keys in KEYED.items() for k in keys]
    items += [(f, payload[f], ref[f], f) for f in SCALAR]
    for name, got, want, field in items:
        if field in tol['relative']:
            err = got / want - 1
            ok = abs(err) <= tol['relative'][field]
        else:
            err = got - want
            ok = abs(err) <= tol['absolute'][field]
        checks[name], errors[name] = ok, err
    note = 'relative errors for %s, absolute errors (field units) otherwise' % sorted(tol['relative'])
    return all(checks.values()), {'checks': checks, 'errors': errors, 'note': note}
