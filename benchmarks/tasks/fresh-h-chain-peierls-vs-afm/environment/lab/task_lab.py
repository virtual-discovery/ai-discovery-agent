"""Lab logic for fresh-h-chain-peierls-vs-afm: new Quantum ESPRESSO SCFs for a periodic hydrogen chain.

validate(request) checks a structured SCF request against the ranges in public.json and returns its
charge; run(request, workdir) writes the pw.x input, runs `mpirun -np CPUS pw.x` in the job directory,
parses the output and returns physical results. Nothing hidden is ever served: the graded values come
from construction-time reference calculations made with this same run() and stored in private/.
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
PSEUDO = SPEC['fixed']['pseudopotential']
CPUS = int(os.environ.get('LAB_MPI_PROCS', '2'))
RY_EV = 13.605693122994
HA_EV = 2 * RY_EV
FIELDS = {'kind', 'bond_lengths_A', 'starting_magnetization', 'k_points', 'smearing', 'degauss_Ry', 'box_A',
          'ecutwfc_Ry'}


def _number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def normalize(request):
    """Return the request with defaults filled in, or raise ValueError with a usable message."""
    if not isinstance(request, dict) or request.get('kind') != 'scf':
        raise ValueError('Use {"kind": "scf", "bond_lengths_A": [...], ...}; see GET /task for the schema')
    unknown = set(request) - FIELDS
    if unknown:
        raise ValueError(f'Unknown fields {sorted(unknown)}; allowed: {sorted(FIELDS)}')
    bonds = request.get('bond_lengths_A')
    if not isinstance(bonds, list) or len(bonds) not in SPEC['atoms_per_cell'] or not all(_number(b) for b in bonds):
        raise ValueError(f'bond_lengths_A must list one bond length (A) per atom in the cell; allowed cell sizes '
                         f'{SPEC["atoms_per_cell"]} atoms')
    lo, hi = SPEC['bond_length_range_A']
    if not all(lo <= b <= hi for b in bonds):
        raise ValueError(f'each bond length must lie in [{lo}, {hi}] A')
    mean = sum(bonds) / len(bonds)
    mlo, mhi = SPEC['mean_spacing_range_A']
    if not mlo - 1e-9 <= mean <= mhi + 1e-9:
        raise ValueError(f'the mean spacing (cell length / atoms) must lie in [{mlo}, {mhi}] A')
    out = {'kind': 'scf', 'bond_lengths_A': [round(float(b), 6) for b in bonds]}
    mag = request.get('starting_magnetization')
    if mag is not None:
        if not isinstance(mag, list) or len(mag) != len(bonds) or not all(_number(m) and -1 <= m <= 1 for m in mag):
            raise ValueError('starting_magnetization must be null (no spin polarization) or one value in [-1, 1] per atom')
        mag = [round(float(m), 4) for m in mag]
    out['starting_magnetization'] = mag
    for key in ('k_points', 'smearing', 'degauss_Ry', 'box_A', 'ecutwfc_Ry'):
        if key not in request:
            raise ValueError(f'{key} is required (no defaults); allowed values {SPEC["allowed"][key]}')
        value = request[key]
        if isinstance(value, bool) or value not in SPEC['allowed'][key]:
            raise ValueError(f'{key} must be one of {SPEC["allowed"][key]}')
        out[key] = value
    if out['k_points'] * len(bonds) > SPEC['max_k_points_times_atoms']:
        raise ValueError(f'k_points x atoms_per_cell may not exceed {SPEC["max_k_points_times_atoms"]}')
    if out['k_points'] * sum(bonds) > SPEC['max_k_points_times_cell_length_A'] + 1e-6:
        raise ValueError(f'k_points x cell length (A) may not exceed {SPEC["max_k_points_times_cell_length_A"]}')
    return out


def charge_of(req):
    """Charge in budget units: proportional to (atoms per cell) x (k points), one unit per 96."""
    return max(1, math.ceil(len(req['bond_lengths_A']) * req['k_points'] / SPEC['charge_unit_k_points_times_atoms']))


def validate(request):
    return charge_of(normalize(request))


def write_input(req, folder):
    bonds, mag = req['bond_lengths_A'], req['starting_magnetization']
    n, length, box = len(bonds), sum(bonds), req['box_A']
    xs = [sum(bonds[:i]) for i in range(n)]
    if mag is None:
        species, labels = ['H'], ['H'] * n
    else:   # one species per distinct starting magnetization
        values = sorted(set(mag), reverse=True)
        species = [f'H{i + 1}' for i in range(len(values))]
        labels = [species[values.index(m)] for m in mag]
    system = [f'ibrav=0, nat={n}, ntyp={len(species)}',
              f"ecutwfc={req['ecutwfc_Ry']}, ecutrho={4 * req['ecutwfc_Ry']}",
              f"occupations='smearing', smearing='{req['smearing']}', degauss={req['degauss_Ry']}"]
    if mag is None:
        system.append('nspin=1')
    else:
        system.append('nspin=2')
        system += [f'starting_magnetization({i + 1})={v}' for i, v in enumerate(sorted(set(mag), reverse=True))]
    lines = ['&CONTROL', "calculation='scf', prefix='hchain', outdir='./work'", f"pseudo_dir='{HERE}'",
             "tprnfor=.true., disk_io='low'", '/', '&SYSTEM', *system, '/',
             '&ELECTRONS', f"conv_thr={SPEC['fixed']['conv_thr_Ry_per_atom'] * n:.1e}, electron_maxstep=300",
             f"mixing_beta={SPEC['fixed']['mixing_beta']}", '/',
             'ATOMIC_SPECIES', *[f'{s} 1.00794 {PSEUDO}' for s in species],
             'CELL_PARAMETERS angstrom', f'{length:.8f} 0 0', f'0 {box:.8f} 0', f'0 0 {box:.8f}',
             'ATOMIC_POSITIONS angstrom', *[f'{labels[i]} {xs[i]:.8f} 0 0' for i in range(n)],
             'K_POINTS automatic', f"{req['k_points']} 1 1 0 0 0", '']
    (folder/'pw.in').write_text('\n'.join(lines))


def band_gap(bands, root):
    """Band gap (eV) on the k mesh when every spin channel holds an integer number of electrons, else 0 (metal)."""
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
    return round(max(0.0, min(bottom) - max(top)), 5)


def parse(folder, req):
    text = (folder/'pw.out').read_text()
    n = len(req['bond_lengths_A'])
    energies = re.findall(r'^!\s+total energy\s+=\s+([-\d.]+)\s+Ry', text, re.M)
    if not energies or 'JOB DONE.' not in text or 'convergence has been achieved' not in text:
        raise ValueError('SCF did not converge within the iteration limit; the job is charged')
    iterations = int(re.findall(r'convergence has been achieved in\s+(\d+)\s+iterations', text)[-1])
    energy_ry = float(energies[-1])
    root = ET.parse(folder/'work/hchain.save/data-file-schema.xml').getroot()
    if abs(2 * float(root.findtext('./output/total_energy/etot')) - energy_ry) > 1e-6:
        raise RuntimeError('XML and printed energies disagree')
    bands = root.find('./output/band_structure')
    fermi = float(bands.findtext('fermi_energy')) * HA_EV
    gap = band_gap(bands, root)
    forces = [float(m.group(1)) * RY_EV / 0.529177210903 for m in
              re.finditer(r'atom\s+\d+ type\s+\d+\s+force =\s+([-\d.]+)\s+[-\d.]+\s+[-\d.]+', text)][-n:]
    result = {'atoms_per_cell': n, 'cell_length_A': round(sum(req['bond_lengths_A']), 6),
              'mean_spacing_A': round(sum(req['bond_lengths_A']) / n, 6),
              'total_energy_Ry_per_cell': energy_ry,
              'fermi_energy_eV': fermi,
              'band_gap_eV': gap,
              'forces_eV_per_A_along_chain': forces,
              'scf_iterations': iterations}
    if req['starting_magnetization'] is not None:
        block = text[text.rfind('Magnetic moment per site'):]
        sites = [float(v) for v in re.findall(r'magn:\s+([-\d.]+)', block)][:n]
        result.update(total_magnetization_muB_per_cell=round(float(root.findtext('./output/magnetization/total')), 5),
                      absolute_magnetization_muB_per_cell=round(float(root.findtext('./output/magnetization/absolute')), 5),
                      site_moments_muB=sites)
    return result


def run(request, workdir):
    req = normalize(request)
    folder = Path(workdir)
    folder.mkdir(parents=True, exist_ok=True)
    write_input(req, folder)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               OMPI_MCA_btl_vader_single_copy_mechanism='none')
    npool = CPUS if req['k_points'] >= 2 * CPUS else 1
    start = time.monotonic()
    with (folder/'pw.out').open('w') as out, (folder/'pw.err').open('w') as err:
        try:
            proc = subprocess.run(['mpirun', '--allow-run-as-root', '--oversubscribe', '-np', str(CPUS), 'pw.x',
                                   '-npool', str(npool), '-in', 'pw.in'], cwd=folder, stdout=out, stderr=err,
                                  env=env, timeout=SPEC['fixed']['timeout_s'])
        except subprocess.TimeoutExpired:
            raise ValueError('pw.x exceeded the time limit; the job is charged') from None
    if proc.returncode:
        raise ValueError('pw.x exited with an error; the job is charged')
    result = parse(folder, req)
    result['wall_seconds'] = round(time.monotonic() - start, 1)
    result['settings'] = {k: req[k] for k in ('k_points', 'smearing', 'degauss_Ry', 'box_A', 'ecutwfc_Ry')}
    result['spin_polarized'] = req['starting_magnetization'] is not None
    shutil.rmtree(folder/'work', ignore_errors=True)   # keep pw.in/pw.out as evidence, drop wavefunctions
    return result


def finite(x):
    return _number(x)


NUMERIC = {'delta_eq': ('1.25', '1.75'), 'dimerization_energy_meV_per_atom': ('1.25', '1.75'),
           'uniform_magnetic_energy_gain_meV_per_atom': ('1.75', '2.0'), 'critical_dimerization': ('1.75', '2.0'),
           'uniform_chain_band_gap_eV': ('1.25',)}
LABELS = ('uniform_nonmagnetic', 'uniform_ferromagnetic', 'uniform_antiferromagnetic', 'dimerized_nonmagnetic',
          'dimerized_magnetic')


def grade(payload, trace, validation):
    private = json.loads((HERE/'private/hidden.json').read_text())
    ref, tol = private['reference'], private['tolerances']
    missing = [f for f in [*NUMERIC, 'ground_state', 'evidence_ids'] if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    for field, keys in NUMERIC.items():
        value = payload[field]
        if not isinstance(value, dict) or set(value) != set(keys) or not all(finite(v) for v in value.values()):
            return False, {'error': f'{field} must map exactly the keys {list(keys)} to finite numbers'}
    states = payload['ground_state']
    if not isinstance(states, dict) or set(states) != {'1.0', '2.0'} or not all(v in LABELS for v in states.values()):
        return False, {'error': f'ground_state must map "1.0" and "2.0" to one of {list(LABELS)}'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    for field, keys in NUMERIC.items():
        for key in keys:
            got, want = payload[field][key], ref[field][key]
            if field in tol['relative']:
                err = got / want - 1
                ok = abs(err) <= tol['relative'][field]
            else:
                err = got - want
                ok = abs(err) <= tol['absolute'][field]
            checks[f'{field}[{key}]'] = ok
            errors[f'{field}[{key}]'] = err
    for key in ('1.0', '2.0'):
        checks[f'ground_state[{key}]'] = states[key] == ref['ground_state'][key]
    return all(checks.values()), {'checks': checks, 'errors': errors}
