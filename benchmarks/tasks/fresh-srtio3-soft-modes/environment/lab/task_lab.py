"""Lab logic for fresh-srtio3-soft-modes: new Quantum ESPRESSO 6.7 calculations for perovskite SrTiO3.

validate(request) checks a structured request against public.json (every field is required, no defaults) and
returns its charge. For `scf` the charge is the canonical size of the request: the number of symmetry-distinct k
points of the mesh (counted here from the space group of the displaced cell, as pw.x reduces the mesh) times a
per-cell cost and a cutoff factor. For `dielectric` (undisplaced five-atom cell only) the charge comes from a table.
run(request, workdir) writes the pw.x input, runs `mpirun -np CPUS pw.x` (and ph.x for `dielectric`) with one thread per
rank, and returns the raw physical results (energy, forces, stress; electronic dielectric tensor and raw Born
charges). pw.in/pw.out (and ph.in/ph.out) stay in the job folder; wavefunctions are deleted.
Nothing hidden is ever served: grade() compares the terminal submission with reference values computed at
construction time through this same run() and stored in private/ (benchmarks/construction/<task id>/).
"""
import itertools
import json
import math
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PUBLIC = json.loads((HERE/'public.json').read_text())
SPEC = PUBLIC['request']
CPUS = int(os.environ.get('LAB_MPI_PROCS', '2'))
BOHR_A = 0.529177210903
RY_EV = 13.605693122994
FIELDS = ('kind', 'cell', 'lattice_A', 'c_A', 'displacements_A', 'ecutwfc_Ry', 'kmesh')
# five-atom cell (fractions of a, a, c): Sr, Ti, O_x (between Ti along x), O_y, O_z
BASIS5 = [('Sr', (0.0, 0.0, 0.0)), ('Ti', (0.5, 0.5, 0.5)), ('O', (0.0, 0.5, 0.5)), ('O', (0.5, 0.0, 0.5)),
          ('O', (0.5, 0.5, 0.0))]
LABELS5 = ['Sr', 'Ti', 'O_x', 'O_y', 'O_z']


# ---------------------------------------------------------------- geometry and symmetry
def geometry(cell_kind, a, c, disp):
    """Cell (rows, Angstrom), species, Cartesian positions (Angstrom) and crystal positions.

    perovskite_5: a1 = (a, 0, 0), a2 = (0, a, 0), a3 = (0, 0, c); atoms Sr, Ti, O_x, O_y, O_z at the ideal sites.
    rotated_10:   a1 = (0, a, c), a2 = (a, 0, c), a3 = (a, a, 0) (the R-point doubling, two formula units);
                  formula unit 1 at the origin, formula unit 2 the same five sites shifted by (a, 0, 0).
    The displacements (Cartesian, Angstrom) are added to the ideal sites in this atom order.
    """
    scale = np.array([a, a, c])
    if cell_kind == 'perovskite_5':
        cell = np.diag(scale)
        shifts = [np.zeros(3)]
    else:
        cell = np.array([[0.0, a, c], [a, 0.0, c], [a, a, 0.0]])
        shifts = [np.zeros(3), np.array([a, 0.0, 0.0])]
    species, cart = [], []
    for s in shifts:
        for el, f in BASIS5:
            species.append(el)
            cart.append(np.array(f) * scale + s)
    cart = np.array(cart) + np.array(disp, float)
    frac = np.linalg.solve(cell.T, cart.T).T
    return cell, species, cart, frac


def atom_labels(cell_kind):
    if cell_kind == 'perovskite_5':
        return list(LABELS5)
    return [lab + '_1' for lab in LABELS5] + [lab + '_2' for lab in LABELS5]


def _lattice_ops(cell, tol=1e-6):
    G = cell @ cell.T
    ops = []
    for entries in itertools.product((-1, 0, 1), repeat=9):
        M = np.array(entries, float).reshape(3, 3)
        if abs(abs(np.linalg.det(M)) - 1) > 1e-9:
            continue
        if np.abs(M.T @ G @ M - G).max() < tol * np.abs(G).max():
            ops.append(M)
    return ops


def _maps(M, t, frac, species, tol):
    for i, x in enumerate(frac):
        y = M @ x + t
        if not any(species[j] == species[i] and np.abs((y - frac[j]) - np.round(y - frac[j])).max() < tol
                   for j in range(len(frac))):
            return False
    return True


def crystal_ops(cell, frac, species, tol=1e-5):
    """Point operations (crystal basis) of the crystal. As in pw.x, fractional translations are allowed unless
    the cell is a supercell (a pure translation maps the crystal onto itself), in which case only operations
    without translation are kept."""
    lattice = _lattice_ops(cell)
    same = [j for j in range(len(frac)) if species[j] == species[0]]
    eye = np.eye(3)
    supercell = any(_maps(eye, frac[j] - frac[0], frac, species, tol) for j in same[1:]
                    if np.abs((frac[j] - frac[0]) - np.round(frac[j] - frac[0])).max() > tol)
    ops = []
    for M in lattice:
        cands = [np.zeros(3)] if supercell else [frac[j] - M @ frac[0] for j in same]
        if any(_maps(M, t, frac, species, tol) for t in cands):
            ops.append(M)
    return ops


def irreducible_k(cell, frac, species, n):
    """Number of k points of the n x n x n Gamma-centred mesh after reduction by the crystal point group and
    time reversal (orbits of the mesh)."""
    ops = crystal_ops(cell, frac, species)
    kops = [np.round(np.linalg.inv(M).T).astype(int) for M in ops]
    seen, orbits = set(), 0
    for idx in itertools.product(range(n), repeat=3):
        k = np.array(idx)
        if tuple(k % n) in seen:
            continue
        orbits += 1
        for M in kops:
            kk = M @ k
            seen.add(tuple(kk % n))
            seen.add(tuple((-kk) % n))
    return orbits, len(ops)


# ---------------------------------------------------------------- request checking and charge
def _number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def normalize(request):
    """Return the canonical request, or raise ValueError with a usable message (invalid requests are free)."""
    if not isinstance(request, dict):
        raise ValueError('JSON object required; see GET /task for the request schema')
    unknown, missing = set(request) - set(FIELDS), [f for f in FIELDS if f not in request]
    if unknown:
        raise ValueError(f'unknown fields {sorted(unknown)}; allowed: {list(FIELDS)}')
    if missing:
        raise ValueError(f'missing fields {missing}: every field is required (no defaults)')
    allowed = SPEC['allowed']
    out = {}
    for key in ('kind', 'cell', 'ecutwfc_Ry'):
        value = request[key]
        if isinstance(value, bool) or value not in allowed[key]:
            raise ValueError(f'{key} must be one of {allowed[key]}')
        out[key] = value
    kms = allowed['kmesh'][out['cell']]
    if isinstance(request['kmesh'], bool) or request['kmesh'] not in kms:
        raise ValueError(f"kmesh for cell {out['cell']} must be one of {kms}")
    out['kmesh'] = request['kmesh']
    for key in ('lattice_A', 'c_A'):
        lo, hi = SPEC[key + '_range']
        if not _number(request[key]) or not lo <= request[key] <= hi:
            raise ValueError(f'{key} must be a number in [{lo}, {hi}] (Angstrom)')
        out[key] = round(float(request[key]), 5)
    if abs(out['c_A'] / out['lattice_A'] - 1) > SPEC['max_abs_c_over_a_minus_1'] + 1e-12:
        raise ValueError(f"c_A / lattice_A must lie within 1 +- {SPEC['max_abs_c_over_a_minus_1']}")
    nat = 5 if out['cell'] == 'perovskite_5' else 10
    d = request['displacements_A']
    if (not isinstance(d, list) or len(d) != nat or not all(isinstance(r, list) and len(r) == 3 for r in d)
            or not all(_number(v) for r in d for v in r)):
        raise ValueError(f"displacements_A must be a list of {nat} Cartesian [dx, dy, dz] (Angstrom), one per atom "
                         f"in the order {atom_labels(out['cell'])}")
    d = np.round(np.array(d, float), 5) + 0.0
    if np.abs(d).max() > SPEC['max_abs_displacement_A']:
        raise ValueError(f"each displacement component must lie within +-{SPEC['max_abs_displacement_A']} A")
    out['displacements_A'] = d.tolist()
    if out['kind'] == 'dielectric':
        if out['cell'] != 'perovskite_5' or np.abs(d).max() > 0:
            raise ValueError('kind "dielectric" is available for the undisplaced perovskite_5 cell only (all displacements 0)')
        table = SPEC['charge']['dielectric']
        if str(out['kmesh']) not in table:
            raise ValueError(f'kmesh for kind "dielectric" must be one of {sorted(int(k) for k in table)}')
    return out


def size_and_charge(req):
    cell, species, _, frac = geometry(req['cell'], req['lattice_A'], req['c_A'], req['displacements_A'])
    nk, nsym = irreducible_k(cell, frac, species, req['kmesh'])
    cost = SPEC['charge']
    if req['kind'] == 'dielectric':
        units = cost['dielectric'][str(req['kmesh'])] * (nk / SPEC['charge']['dielectric_cubic_k_points'][str(req['kmesh'])])
        units = units * cost['cutoff_factor'][str(req['ecutwfc_Ry'])]
    else:
        units = nk * cost['scf_k_point_weight'][req['cell']] * cost['cutoff_factor'][str(req['ecutwfc_Ry'])]
    return nk, nsym, max(1, math.ceil(units - 1e-9))


def validate(request):
    req = normalize(request)
    nk, nsym, charge = size_and_charge(req)
    if charge > SPEC['charge']['max_units_per_request']:
        raise ValueError(f"request too large: {charge} units (maximum {SPEC['charge']['max_units_per_request']} "
                         f"per request); use a smaller mesh or cutoff, or a more symmetric displacement pattern")
    return charge


# ---------------------------------------------------------------- pw.x / ph.x input, execution, parsing
def pw_input(req, cell, species, frac):
    fx = SPEC['fixed']
    nat = len(species)
    lines = ['&CONTROL', "calculation='scf', prefix='sto', outdir='./work', pseudo_dir='%s'" % HERE,
             "tprnfor=.true., tstress=.true., disk_io='low'", '/',
             '&SYSTEM', f'ibrav=0, nat={nat}, ntyp=3',
             f"ecutwfc={req['ecutwfc_Ry']}, ecutrho={fx['dual'] * req['ecutwfc_Ry']}, occupations='fixed'", '/',
             '&ELECTRONS', f"conv_thr={fx['conv_thr_Ry_per_atom'] * nat:.1e}, mixing_beta={fx['mixing_beta']}, "
             f"electron_maxstep={fx['electron_maxstep']}", '/', 'ATOMIC_SPECIES']
    for el in ('Sr', 'Ti', 'O'):
        lines.append(f"{el} {fx['mass_amu'][el]} {SPEC['pseudopotentials'][el]}")
    lines += ['CELL_PARAMETERS angstrom'] + ['  '.join(f'{v:.12f}' for v in row) for row in cell]
    lines += ['ATOMIC_POSITIONS crystal'] + [f'{el} ' + '  '.join(f'{v:.12f}' for v in row)
                                             for el, row in zip(species, frac)]
    n = req['kmesh']
    lines += ['K_POINTS automatic', f'{n} {n} {n} 0 0 0', '']
    return '\n'.join(lines)


def ph_input(req):
    return ('dielectric response\n&INPUTPH\n'
            f" prefix='sto', outdir='./work', tr2_ph={SPEC['fixed']['tr2_ph']}, epsil=.true., trans=.false.\n"
            '/\n0.0 0.0 0.0\n')


def _run(args, folder, infile, outfile, timeout):
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               OMPI_MCA_btl_vader_single_copy_mechanism='none')
    with (folder/outfile).open('w') as out, (folder/(outfile + '.err')).open('w') as err:
        try:
            proc = subprocess.run(args + ['-in', infile], cwd=folder, stdout=out, stderr=err, env=env, timeout=timeout)
        except subprocess.TimeoutExpired:
            raise ValueError(f'{outfile}: exceeded the time limit; the job is charged') from None
    text = (folder/outfile).read_text(errors='ignore')
    if proc.returncode or 'JOB DONE.' not in text:
        raise ValueError(f'{outfile}: the program did not finish; the job is charged')
    return text


def _mpi(program):
    return ['mpirun', '--allow-run-as-root', '--oversubscribe', '-np', str(CPUS), program, '-npool', str(CPUS)]


def parse_scf(text, nat):
    if 'convergence NOT achieved' in text or 'convergence has been achieved' not in text:
        raise ValueError('SCF did not converge; the job is charged (see pw.out)')
    energy = float(re.findall(r'^!\s+total energy\s+=\s+([-\d.]+)\s+Ry', text, re.M)[-1])
    block = text[text.rfind('total   stress'):].split('\n')
    stress_kbar = [[float(x) for x in block[i].split()[3:6]] for i in (1, 2, 3)]
    chunk = text.split('Forces acting on atoms')[-1]
    forces = [[float(x) for x in m.groups()] for m in
              re.finditer(r'atom\s+\d+ type\s+\d+\s+force =\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', chunk)][:nat]
    nk = int(re.findall(r'number of k points=\s+(\d+)', text)[-1])
    nsym = re.findall(r'(\d+) Sym\. Ops\.', text)
    homo = re.findall(r'highest occupied level \(ev\):\s+([-\d.]+)', text)
    return {'total_energy_Ry': energy,
            'stress_GPa': [[round(v / 10.0, 6) for v in row] for row in stress_kbar],
            'pressure_GPa': round(sum(stress_kbar[i][i] for i in range(3)) / 30.0, 6),
            'forces_Ry_per_bohr': forces,
            'pw_irreducible_k_points': nk, 'pw_symmetry_operations': int(nsym[-1]) if nsym else None,
            'highest_occupied_level_eV': float(homo[-1]) if homo else None,
            'scf_iterations': int(re.findall(r'convergence has been achieved in\s+(\d+)\s+iterations', text)[-1])}


def parse_ph(text, nat):
    """Electronic dielectric tensor and Born effective charges (dF/dE, no sum rule) printed by ph.x."""
    lines = text.split('\n')
    k = max(i for i, line in enumerate(lines) if 'Dielectric constant in cartesian axis' in line)
    eps = [[float(x) for x in lines[k + 2 + r].replace('(', ' ').replace(')', ' ').split()] for r in range(3)]
    k = max(i for i, line in enumerate(lines) if 'Effective charges (d Force / dE) in cartesian axis' in line)
    z = []
    for atom in range(nat):
        start = k + 2 + atom * 4
        z.append([[float(x) for x in lines[start + 1 + r].replace('(', ' ').replace(')', ' ').split()[1:]]
                  for r in range(3)])
    return eps, z


def run(request, workdir):
    req = normalize(request)
    nk_expected, nsym_expected, _ = size_and_charge(req)
    folder = Path(workdir)
    folder.mkdir(parents=True, exist_ok=True)
    timeout = SPEC['fixed']['timeout_s']
    cell, species, cart, frac = geometry(req['cell'], req['lattice_A'], req['c_A'], req['displacements_A'])
    nat = len(species)
    start = time.monotonic()
    res = {}
    try:
        (folder/'pw.in').write_text(pw_input(req, cell, species, frac))
        res.update(parse_scf(_run(_mpi('pw.x'), folder, 'pw.in', 'pw.out', timeout), nat))
        res['pw_seconds'] = round(time.monotonic() - start, 1)
        if req['kind'] == 'dielectric':
            t1 = time.monotonic()
            (folder/'ph.in').write_text(ph_input(req))
            eps, z = parse_ph(_run(_mpi('ph.x'), folder, 'ph.in', 'ph.out', timeout), nat)
            res['dielectric_tensor_eps_inf'] = eps
            res['born_charges_raw'] = {lab: zz for lab, zz in zip(atom_labels(req['cell']), z)}
            res['ph_seconds'] = round(time.monotonic() - t1, 1)
    finally:      # keep inputs and outputs as evidence, never the wavefunctions or scratch
        shutil.rmtree(folder/'work', ignore_errors=True)
        for p in folder.glob('_ph*'):
            shutil.rmtree(p, ignore_errors=True)
    res.update({'atoms': atom_labels(req['cell']), 'species': species,
                'cell_vectors_A': [[round(v, 10) for v in row] for row in cell],
                'cell_volume_A3': round(float(abs(np.linalg.det(cell))), 8),
                'positions_cartesian_A': [[round(v, 8) for v in row] for row in cart],
                'atomic_masses_amu': SPEC['fixed']['mass_amu'],
                'settings': req, 'charged_k_points': nk_expected, 'charged_symmetry_operations': nsym_expected,
                'wall_seconds': round(time.monotonic() - start, 1)})
    return res


# ---------------------------------------------------------------- grading
NUMERIC = ('a0_A', 'gamma_soft_mode_cm-1', 'gamma_LO_highest_cm-1', 'eps_inf', 'born_charge_Ti',
           'polar_well_depth_meV_per_fu', 'R_mode_cm-1', 'tilt_angle_deg', 'tilt_well_depth_meV_per_fu',
           'strained_c_A', 'strained_inplane_polar_mode_cm-1')
CHOICES = {'larger_energy_gain': ('polar', 'tilt')}


def grade(payload, trace, validation):
    private = json.loads((HERE/'private/hidden.json').read_text())
    ref, tol = private['reference'], private['tolerances']
    missing = [f for f in NUMERIC + tuple(CHOICES) + ('evidence_ids',) if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    bad = [f for f in NUMERIC if not _number(payload[f])]
    if bad:
        return False, {'error': 'fields must be finite numbers: %s' % bad}
    bad = [f for f, options in CHOICES.items() if payload[f] not in options]
    if bad:
        return False, {'error': '%s must be one of %s' % (bad, [CHOICES[f] for f in bad])}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    for f in NUMERIC:
        err = payload[f] - ref[f]
        if f in tol['relative']:
            err = err / abs(ref[f])
            checks[f] = abs(err) <= tol['relative'][f]
        else:
            checks[f] = abs(err) <= tol['absolute'][f]
        errors[f] = err
    for f in CHOICES:
        checks[f] = payload[f] == ref[f]
    return all(checks.values()), {'checks': checks, 'errors': errors,
                                  'note': 'relative errors for %s, absolute otherwise' % sorted(tol['relative'])}
