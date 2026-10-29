"""Lab logic for fresh-ge-internal-strain: new Quantum ESPRESSO 6.7 pw.x calculations for diamond Ge.

validate(request) checks a structured request against public.json and returns its charge, computed from the
canonical size of the request (the number of symmetry-distinct k points pw.x will use, counted here exactly as
pw.x counts them, times a cost factor for the pseudopotential/cutoff and for structural relaxation).
run(request, workdir) writes the pw.x input, runs `mpirun -np CPUS pw.x` (threads pinned to 1), parses the
output and returns raw physical results; pw.in/pw.out are kept as evidence. Nothing hidden is ever served:
grade() compares the terminal submission with reference values computed at construction time through this
same run() and stored in private/ (see benchmarks/construction/fresh-ge-internal-strain/).
"""
import itertools
import json
import math
import os
import re
import shutil
import subprocess
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PUBLIC = json.loads((HERE/'public.json').read_text())
SPEC = PUBLIC['request']
CPUS = int(os.environ.get('LAB_MPI_PROCS', '2'))
BOHR_A = 0.529177210903
RY_EV = 13.605693122994
FIELDS = ('kind', 'pseudopotential', 'lattice_A', 'strain', 'displacement_A', 'ecutwfc_Ry', 'kmesh', 'kshift',
          'occupations', 'degauss_Ry', 'forc_conv_thr_Ry_per_bohr')


# ---------------------------------------------------------------- symmetry and canonical request size
def _cubic_ops():
    ops = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            R = np.zeros((3, 3))
            for i, (p, s) in enumerate(zip(perm, signs)):
                R[i, p] = s
            ops.append(R)
    return ops


OPS = _cubic_ops()
# operations of the site group Td (signed permutations with an even number of sign changes) keep each sublattice;
# the other 24 operations of the diamond space group exchange the two sublattices
TD = [round(float(np.prod(R.sum(axis=1)))) == 1 for R in OPS]


def symmetry(eps, u, tol=1e-9):
    """Point operations of the strained lattice and of the strained, displaced crystal (Cartesian matrices)."""
    lattice = [R for R in OPS if np.abs(R @ eps @ R.T - eps).max() < tol]
    crystal = [R for R, td in zip(OPS, TD) if np.abs(R @ eps @ R.T - eps).max() < tol
               and np.abs(R @ u - (u if td else -u)).max() < tol]
    return lattice, crystal


def cell_and_positions(a, eps, u):
    F = np.eye(3) + eps
    cell = 0.5 * a * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], float) @ F.T      # rows a1, a2, a3 (Angstrom)
    cart = np.array([-u / 2, F @ (0.25 * a * np.ones(3)) + u / 2])
    return cell, cart, np.linalg.solve(cell.T, cart.T).T


def irreducible_k(cell, eps, u, n, shift):
    """Number of k points pw.x keeps for an n x n x n Monkhorst-Pack mesh (kpoint_grid + irreducible_BZ)."""
    B = np.linalg.inv(cell).T
    lattice, crystal = symmetry(eps, u)
    frac = lambda R: np.round(np.linalg.inv(B.T) @ R @ B.T)
    Lf, Gf = [frac(R) for R in lattice], [frac(R) for R in crystal]
    grid = (np.array(list(itertools.product(range(n), repeat=3)), float) + shift / 2) / n
    key = lambda k: tuple(np.round((k - np.floor(k + 1e-9)) * 2 * n).astype(int) % (2 * n))
    index = {key(k): i for i, k in enumerate(grid)}
    parent = list(range(len(grid)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, k in enumerate(grid):
        for M in Lf:
            for kk in (M @ k, -(M @ k)):
                j = index.get(key(kk))
                if j is not None:
                    x, y = find(i), find(j)
                    if x != y:
                        parent[max(x, y)] = min(x, y)
    reps = sorted({find(i) for i in range(len(grid))})
    if len(crystal) == len(lattice):
        return len(reps), len(crystal)
    total = 0
    for r in reps:
        star = {}
        for M in Lf:
            star.setdefault(key(M @ grid[r]), M @ grid[r])
        seen = set()
        for kkey, kk in star.items():
            if kkey not in seen:
                total += 1
                for M in Gf:
                    seen.update((key(M @ kk), key(-(M @ kk))))
    return total, len(crystal)


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
        raise ValueError(f'missing fields {missing}: every field is required (no defaults; use null where allowed)')
    allowed = SPEC['allowed']
    out = {}
    for key in ('kind', 'pseudopotential', 'ecutwfc_Ry', 'kmesh', 'kshift', 'occupations'):
        value = request[key]
        if isinstance(value, bool) or value not in allowed[key]:
            raise ValueError(f'{key} must be one of {allowed[key]}')
        out[key] = value
    lo, hi = SPEC['lattice_A_range']
    if not _number(request['lattice_A']) or not lo <= request['lattice_A'] <= hi:
        raise ValueError(f'lattice_A must be a number in [{lo}, {hi}]')
    out['lattice_A'] = round(float(request['lattice_A']), 6)
    eps = request['strain']
    if (not isinstance(eps, list) or len(eps) != 3 or not all(isinstance(r, list) and len(r) == 3 for r in eps)
            or not all(_number(v) for r in eps for v in r)):
        raise ValueError('strain must be a 3x3 list of numbers (symmetric strain tensor, Cartesian cube axes)')
    eps = np.round(np.array(eps, float), 6)
    if np.abs(eps - eps.T).max() > 1e-9:
        raise ValueError('strain must be symmetric (eps_ij = eps_ji)')
    if np.abs(eps).max() > SPEC['max_abs_strain_component']:
        raise ValueError(f"each strain component must lie within +-{SPEC['max_abs_strain_component']}")
    if np.any((np.abs(eps) > 0) & (np.abs(eps) < SPEC['min_nonzero_strain_component'])):
        raise ValueError(f"nonzero strain components must have magnitude >= {SPEC['min_nonzero_strain_component']}")
    out['strain'] = eps.tolist()
    u = request['displacement_A']
    if not isinstance(u, list) or len(u) != 3 or not all(_number(v) for v in u):
        raise ValueError('displacement_A must be a list of three numbers (Angstrom)')
    u = np.round(np.array(u, float), 6)
    if np.abs(u).max() > SPEC['max_abs_displacement_component_A']:
        raise ValueError(f"each displacement component must lie within +-{SPEC['max_abs_displacement_component_A']} A")
    if np.any((np.abs(u) > 0) & (np.abs(u) < SPEC['min_nonzero_displacement_component_A'])):
        raise ValueError(f"nonzero displacement components must have magnitude >= {SPEC['min_nonzero_displacement_component_A']} A")
    out['displacement_A'] = u.tolist()
    if out['occupations'] == 'fixed':
        if request['degauss_Ry'] is not None:
            raise ValueError('degauss_Ry must be null for fixed occupations')
        out['degauss_Ry'] = None
    else:
        if isinstance(request['degauss_Ry'], bool) or request['degauss_Ry'] not in allowed['degauss_Ry']:
            raise ValueError(f"degauss_Ry must be one of {allowed['degauss_Ry']} for smearing")
        out['degauss_Ry'] = request['degauss_Ry']
    thr = request['forc_conv_thr_Ry_per_bohr']
    if out['kind'] == 'scf':
        if thr is not None:
            raise ValueError('forc_conv_thr_Ry_per_bohr must be null for kind scf')
        out['forc_conv_thr_Ry_per_bohr'] = None
    else:
        if isinstance(thr, bool) or thr not in allowed['forc_conv_thr_Ry_per_bohr']:
            raise ValueError(f"forc_conv_thr_Ry_per_bohr must be one of {allowed['forc_conv_thr_Ry_per_bohr']} for relax")
        out['forc_conv_thr_Ry_per_bohr'] = thr
    return out


def size_and_charge(req):
    eps, u = np.array(req['strain']), np.array(req['displacement_A'])
    cell, _, _ = cell_and_positions(req['lattice_A'], eps, u)
    nk, nsym = irreducible_k(cell, eps, u, req['kmesh'], req['kshift'])
    cost = SPEC['charge']
    units = nk * cost['factor'][req['pseudopotential']][str(req['ecutwfc_Ry'])] / cost['k_points_per_unit']
    if req['kind'] == 'relax':
        units *= cost['relax_multiplier']
    return nk, nsym, max(1, math.ceil(units - 1e-9))


def validate(request):
    req = normalize(request)
    nk, nsym, charge = size_and_charge(req)
    if charge > SPEC['charge']['max_units_per_request']:
        raise ValueError(f"request too large: {nk} symmetry-distinct k points -> {charge} units; the maximum per "
                         f"request is {SPEC['charge']['max_units_per_request']} (use a smaller mesh or cutoff, or a "
                         f"more symmetric cell)")
    return charge


# ---------------------------------------------------------------- pw.x input, execution, parsing
def write_input(req, folder):
    upf, zval, mass = req['pseudopotential'], SPEC['pseudopotentials'][req['pseudopotential']]['z_valence'], SPEC['fixed']['mass_amu']
    eps, u = np.array(req['strain']), np.array(req['displacement_A'])
    cell, _, frac = cell_and_positions(req['lattice_A'], eps, u)
    nbnd = zval + SPEC['fixed']['empty_bands']          # 2 atoms x zval electrons -> zval occupied bands
    system = ['ibrav=0, nat=2, ntyp=1', f"ecutwfc={req['ecutwfc_Ry']}, ecutrho={SPEC['fixed']['dual'] * req['ecutwfc_Ry']}",
              f'nbnd={nbnd}']
    if req['occupations'] == 'fixed':
        system.append("occupations='fixed'")
    else:
        system.append(f"occupations='smearing', smearing='mv', degauss={req['degauss_Ry']}")
    control = [f"calculation='{req['kind']}', prefix='ge', outdir='./work'", f"pseudo_dir='{HERE}'",
               "tprnfor=.true., tstress=.true., disk_io='low'"]
    if req['kind'] == 'relax':
        control.append(f"forc_conv_thr={req['forc_conv_thr_Ry_per_bohr']}, "
                       f"etot_conv_thr={SPEC['fixed']['relax_etot_conv_thr_Ry']}, nstep={SPEC['fixed']['relax_max_steps']}")
    lines = ['&CONTROL', *control, '/', '&SYSTEM', *system, '/', '&ELECTRONS',
             f"conv_thr={SPEC['fixed']['conv_thr_Ry']}, mixing_beta={SPEC['fixed']['mixing_beta']}, electron_maxstep=200", '/']
    if req['kind'] == 'relax':
        lines += ['&IONS', "ion_dynamics='bfgs'", '/']
    lines += ['ATOMIC_SPECIES', f'Ge {mass} {upf}', 'CELL_PARAMETERS angstrom']
    lines += ['  '.join(f'{v:.12f}' for v in row) for row in cell]
    lines += ['ATOMIC_POSITIONS crystal'] + ['Ge ' + '  '.join(f'{v:.12f}' for v in row) for row in frac]
    n, s = req['kmesh'], req['kshift']
    lines += ['K_POINTS automatic', f'{n} {n} {n} {s} {s} {s}', '']
    (folder/'pw.in').write_text('\n'.join(lines))
    return cell, frac


def parse(folder, req, cell):
    text = (folder/'pw.out').read_text(errors='ignore')
    if 'JOB DONE.' not in text:
        raise ValueError('pw.x did not finish; the job is charged (see pw.out)')
    if 'convergence NOT achieved' in text:
        raise ValueError('SCF did not converge within 200 iterations; the job is charged')
    energies = [float(x) for x in re.findall(r'^!\s+total energy\s+=\s+([-\d.]+)\s+Ry', text, re.M)]
    block = text[text.rfind('total   stress'):].split('\n')
    stress_kbar = [[float(x) for x in block[i].split()[3:6]] for i in (1, 2, 3)]
    forces = [[[float(x) for x in m.groups()] for m in
               re.finditer(r'atom\s+\d+ type\s+\d+\s+force =\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', chunk)]
              for chunk in text.split('Forces acting on atoms')[1:]]
    forces = [f[:2] for f in forces]
    nk = int(re.findall(r'number of k points=\s+(\d+)', text)[-1])
    nsym = re.findall(r'(\d+) Sym\. Ops\.', text)
    iterations = [int(x) for x in re.findall(r'convergence has been achieved in\s+(\d+)\s+iterations', text)]
    root = ET.parse(folder/'work/ge.save/data-file-schema.xml').getroot()
    out = root.find('output')
    pos_bohr = np.array([[float(x) for x in at.text.split()] for at in out.find('atomic_structure/atomic_positions')])
    frac_final = np.linalg.solve(cell.T, (pos_bohr * BOHR_A).T).T
    bands = out.find('band_structure')
    res = {'total_energy_Ry': energies[-1],
           'stress_GPa': [[round(v / 10.0, 5) for v in row] for row in stress_kbar],
           'pressure_GPa': round(sum(stress_kbar[i][i] for i in range(3)) / 30.0, 5),
           'forces_Ry_per_bohr': forces[-1],
           'irreducible_k_points': nk, 'symmetry_operations': int(nsym[-1]) if nsym else None,
           'scf_iterations': iterations[-1] if req['kind'] == 'scf' else iterations}
    if req['occupations'] == 'smearing':
        res['fermi_energy_eV'] = round(float(bands.findtext('fermi_energy')) * 2 * RY_EV, 5)
    else:
        res['highest_occupied_level_eV'] = round(float(bands.findtext('highestOccupiedLevel')) * 2 * RY_EV, 5)
        lu = bands.findtext('lowestUnoccupiedLevel')
        if lu is not None:
            res['lowest_unoccupied_level_eV'] = round(float(lu) * 2 * RY_EV, 5)
    best = None
    for ks in bands.findall('ks_energies'):
        k = np.array([float(x) for x in ks.findtext('k_point').split()])        # Cartesian, units 2 pi / alat
        if best is None or np.linalg.norm(k) < best[0] - 1e-9:
            best = (np.linalg.norm(k), k, [round(float(x) * 2 * RY_EV, 4) for x in ks.findtext('eigenvalues').split()])
    alat_A = float(out.find('atomic_structure').get('alat')) * BOHR_A
    res['mesh_point_nearest_gamma'] = {'k_cartesian_inverse_A': [round(v * 2 * math.pi / alat_A, 6) for v in best[1]],
                                       'band_energies_eV': best[2]}
    if req['kind'] == 'relax':
        steps = [{'total_energy_Ry': e, 'max_force_component_Ry_per_bohr': max(abs(v) for atom in f for v in atom)}
                 for e, f in zip(energies, forces)]
        res['ionic_steps'] = steps
        res['relax_converged'] = 'bfgs converged' in text
        res['final_positions_crystal'] = [[round(v, 10) for v in row] for row in frac_final]
    return res


def run(request, workdir):
    req = normalize(request)
    nk_expected, nsym_expected, _ = size_and_charge(req)
    folder = Path(workdir)
    folder.mkdir(parents=True, exist_ok=True)
    cell, frac = write_input(req, folder)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               OMPI_MCA_btl_vader_single_copy_mechanism='none')
    start = time.monotonic()
    with (folder/'pw.out').open('w') as out, (folder/'pw.err').open('w') as err:
        try:
            proc = subprocess.run(['mpirun', '--allow-run-as-root', '--oversubscribe', '-np', str(CPUS), 'pw.x',
                                   '-npool', str(CPUS), '-in', 'pw.in'], cwd=folder, stdout=out, stderr=err, env=env,
                                  timeout=SPEC['fixed']['timeout_s'])
        except subprocess.TimeoutExpired:
            raise ValueError('pw.x exceeded the time limit; the job is charged') from None
    if proc.returncode:
        raise ValueError('pw.x exited with an error; the job is charged (see pw.out)')
    result = parse(folder, req, cell)
    result['wall_seconds'] = round(time.monotonic() - start, 1)
    result['cell_vectors_A'] = [[round(v, 10) for v in row] for row in cell]
    result['input_positions_crystal'] = [[round(v, 10) for v in row] for row in frac]
    result['settings'] = req
    result['charged_k_points'] = nk_expected
    shutil.rmtree(folder/'work', ignore_errors=True)        # keep pw.in/pw.out as evidence, drop wavefunctions
    return result


# ---------------------------------------------------------------- grading
def grade(payload, trace, validation):
    private = json.loads((HERE/'private/hidden.json').read_text())
    ref, tol = private['reference'], private['tolerances']
    fields = list(ref)
    missing = [f for f in fields + ['evidence_ids'] if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    bad = [f for f in fields if not _number(payload[f])]
    if bad:
        return False, {'error': f'fields must be finite numbers: {bad}'}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    for field in fields:
        if field in tol['relative']:
            err = payload[field] / ref[field] - 1
            ok = abs(err) <= tol['relative'][field]
        else:
            err = payload[field] - ref[field]
            ok = abs(err) <= tol['absolute'][field]
        checks[field] = ok
        errors[field] = err
    return all(checks.values()), {'checks': checks, 'errors': errors,
                                  'note': 'relative errors for %s, absolute otherwise' % sorted(tol['relative'])}
