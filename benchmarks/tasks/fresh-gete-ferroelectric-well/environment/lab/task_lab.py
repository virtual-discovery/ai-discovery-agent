"""Lab logic for fresh-gete-ferroelectric-well: new Quantum ESPRESSO 6.7 pw.x calculations for rhombohedral GeTe.

validate(request) checks a structured request against public.json (no defaults) and returns its charge, computed
from the canonical size of the request: the number of symmetry-distinct k points pw.x will use (counted here as
pw.x counts them), times factors for the pseudopotential set (spin-orbit) and the cutoff.
run(request, workdir) writes the pw.x input, runs `mpirun -np CPUS pw.x` (threads pinned to 1), parses the output
and returns raw physical results; pw.in/pw.out are kept as evidence and the pw.x scratch is always deleted.
Nothing hidden is ever served: grade() compares the terminal submission with reference values computed at
construction time through this same run() and stored in private/ (benchmarks/construction/fresh-gete-ferroelectric-well/).
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
FIXED = SPEC['fixed']
CPUS = int(os.environ.get('LAB_MPI_PROCS', '2'))
BOHR_A = 0.529177210903
RY_EV = 13.605693122994
RY_BOHR3_GPA = 14710.507848260711         # 1 Ry/bohr^3 in GPa
FIELDS = ('kind', 'pseudopotential_set', 'a_rh_A', 'alpha_deg', 'x_Te', 'ecutwfc_Ry', 'kmesh', 'kshift', 'occupations',
          'degauss_Ry')


# ---------------------------------------------------------------- geometry and symmetry
def cell_vectors(a, alpha_deg):
    """Rows a1, a2, a3 (Angstrom): a_i = a (p, q, q) and cyclic, threefold axis along Cartesian [111], so that at
    alpha = 60 deg the cell is the fcc primitive cell with its cube axes along x, y, z (pw.x then finds all 48
    operations of rocksalt)."""
    c = math.cos(math.radians(alpha_deg))
    q = (math.sqrt(1 + 2 * c) + math.sqrt(1 - c)) / 3
    p = q - math.sqrt(1 - c)
    return a * np.array([[p, q, q], [q, p, q], [q, q, p]])


def volume(a, alpha_deg):
    c = math.cos(math.radians(alpha_deg))
    return a ** 3 * math.sqrt(1 - 3 * c * c + 2 * c ** 3)


_CANDIDATES = np.array([m for m in itertools.product((-1, 0, 1), repeat=9)], float).reshape(-1, 3, 3)
_CANDIDATES = _CANDIDATES[np.abs(np.abs(np.linalg.det(_CANDIDATES)) - 1) < 1e-9]


def point_groups(cell, x, tol=1e-6):
    """Integer matrices M (crystal coordinates, r -> M r) of the lattice point group and of the crystal point group."""
    G = cell @ cell.T
    diff = np.einsum('nji,jk,nkl->nil', _CANDIDATES, G, _CANDIDATES) - G
    lattice = _CANDIDATES[np.abs(diff).max(axis=(1, 2)) < tol * np.abs(G).max()]
    tau = np.array([x, x, x])
    moved = lattice @ tau
    crystal = lattice[np.abs((moved - tau) - np.round(moved - tau)).max(axis=1) < 1e-7]
    return lattice, crystal


def irreducible_k(cell, x, n, shift):
    """Number of k points pw.x keeps for an n x n x n Monkhorst-Pack mesh (kpoint_grid + irreducible_BZ).

    kpoint_grid merges mesh points related by a lattice operation or by time reversal (only images that fall on the
    mesh count); irreducible_BZ then splits the lattice star of each kept point into orbits of the crystal group.
    k points are handled as integer coordinates on the (2n)^3 grid of half mesh steps."""
    lattice, crystal = point_groups(cell, x)
    m = 2 * n
    ops_l = np.array([M.T for M in lattice] + [-M.T for M in lattice]).astype(int)
    ops_c = np.array([M.T for M in crystal] + [-M.T for M in crystal]).astype(int)
    grid = 2 * np.array(list(itertools.product(range(n), repeat=3))) + shift          # integer coordinates
    code = lambda v: ((v[..., 0] % m) * m + v[..., 1] % m) * m + v[..., 2] % m        # noqa: E731
    lookup = -np.ones(m ** 3, int)
    lookup[code(grid)] = np.arange(len(grid))
    images = np.einsum('oij,pj->opi', ops_l, grid)                                      # (ops, points, 3)
    idx = lookup[code(images)]
    rep = np.where(idx >= 0, idx, len(grid)).min(axis=0)
    reps = np.unique(rep)
    if len(crystal) == len(lattice):
        return len(reps), len(crystal)
    total = 0
    for r in reps:
        star = np.unique(code(images[:len(lattice), r]))
        pts = np.stack([star // (m * m), (star // m) % m, star % m], axis=1)
        orbit = code(np.einsum('oij,pj->opi', ops_c, pts)).min(axis=0)
        total += len(np.unique(orbit))
    return total, len(crystal)


# ---------------------------------------------------------------- request checking and charge
def _number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def normalize(request):
    """Return the canonical request, or raise ValueError with a usable message (invalid requests are free)."""
    if not isinstance(request, dict):
        raise ValueError('JSON object required; GET /task (python3 /app/lab_client.py task) lists the schema')
    unknown, missing = set(request) - set(FIELDS), [f for f in FIELDS if f not in request]
    if unknown:
        raise ValueError(f'unknown fields {sorted(unknown)}; allowed: {list(FIELDS)}')
    if missing:
        raise ValueError(f'missing fields {missing}: every field is required (no defaults; use null where allowed)')
    allowed = SPEC['allowed']
    out = {}
    for key in ('kind', 'pseudopotential_set', 'ecutwfc_Ry', 'kmesh', 'kshift', 'occupations'):
        value = request[key]
        if isinstance(value, bool) or value not in allowed[key]:
            raise ValueError(f'{key} must be one of {allowed[key]}')
        out[key] = value
    for key, digits in (('a_rh_A', 6), ('alpha_deg', 5), ('x_Te', 6)):
        lo, hi = SPEC[key + '_range']
        if not _number(request[key]) or not lo <= request[key] <= hi:
            raise ValueError(f'{key} must be a number in [{lo}, {hi}]')
        out[key] = round(float(request[key]), digits)
    if out['occupations'] == 'fixed':
        if request['degauss_Ry'] is not None:
            raise ValueError('degauss_Ry must be null for fixed occupations')
        out['degauss_Ry'] = None
    else:
        if isinstance(request['degauss_Ry'], bool) or request['degauss_Ry'] not in allowed['degauss_Ry']:
            raise ValueError(f"degauss_Ry must be one of {allowed['degauss_Ry']} for smearing")
        out['degauss_Ry'] = request['degauss_Ry']
    return {k: out[k] for k in FIELDS}


def size_and_charge(req):
    cell = cell_vectors(req['a_rh_A'], req['alpha_deg'])
    nk, nsym = irreducible_k(cell, req['x_Te'], req['kmesh'], req['kshift'])
    rule = SPEC['charge']
    units = (nk * rule['set_factor'][req['pseudopotential_set']] * rule['ecut_factor'][str(req['ecutwfc_Ry'])]
             / rule['k_points_per_unit'])
    return nk, nsym, max(1, math.ceil(units - 1e-9))


def validate(request):
    req = normalize(request)
    nk, nsym, charge = size_and_charge(req)
    cap = SPEC['charge']['max_units_per_request']
    if charge > cap:
        raise ValueError(f'request too large: {nk} symmetry-distinct k points -> {charge} units; the maximum per request '
                         f'is {cap} (use a smaller mesh or cutoff, or a more symmetric structure)')
    return charge


# ---------------------------------------------------------------- pw.x input, execution, parsing
def write_input(req, folder):
    pp = SPEC['pseudopotential_sets'][req['pseudopotential_set']]
    soc = pp['spin_orbit']
    masses = FIXED['masses_amu']
    nelec = 10 if req['pseudopotential_set'] != 'sssp' else 20
    nbnd = (nelec if soc else nelec // 2) + FIXED['empty_bands'] * (2 if soc else 1)
    system = ['ibrav=0, nat=2, ntyp=2',
              f"ecutwfc={req['ecutwfc_Ry']}, ecutrho={FIXED['dual'] * req['ecutwfc_Ry']}", f'nbnd={nbnd}']
    if req['occupations'] == 'fixed':
        system.append("occupations='fixed'")
    else:
        system.append(f"occupations='smearing', smearing='mv', degauss={req['degauss_Ry']}")
    if soc:
        system.append('noncolin=.true., lspinorb=.true.')
    control = [f"calculation='{req['kind']}', prefix='gete', outdir='./work'", f"pseudo_dir='{HERE}'",
               "tprnfor=.true., tstress=.true., disk_io='low'"]
    lines = ['&CONTROL', *control, '/', '&SYSTEM', *system, '/', '&ELECTRONS',
             f"conv_thr={FIXED['conv_thr_Ry']}, mixing_beta={FIXED['mixing_beta']}, electron_maxstep=200", '/']
    lines += ['ATOMIC_SPECIES', f"Ge {masses['Ge']} {pp['Ge']}", f"Te {masses['Te']} {pp['Te']}"]
    lines += ['CELL_PARAMETERS angstrom'] + ['  '.join(f'{v:.12f}' for v in row)
                                             for row in cell_vectors(req['a_rh_A'], req['alpha_deg'])]
    lines += ['ATOMIC_POSITIONS crystal', 'Ge 0.0 0.0 0.0', f"Te {req['x_Te']:.10f} {req['x_Te']:.10f} {req['x_Te']:.10f}"]
    n, s = req['kmesh'], req['kshift']
    lines += ['K_POINTS automatic', f'{n} {n} {n} {s} {s} {s}', '']
    (folder/'pw.in').write_text('\n'.join(lines))
    return nelec


def _cell_params(cell):
    a = np.linalg.norm(cell, axis=1)
    cosines = [cell[i] @ cell[j] / (a[i] * a[j]) for i, j in ((0, 1), (1, 2), (0, 2))]
    return float(a.mean()), float(np.degrees(np.arccos(np.mean(cosines)))), float(abs(np.linalg.det(cell)))


def parse(folder, req, nelec):
    text = (folder/'pw.out').read_text(errors='ignore')
    if 'JOB DONE.' not in text:
        raise ValueError('pw.x did not finish; the job is charged (see pw.out)')
    if 'convergence NOT achieved' in text:
        raise ValueError('SCF did not converge within 200 iterations; the job is charged')
    soc = SPEC['pseudopotential_sets'][req['pseudopotential_set']]['spin_orbit']
    energies = [float(x) for x in re.findall(r'^!\s+total energy\s+=\s+([-\d.]+)\s+Ry', text, re.M)]
    stresses = []
    for m in re.finditer(r'total   stress .*P=\s*([-\d.]+)\s*\n(.*)\n(.*)\n(.*)\n', text):
        stresses.append(([[float(v) / 10 for v in m.group(i).split()[3:6]] for i in (2, 3, 4)], float(m.group(1)) / 10))
    forces = [[[float(v) for v in f.groups()] for f in
               re.finditer(r'atom\s+\d+ type\s+\d+\s+force =\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', chunk)][:2]
              for chunk in text.split('Forces acting on atoms')[1:]]
    nk = int(re.findall(r'number of k points=\s+(\d+)', text)[-1])
    nsym = re.findall(r'(\d+) Sym\. Ops\.', text)
    root = ET.parse(folder/'work/gete.save/data-file-schema.xml').getroot()
    out = root.find('output')
    cell = np.array([[float(v) for v in out.find(f'atomic_structure/cell/a{i}').text.split()] for i in (1, 2, 3)]) * BOHR_A
    pos = np.array([[float(v) for v in at.text.split()] for at in out.find('atomic_structure/atomic_positions')]) * BOHR_A
    frac = np.linalg.solve(cell.T, pos.T).T
    rel = frac[1] - frac[0]
    rel = rel - np.floor(rel)
    a, alpha, vol = _cell_params(cell)
    res = {'total_energy_Ry': energies[-1], 'stress_GPa': [[round(v, 4) for v in row] for row in stresses[-1][0]],
           'pressure_GPa': round(stresses[-1][1], 4), 'forces_Ry_per_bohr': forces[-1],
           'final_structure': {'a_rh_A': round(a, 6), 'alpha_deg': round(alpha, 5), 'x_Te': round(float(rel.mean()), 7),
                               'volume_A3': round(vol, 5)},
           'irreducible_k_points': nk, 'symmetry_operations': int(nsym[-1]) if nsym else None}
    bands = out.find('band_structure')
    if req['occupations'] == 'smearing':
        res['fermi_energy_eV'] = round(float(bands.findtext('fermi_energy')) * 2 * RY_EV, 5)
        ts = re.findall(r'smearing contrib\. \(-TS\)\s+=\s+([-\d.]+)\s+Ry', text)
        res['smearing_contribution_minus_TS_Ry'] = float(ts[-1]) if ts else None
    eig, kpts = [], []
    for ks in bands.findall('ks_energies'):
        kpts.append([float(v) for v in ks.findtext('k_point').split()])
        eig.append([float(v) * 2 * RY_EV for v in ks.findtext('eigenvalues').split()])
    eig, kpts = np.array(eig), np.array(kpts)
    alat = float(out.find('atomic_structure').get('alat')) * BOHR_A
    kcryst = kpts @ (cell / alat).T            # k . a_i in units of 2 pi
    nocc = nelec if soc else nelec // 2
    top, bottom = eig[:, nocc - 1], eig[:, nocc]
    res['band_edges_on_mesh'] = {
        'occupied_bands_per_k_insulator_filling': nocc,
        'valence_top_eV': round(float(top.max()), 4), 'valence_top_k_crystal': [round(v, 4) for v in kcryst[top.argmax()]],
        'conduction_bottom_eV': round(float(bottom.min()), 4),
        'conduction_bottom_k_crystal': [round(v, 4) for v in kcryst[bottom.argmin()]],
        'indirect_gap_eV': round(float(bottom.min() - top.max()), 4),
        'smallest_direct_gap_eV': round(float((bottom - top).min()), 4),
        'smallest_direct_gap_k_crystal': [round(v, 4) for v in kcryst[(bottom - top).argmin()]]}
    return res


def run(request, workdir):
    req = normalize(request)
    nk_expected, nsym_expected, _ = size_and_charge(req)
    folder = Path(workdir)
    folder.mkdir(parents=True, exist_ok=True)
    nelec = write_input(req, folder)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
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
            raise ValueError('pw.x exited with an error; the job is charged (see pw.out)')
        result = parse(folder, req, nelec)
    finally:        # keep pw.in/pw.out as evidence; always drop wavefunctions and scratch
        shutil.rmtree(folder/'work', ignore_errors=True)
    result['wall_seconds'] = round(time.monotonic() - start, 1)
    result['charged_k_points'] = nk_expected
    result['settings'] = req
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
