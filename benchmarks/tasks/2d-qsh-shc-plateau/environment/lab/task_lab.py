"""Lab logic for 2d-qsh-shc-plateau: a numpy evaluator of archived Wannier Hamiltonians plus archived SHC files.

Computing requests evaluate the archived Wannier90 Hamiltonians (aiida_hr.dat of eight monolayer insulators of
Materials Cloud record rfcdy-wcy32, Wannier_Hamiltonians/insulators/, stored as data/<material>.npz with Wannier90's
Wigner-Seitz image selection applied) with qsh_tb.py: band energies, band edges on a mesh, hybrid Wannier charge
centres, ribbon spectra and the Kubo spin Hall conductivity. `archived_shc` returns the authors' archived Wannier90
SHC Fermi scans of Top_materials/ (relaxed structures). Graded values are held in private/hidden.json.
"""
import json
import math
import subprocess
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import qsh_tb  # noqa: E402

PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
META = json.loads((HERE/'data/meta.json').read_text())
MATERIALS = tuple(META['materials'])
ARCHIVED = json.loads((HERE/'data/archived_shc.json').read_text())
NPROC = 2
WCC_UNIT = (6000, 24000)      # len(k) * n_loop up to which a wcc request costs 1 or 2 units (3 above)
MAX_K = 160000
K_PER_UNIT = 40000
SPIN = {'x': 0, 'y': 1, 'z': 2}
SOURCE_TB = ('new evaluation of the archived Wannier90 Hamiltonian (aiida_hr.dat with use_ws_distance image selection); '
             'energies in eV on the scale of the archived calculation; Wannier bands numbered from 1')
SOURCE_ARCHIVE = 'retrieved author-archived result (no new calculation)'
_MODELS = {}


def model(name):
    if name not in _MODELS:
        _MODELS[name] = qsh_tb.Model(HERE/'data'/f'{name}.npz', META['materials'][name]['cell_A'])
    return _MODELS[name]


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def integer(x):
    return type(x) is int


def check_material(request):
    if request.get('material') not in MATERIALS:
        raise ValueError('material must be one of ' + ', '.join(MATERIALS))
    return META['materials'][request['material']]['num_wann']


def kpoints(value, limit):
    if (not isinstance(value, list) or not 1 <= len(value) <= limit or
            not all(isinstance(p, list) and len(p) == 2 and all(number(x) and abs(x) <= 10 for x in p) for p in value)):
        raise ValueError(f'k_frac must be a list of 1 to {limit} points [k1, k2] (reduced coordinates)')
    return np.c_[np.array(value, float), np.zeros(len(value))]


def shc_size(request):
    n1, n2 = request['mesh']
    total = n1 * n2
    if 'refine' in request:
        r = request['refine']
        total += r['factor'] ** 2 * r['max_points']
    return total


def validate(request):
    kind = request.get('kind')
    if kind == 'bands':
        if not {'kind', 'material', 'k_frac'} <= set(request) <= {'kind', 'material', 'k_frac', 'spin'}:
            raise ValueError('Use {"kind": "bands", "material": ..., "k_frac": [[k1, k2], ...] (<= 1000), "spin": false}')
        check_material(request)
        kpoints(request['k_frac'], 1000)
        if not isinstance(request.get('spin', False), bool):
            raise ValueError('spin must be true or false')
        return 1
    if kind == 'band_edges':
        if set(request) != {'kind', 'material', 'band', 'mesh'}:
            raise ValueError('Use {"kind": "band_edges", "material": ..., "band": n, "mesh": N (12-300)}')
        nw = check_material(request)
        if not integer(request['band']) or not 1 <= request['band'] <= nw - 1:
            raise ValueError(f'band must be an integer from 1 to {nw - 1}')
        if not integer(request['mesh']) or not 12 <= request['mesh'] <= 300:
            raise ValueError('mesh must be an integer from 12 to 300')
        return 1
    if kind == 'wcc':
        if set(request) != {'kind', 'material', 'bands', 'loop_axis', 'k', 'n_loop'}:
            raise ValueError('Use {"kind": "wcc", "material": ..., "bands": n, "loop_axis": 0 or 1, "k": [...] (<= 201), '
                             '"n_loop": 20-300}')
        nw = check_material(request)
        if not integer(request['bands']) or not 1 <= request['bands'] <= nw - 1:
            raise ValueError(f'bands must be an integer from 1 to {nw - 1}')
        if request['loop_axis'] not in (0, 1) or isinstance(request['loop_axis'], bool):
            raise ValueError('loop_axis must be 0 or 1')
        k = request['k']
        if not isinstance(k, list) or not 1 <= len(k) <= 201 or not all(number(x) and abs(x) <= 10 for x in k):
            raise ValueError('k must be a list of 1 to 201 reduced coordinates')
        if not integer(request['n_loop']) or not 20 <= request['n_loop'] <= 300:
            raise ValueError('n_loop must be an integer from 20 to 300')
        size = len(k) * request['n_loop']
        return 1 if size <= WCC_UNIT[0] else (2 if size <= WCC_UNIT[1] else 3)
    if kind == 'shc':
        keys = {'kind', 'material', 'spin', 'mesh', 'eta_eV', 'position', 'mu_eV'}
        if not keys <= set(request) <= keys | {'refine'}:
            raise ValueError('Use {"kind": "shc", "material": ..., "spin": ["z"], "mesh": [N1, N2], "eta_eV": 0.001, '
                             '"position": "lattice" or "wannier_centres", "mu_eV": [...], "refine": {...} (optional)}')
        check_material(request)
        sp = request['spin']
        if not isinstance(sp, list) or not sp or len(set(sp)) != len(sp) or not all(s in SPIN for s in sp):
            raise ValueError('spin must be a non-empty list of distinct axes from "x", "y", "z"')
        mesh = request['mesh']
        if not isinstance(mesh, list) or len(mesh) != 2 or not all(integer(n) and 12 <= n <= 300 for n in mesh):
            raise ValueError('mesh must be [N1, N2] with integers from 12 to 300')
        if not number(request['eta_eV']) or not 1e-4 <= request['eta_eV'] <= 0.05:
            raise ValueError('eta_eV must be between 0.0001 and 0.05')
        if request['position'] not in ('lattice', 'wannier_centres'):
            raise ValueError('position must be "lattice" or "wannier_centres"')
        mu = request['mu_eV']
        if not isinstance(mu, list) or not 1 <= len(mu) <= 41 or not all(number(x) and abs(x) < 50 for x in mu):
            raise ValueError('mu_eV must be a list of 1 to 41 energies (eV)')
        if 'refine' in request:
            r = request['refine']
            if (not isinstance(r, dict) or set(r) != {'factor', 'threshold_A2', 'max_points'} or
                    not integer(r['factor']) or not 2 <= r['factor'] <= 16 or not number(r['threshold_A2']) or
                    r['threshold_A2'] <= 0 or not integer(r['max_points']) or not 1 <= r['max_points'] <= 5000):
                raise ValueError('refine must be {"factor": 2-16, "threshold_A2": > 0, "max_points": 1-5000}')
        total = shc_size(request)
        if total > MAX_K:
            raise ValueError(f'too many k-points: N1*N2 + factor^2*max_points = {total} > {MAX_K}')
        return math.ceil(total / K_PER_UNIT) + len(sp) - 1
    if kind == 'ribbon':
        if set(request) != {'kind', 'material', 'finite_axis', 'width', 'k', 'window_eV'}:
            raise ValueError('Use {"kind": "ribbon", "material": ..., "finite_axis": 0 or 1, "width": cells, '
                             '"k": [...] (<= 31), "window_eV": [lo, hi]}')
        nw = check_material(request)
        if request['finite_axis'] not in (0, 1) or isinstance(request['finite_axis'], bool):
            raise ValueError('finite_axis must be 0 or 1')
        if not integer(request['width']) or not 4 <= request['width'] or request['width'] * nw > 1300:
            raise ValueError(f'width must be an integer >= 4 with width * num_wann <= 1300 (num_wann = {nw})')
        k = request['k']
        if not isinstance(k, list) or not 1 <= len(k) <= 31 or not all(number(x) and abs(x) <= 10 for x in k):
            raise ValueError('k must be a list of 1 to 31 reduced coordinates')
        w = request['window_eV']
        if not isinstance(w, list) or len(w) != 2 or not all(number(x) for x in w) or not 0 < w[1] - w[0] <= 0.6:
            raise ValueError('window_eV must be [lo, hi] with 0 < hi - lo <= 0.6 eV')
        return 2
    if kind == 'archived_shc':
        if set(request) != {'kind', 'folder'} or request['folder'] not in ARCHIVED:
            raise ValueError('Use {"kind": "archived_shc", "folder": one of ' + ', '.join(ARCHIVED) + '}')
        return 1
    raise ValueError('kind must be one of bands, band_edges, wcc, shc, ribbon, archived_shc')


def run_shc(request, workdir):
    name = request['material']
    m = model(name)
    gammas = [SPIN[s] for s in request['spin']]
    mus = [float(x) for x in request['mu_eV']]
    n1, n2 = request['mesh']
    refine = request.get('refine')
    spec = {'npz': str(HERE/'data'/f'{name}.npz'), 'cell_A': META['materials'][name]['cell_A'], 'gammas': gammas,
            'eta': float(request['eta_eV']), 'mus': mus, 'convention': request['position'], 'mesh': [n1, n2],
            'refine': None if refine is None else [refine['factor'], refine['threshold_A2'], refine['max_points']]}
    (workdir/'request.json').write_text(json.dumps(spec))

    def workers(phase):
        procs = [subprocess.Popen([sys.executable, str(HERE/'shc_worker.py'), str(workdir), phase, str(i), str(NPROC)],
                                  stderr=subprocess.PIPE) for i in range(NPROC)]
        for p in procs:
            _, err = p.communicate()
            if p.returncode:
                raise RuntimeError('worker failed: ' + err.decode()[-300:])
    workers('coarse')
    parts = [np.load(workdir/f'coarse-{i}.npy') for i in range(NPROC)]
    c = np.zeros((n2, n1, len(gammas), len(mus)))
    for i in range(NPROC):
        c[i::NPROC] = parts[i]
    c = c.reshape(n2 * n1, len(gammas), len(mus))
    nref = 0
    if spec['refine'] is not None:
        cand = qsh_tb.select_refine(c, spec['refine'])
        if len(cand):
            kpts = np.c_[(cand % n1) / n1, (cand // n1) / n2]
            np.save(workdir/'candidates.npy', kpts)
            workers('refine')
            for i in range(NPROC):
                c[cand[i::NPROC]] = np.load(workdir/f'refine-{i}.npy')
        nref = int(len(cand))
    sheet = np.pi * c.sum(0) / (m.area * n1 * n2)
    w90 = sheet * qsh_tb.E2_OVER_H / (m.c * 1e-8)
    for f in workdir.glob('*.npy'):
        f.unlink()
    return {'material': name, 'spin': request['spin'], 'mesh': [n1, n2], 'eta_eV': request['eta_eV'],
            'position': request['position'], 'mu_eV': mus,
            'sigma_xy_hbar_over_e_S_per_cm': {s: [float('%.6g' % v) for v in w90[i]] for i, s in enumerate(request['spin'])},
            'refined_points': nref, 'kpoints_evaluated': int(n1 * n2 + nref * (refine['factor'] ** 2 if refine else 0)),
            'units': 'Wannier90 shc units, (hbar/e) S/cm, with the volume of the archived cell (as aiida-shc-fermiscan.dat)',
            'source': SOURCE_TB}


def run(request, workdir):
    kind = request['kind']
    if kind == 'archived_shc':
        return {**ARCHIVED[request['folder']], 'source': SOURCE_ARCHIVE}
    if kind == 'shc':
        return run_shc(request, workdir)
    name = request['material']
    m = model(name)
    out = {'material': name, 'num_wann': m.nw, 'source': SOURCE_TB}
    if kind == 'bands':
        k = kpoints(request['k_frac'], 1000)
        if request.get('spin'):
            e, s = m.spin_expectation(k)
            out.update(spin_expectation=np.round(s, 4).tolist(),
                       spin_note='<sigma_x>, <sigma_y>, <sigma_z> of each band (projected-spinor Pauli operator)')
        else:
            e = m.bands(k)
        out.update(k_frac=k[:, :2].round(8).tolist(), energies_eV=np.round(e, 6).tolist())
        return out
    if kind == 'band_edges':
        n, N = request['band'], request['mesh']
        be = m.band_edges(n, (N, N))
        out.update(band=n, mesh=N, max_of_band_n={'energy_eV': round(be['vbm'][0], 6), 'k_frac': be['vbm'][1]},
                   min_of_band_n_plus_1={'energy_eV': round(be['cbm'][0], 6), 'k_frac': be['cbm'][1]},
                   smallest_direct_gap={'gap_eV': round(be['direct'][0], 6), 'k_frac': be['direct'][1]},
                   note='extrema over the Gamma-centred N x N mesh of reduced coordinates (no refinement)')
        return out
    if kind == 'wcc':
        w = m.wcc(request['bands'], [float(x) for x in request['k']], request['n_loop'], request['loop_axis'])
        out.update(bands=request['bands'], loop_axis=request['loop_axis'], k=request['k'], n_loop=request['n_loop'],
                   wcc=np.round(w, 6).tolist(),
                   note='eigenphases / 2 pi of the discrete Wilson loop of bands 1..n along reciprocal axis loop_axis, '
                        'at each reduced k along the other axis; centres in fractions of the lattice vector, in [0, 1)')
        return out
    if kind == 'ribbon':
        res = m.ribbon([float(x) for x in request['k']], request['width'], request['finite_axis'], tuple(request['window_eV']))
        out.update(finite_axis=request['finite_axis'], width=request['width'], k=request['k'], window_eV=request['window_eV'],
                   states=[{'energies_eV': np.round(e, 6).tolist(), 'weight_first_3_cells': np.round(l, 4).tolist(),
                            'weight_last_3_cells': np.round(r, 4).tolist()} for e, l, r in res],
                   note='ribbon finite along lattice vector finite_axis (open ends, no relaxation or passivation), '
                        'periodic along the other lattice vector with reduced momentum k')
        return out
    raise ValueError('unknown kind')


def grade(payload, trace, validation):
    t, th = PRIVATE['targets'], PRIVATE['thresholds']
    fields = ['z2', 'plateau_sigma_z', 'spin_direction', 'evidence_ids']
    missing = [f for f in fields if f not in payload]
    if missing:
        return False, {'missing_fields': missing}
    z2, pl, sd = payload['z2'], payload['plateau_sigma_z'], payload['spin_direction']
    if not isinstance(z2, dict) or set(z2) != set(MATERIALS) or not all(v in (0, 1) and not isinstance(v, bool) for v in z2.values()):
        return False, {'error': 'z2 must give 0 or 1 for every material'}
    if not isinstance(pl, dict) or set(pl) != set(MATERIALS) or not all(number(v) for v in pl.values()):
        return False, {'error': 'plateau_sigma_z must give a finite number for every material'}
    want = set(t['spin_direction'])
    if (not isinstance(sd, dict) or set(sd) != want or
            not all(isinstance(v, dict) and set(v) == {'theta_deg', 'phi_deg'} and all(number(x) for x in v.values())
                    for v in sd.values())):
        return False, {'error': 'spin_direction must give {"theta_deg": x, "phi_deg": x} for ' + ', '.join(sorted(want))}
    executed = {r['id'] for r in trace if r['status'] == 'complete'}
    ids = payload['evidence_ids']
    checks = {'evidence_executed': isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids)
              and set(ids) <= executed}
    errors = {}
    for mname in MATERIALS:
        checks[f'z2:{mname}'] = z2[mname] == t['z2'][mname]
        ref = t['plateau_sigma_z'][mname]
        errors[f'plateau:{mname}'] = pl[mname] - ref
        checks[f'plateau:{mname}'] = abs(errors[f'plateau:{mname}']) <= max(th['plateau_rel'] * abs(ref), th['plateau_abs_floor'])
    for mname, ref in t['spin_direction'].items():
        dth = sd[mname]['theta_deg'] - ref['theta_deg']
        dph = (sd[mname]['phi_deg'] - ref['phi_deg'] + 180.0) % 360.0 - 180.0
        errors[f'theta:{mname}'], errors[f'phi:{mname}'] = dth, dph
        checks[f'theta:{mname}'] = abs(dth) <= th['theta_deg']
        checks[f'phi:{mname}'] = abs(dph) <= th['phi_deg']
    return all(checks.values()), {'checks': checks, 'errors': errors,
                                  'error_units': 'plateau in (hbar/e)(e^2/h); angles in degrees'}
