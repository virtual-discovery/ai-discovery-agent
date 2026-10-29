"""Lab logic for tan-pressure-kappa-channels: a phonon Boltzmann-transport service for theta-TaN under pressure.

The lab evaluates harmonic phonons, three-phonon and isotope scattering, and the RTA or iterative BTE
(tanbte.py, ShengBTE conventions) from the authors' archived force constants. Four-phonon rates are lookups
precomputed at construction from the archived fourth-order constants (FourPhonon conventions, Monte Carlo
sampling of the 16x16x16 grid); phonon-electron rates are the archived BTE.sr_phel_metal tables. At 40 and 80 GPa
the fourth-order constants and phonon-electron tables are withheld, and four-phonon rates at 60 GPa / 600 K are
withheld; reference conductivities computed with the same code path grade the submission.
"""
import json
import math
from pathlib import Path
import threading
import numpy as np
import tanbte as tb

HERE = Path(__file__).parent
DATA = HERE/'data'
META = json.loads((DATA/'meta.json').read_text())
PUBLIC = json.loads((HERE/'public.json').read_text())
PRIVATE = json.loads((HERE/'private/hidden.json').read_text())
PRESSURES = META['pressures_GPa']
MESHES = META['meshes']
_LOCK = threading.Lock()
_CRYSTAL, _GRID, _THREE = {}, {}, {}


def crystal(p):
    if p not in _CRYSTAL:
        with np.load(DATA/('P%d.npz' % p)) as z:
            poscar = {'A': z['A'], 'species': ['Ta', 'N'], 'frac': z['frac']}
            fc3 = {'R': z['fc3_R'], 'atoms': z['fc3_atoms'], 'fc': z['fc3_fc']}
            _CRYSTAL[p] = tb.Crystal(poscar, z['fc2'], fc3=fc3)
    return _CRYSTAL[p]


def grid(p, mesh):
    key = (p, mesh)
    if key not in _GRID:
        _GRID[key] = tb.Grid(crystal(p), (mesh,) * 3)
    return _GRID[key]


def three(p, mesh, sb, T):
    key = (p, mesh, sb, T)
    if key not in _THREE:
        g = grid(p, mesh)
        rp, rm, procs = g.three_phonon(sb, T)
        iso = g.isotope_rates(sb)
        if len(_THREE) > 6:
            _THREE.pop(next(iter(_THREE)))
        _THREE[key] = (rp, rm, iso, procs)
    return _THREE[key]


def fourph(p, T):
    f = DATA/('fourph_P%d_T%d.npz' % (p, T))
    if not f.exists():
        return None
    with np.load(f) as z:
        return {k: z[k] for k in ('pp', 'pm', 'mm')}


def phel(p):
    f = DATA/'archive'/str(p)/'BTE.sr_phel_metal'
    if not f.exists():
        return None
    return tb.read_phel(f, grid(p, 16))[1]


def sig(a, n=6):
    a = np.asarray(a, float)
    out = np.zeros_like(a)
    nz = a != 0
    mag = np.floor(np.log10(np.abs(a[nz])))
    out[nz] = np.round(a[nz] * 10.0 ** (n - 1 - mag)) / 10.0 ** (n - 1 - mag)
    return out.tolist()


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def check_common(r, keys):
    extra = set(r) - set(keys) - {'kind'}
    if extra:
        raise ValueError('unknown fields %s' % sorted(extra))
    if r.get('pressure_GPa') not in PRESSURES:
        raise ValueError('pressure_GPa must be one of %s' % PRESSURES)
    if 'mesh' in keys and r.get('mesh', 16) not in MESHES:
        raise ValueError('mesh must be one of %s (n x n x n, Gamma-centred)' % MESHES)
    if 'T_K' in keys and not (number(r.get('T_K')) and 50 <= r['T_K'] <= 1200):
        raise ValueError('T_K must be a number between 50 and 1200')
    if 'scalebroad' in keys and not (number(r.get('scalebroad')) and 0.05 <= r['scalebroad'] <= 2.0):
        raise ValueError('scalebroad must be a number between 0.05 and 2')


def validate(r):
    kind = r.get('kind')
    costs = {k: v['cost'] for k, v in PUBLIC['requests'].items()}
    if kind == 'archived_file':
        check_common(r, ['pressure_GPa', 'file'])
        files = META['archived_files'][str(r['pressure_GPa'])]
        if r.get('file') not in files:
            raise ValueError('file must be one of %s at this pressure (others are withheld or not served)' % files)
        return costs[kind]
    if kind == 'harmonic':
        check_common(r, ['pressure_GPa', 'mesh'])
        return costs[kind]
    if kind == 'rates':
        check_common(r, ['pressure_GPa', 'mesh', 'T_K', 'scalebroad'])
        return costs[kind]
    if kind == 'kappa':
        check_common(r, ['pressure_GPa', 'mesh', 'T_K', 'scalebroad', 'channels', 'solver', 'extra_rates_ps'])
        ch = r.get('channels')
        if not isinstance(ch, list) or not ch or not set(ch) <= {'3ph', '4ph', 'isotope', 'phel'} \
                or len(set(ch)) != len(ch):
            raise ValueError('channels must be a non-empty list drawn from "3ph", "4ph", "isotope", "phel"')
        if r.get('solver') not in ('RTA', 'iterative'):
            raise ValueError('solver must be "RTA" or "iterative"')
        if r['solver'] == 'iterative' and '3ph' not in ch:
            raise ValueError('the iterative solver iterates three-phonon scattering: include "3ph"')
        p, mesh, T = r['pressure_GPa'], r.get('mesh', 16), r['T_K']
        if '4ph' in ch:
            if mesh != 16 or abs(r['scalebroad'] - 0.1) > 1e-12:
                raise ValueError('four-phonon rates are precomputed lookups on the archived protocol only '
                                 '(mesh 16, scalebroad 0.1)')
            if not (DATA/('fourph_P%d_T%d.npz' % (p, T))).exists():
                raise ValueError('four-phonon rates are not available at %s GPa and %s K; available (GPa, K): %s'
                                 % (p, T, META['fourph_available']))
        if 'phel' in ch:
            if mesh != 16:
                raise ValueError('the archived phonon-electron rates exist on the 16x16x16 grid only')
            if 'BTE.sr_phel_metal' not in META['archived_files'][str(p)]:
                raise ValueError('the phonon-electron rates at this pressure are withheld')
        if 'extra_rates_ps' in r:
            x = r['extra_rates_ps']
            if mesh != 16 or not isinstance(x, list) or len(x) != META['nirr16'] \
                    or not all(isinstance(row, list) and len(row) == 6 and all(number(v) and v >= 0 for v in row)
                               for row in x):
                raise ValueError('extra_rates_ps must be a %d x 6 list of non-negative rates (ps^-1) on the '
                                 'catalog irreducible q-points of the 16x16x16 grid (mesh 16 only)' % META['nirr16'])
        return costs[kind] * (2 if mesh > 16 else 1)
    raise ValueError('kind must be "archived_file", "harmonic", "rates" or "kappa"')


def run(r, workdir):
    with _LOCK:
        return _run(r)


def _run(r):
    kind, p = r['kind'], r['pressure_GPa']
    if kind == 'archived_file':
        text = (DATA/'archive'/str(p)/r['file']).read_text()
        return {'pressure_GPa': p, 'file': r['file'], 'text': text,
                'source': 'author-archived file, served verbatim (Zenodo 10564835)'}
    mesh = r.get('mesh', 16)
    g = grid(p, mesh)
    if kind == 'harmonic':
        wi, vi = g.w[g.irr], g.v[g.irr]
        wall = g.w[g.valid]
        hist, edges = np.histogram(wall, bins=80, range=(0, 1.02 * wall.max()))
        ac, op = g.w[:, :3].max(), g.w[:, 3:].min()
        return {'pressure_GPa': p, 'mesh': [mesh] * 3, 'lattice_A': sig(g.cr.A, 8),
                'irreducible_q_reduced': sig(g.q[g.irr], 8), 'multiplicity': g.mult.tolist(),
                'omega_rad_ps': sig(wi), 'group_velocity_km_s': sig(vi),
                'dos_bin_edges_rad_ps': sig(edges), 'dos_counts': hist.tolist(),
                'acoustic_max_rad_ps': float(ac), 'optical_min_rad_ps': float(op),
                'note': 'branches sorted by frequency at each q; velocities from the analytic derivative of D(q)'}
    T, sb = float(r['T_K']), float(r['scalebroad'])
    rp, rm, iso, procs = three(p, mesh, sb, T)
    if kind == 'rates':
        out = {'pressure_GPa': p, 'mesh': [mesh] * 3, 'T_K': T, 'scalebroad': sb,
               'irreducible_q_reduced': sig(g.q[g.irr], 8), 'multiplicity': g.mult.tolist(),
               'omega_rad_ps': sig(g.w[g.irr]),
               'rate_3ph_absorption_ps': sig(rp, 5), 'rate_3ph_emission_ps': sig(rm, 5),
               'rate_isotope_ps': sig(iso, 5)}
        f4 = fourph(p, int(round(T))) if (mesh == 16 and abs(sb - 0.1) < 1e-12 and float(T).is_integer()) else None
        if f4 is not None:
            out.update(rate_4ph_pp_ps=sig(f4['pp'], 5), rate_4ph_pm_ps=sig(f4['pm'], 5), rate_4ph_mm_ps=sig(f4['mm'], 5))
        else:
            out['rate_4ph'] = 'not available for this request (lookups: mesh 16, scalebroad 0.1, (GPa, K) in %s)' \
                % META['fourph_available']
        out['units'] = 'angular frequency rad/ps; scattering rates 1/tau in ps^-1'
        return out
    ch = r['channels']
    rate = np.zeros_like(rp)
    if '3ph' in ch:
        rate += rp + rm
    if 'isotope' in ch:
        rate += iso
    if '4ph' in ch:
        f4 = fourph(p, int(round(T)))
        rate += f4['pp'] + f4['pm'] + f4['mm']
    if 'phel' in ch:
        rate += phel(p)
    if 'extra_rates_ps' in r:
        rate += np.array(r['extra_rates_ps'], float)
    k, info = g.kappa(T, rate, procs if r['solver'] == 'iterative' else None, iterate=r['solver'] == 'iterative')
    out = {'pressure_GPa': p, 'mesh': [mesh] * 3, 'T_K': T, 'scalebroad': sb, 'channels': ch,
           'solver': r['solver'], 'extra_rates_supplied': 'extra_rates_ps' in r,
           'kappa_W_mK': sig(k, 6), 'zero_rate_modes_on_full_grid': info['zero_rate_modes']}
    if r['solver'] == 'iterative':
        out.update(iterations=info['iterations'], converged=info['converged'],
                   kappa_RTA_W_mK=sig(info['kappa_rta'], 6),
                   iteration_history_diag_W_mK=[sig(h, 5) for h in info['history_diag'][:: max(1, len(info['history_diag']) // 25)]])
    return out


def grade(payload, trace, validation):
    ref, tol = PRIVATE['reference'], PRIVATE['tolerance']
    fields = list(ref)
    vals = {f: payload.get(f) for f in fields}
    bad = [f for f, v in vals.items() if not number(v)]
    if bad or 'evidence_ids' not in payload:
        return False, {'error': 'missing or non-finite fields', 'fields': bad + ([] if 'evidence_ids' in payload else ['evidence_ids'])}
    executed = {t['id'] for t in trace if t['status'] == 'complete'}
    ids = payload['evidence_ids']
    evidence = isinstance(ids, list) and len(ids) > 0 and all(type(i) is int for i in ids) and set(ids) <= executed
    checks = {'evidence_executed': evidence}
    errors = {}
    for f in fields:
        if tol[f]['type'] == 'rel':
            e = vals[f] / ref[f] - 1
        else:
            e = vals[f] - ref[f]
        errors[f] = round(e, 4)
        checks[f] = abs(e) <= tol[f]['value']
    return all(checks.values()), {'checks': checks, 'errors': errors}
