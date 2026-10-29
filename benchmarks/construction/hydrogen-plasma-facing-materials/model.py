"""Shared analysis of the archived QE runs (construction side). Pure functions of pool.json records."""
import json
from pathlib import Path
import numpy as np

RY = 13.605693122994          # eV per Ry (CODATA 2018), as used in the reference
CACHE = Path.home()/'.cache/ai-theorist/hydrogen-plasma-facing-materials'


def load_pool(cache=CACHE):
    return json.loads((Path(cache)/'pool.json').read_text())


def energy(rec):
    """Final energy of a run: 'Final energy' of a converged relaxation, else the last SCF total energy."""
    o = rec['output']
    return o['final_energy_Ry'] if o['final_energy_Ry'] is not None else o['last_total_energy_Ry']


def converged(rec):
    msg = rec['output']['relaxation_message']
    return msg is None or msg.startswith('bfgs converged')


def final_positions(rec):
    return rec['output'].get('positions_crystal_final') or rec['positions_crystal_initial']


def h_environment(rec, initial=False, n=8):
    """Sorted (distance, element) of the n host atoms nearest to the interstitial H (minimum image)."""
    pos = rec['positions_crystal_initial'] if initial else final_positions(rec)
    cell = np.array(rec['cell_A'])
    sp = [p[0] for p in pos]
    x = np.array([p[1:] for p in pos])
    iH = sp.index('H')
    d = x - x[iH]
    d -= np.round(d)
    shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    r = np.linalg.norm((d[:, None, :] + shifts[None]) @ cell, axis=2).min(axis=1)
    r[iH] = np.inf
    order = np.argsort(r)[:n]
    return [(float(r[k]), sp[k]) for k in order]


def hife(pool, mat, site, ref='H2'):
    """H interstitial formation energy (eV) of one archived relaxation."""
    e_h = energy(pool['H-IFE/H2']) / 2 if ref == 'H2' else energy(pool['H-IFE/H'])
    return (energy(pool[f'H-IFE/{mat}/site-{site}']) - energy(pool[f'H-IFE/{mat}/supercell']) - e_h) * RY


def site_ids(pool, mat):
    return sorted((k.split('site-')[1] for k in pool if k.startswith(f'H-IFE/{mat}/site-')), key=int)


def distinct_sites(pool, mat, e_tol=0.005, r_tol=0.05, ref="H2"):
    """Cluster the archived relaxations of one host into distinct relaxed configurations.

    Two runs are the same site when their formation energies agree within e_tol (eV) and the distances of
    the eight nearest host atoms (with species) agree within r_tol (A); r_tol=None clusters by energy only. Returns clusters sorted by energy;
    each cluster: dict(sites, energies, rep_energy (lowest converged), env).
    """
    runs = []
    for s in site_ids(pool, mat):
        rec = pool[f'H-IFE/{mat}/site-{s}']
        runs.append(dict(site=s, e=hife(pool, mat, s, ref), env=h_environment(rec), conv=converged(rec)))
    clusters = []
    for r in sorted(runs, key=lambda r: r['e']):
        for c in clusters:
            same_env = r_tol is None or all(a[1] == b[1] and abs(a[0] - b[0]) < r_tol for a, b in zip(c['env'], r['env']))
            if abs(c['energies'][0] - r['e']) < e_tol and same_env:
                c['sites'].append(r['site'])
                c['energies'].append(r['e'])
                c['conv'].append(r['conv'])
                break
        else:
            clusters.append(dict(sites=[r['site']], energies=[r['e']], conv=[r['conv']], env=r['env']))
    for c in clusters:
        conv_e = [e for e, ok in zip(c['energies'], c['conv']) if ok] or c['energies']
        c['rep_energy'] = min(conv_e)
        c['rep_site'] = c['sites'][c['energies'].index(c['rep_energy'])]
    return sorted(clusters, key=lambda c: c['rep_energy'])


def slab_area(rec):
    c = np.array(rec['cell_A'])
    return float(np.linalg.norm(np.cross(c[0], c[1])))


def surface_energy(pool, mat, label):
    """Eq. (11) of the paper: (E_slab - N_slab/N_bulk E_bulk) / 2A, eV/A^2, relaxed slab and relaxed bulk cell."""
    slab, uc = pool[f'SBE/{mat}/slab/{label}'], pool[f'SBE/{mat}/unitcell']
    return (energy(slab) - slab['input']['nat'] / uc['input']['nat'] * energy(uc)) * RY / (2 * slab_area(slab))


def terminations(pool, mat):
    return sorted({k.split('/')[3] for k in pool if k.startswith(f'SBE/{mat}/slab_supercell/')})


def vacancies(pool, mat, term):
    pre = f'SBE/{mat}/vacancy/{term}/'
    return sorted(k[len(pre):] for k in pool if k.startswith(pre))


def sbe(pool, mat, term, label, atom_id=None):
    el = label.split('_')[0]
    atom = pool[atom_id or f'SBE/{mat}/atom/{el}']
    return (energy(pool[f'SBE/{mat}/vacancy/{term}/{label}']) + energy(atom)
            - energy(pool[f'SBE/{mat}/slab_supercell/{term}'])) * RY


def cohesive(pool, mat, bulk_id=None, atom_ids=None):
    """[sum of isolated-atom energies - bulk energy] / N, eV/atom."""
    bulk = pool[bulk_id or f'SBE/{mat}/unitcell']
    counts = bulk['input']['species_counts']
    atoms = sum(n * energy(pool[(atom_ids or {}).get(el, f'SBE/{mat}/atom/{el}')]) for el, n in counts.items())
    return (atoms - energy(bulk)) / bulk['input']['nat'] * RY
