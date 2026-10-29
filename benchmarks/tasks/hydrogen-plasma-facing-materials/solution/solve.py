"""Reference solution for hydrogen-plasma-facing-materials; runs as the Harbor oracle and in check_task."""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CAT = json.loads((APP/'catalog.json').read_text())
RY = 13.605693122994
EVIDENCE, CACHE = [], {}


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(cid):
    if cid not in CACHE:
        r = lab('evaluate', {'kind': 'calculation', 'id': cid})
        while r['status'] == 'running':
            time.sleep(1)
            r = lab('job/%d' % r['id'])
        if r['status'] != 'complete':
            raise RuntimeError(r)
        EVIDENCE.append(r['id'])
        CACHE[cid] = r['result']
    return CACHE[cid]


def energy(run):
    o = run['output']
    return o['final_energy_Ry'] if o['final_energy_Ry'] is not None else o['last_total_energy_Ry']


def cutoffs(run):
    return run['input']['ecutwfc_Ry'], run['input']['ecutrho_Ry']


def environment(run, initial=False, n=8):
    pos = run['positions_crystal_initial'] if initial else (run['output'].get('positions_crystal_final')
                                                            or run['positions_crystal_initial'])
    cell = np.array(run['cell_A'])
    sp = [p[0] for p in pos]
    x = np.array([p[1:] for p in pos], float)
    iH = sp.index('H')
    d = x - x[iH]
    d -= np.round(d)
    shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    r = np.linalg.norm((d[:, None, :] + shifts[None]) @ cell, axis=2).min(axis=1)
    r[iH] = np.inf
    return [(float(r[k]), sp[k]) for k in np.argsort(r)[:n]]


def hydrogen(mat, e_half_h2):
    info = CAT['materials'][mat]['hydrogen_interstitial']
    bulk = buy(info['pristine_supercell'])
    runs = []
    for s in info['interstitial_runs']:
        run = buy(s['id'])
        assert cutoffs(run) == cutoffs(bulk) and run['input']['k_points'] == bulk['input']['k_points']
        assert np.allclose(run['cell_A'], bulk['cell_A'])
        msg = run['output']['relaxation_message'] or ''
        runs.append({'site': s['generation_number'], 'run': run, 'e': (energy(run) - energy(bulk) - e_half_h2) * RY,
                     'conv': msg.startswith('bfgs converged'), 'env': environment(run)})
    # The archive calls every run an "inequivalent site", but many relax into the same configuration:
    # same energy (to meV) and same coordination shell -> one site.
    clusters = []
    for r in sorted(runs, key=lambda r: r['e']):
        for c in clusters:
            same = all(a[1] == b[1] and abs(a[0] - b[0]) < 0.05 for a, b in zip(c[0]['env'], r['env']))
            if abs(c[0]['e'] - r['e']) < 0.005 and same:
                c.append(r)
                break
        else:
            clusters.append([r])
    reps = []
    for c in clusters:
        good = [r for r in c if r['conv']] or c
        reps.append(min(good, key=lambda r: r['e']))
    reps.sort(key=lambda r: r['e'])
    best = reps[0]
    d, el = best['env'][0]
    return {'E_f': best['e'], 'sites': [r['e'] for r in reps], 'nearest': (el, d), 'clusters': [[r['site'] for r in c] for c in clusters],
            'unconverged': [r['site'] for r in runs if not r['conv']], 'env': best['env'][:4]}


def surface_binding(mat, elements):
    sb = CAT['materials'][mat]['surface_binding']
    supers = {s['termination']: buy(s['id']) for s in sb['slab_supercells']}
    first = next(iter(supers.values()))
    for s in supers.values():   # same composition and cell: the energy difference is 2A * (difference of gamma)
        assert s['input']['species_counts'] == first['input']['species_counts'] and np.allclose(s['cell_A'], first['cell_A'])
    term = min(supers, key=lambda t: energy(supers[t]))
    out, notes = {}, []
    for el in elements:
        atom = buy(sb['isolated_atoms'][el])
        if cutoffs(atom) != cutoffs(supers[term]):
            notes.append(f'{mat}: the archived {el} atom uses {cutoffs(atom)} Ry but the slabs use {cutoffs(supers[term])} Ry')
        vals = []
        for v in sb['surface_vacancies']:
            if v['slab_supercell'].endswith('/' + term) and v['removed_element'] == el:
                vac = buy(v['id'])
                assert cutoffs(vac) == cutoffs(supers[term])
                vals.append((energy(vac) + energy(atom) - energy(supers[term])) * RY)
        out[el] = min(vals)
    return term, out, notes


def cohesive(mat):
    sb = CAT['materials'][mat]['surface_binding']
    bulk = buy(sb['bulk_unit_cell'])
    total = 0.0
    for el, n in bulk['input']['species_counts'].items():
        atom = buy(sb['isolated_atoms'][el])
        assert cutoffs(atom) == cutoffs(bulk), (mat, el)
        total += n * energy(atom)
    return (total - energy(bulk)) / bulk['input']['nat'] * RY


def main():
    graded = CAT['graded_materials']
    h2 = buy(CAT['hydrogen_references']['H2_molecule'])
    e_half = energy(h2) / 2
    H = {m: hydrogen(m, e_half) for m in graded['hydrogen']}
    S, terms, notes = {}, {}, []
    for m, els in graded['surface_binding'].items():
        terms[m], S[m], n = surface_binding(m, els)
        notes += n
    C = {m: cohesive(m) for m in graded['cohesive_energy']}
    submission = {
        'h_formation_energy_eV': {m: round(h['E_f'], 4) for m, h in H.items()},
        'distinct_site_energies_eV': {m: [round(e, 4) for e in h['sites']] for m, h in H.items()},
        'preferred_site': {m: {'nearest_host_element': h['nearest'][0], 'nearest_host_distance_A': round(h['nearest'][1], 4)}
                           for m, h in H.items()},
        'surface_binding_energy_eV': {m: {el: round(v, 4) for el, v in d.items()} for m, d in S.items()},
        'cohesive_energy_eV_per_atom': {m: round(v, 4) for m, v in C.items()},
        'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    lines = [f"- {m}: E_f = {h['E_f']:.3f} eV; distinct sites {[round(e, 3) for e in h['sites']]} from runs "
             f"{h['clusters']}; preferred site nearest host {h['nearest'][0]} at {h['nearest'][1]:.3f} A"
             + (f"; unconverged runs {h['unconverged']}" if h['unconverged'] else '') for m, h in H.items()]
    (APP/'report.md').write_text(f"""# Plasma-facing candidates: hydrogen uptake and surface binding (reference solution)

## Hydrogen interstitials
E_f = E(supercell + H) - E(supercell) - 1/2 E(H2), from the archived relaxed runs (PBEsol, 60/480 Ry, same supercell,
cell and k-mesh for host and host + H), no zero-point energy. The isolated-atom reservoir would shift every value by
-2.312 eV (the H2 binding energy per H) and is not the equilibrium reference for dissolution from H2 gas.
The archive labels every run an "inequivalent site", but the relaxed outputs show that many runs end in the same
configuration (energies equal to < 1 meV, identical coordination shells). Runs were merged only when both agree;
distinct configurations closer than 0.06 eV (TiN, HfB2) differ in coordination and stay separate. The graphite run 39
stopped at the BFGS step limit inside the basin of run 17 and is counted with it.
{chr(10).join(lines)}
Site competition is not textbook: in hcp Re the octahedral site beats the tetrahedral one by 0.56 eV; in TaC the H
binds to C (C-H 1.15 A) rather than sitting in the Ta4 void (0.34 eV higher); in TiN the Ti4 tetrahedral site wins by
only 0.05 eV over the N-H configuration. Gaps of ~0.05 eV (TiN, HfB2) are within typical site-to-site differences of the
H zero-point energy, so the site order there is not robust; the archive has no vibrational data, and supercell-size
convergence (single cell per host, fixed lattice, H-cell pressures up to ~20 kbar) cannot be checked either.
Of the graded hosts only Re dissolves H exothermically from H2 at this level (and only by 0.08 eV, inside the missing
zero-point correction); graphite, the carbides, the nitride and the boride are clearly endothermic.

## Surface binding energies
E_sb = E(atom) + E(slab with the atom removed, unrelaxed) - E(slab), on the termination of the authors' facet with the
lowest surface energy (terminations: {terms}), smallest value per element. Isolated atoms must use the same cutoff as
the slab: the W atom of the W set is at 30/240 Ry and that of the W-compound sets at 45/360 Ry, 0.13 eV apart, so reusing
one W atom across materials shifts W SBEs by 0.13 eV (the published compound table appears to do this). Notes: {notes}.
The TiB2 Ti atom is a byte-identical copy of the TiN Ti atom (60/480 Ry) while the TiB2 slabs use 35/280 Ry; the B
atoms (35/50/60 Ry) and the TiB2 cohesive energy computed at 60 Ry (H-IFE supercell) versus with the mixed set differ by
only a few meV, so the Ti value is reported with that caveat.
SBE values: {S}.

## Cohesive energies and the sublimation proxy
Cohesive energies use the relaxed SBE unit cell and atoms at the same cutoff: {C}. Mixing the 60/480 Ry H-IFE
supercell with the lower-cutoff atoms would shift W by +0.21 eV and Re by +0.06 eV (WC by +0.01 eV).
For the bcc/hcp metals the flat-surface SBE is 1.2-1.35 times the cohesive energy (unrelaxed vacancy in a complete
layer), and in the compounds the element- and termination-specific SBEs spread from ~5 to ~17 eV around cohesive
energies of 8-9 eV, so the cohesive (sublimation) energy is not a reliable SBE proxy; termination and the removed atom
matter more.

## What these data do not support
No zero-point or vibrational free energies, hence no site order or solubility ranking at operating temperature for
near-degenerate sites; no supercell-size, spin (all bulk runs are non-magnetic, including Cr) or cell-relaxation
checks for the interstitials; relaxed surface vacancies are not available (SBEs are unrelaxed upper bounds).
""")
    lab('submit', submission)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
