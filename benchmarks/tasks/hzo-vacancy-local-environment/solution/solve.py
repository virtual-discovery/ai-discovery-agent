"""Reference solution for hzo-vacancy-local-environment; runs as the Harbor oracle and in check_task."""
import json
import math
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
KT = 8.617333262e-5 * 700.0
MU_POOR = CATALOG['oxygen_chemical_potential_O_poor_eV']
CLASSES = [(4, k) for k in range(5)] + [(3, k) for k in range(4)]     # (first-shell coordination, n_Zr)
# Purchases at x = 0.25 and 0.75: matched environments spanning both sublattices and the Zr counts present,
# weighted towards the classes that dominate the population; first k purchasable sites of each class.
PLAN = {'x0.25': [((4, 0), 1), ((4, 1), 2), ((4, 2), 2), ((4, 3), 1), ((3, 0), 1), ((3, 1), 1), ((3, 2), 1)],
        'x0.75': [((4, 4), 3), ((4, 3), 2), ((4, 2), 1), ((3, 1), 1), ((3, 2), 1), ((3, 3), 1)]}
LDA_SITES = [33, 38, 51, 68, 85, 41, 46, 58, 75, 92]   # labels of CASTEP cells are only known after purchase


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def parse_cif(path):
    lat, sp, frac = {}, [], []
    for line in Path(path).read_text().splitlines():
        s = line.split()
        if line.startswith('_cell_length_'):
            lat[s[0]] = float(s[1])
        elif len(s) == 7 and s[0] in ('Hf', 'Zr', 'O'):
            sp.append(s[0])
            frac.append([float(v) for v in s[3:6]])
    return np.diag([lat['_cell_length_a'], lat['_cell_length_b'], lat['_cell_length_c']]), sp, np.array(frac)


def environments(cell, sp, frac):
    cat = [i for i, s in enumerate(sp) if s != 'O']
    out = {}
    for i, s in enumerate(sp):
        if s == 'O':
            d = frac[cat] - frac[i]
            d -= np.round(d)
            r = np.linalg.norm(d @ cell, axis=1)
            nb = [cat[j] for j in np.where(r <= 2.3)[0]]
            out[i] = (len(nb), sum(sp[j] == 'Zr' for j in nb))
    return out


def features(cls, x):
    n, z = cls
    o3 = float(n == 3)
    return [1.0, o3, o3 * z, (1 - o3) * z, x]


def binom(x, cls):
    n, k = cls
    return 0.5 * math.comb(n, k) * x ** k * (1 - x) ** (n - k)


def main():
    env = {c: environments(*parse_cif(APP/v['cif'])) for c, v in CATALOG['supercells'].items()}
    xval = {c: v['x'] for c, v in CATALOG['supercells'].items()}
    evidence, buy = [], {}

    def purchase(req):
        r = lab('evaluate', req)
        while r['status'] == 'running':
            time.sleep(1)
            r = lab('job/%d' % r['id'])
        if r['status'] != 'complete':
            raise RuntimeError(r)
        evidence.append(r['id'])
        return r['result']

    # 1. Pin the oxygen chemical potential of the released x = 0.50 table with one vacancy + pristine purchase.
    table = list(__import__('csv').DictReader(open(APP/'x050_site_table.csv')))
    e_bulk = {c: purchase({'kind': 'pristine', 'composition': c})['total_energy_eV'] for c in ('x0.50', 'x0.25', 'x0.75')}
    probe = int(table[0]['O_index'])
    e_probe = purchase({'kind': 'vacancy', 'composition': 'x0.50', 'site': probe})['total_energy_eV']
    mu_table = float(table[0]['E_vac']) - (e_probe - e_bulk['x0.50'])
    rows = [(0.5, env['x0.50'][int(r['O_index'])], float(r['E_vac']) - mu_table + MU_POOR) for r in table]
    # 2. Matched environments at the other compositions.
    for c, plan in PLAN.items():
        avail = CATALOG['vacancy_sites_available'][c]
        for cls, k in plan:
            for s in [s for s in avail if env[c][s] == cls][:k]:
                e = purchase({'kind': 'vacancy', 'composition': c, 'site': s})['total_energy_eV'] - e_bulk[c] + MU_POOR
                buy[(c, s)] = e
                rows.append((xval[c], env[c][s], e))
    # 3. LDA (CASTEP) vacancies: relative energies within one functional, with the authors' labels of each calculation.
    lda = [purchase({'kind': 'castep_vacancy', 'functional': 'LDA', 'site': s}) for s in LDA_SITES]
    mismatch = sum((4 if r['vacancy_type'] == 'O4' else 3, r['n_Zr']) != env['x0.50'][r['site']] for r in lda)
    A = np.array([[1.0, float(r['vacancy_type'] == 'O3'), r['n_Zr']] for r in lda])
    _, d_lda, b_lda = np.linalg.lstsq(A, np.array([r['final_enthalpy_eV'] for r in lda]), rcond=None)[0]
    # 4. Local model with sublattice offset, per-sublattice Zr slopes and a linear composition term.
    A = np.array([features(cls, x) for x, cls, e in rows])
    y = np.array([e for x, cls, e in rows])
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    resid = y - A @ coef
    sigma = float(np.sqrt(np.mean(resid ** 2)))
    energy = lambda c, cls: float(np.dot(coef, features(cls, xval[c])))
    pred = {p['id']: energy(p['composition'], env[p['composition']][p['site']]) for p in CATALOG['sites_to_predict']}
    # 5. Equilibrium populations in a macroscopic random solid solution (binomial environment probabilities x
    #    Boltzmann factors); Gaussian within-class scatter raises <exp(-E/kT)> by exp(sigma^2 / 2 kT^2).
    dist, frac, o3 = {}, {}, {}
    for c in ('x0.25', 'x0.75'):
        w = {cls: binom(xval[c], cls) * math.exp(-energy(c, cls) / KT) for cls in CLASSES}
        z = sum(w.values())
        dist[c] = [sum(v for cls, v in w.items() if cls[1] == k) / z for k in range(5)]
        frac[c] = z * math.exp(sigma ** 2 / (2 * KT ** 2))
        o3[c] = sum(v for cls, v in w.items() if cls[0] == 3) / z
    cell_counts = {c: sum(1 for s, cls in env[c].items() if cls == (4, 4)) for c in env}
    submission = {'predicted_formation_energies_eV': pred, 'vacancy_Zr_neighbour_distribution': dist,
                  'vacancy_site_fraction': frac, 'lda_common_slope_model': {'offset_O3_eV': float(d_lda), 'slope_per_Zr_eV': float(b_lda)},
                  'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# Oxygen vacancies in Hf1-xZrxO2 (reference solution)

Conventions: E_f = E(vacancy cell) - E(pristine cell, same composition) + mu_O, mu_O = {MU_POOR} eV (O-poor);
first-shell cations within 2.3 A of the oxygen in the staged structures.

The released x = 0.50 table carries O-rich energies: one purchased vacancy and pristine cell give
mu_O(table) = {mu_table:.2f} eV, so the table was shifted by {MU_POOR - mu_table:.2f} eV to the O-poor limit.

What controls E_f. Two oxygen sublattices: threefold (O3; bonds 2.03-2.13 A) and fourfold (O4; 2.11-2.25 A).
Model fitted to the 64 x = 0.50 sites plus {len(buy)} matched sites at x = 0.25 and 0.75:
E_f = {coef[0]:.3f} + {coef[1]:.3f} [O3] + ({coef[2]:.3f} n_Zr on O3) + ({coef[3]:.3f} n_Zr on O4) + ({coef[4]:.3f} x) eV,
residual {sigma:.3f} eV. O3 vacancies cost ~{coef[1]:.2f} eV more at equal Zr count (shorter, stronger bonds); each Zr
neighbour lowers E_f by ~0.08-0.10 eV (weaker Zr-O bonding). At fixed first-shell environment E_f falls by
~{-coef[4] * 0.25:.3f} eV per 0.25 increase of x, a composition effect that pooled, composition-blind fits misassign to
the Zr-count slope because each cell's environment mix follows its composition. All cells share one lattice, so
this shift cannot be separated from a fixed-volume strain or electrostatic effect with these data. The geometric
descriptors of the ideal cells (M-O distances, M-O-M angle groups) are the same for every site of a sublattice:
they reproduce the sublattice split but cannot explain the Zr dependence.

Populations (700 K, dilute, random solid solution): environment probabilities are binomial, not those of the
96-atom cells (the x = 0.25 cell has {cell_counts['x0.25']} 4-Zr fourfold site of 64; a random alloy has 0.2 %), times
Boltzmann factors. Vacancy n_Zr distribution: x = 0.25 {[round(v, 3) for v in dist['x0.25']]},
x = 0.75 {[round(v, 3) for v in dist['x0.75']]}; O3 share {o3['x0.25']:.4f} / {o3['x0.75']:.4f}; site fraction
{frac['x0.25']:.3g} / {frac['x0.75']:.3g}. Rare Zr-rich sites hold a disproportionate share at low x.

Functional check (CASTEP-LDA, relative energies only; no pristine or reference energies exist). The authors' labels
of the LDA calculations differ from the x = 0.50 VASP cell's environments at {mismatch} of {len(lda)} sites, so the
CASTEP cells are separate decorations and were classified by their own labels. LDA: offset {d_lda:.3f} eV,
slope {b_lda:.3f} eV per Zr (common-slope description). The Zr softening is robust to the functional; the
O3-O4 offset is not (PBE ~0.25 eV), so the small O3 share of vacancies is functional-sensitive.

Not supported: charged vacancies (neutral only), more than one decoration per composition, a mechanism for the
composition shift, and supercell-size convergence (the 3x3x3 data were not purchased here).
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
