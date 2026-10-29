"""Reference solution for fresh-li-bcc-soft-shear; runs as the Harbor oracle and in check_task.

Workflow: equation of state on a mesh where pressures are converged -> a0 and B; k-mesh/smearing ladder of the
shear response itself (symmetric +-1 % volume-conserving tetragonal strains: the Bain path has a large
third-order term, so one-sided strains are biased) -> converged C'; symmetric +-1 % trigonal [111] shear -> C44;
the compressed cell with the same strains, using stress-strain slopes (energy curvatures of volume-changing
strains include the -P dV work).
"""
import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
RY_EV = 13.605693122994
EVA3_GPA = 160.21766208


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def compute(request, log):
    """Submit one calculation and poll it to completion."""
    record = lab('evaluate', request)
    while record['status'] == 'running':
        time.sleep(5)
        record = lab('job/%d' % record['id'])
    if record['status'] != 'complete':
        raise RuntimeError('calculation %s failed: %s' % (request, record.get('error')))
    log.append(record['id'])
    print(json.dumps({'id': record['id'], 'request': request, 'charge': record['charge'],
                      'cumulative': record['cumulative_charge'], 'wall_s': record['result']['wall_seconds']}), flush=True)
    return record['result']


# ---------------------------------------------------------------- analysis helpers (pure numpy)

def cauchy(result):
    """Cauchy stress with tension positive (the lab reports Quantum ESPRESSO's sign, positive = compressive)."""
    return -np.array(result['stress_GPa'])


def tetragonal_strain(d):
    """Volume-conserving tetragonal strain (d, d, e3) with (1 + d)^2 (1 + e3) = 1."""
    return [d, d, round(1.0 / (1.0 + d) ** 2 - 1.0, 8), 0, 0, 0]


def trigonal_strain(e):
    return [0, 0, 0, e, e, e]


def cprime_stress(plus, minus):
    """C' (or B' under pressure) from the symmetric difference of the stress anisotropy sigma_xx - sigma_zz."""
    aniso = lambda r: cauchy(r)[0][0] - cauchy(r)[2][2]
    lever = lambda r: r['strain_voigt'][0] - r['strain_voigt'][2]
    return float((aniso(plus) - aniso(minus)) / (2.0 * (lever(plus) - lever(minus))))


def cprime_one_sided(strained):
    """One-sided C' from a single strained cell (biased by the third-order term)."""
    s, e = cauchy(strained), strained['strain_voigt']
    return float((s[0][0] - s[2][2]) / (2.0 * (e[0] - e[2])))


def shear_stress(result):
    s = cauchy(result)
    return float(np.mean([s[1][2], s[0][2], s[0][1]]))


def c44_stress(plus, minus):
    """C44 (B44 under pressure) from the symmetric difference of the trigonal shear stress (engineering e4)."""
    return float((shear_stress(plus) - shear_stress(minus)) / (plus['strain_voigt'][3] - minus['strain_voigt'][3]))


def deviatoric_norm(voigt):
    e = np.array(voigt[:3], float)
    return float((e ** 2).sum() - e.sum() ** 2 / 3.0)


def cprime_energy(unstrained, plus, minus):
    """C' from the symmetric second difference of F for volume-conserving tetragonal strains:
    dF/V = C' (sum e_i^2 - (sum e_i)^2 / 3) to second order (exact volume conservation removes the P dV term)."""
    v = unstrained['cell_volume_A3']
    d2 = plus['free_energy_Ry'] + minus['free_energy_Ry'] - 2.0 * unstrained['free_energy_Ry']
    return float(d2 * RY_EV / v * EVA3_GPA / (deviatoric_norm(plus['strain_voigt']) + deviatoric_norm(minus['strain_voigt'])))


def c44_energy_linear(unstrained, plus, minus):
    """Energy curvature of the (volume-changing) trigonal shear: dF/V = (3/2) X e^2 with X = B44 + P/2."""
    e = plus['strain_voigt'][3]
    v = unstrained['cell_volume_A3']
    d2 = plus['free_energy_Ry'] + minus['free_energy_Ry'] - 2.0 * unstrained['free_energy_Ry']
    return float(d2 * RY_EV / v * EVA3_GPA / (3.0 * e * e))


def birch_murnaghan(lattice, energy_ry):
    """Third-order Birch-Murnaghan fit of F(V), V = a^3/2 (bcc); returns a0 (A) and B0 (GPa)."""
    v = np.array(lattice, float) ** 3 / 2.0
    x = v ** (-2.0 / 3.0)
    c = np.polyfit(x, energy_ry, 3)
    roots = np.roots(np.polyder(c))
    x0 = [r.real for r in roots if abs(r.imag) < 1e-12 and np.polyval(np.polyder(c, 2), r.real) > 0
          and x.min() * 0.9 < r.real < x.max() * 1.1][0]
    v0 = x0 ** -1.5
    d2 = np.polyval(np.polyder(c, 2), x0) * (2.0 / 3.0 * v0 ** (-5.0 / 3.0)) ** 2 \
        + np.polyval(np.polyder(c, 1), x0) * (10.0 / 9.0 * v0 ** (-8.0 / 3.0))
    return float((2.0 * v0) ** (1.0 / 3.0)), float(v0 * d2 * RY_EV * EVA3_GPA)


def main():
    evidence = []
    s_main = 0.01
    # 1. Equation of state. Pressures and energies converge on modest meshes (checked against k48 in construction).
    grid = [3.39, 3.41, 3.43, 3.45, 3.47]
    eos = [compute({'lattice_A': a, 'kmesh': 32, 'smearing_Ry': s_main}, evidence) for a in grid]
    a0, b0 = birch_murnaghan(grid, [r['free_energy_Ry'] for r in eos])
    a0 = round(a0, 4)
    p_eos = [r['pressure_GPa'] for r in eos]
    # 2. Shear-response ladder at a0: symmetric +-1 % volume-conserving tetragonal strains, stress route.
    def tet_pair(a, k, s):
        return [compute({'lattice_A': a, 'strain_voigt': tetragonal_strain(d), 'kmesh': k, 'smearing_Ry': s}, evidence)
                for d in (0.01, -0.01)]
    ladder = {}
    for k, s in [(40, 0.01), (48, 0.01), (48, 0.005)]:
        ladder[(k, s)] = tet_pair(a0, k, s)
    cp = {key: cprime_stress(*pair) for key, pair in ladder.items()}
    one_sided = {key: (cprime_one_sided(pair[0]), cprime_one_sided(pair[1])) for key, pair in ladder.items()}
    c_prime = float(np.mean([cp[(48, 0.01)], cp[(48, 0.005)]]))
    # 3. C44 from the trigonal [111] shear (12 symmetry operations: cheaper than a single e6 shear).
    s_c44 = 0.005
    tri = [compute({'lattice_A': a0, 'strain_voigt': trigonal_strain(e), 'kmesh': 48, 'smearing_Ry': s_c44}, evidence)
           for e in (0.01, -0.01)]
    c44 = c44_stress(*tri)
    # 4. Compressed cell a = 3.30 A: pressure, and stress-strain slopes about the pressurized cell.
    ac = 3.30
    cub_c = compute({'lattice_A': ac, 'kmesh': 48, 'smearing_Ry': s_main}, evidence)
    tet_c = tet_pair(ac, 48, s_main)
    tri_c = [compute({'lattice_A': ac, 'strain_voigt': trigonal_strain(e), 'kmesh': 48, 'smearing_Ry': s_main},
                     evidence) for e in (0.01, -0.01)]
    p_c = cub_c['pressure_GPa']
    bp_c, b44_c = cprime_stress(*tet_c), c44_stress(*tri_c)
    cp_energy_c = cprime_energy(cub_c, *tet_c)
    c44_energy_c = c44_energy_linear(cub_c, *tri_c)
    submission = {'a0_A': a0, 'bulk_modulus_GPa': b0, 'Cprime_GPa': c_prime, 'C44_GPa': c44,
                  'zener_anisotropy': c44 / c_prime, 'mechanically_stable_at_a0': bool(c_prime > 0 and c44 > 0 and b0 > 0),
                  'compressed_pressure_GPa': p_c, 'compressed_Bprime_GPa': bp_c, 'compressed_B44_GPa': b44_c,
                  'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    c11, c12 = b0 + 4 * c_prime / 3, b0 - 2 * c_prime / 3
    fmt = lambda d: ', '.join(f'k{k}/{s} Ry: {v:.3f}' for (k, s), v in d.items())
    (APP/'report.md').write_text(f"""# bcc Li elastic constants in PBE (reference solution)

**Equilibrium.** Five SCF free energies (k32, MV 0.01 Ry, 40 Ry) between 3.39 and 3.47 A, Birch-Murnaghan fit:
a0 = {a0:.4f} A, B = {b0:.2f} GPa (pressures {', '.join(f'{p:.3f}' for p in p_eos)} GPa change sign at the same
a0). Total energies and pressures are converged already on coarse meshes (in the construction ladder k16 with 0.02 Ry
gives the same pressure to 0.005 GPa), so they say nothing about the convergence of the shear response. All strains are applied about a0,
not about the experimental lattice parameter (3.48 A at 78 K), which is under about -0.5 GPa in this protocol.

**C' (tetragonal shear).** The energy of a 1 % tetragonal strain is only ~0.1 meV/atom (stress anisotropy
~0.05 GPa), set by Fermi-surface states near the N points. I use volume-conserving tetragonal strains
(d, d, (1+d)^-2 - 1) and the stress anisotropy sigma_xx - sigma_zz (dsigma = 2 C' d(e1 - e3)). The Bain path is
strongly asymmetric: one-sided values at +1 % / -1 % are {fmt({k: v[0] for k, v in one_sided.items()})} /
{fmt({k: v[1] for k, v in one_sided.items()})} GPa, so only symmetric +-1 % differences are used (the odd term
cancels). Symmetric C': {fmt(cp)}. In the construction ladder the coarse level whose total energy is converged (k16,
0.02 Ry) gives ~1.5 GPa, small smearing on coarse meshes gives noisy values, even negative ones, and 0.02 Ry is
k-converged by k32 but about 10 % low. I take C' = {c_prime:.3f} GPa from k48 at 0.01 and 0.005 Ry, which agree.

**C44.** Trigonal [111] shear e4 = e5 = e6 = +-1 % (12 symmetry operations, much cheaper than a single e6 shear),
shear stress sigma_4 = C44 e4 with engineering strain (dF/V = (3/2) C44 e^2 for the three components):
C44 = {c44:.3f} GPa (k48, {s_c44} Ry; C44 also rises by a few tenths of a GPa from 0.02 to 0.005 Ry). Zener anisotropy A = C44/C' = {c44 / c_prime:.2f}; C11 = {c11:.2f},
C12 = {c12:.2f} GPa from B and C'.

**Stability.** Born criteria C11 - C12 > 0, C44 > 0, C11 + 2C12 > 0 are all met: bcc Li is mechanically
(elastically) stable at T = 0 in this protocol, but only barely along the Bain path. This does not decide the
thermodynamic competition with fcc, hcp and 9R, which are within meV/atom and need total-energy and zero-point
calculations not done here. The small C' is the near-cancellation of the electrostatic (Ewald) stiffness by the
band-structure (Fermi-surface) term of a nearly free-electron metal whose Fermi surface approaches the N-point
zone boundary; C12 - C44 = {c12 - c44:.2f} GPa != 0 (Cauchy violation) reflects non-central, volume-dependent
metallic bonding.

**Compression (a = 3.30 A).** P = {p_c:.3f} GPa. With a pre-stress the energy curvature of a strain that
changes the volume contains the -P dV work: the trigonal energy curvature gives {c44_energy_c:.3f} GPa = B44 + P/2,
not B44. Stress-strain slopes about the compressed cell give the stability coefficients directly:
B' = {bp_c:.3f} GPa (volume-conserving tetragonal energy curvature, which has no P dV term: {cp_energy_c:.3f} GPa) and
B44 = {b44_c:.3f} GPa. Under compression C44 and B stiffen, but the tetragonal shear softens
(B' = {bp_c:.3f} < C' = {c_prime:.3f} GPa): the Bain path gets softer, consistent with the bcc -> close-packed
transitions of Li under pressure. Scaling C' with B, or reusing the a0 value, would get the sign of the change wrong.
The crystal is still elastically stable at 2 GPa.

**Experiment and limits.** Experiment at 78 K gives C' ~ 1.2 GPa, C44 ~ 10.8 GPa, A ~ 9. The converged static PBE C'
is lower and A higher. Static-lattice values omit zero-point motion and thermal expansion (both expand the lattice;
light Li has a large zero-point energy) and anharmonic renormalisation of the soft shear, and smearing-converged C'
depends on Fermi-surface details that PBE need not get right; the comparison is therefore qualitative. These are harmonic, static-lattice PBE numbers
for one pseudopotential, local to the two volumes studied; they do not establish the martensitic transformation
path or temperature.
""")
    lab('submit', submission)
    for _ in range(120):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
