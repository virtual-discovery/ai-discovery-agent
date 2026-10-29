"""Reference solution for fresh-inp-polar-phonons; runs as the Harbor oracle and in check_task.

Workflow: equation of state on a converged shifted mesh -> a0; k-point ladder of the electric-field response
at a0 (Gamma-centred meshes are avoided: they contain the Gamma point, where the PBE gap is tiny) -> converged
eps_inf and Z*; phonons at Gamma on a coarser mesh (the TO frequency converges fast); LO from the nonanalytic
term with the converged tensors; centred volume differences; the strained cell with its own tensors.
"""
import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
RY_CM = 109737.31568160          # 1 Ry in cm^-1
AMU_RY = 911.44424310865645      # 1 amu in Rydberg mass units (m_e / 2)
BOHR_A = 0.529177210903
RY_EV = 13.605693122994


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

def asr_charges(result):
    """Born charge tensors with the charge-neutrality (acoustic) sum rule imposed by removing the mean."""
    z = np.array([result['born_charges_raw']['In'], result['born_charges_raw']['P']])
    return z - z.mean(axis=0)


def asr_force_constants(result):
    phi = np.array(result['force_constants_gamma_Ry_per_bohr2'])
    fixed = phi.copy()
    for i in range(2):
        for a in range(3):
            for b in range(3):
                fixed[3 * i + a, 3 * i + b] -= sum(phi[3 * i + a, 3 * j + b] for j in range(2))
    return 0.5 * (fixed + fixed.T)


def modes(result, q=None, zstar=None, eps=None):
    """Zone-centre frequencies (cm^-1) and eigenvectors; with q, adds the nonanalytic (macroscopic-field) term
    4 pi e^2/Omega (q.Z_i)(q.Z_j)/(q.eps.q) (Rydberg units, e^2 = 2)."""
    phi = asr_force_constants(result)
    if q is not None:
        zstar = asr_charges(result) if zstar is None else np.array(zstar)
        eps = np.array(result['dielectric_tensor_eps_inf']) if eps is None else np.array(eps)
        q = np.array(q, float) / np.linalg.norm(q)
        qz = np.concatenate([q @ zstar[i] for i in range(2)])
        omega = result['cell_volume_A3'] / BOHR_A ** 3
        phi = phi + 4 * np.pi * 2.0 / omega * np.outer(qz, qz) / (q @ eps @ q)
    m = np.repeat([result['atomic_masses_amu']['In'] * AMU_RY, result['atomic_masses_amu']['P'] * AMU_RY], 3)
    w2, vec = np.linalg.eigh(phi / np.sqrt(np.outer(m, m)))
    return np.sign(w2) * np.sqrt(np.abs(w2)) * RY_CM, vec


def to_lo(result, zstar=None, eps=None):
    """Cubic cell: (omega_TO, omega_LO) with the LO along z."""
    w_to = modes(result)[0][-1]
    w_lo = modes(result, (0, 0, 1), zstar, eps)[0][-1]
    return float(w_to), float(w_lo)


def strained_lo(result):
    """Tetragonal cell: the LO for q along z and along x, each the highest mode including its own
    nonanalytic term (the field couples only to the TO component polarized along q)."""
    return float(modes(result, (0, 0, 1))[0][-1]), float(modes(result, (1, 0, 0))[0][-1])


def birch_murnaghan_a0(lattice, energy_ry):
    """Third-order Birch-Murnaghan fit of E(V) (V = a^3/4); returns a0 and B0 (GPa)."""
    v = np.array(lattice) ** 3 / 4.0
    x = v ** (-2.0 / 3.0)
    c = np.polyfit(x, energy_ry, 3)
    roots = np.roots(np.polyder(c))
    x0 = [r.real for r in roots if abs(r.imag) < 1e-12 and np.polyval(np.polyder(c, 2), r.real) > 0
          and x.min() * 0.9 < r.real < x.max() * 1.1][0]
    v0 = x0 ** -1.5
    d2 = np.polyval(np.polyder(c, 2), x0) * (2.0 / 3.0 * v0 ** (-5.0 / 3.0)) ** 2 \
        + np.polyval(np.polyder(c, 1), x0) * (10.0 / 9.0 * v0 ** (-8.0 / 3.0))
    b0 = v0 * d2 * RY_EV * 160.21766208      # V d2E/dV2 in Ry/A^3 -> eV/A^3 -> GPa
    return float((4 * v0) ** (1.0 / 3.0)), float(b0)


def main():
    evidence = []
    # 1. Equation of state (SCF energies converge on coarse meshes; 6x6x6 shifted is ample).
    grid = [5.88, 5.92, 5.96, 6.00, 6.04]
    eos = [compute({'kind': 'scf', 'lattice_A': a, 'kmesh': 6, 'kshift': 1}, evidence) for a in grid]
    a0, b0 = birch_murnaghan_a0(grid, [r['total_energy_Ry'] for r in eos])
    a0 = round(a0, 4)
    # 2. Band edges on a Gamma-centred mesh: the direct PBE gap at Gamma (shifted meshes never sample Gamma).
    gap_gamma = compute({'kind': 'scf', 'lattice_A': a0, 'kmesh': 4, 'kshift': 0}, evidence)['band_gap_on_mesh_eV']
    # 3. k ladder of the electric-field response at a0 (shifted meshes only).
    ladder = {k: compute({'kind': 'dielectric', 'lattice_A': a0, 'kmesh': k, 'kshift': 1}, evidence) for k in (4, 6, 8)}
    eps_k = {k: r['dielectric_tensor_eps_inf'][0][0] for k, r in ladder.items()}
    z_k = {k: asr_charges(r)[0][0][0] for k, r in ladder.items()}
    converged = abs(eps_k[8] - eps_k[6]) < 0.05 and abs(z_k[8] - z_k[6]) < 0.01
    eps_inf, z_in = eps_k[8], z_k[8]
    # 4. Phonons at a0 (6x6x6 shifted; the TO converges on coarse meshes) and the LO with the converged tensors.
    ph0 = compute({'kind': 'dfpt', 'lattice_A': a0, 'kmesh': 6, 'kshift': 1}, evidence)
    zt = np.array([np.eye(3) * z_in, -np.eye(3) * z_in])
    w_to, w_lo = to_lo(ph0, zt, np.eye(3) * eps_inf)
    eps_0 = eps_inf * (w_lo / w_to) ** 2
    # 5. Centred hydrostatic differences (+-1% in a, about +-3% in V) with one mesh (6x6x6 shifted DFPT)
    #    at all three volumes, so that the residual k-point error cancels in the differences.
    side = {}
    for s in (-1, 1):
        side[s] = compute({'kind': 'dfpt', 'lattice_A': round(a0 * (1 + 0.01 * s), 4), 'kmesh': 6, 'kshift': 1}, evidence)
    dlnv = np.log(side[1]['cell_volume_A3'] / side[-1]['cell_volume_A3'])

    def q(r):
        to, lo = to_lo(r)
        return {'Z': asr_charges(r)[0][0][0], 'eps': r['dielectric_tensor_eps_inf'][0][0], 'TO': to, 'S': lo - to}
    lo_side, hi_side = q(side[-1]), q(side[1])
    deriv = {k: float(np.log(hi_side[k] / lo_side[k]) / dlnv) for k in ('Z', 'eps', 'TO', 'S')}
    # 6. [001]-strained cell with its own tensors.
    st = compute({'kind': 'dfpt', 'lattice_A': a0, 'strain_zz': -0.02, 'kmesh': 6, 'kshift': 1}, evidence)
    lo_z, lo_x = strained_lo(st)
    w_st = sorted(modes(st)[0][3:])
    submission = {'a0_A': a0, 'omega_TO_cm-1': w_to, 'omega_LO_cm-1': w_lo, 'born_charge_In': float(z_in),
                  'eps_inf': float(eps_inf), 'eps_0': float(eps_0),
                  'dlnZ_dlnV': deriv['Z'], 'dln_eps_inf_dlnV': deriv['eps'], 'gruneisen_TO': -deriv['TO'],
                  'LOTO_splitting_under_compression': 'grows' if deriv['S'] < 0 else 'shrinks',
                  'strained_omega_LO_q_z_cm-1': lo_z, 'strained_omega_LO_q_x_cm-1': lo_x, 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    split2 = 2 * deriv['Z'] - deriv['eps'] - 1
    eps_exp = 9.61
    lo_exp_eps = float(np.sqrt(w_to ** 2 + (w_lo ** 2 - w_to ** 2) * eps_inf / eps_exp))
    (APP/'report.md').write_text(f"""# InP polar phonons in DFPT-PBE (reference solution)

**Equilibrium.** Five SCF energies (6x6x6 shifted mesh, 50 Ry) fitted with a third-order Birch-Murnaghan form give
a0 = {a0:.4f} A and B0 = {b0:.1f} GPa (static lattice; the experimental 5.869 A is ~1.5% smaller, the usual PBE
overestimate). All graded quantities are evaluated at this a0, not at the experimental lattice parameter.

**k-point convergence.** Total energies and the TO frequency converge on coarse meshes, but eps_inf and Z* do not.
Gamma-centred meshes contain the Gamma point, where the PBE gap is small ({gap_gamma:.3f} eV on the 4x4x4 Gamma-centred
mesh at a0, versus >1 eV on the shifted meshes, which never sample Gamma), and overweight it: eps_inf comes out at
tens to ~10^2 and decays only slowly with N, with a large violation of charge neutrality. Shifted meshes:
eps_inf = {eps_k[4]:.3f}, {eps_k[6]:.3f}, {eps_k[8]:.3f} and Z*(In, sum rule imposed) = {z_k[4]:.4f}, {z_k[6]:.4f},
{z_k[8]:.4f} for 4, 6, 8 (converged: {converged}). The 4x4x4 values are wrong by ~0.5 in eps_inf and ~0.05 in Z*,
yet the LO frequency barely moves because the errors cancel in Z*^2/eps_inf. I use the 8x8x8 tensors.

**LO and eps_0.** The analytic (ph.x) Gamma frequency is the TO, {w_to:.2f} cm^-1. The LO adds the macroscopic-field
term 4 pi e^2 (q.Z*)^2/(Omega mu q.eps_inf.q) (clamped-ion eps_inf, not eps_0), giving {w_lo:.2f} cm^-1 (splitting
{w_lo - w_to:.2f} cm^-1). Lyddane-Sachs-Teller: eps_0 = eps_inf (w_LO/w_TO)^2 = {eps_0:.3f}. Sum rules were imposed on
Z* (mean removed) and on the force constants (acoustic modes at zero).

**Compression.** Centred differences at a0(1 +- 1%) on one mesh: dlnZ*/dlnV = {deriv['Z']:.3f}, dln eps_inf/dlnV =
{deriv['eps']:.3f}, gamma_TO = {-deriv['TO']:.3f}, dln(w_LO - w_TO)/dlnV = {deriv['S']:.3f}. Since
w_LO^2 - w_TO^2 ~ Z*^2/(Omega eps_inf mu), dln(w_LO^2 - w_TO^2)/dlnV = 2 dlnZ* - dln eps_inf - 1 = {split2:.3f}:
Z* falls under compression (covalency increases), which shrinks the splitting; eps_inf falls (the gap opens), which
enlarges it; the 1/Omega factor enlarges it. The TO stiffens (gamma_TO > 0), which reduces w_LO - w_TO for a given
w_LO^2 - w_TO^2. Net: the splitting {'grows' if deriv['S'] < 0 else 'shrinks'} under compression. A Z*-only argument
(2 dlnZ* - 1 = {2 * deriv['Z'] - 1:.2f} > 0) would predict the opposite; the fall of eps_inf dominates. The net
d ln(w_LO - w_TO)/d ln V is a near-cancellation and is sensitive to the k mesh (it drifts with the mesh because the
expanded cell's smaller gap converges more slowly), so only its sign is robust.

**Experiment.** Experiment (room temperature) gives TO/LO ~304/345 cm^-1 and eps_inf ~9.6. This protocol's splitting
is {w_lo - w_to:.1f} cm^-1. Replacing the PBE eps_inf by 9.61 while keeping Z*, the TO and Omega gives an LO of
{lo_exp_eps:.1f} cm^-1 (splitting {lo_exp_eps - w_to:.1f} cm^-1), so the PBE overestimate of eps_inf (gap
underestimate) accounts for most of the too-small splitting; the soft TO is mostly the expanded PBE volume.

**Strain.** At strain_zz = -0.02 the TO triplet splits into a z-polarized singlet and an xy doublet
({w_st[0]:.2f}, {w_st[1]:.2f}, {w_st[2]:.2f} cm^-1), and Z* and eps_inf become tetragonal tensors
(Z*_zz = {asr_charges(st)[0][2][2]:.4f}, Z*_xx = {asr_charges(st)[0][0][0]:.4f}, eps_zz =
{st['dielectric_tensor_eps_inf'][2][2]:.3f}, eps_xx = {st['dielectric_tensor_eps_inf'][0][0]:.3f}). The macroscopic field
couples only to the TO component along q, so q || z dresses the singlet (LO {lo_z:.2f} cm^-1) and q || x dresses one
doublet component (LO {lo_x:.2f} cm^-1).

**Limits.** These are static-lattice, harmonic, PBE numbers at one pseudopotential set; they do not include
anharmonic or zero-point shifts, the true (larger) gap, or the zincblende-rocksalt transition at higher pressure.
Volume derivatives are local (a0 +- 1%); the trend is not established far from a0.
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
