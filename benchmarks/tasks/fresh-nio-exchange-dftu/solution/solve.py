"""Reference solution for fresh-nio-exchange-dftu; runs as the Harbor oracle and in check_task.

It requests new pw.x calculations from the lab (about 35-45 minutes of compute). The estimator functions at
module level are also used by the construction script on the reference calculations.
"""
import json
import os
import subprocess
import time
from pathlib import Path

APP = Path(os.environ.get('APP_DIR', '/app'))
RY_MEV = 13605.693122994
KB_MEV = 0.08617333262
S = 1.0
STANDARD = dict(kshift=0, ecutwfc_Ry=50, smearing='gaussian', degauss_Ry=0.005)
CELL = {'pr': 'primitive_1fu', 'rh': 'rhombohedral_2fu', 'te': 'tetragonal_2fu'}
START = {('pr', 'fm'): [0.5], ('rh', 'fm'): [0.5, 0.5], ('rh', 'af2'): [0.5, -0.5],
         ('te', 'fm'): [0.5, 0.5], ('te', 'af1'): [0.5, -0.5]}


# ---------- estimators (pure functions of lab results) ----------
def energy(res):
    """Zero-smearing estimate of the energy per formula unit (per Ni), meV.

    For Gaussian (and Fermi-Dirac) smearing E(sigma -> 0) ~ (E + F)/2 = F - (-TS)/2; Marzari-Vanderbilt F is
    already accurate to third order. Insulators have TS = 0 and all three coincide.
    """
    free, ts = res['free_energy_Ry_per_cell'], res['smearing_contribution_minus_TS_Ry_per_cell']
    e0 = free if res['request']['smearing'] == 'mv' else free - ts / 2
    return e0 / res['formula_units_per_cell'] * RY_MEV


def exchange(fm_rh, af2_rh, fm_te, af1_te):
    """J1, J2 (meV) for E = E0 + sum_{i<j} J_ij s_i s_j, s = +-1 (S = 1), each pair once.

    Per Ni: E_FM = E0 + 6 J1 + 3 J2, E_AFM-II = E0 - 3 J2, E_AFM-I = E0 - 2 J1 + 3 J2. Each difference is taken
    between two configurations of the same cell and k mesh: J1 = (E_FM - E_AFM-I)/8 from the tetragonal cell,
    J1 + J2 = (E_FM - E_AFM-II)/6 from the rhombohedral cell.
    """
    j1 = (energy(fm_te) - energy(af1_te)) / 8
    j2 = (energy(fm_rh) - energy(af2_rh)) / 6 - j1
    return j1, j2


def neel_mean_field(j2):
    """Molecular-field T_N (K) of AFM-II for the quantum spin-S model: k T_N = S(S+1)/3 * 6 J2 (J1 cancels)."""
    return S * (S + 1) / 3 * 6 * j2 / KB_MEV


def gap_closing(points):
    """points: (U, FM band gap) pairs. U where the FM gap closes: straight line through the two insulating
    points of smallest U, kept inside the bracket set by the largest metallic U below them."""
    ins = sorted((u, g) for u, g in points if g > 0)
    met = [u for u, g in points if g == 0 and u < ins[0][0]]
    (u1, g1), (u2, g2) = ins[0], ins[1]
    uc = u1 - g1 * (u2 - u1) / (g2 - g1)
    if met and not max(met) <= uc <= u1:
        uc = 0.5 * (max(met) + u1)
    return uc


# ---------- lab access ----------
def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


EVIDENCE = []


def scf(U, cell, state, k, mag=None, **settings):
    mesh = [k] * 3 if isinstance(k, int) else list(k)
    request = {'kind': 'scf', 'cell': CELL[cell], 'starting_magnetization': mag or START[(cell, state)],
               'hubbard_U_eV': U, 'kpoints': mesh, **dict(STANDARD, **settings)}
    job = lab('evaluate', request)
    while job['status'] == 'running':
        time.sleep(5)
        job = lab('job/%d' % job['id'])
    if job['status'] != 'complete':
        raise RuntimeError('calculation failed: %s' % job.get('error'))
    EVIDENCE.append(job['id'])
    return job['result']


def main():
    R = {}
    # U = 4 eV (smallest gaps among the graded U): k convergence of each same-cell difference.
    for k in (4, 6):
        R[f'rh4_k{k}'] = (scf(4.0, 'rh', 'fm', k), scf(4.0, 'rh', 'af2', k))
    for mesh in ((8, 8, 6), (10, 10, 8)):
        R[f'te4_k{mesh[0]}'] = (scf(4.0, 'te', 'fm', mesh), scf(4.0, 'te', 'af1', mesh))
    # U = 8 eV with the meshes converged at U = 4 eV (larger gaps converge faster).
    R['rh8'] = (scf(8.0, 'rh', 'fm', 6), scf(8.0, 'rh', 'af2', 6))
    R['te8'] = (scf(8.0, 'te', 'fm', (8, 8, 6)), scf(8.0, 'te', 'af1', (8, 8, 6)))
    # Plain PBE: AFM-II insulator in its cell; the metallic FM in the cheap primitive cell with dense meshes.
    af2_0 = scf(0.0, 'rh', 'af2', 6)
    fm0_k12, fm0_k16 = scf(0.0, 'pr', 'fm', 12), scf(0.0, 'pr', 'fm', 16)
    fm0_m1 = scf(0.0, 'pr', 'fm', 12, mag=[1.0])
    # Where the FM state turns insulating (primitive cell, gap on a 12^3 mesh).
    gaps = [(u, scf(u, 'pr', 'fm', 12)['band_gap_eV']) for u in (2.0, 2.2, 2.5)]

    j4 = exchange(*R['rh4_k6'], *R['te4_k10'])
    j8 = exchange(*R['rh8'], *R['te8'])
    de0 = energy(fm0_k16) - energy(af2_0)
    uc = gap_closing(gaps + [(4.0, R['rh4_k6'][0]['band_gap_eV'])])
    submission = {'J1_meV': {'4.0': j4[0], '8.0': j8[0]}, 'J2_meV': {'4.0': j4[1], '8.0': j8[1]},
                  'neel_temperature_mean_field_K': {'4.0': neel_mean_field(j4[1]), '8.0': neel_mean_field(j8[1])},
                  'pbe_fm_minus_afm2_meV_per_Ni': de0, 'fm_gap_closing_U_eV': uc,
                  'evidence_ids': list(EVIDENCE)}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')

    def d(pair):
        return energy(pair[0]) - energy(pair[1])
    m = lambda res: res['Ni_moments_muB'][0]   # noqa: E731
    (APP/'report.md').write_text(f"""# NiO exchange constants from PBE+U total energies (reference solution)

**Mapping.** E = E0 + sum over Ni pairs (each once) J_ij s_i s_j with s = +-1 (S = 1, S^2 = 1). Per Ni the
fcc Ni sublattice gives E_FM = E0 + 6J1 + 3J2, E_AFM-II = E0 - 3J2 ([111] planes: 6 of 12 nearest neighbours
antiparallel, all 6 next-nearest antiparallel) and E_AFM-I = E0 - 2J1 + 3J2 ([001] planes). FM and AFM-II
alone only give J1 + J2; AFM-I is needed for J1. AFM-II needs the rhombohedral 2-f.u. cell (spins +-), AFM-I
the tetragonal 2-f.u. cell. Energies are per formula unit (the lab reports Ry per cell).

**Same-cell differences.** Each difference is taken between configurations computed in the same cell with the
same k mesh: J1 = (E_FM - E_AFM-I)/8 in the tetragonal cell, J1 + J2 = (E_FM - E_AFM-II)/6 in the rhombohedral
cell. Absolute energies per f.u. from different cells are not interchangeable at moderate meshes (the k error
of the FM energy is several meV per f.u. on a 4^3 rhombohedral mesh while J1 is ~1-2 meV x 8).

**Convergence (U = 4 eV, the smaller gaps).** E_FM - E_AFM-II per Ni: {d(R['rh4_k4']):.2f} (4^3) and
{d(R['rh4_k6']):.2f} meV (6^3). E_AFM-I - E_FM: {-d(R['te4_k8']):.2f} (8x8x6) and {-d(R['te4_k10']):.2f} meV
(10x10x8). Energy differences are insensitive to the cutoff (50 Ry used; 40-75 Ry change them by about 1%),
and a small Gaussian smearing (0.005 Ry) leaves the gapped states untouched; metal-style smearing (0.02 Ry
Fermi-Dirac) would close the 0.77 eV FM gap at U = 4 eV.

**Results.** U = 4 eV: J1 = {j4[0]:.3f}, J2 = {j4[1]:.3f} meV. U = 8 eV: J1 = {j8[0]:.3f}, J2 = {j8[1]:.3f} meV.
J2 (180 deg Ni-O-Ni superexchange) dominates and falls roughly as 1/U (ratio {j4[1] / j8[1]:.2f} for a
doubling of U), as expected for superexchange ~ t^2/(effective charge-transfer/Hubbard gap), with the gap
growing with U. J1 is small and ferromagnetic (90 deg path, Goodenough-Kanamori). U is a protocol choice, not a
fitted truth, so the J's are quoted at fixed U.

**Neel temperature.** Molecular field for AFM-II: k T_N = S(S+1)/3 * 6 J2 (J1 cancels because each Ni has six
parallel and six antiparallel nearest neighbours): {neel_mean_field(j4[1]):.0f} K (U = 4 eV) and
{neel_mean_field(j8[1]):.0f} K (U = 8 eV). Mean field overestimates T_N (quantum Monte Carlo/RPA reduce it by
~20-30% for the fcc lattice), and the exchange striction (rhombohedral distortion) below T_N is not modelled.

**Small U.** Plain PBE (U = 0): AFM-II is an insulator with moment {m(af2_0):.2f} muB; FM is a metal whose
moment depends on the k mesh and smearing ({m(fm0_k12):.3f} muB at 12^3, {m(fm0_k16):.3f} at 16^3), so its energy
must be converged in the cheap primitive cell with small smearing (zero-smearing estimate (E + F)/2). Starting
moments of 0.5 and 1.0 give the same state ({m(fm0_m1):.3f} muB); much smaller starting moments can land in a
low-moment metastable FM solution, so the lowest-energy FM state must be checked. E_FM - E_AFM-II = {de0:.1f} meV per Ni at U = 0. The FM gap opens only at
U = {uc:.2f} eV (FM gaps on a 12^3 mesh: {gaps}); below it the FM reference is (half-)metallic and the
moments are not rigid, so a Heisenberg fit through U = 0 is not meaningful.

**Limits.** Fixed cubic a = 4.17 A, collinear spins, Dudarev U with non-orthogonalized atomic projectors and
these pseudopotentials: the J's are protocol-specific (linear-response U or other projectors would give other
numbers); no spin-orbit coupling, no lattice relaxation, no magnon or T_N beyond mean field.
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
