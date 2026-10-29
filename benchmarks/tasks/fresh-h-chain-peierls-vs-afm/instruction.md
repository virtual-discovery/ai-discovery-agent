# Hydrogen chain: Peierls dimerization versus magnetism in PBE

An infinite chain of hydrogen atoms with mean spacing R has one electron per atom, so its band is half filled. Within spin-polarized PBE, does the chain lower its energy by bond dimerization (Peierls distortion, H₂ formation), by spin-symmetry breaking (magnetic order), or by both? How do the two instabilities compete as R grows from 1.0 to 2.0 Å, and over what range of distortion do they coexist?

## Calculation service

The lab runs new Quantum ESPRESSO 6.7 `pw.x` SCF calculations (PBE, SG15 ONCV norm-conserving hydrogen pseudopotential, density cutoff 4 × the wavefunction cutoff). The chain lies along x in a square box; the cell holds 1, 2 or 4 atoms, and you give the bond lengths in order (atom i sits at the sum of the preceding bond lengths). `GET /task` (`python3 /app/lab_client.py task`) returns the full request schema and the allowed values.

```sh
python3 /app/lab_client.py evaluate '{"kind": "scf", "bond_lengths_A": [1.4, 1.6], "starting_magnetization": null,
    "k_points": 48, "smearing": "mv", "degauss_Ry": 0.005, "box_A": 10, "ecutwfc_Ry": 80}'
python3 /app/lab_client.py job/1
python3 /app/lab_client.py status
```

- **`starting_magnetization`:** `null` for a non-spin-polarized calculation, or one value in [−1, 1] per atom for a collinear spin-polarized one.
- **Numerical settings** (all required): `k_points` along the chain (Γ-centred; k-points × atoms ≤ 384 and k-points × cell length ≤ 480 Å), `smearing` (`mv` Marzari–Vanderbilt, `gaussian`, `fd` Fermi–Dirac), `degauss_Ry`, `box_A` (transverse size) and `ecutwfc_Ry`.
- **Results:** total energy in Ry per cell, Fermi energy, band gap (eV, 0 for a metal), force on each atom along the chain, and for spin-polarized runs the total and absolute magnetization per cell and the moment on each site.
- **Cost:** one unit per 96 (k-points × atoms per cell), rounded up, minimum 1: a 2-atom cell costs 1 unit at 48 k-points, 2 at 96 and 4 at 192. The budget is 36 units.
- **Timing:** one calculation runs at a time and takes 0.5–6 minutes (longer close to a magnetic instability, where the SCF converges slowly). `evaluate` returns within 20 s with a job id; poll `job/<id>` until it is `complete` or `failed` and never resubmit a running job. Failed jobs stay charged.

## Conventions

- In a 2-atom cell the bonds alternate between R(1 − δ) and R(1 + δ); δ is the dimerization.
- Energies are per atom, in meV.
- Every graded quantity refers to the zero-smearing (zero electronic temperature) limit with converged k-point sampling, cutoff and box, for collinear spins.
- A *uniform* chain has δ = 0. The *magnetic energy gain* of the uniform chain is E(non-spin-polarized) − E(lowest-energy spin-polarized solution) at the same R; it is 0 if no magnetic solution exists.
- The *dimerization energy* is E(lowest-energy uniform chain, any collinear spin arrangement) − E(equilibrium dimerized chain) at the same R; δ_eq is the dimerization of the energy minimum.
- At fixed R, the *critical dimerization* δ_c is the bond alternation above which the lowest-energy spin-polarized solution has zero site moments; report 0 if that is already the case at δ = 0.

## Deliverables

Save `/app/submission.json`:

```json
{"delta_eq": {"1.25": 0.0, "1.75": 0.0},
 "dimerization_energy_meV_per_atom": {"1.25": 0.0, "1.75": 0.0},
 "uniform_magnetic_energy_gain_meV_per_atom": {"1.75": 0.0, "2.0": 0.0},
 "critical_dimerization": {"1.75": 0.0, "2.0": 0.0},
 "uniform_chain_band_gap_eV": {"1.25": 0.0},
 "ground_state": {"1.0": "uniform_nonmagnetic", "2.0": "uniform_nonmagnetic"},
 "evidence_ids": [1, 2]}
```

Values are placeholders. Keys are R in Å.

- **`uniform_chain_band_gap_eV`:** the band gap of the lowest-energy uniform chain at R = 1.25 Å.
- **`ground_state`:** the lowest-energy state over all δ and collinear spin arrangements at that R, one of `uniform_nonmagnetic`, `uniform_ferromagnetic`, `uniform_antiferromagnetic`, `dimerized_nonmagnetic`, `dimerized_magnetic`.
- **`evidence_ids`:** the calculations your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- which instability wins at each R and why, and what the coexistence of dimerization and magnetism looks like;
- how the magnetic solutions were found and checked;
- how you established numerical convergence, and which quantities needed it;
- the physical meaning of a broken-symmetry magnetic solution in a one-dimensional chain;
- which conclusions these calculations do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks every field of the submission; physical reasoning is reviewed separately.
