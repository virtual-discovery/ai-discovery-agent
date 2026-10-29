# NiO: Heisenberg exchange from PBE+U total energies and its U dependence

NiO is the textbook superexchange antiferromagnet: Ni²⁺ (S = 1) ions on an fcc lattice, ordered in ferromagnetic (111) planes stacked antiferromagnetically (AFM-II). Using new spin-polarized PBE+U calculations, determine the nearest-neighbour (90° Ni–O–Ni) and next-nearest-neighbour (180° Ni–O–Ni) exchange constants J1 and J2, how they change with the Hubbard U, the mean-field Néel temperature they imply, and how the picture changes in the plain-PBE limit.

## Calculation service

The lab runs new Quantum ESPRESSO 6.7 `pw.x` SCF calculations of rocksalt NiO at a fixed cubic lattice parameter a = 4.17 Å: collinear spin-polarized PBE+U with the simplified (Dudarev) U on Ni 3d and non-orthogonalized atomic projectors, the SSSP 1.3 precision pseudopotentials `ni_pbe_v1.4.uspp.F.UPF` (GBRV) and `O.pbe-n-kjpaw_psl.0.1.UPF` (PSlibrary), density cutoff 8 × the wavefunction cutoff. `GET /task` (`python3 /app/lab_client.py task`) returns the full request schema: the four cells offered (vectors and atomic positions), the allowed values and the charge rule.

```sh
python3 /app/lab_client.py evaluate '{"kind": "scf", "cell": "rhombohedral_2fu", "starting_magnetization": [0.5, -0.5],
    "hubbard_U_eV": 6.0, "kpoints": [4, 4, 4], "kshift": 0, "ecutwfc_Ry": 50, "smearing": "gaussian", "degauss_Ry": 0.005}'
python3 /app/lab_client.py job/1
python3 /app/lab_client.py status
```

- **Cells:** `primitive_1fu` (1 Ni), `rhombohedral_2fu` and `tetragonal_2fu` (2 Ni each), `cubic_4fu` (4 Ni). `starting_magnetization` gives one value in [−1, 1] per Ni in the order listed for the cell (null = non-spin-polarized).
- **Settings (all required):** `hubbard_U_eV` 0–9 eV in steps of 0.1 (0 is plain PBE), the Monkhorst–Pack mesh `kpoints` [n1, n2, n3] along the cell vectors and `kshift`, `ecutwfc_Ry` (40, 50, 60 or 75), `smearing` (`gaussian`, `mv`, `fd`) and `degauss_Ry`.
- **Results:** free energy F = E − TS as pw.x prints it, the smearing term −TS and the internal energy E (Ry per cell); Fermi energy and band gap on the mesh (eV, 0 for a metal); Ni and O moments, total and absolute magnetization; Ni 3d occupation-matrix traces and eigenvalues per spin; SCF iterations and wall time.
- **Cost:** ceil(w × n1·n2·n3 / 1000 × f) units, at least 1, with cell weights w = 0.4 (primitive), 10 (rhombohedral), 2.5 (tetragonal), 20 (cubic) and cutoff factors f = 0.72, 1, 1.31, 1.84 for 40, 50, 60, 75 Ry; for example a rhombohedral [4,4,4] mesh at 50 Ry costs 1 unit and [6,6,6] costs 3. The budget is 36 units.
- **Timing:** one calculation runs at a time and takes 10 s to about 8 minutes. `evaluate` returns within 20 s with a job id; poll `job/<id>` until it is `complete` or `failed`, and never resubmit a running job. Failed jobs (for example an SCF that does not converge in 200 iterations) stay charged.

## Conventions

- **Exchange model:** E = E0 + Σ_{i<j} J_ij s_i s_j for collinear configurations, with s_i = ±1 the direction of the Ni spin (S = 1, so S·S = S² = 1 in the energy mapping), every Ni pair counted once, J > 0 antiferromagnetic. J1 couples the 12 nearest Ni neighbours (a/√2), J2 the 6 next-nearest (a); longer-range couplings are neglected. Energies in meV.
- **Graded protocol:** every graded number refers to the lab's fixed protocol at the stated U, converged with respect to k-point sampling and cutoff and in the zero-smearing (zero electronic temperature) limit.
- **Néel temperature:** the molecular-field (Weiss mean-field) Néel temperature of AFM-II order for the quantum spin-1 Hamiltonian H = Σ_{i<j} J_ij **S**_i·**S**_j with your J1 and J2, in K.
- **Plain PBE:** at U = 0, E(FM) − E(AFM-II) per Ni (meV), for the lowest-energy ferromagnetic state.
- **FM gap closing:** the value of U (eV) below which the ferromagnetic state has no band gap.

## Deliverables

Save `/app/submission.json`:

```json
{"J1_meV": {"4.0": 0.0, "8.0": 0.0},
 "J2_meV": {"4.0": 0.0, "8.0": 0.0},
 "neel_temperature_mean_field_K": {"4.0": 0.0, "8.0": 0.0},
 "pbe_fm_minus_afm2_meV_per_Ni": 0.0,
 "fm_gap_closing_U_eV": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders. Keys are U in eV.

- **`evidence_ids`:** the calculations your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the energy mapping, which magnetic configurations and cells you used and why;
- how you established numerical convergence, and which quantities needed it;
- how J1 and J2 depend on U, and the physical mechanism behind it;
- the mean-field Néel temperature and its limitations;
- what changes at small U and whether a Heisenberg description still applies there;
- which conclusions these calculations do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks every field of the submission; physical reasoning is reviewed separately.
