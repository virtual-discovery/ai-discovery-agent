# Germanium under pressure: internal strain, clamped- versus relaxed-ion shear, and the Raman mode (PBE)

In the diamond structure neither atom of the primitive cell sits at a centre of inversion, so a shear strain exerts forces on the two sublattices and they shift relative to each other (internal strain). Using new Quantum ESPRESSO calculations, determine for germanium compressed hydrostatically to 6 GPa (diamond-anvil-cell conditions, below the β-tin transition) how much of the shear stiffness this internal relaxation removes in PBE, what Kleinman parameter it implies, and whether the softening is consistent with the zone-centre optical (Raman) phonon of the same protocol. Determine also the protocol's lattice parameters at zero pressure and at 6 GPa, its direct gap at Γ at zero pressure, the other two cubic elastic constants at 6 GPa, and the volume dependence of the Raman frequency there.

## Calculation service

Every request runs new `pw.x` calculations for the two-atom primitive cell of diamond Ge (PBE). One calculation runs at a time and takes from a few seconds to about ten minutes; `evaluate` returns within 20 s with a job id while it runs, and you poll `job/<id>` until it is `complete` or `failed` (never resubmit a running job).

```sh
python3 /app/lab_client.py task          # full request schema, allowed values, cost rule and budget
python3 /app/lab_client.py evaluate '{"kind": "scf", "pseudopotential": "ge_pbe_v1.4.uspp.F.UPF", "lattice_A": 5.70,
    "strain": [[0,0,0],[0,0,0],[0,0,0]], "displacement_A": [0,0,0], "ecutwfc_Ry": 45, "kmesh": 6, "kshift": 1,
    "occupations": "fixed", "degauss_Ry": null, "forc_conv_thr_Ry_per_bohr": null}'
python3 /app/lab_client.py job/1
python3 /app/lab_client.py status
```

Every field is required (no defaults):

- **`kind`:** `"scf"`, or `"relax"` (BFGS relaxation of both atomic positions at fixed cell, starting from the given displacement).
- **Geometry:** `lattice_A` (cubic lattice parameter, 5.4–6.0 Å); `strain`, a symmetric 3×3 strain tensor ε in the cube axes (cell vectors and ideal atomic sites are mapped r → (1 + ε) r; components within ±0.03); `displacement_A`, the displacement u of atom 2 relative to atom 1 beyond the homogeneous strain (Cartesian Å; atom 1 is placed at −u/2 and atom 2 at (1 + ε)(a/4)(1,1,1) + u/2).
- **Numerics:** `pseudopotential` (two Ge files: GBRV ultrasoft with the 3d shell in valence, or PSlibrary ultrasoft with the 3d shell frozen in the core), `ecutwfc_Ry` (30, 35, 45 or 60; density cutoff 8×), `kmesh` (n×n×n Monkhorst–Pack, n = 4–16) and `kshift` (0 contains Γ, 1 is offset by half a step), `occupations` (`"fixed"` or `"smearing"`, Marzari–Vanderbilt) with `degauss_Ry` (null when fixed), and `forc_conv_thr_Ry_per_bohr` for relaxations (null for scf).
- **Results:** total energy (Ry per cell), stress tensor and pressure (GPa, in the pw.x sign convention: the negative of the Cauchy stress), forces on both atoms (Ry/bohr), band information (highest occupied and lowest unoccupied levels or the Fermi energy, and all band energies at the mesh point nearest Γ), the cell, the number of k points and symmetry operations, and for relaxations the ionic steps and final positions (crystal coordinates).
- **Cost:** proportional to the number of symmetry-distinct k points of the request (lower-symmetry cells need more), times a factor for the pseudopotential and cutoff and ×3 for a relaxation, with a minimum of 1 unit and a maximum of 15 per request; `task` gives the exact rule. The budget is 32 units. Failed jobs stay charged.

## Conventions

- Report static-lattice values of this protocol: PBE with the 3d-in-valence pseudopotential `ge_pbe_v1.4.uspp.F.UPF`, converged with respect to the cutoff, k-point sampling and occupations the lab exposes.
- `a0_A`: lattice parameter at zero pressure; `a_6GPa_A`: lattice parameter at a hydrostatic pressure of 6.0 GPa.
- `gamma_gap_at_a0_eV`: E(Γ2′) − E(Γ25′) at a0, where Γ2′ is the s-like antibonding singlet and Γ25′ the p-like bonding triplet (the valence-band top of an ordinary diamond-structure semiconductor); signed.
- Elastic quantities refer to the crystal at a_6GPa_A, in the small-strain limit: the coefficients C = dσ/dε of the Cauchy stress (tensile positive) with respect to a strain applied to the compressed cell, in GPa, Voigt notation with engineering shear strain (ε4 = 2ε_yz). C11 and C12 need no internal relaxation in this structure.
- `C44_clamped_ion_GPa`: homogeneous shear with the atoms fixed at their strained ideal sites; `C44_relaxed_ion_GPa`: the atomic positions relaxed.
- `kleinman_zeta`: for a pure shear ε4, the relaxed relative displacement is u_x = −ζ a ε4 / 4 (u as defined above, a = a_6GPa_A; ζ = 0 is a homogeneous deformation, ζ = 1 keeps the nearest-neighbour bond lengths fixed to first order).
- `raman_frequency_cm-1`: harmonic zone-centre optical frequency at a_6GPa_A, in cm⁻¹; `raman_gruneisen`: γ = −d ln ω / d ln V at a_6GPa_A under hydrostatic strain.

## Deliverables

Save `/app/submission.json`:

```json
{"a0_A": 0.0, "gamma_gap_at_a0_eV": 0.0, "a_6GPa_A": 0.0,
 "C11_GPa": 0.0, "C12_GPa": 0.0, "C44_clamped_ion_GPa": 0.0, "C44_relaxed_ion_GPa": 0.0, "kleinman_zeta": 0.0,
 "raman_frequency_cm-1": 0.0, "raman_gruneisen": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders. `evidence_ids` lists the calculations your numbers rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how the lattice parameters and the numerical settings were chosen, and the evidence that each graded quantity is converged;
- the electronic structure near the Fermi level at zero pressure and at 6 GPa;
- how the clamped- and relaxed-ion shear constants and ζ were obtained, and whether they are consistent with the zone-centre force constant (state the relation you used and its numerical check);
- how the Raman frequency and its Grüneisen parameter were obtained;
- how the protocol compares with experiment, and which conclusions these calculations do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks every field of the submission; physical reasoning is reviewed separately.
