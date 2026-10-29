# SrTiO3: polar and antiferrodistortive soft modes in PBE

Cubic perovskite SrTiO3 has two candidate lattice instabilities: a polar (ferroelectric) mode at the zone centre Γ and an antiferrodistortive rotation of the TiO6 octahedra at the zone-corner point R. Using new Quantum ESPRESSO calculations, determine what this computational protocol (PBE, the lab's pseudopotentials, static lattice) predicts at its own equilibrium lattice parameter: which of the two instabilities exist, how strong each is (harmonic frequencies and energy wells), which gives the larger energy gain, the dielectric screening and the LO frequencies at Γ, and how 1 % biaxial in-plane tension changes the in-plane polar mode.

## Calculation service

Every request runs new calculations: `pw.x` (self-consistent energy, forces and stress) and, for `"dielectric"`, `ph.x` (electric-field response at q = 0). One calculation runs at a time and takes from under a minute to about ten minutes; `evaluate` returns within 20 s with a job id while it runs, and you poll `job/<id>` until it is `complete` or `failed` (never resubmit a running job).

```sh
python3 /app/lab_client.py task          # full request schema, allowed values, cost rule and budget
python3 /app/lab_client.py evaluate '{"kind": "scf", "cell": "perovskite_5", "lattice_A": 3.95, "c_A": 3.95,
    "displacements_A": [[0,0,0],[0,0,0.01],[0,0,0],[0,0,0],[0,0,0]], "ecutwfc_Ry": 60, "kmesh": 6}'
python3 /app/lab_client.py job/1
python3 /app/lab_client.py status
```

Every field is required (no defaults):

- **`kind`:** `"scf"`, or `"dielectric"` (scf plus the electronic dielectric tensor ε∞ and the Born effective charges as computed, dF/dE with no sum rule imposed; undisplaced `perovskite_5` cell only).
- **`cell`:** `"perovskite_5"`, the five-atom cell a1 = (a, 0, 0), a2 = (0, a, 0), a3 = (0, 0, c) with atoms Sr (0, 0, 0), Ti (½, ½, ½), O_x (0, ½, ½), O_y (½, 0, ½), O_z (½, ½, 0) in fractions of (a, a, c) (O_x lies between Ti neighbours along x, and so on); or `"rotated_10"`, the cell doubled at R, a1 = (0, a, c), a2 = (a, 0, c), a3 = (a, a, 0), holding two formula units: atoms 1–5 at the five-atom sites, atoms 6–10 at the same sites shifted by (a, 0, 0).
- **Geometry:** `lattice_A` (a) and `c_A` (c, with |c/a − 1| ≤ 0.05), and `displacements_A`, one Cartesian displacement [dx, dy, dz] in Å per atom in the cell's atom order (zeros for the ideal structure; components within ±0.4 Å).
- **Numerics:** `ecutwfc_Ry` (45, 60, 75 or 90; density cutoff 8×) and `kmesh` (n of an n×n×n Γ-centred mesh in the reciprocal basis of the chosen cell: 4, 6 or 8 for `perovskite_5`, 2, 3, 4 or 6 for `rotated_10`).
- **Results:** total energy (Ry per cell), forces on every atom, stress, and for `"dielectric"` the ε∞ tensor and raw Born charge tensors; also the masses used.
- **Cost:** for `"scf"`, proportional to the number of symmetry-distinct k points of the displaced cell (lower-symmetry displacements cost more) with a per-cell weight (`rotated_10` about 5× per k point) and a cutoff factor; `"dielectric"` 4, 5 or 9 units for meshes 4, 6 or 8 at 60 Ry. Minimum 1 unit, at most 10 per request; `task` gives the exact rule and each result reports the k points charged. The budget is 80 units. Failed jobs stay charged.

## Conventions

- **Protocol:** PBE with the lab's pseudopotentials, static lattice, converged with respect to the cutoff and k-point sampling the lab exposes. All quantities are at the protocol's own zero-pressure lattice parameter a0 (cubic cell) unless stated.
- **Frequencies** in cm⁻¹; an imaginary frequency (instability) is reported as a negative number.
- `gamma_soft_mode_cm-1`: the lowest polar optical mode at Γ (transverse, harmonic). `gamma_LO_highest_cm-1`: the highest longitudinal optical frequency for q → 0. `eps_inf`: the electronic dielectric constant; `born_charge_Ti`: the Born effective charge of Ti after imposing charge neutrality.
- `polar_well_depth_meV_per_fu`: the energy lowering (positive number, meV per formula unit) at the minimum of the energy along the eigendisplacement of the Γ soft mode polarized along z (all atoms displaced in the fixed proportions of that harmonic mode, Sr held fixed), five-atom cell at fixed cubic a0.
- **Octahedral rotation:** the a0a0c⁻ pattern, rotation about z alternating in sign between neighbouring octahedra in all three directions. The rotation angle θ is defined by the displacement u of the O_x and O_y atoms perpendicular to their Ti–O–Ti line: tan θ = 2u/a. `tilt_angle_deg` and `tilt_well_depth_meV_per_fu` (positive, meV per formula unit) describe the energy minimum along this pattern at fixed cubic cell a0; `R_mode_cm-1` is the harmonic frequency of this mode at θ → 0.
- `larger_energy_gain`: `"polar"` or `"tilt"`, whichever of the two well depths is larger.
- **Strained state:** a = b = a0 (1 + 0.01) (biaxial in-plane tension), c relaxed to zero out-of-plane stress, atoms at the ideal sites. `strained_c_A` is that c; `strained_inplane_polar_mode_cm-1` the lowest polar mode polarized along x in this cell.

## Deliverables

Save `/app/submission.json`:

```json
{"a0_A": 0.0,
 "gamma_soft_mode_cm-1": 0.0, "gamma_LO_highest_cm-1": 0.0, "eps_inf": 0.0, "born_charge_Ti": 0.0,
 "polar_well_depth_meV_per_fu": 0.0,
 "R_mode_cm-1": 0.0, "tilt_angle_deg": 0.0, "tilt_well_depth_meV_per_fu": 0.0,
 "larger_energy_gain": "polar or tilt",
 "strained_c_A": 0.0, "strained_inplane_polar_mode_cm-1": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders. `evidence_ids` lists the calculations your numbers rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how a0 and the numerical settings were chosen, and the evidence that each graded quantity is converged;
- how the force constants, the soft-mode frequency and eigenvector, and the LO frequencies were obtained from the raw outputs;
- how the energy wells were mapped and referenced, and how the two instabilities compare and interact;
- how the strained state was found and what the strain does to the polar instability;
- how the protocol's predictions compare with real SrTiO3 (quantum paraelectric, 105 K antiferrodistortive transition), and which conclusions these calculations do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks every field of the submission; physical reasoning is reviewed separately.
