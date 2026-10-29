# Monolayer MoS2 under strain: where does the band gap stop being direct? (PBE with spin–orbit coupling)

Unstrained monolayer 1H-MoS2 has a direct band gap at K. Using new Quantum ESPRESSO calculations, determine over which range of biaxial in-plane strain the gap of this protocol stays direct, which band edge leaves K at each end of that range and where it goes, how the K-valley band edges move with biaxial strain on an absolute energy scale, and what the gap is when the monolayer is compressed by 4 % along a zigzag direction under uniaxial stress. Determine also the protocol's equilibrium lattice parameter.

## Calculation service

Every request runs new calculations for the three-atom cell of the monolayer: `pw.x` (scf, or a relaxation of the atomic positions at fixed cell), then `pp.x` for the planar-averaged electrostatic potential, then `pw.x` band energies at k points you choose. One calculation runs at a time and takes from about one to ten minutes; `evaluate` returns within 20 s with a job id while it runs, and you poll `job/<id>` until it is `complete` or `failed` (never resubmit a running job).

```sh
python3 /app/lab_client.py task          # full request schema, allowed values, cost rule and budget
python3 /app/lab_client.py evaluate '{"kind": "scf", "relativistic": "scalar", "lattice_A": 3.20,
    "strain": [[0, 0], [0, 0]], "sulfur_height_A": 1.56, "sulfur_shift_A": [0, 0], "cell_height_A": 20,
    "ecutwfc_Ry": 45, "kmesh": 12, "band_kpoints": [[0, 0], [0.333333, 0.333333]]}'
python3 /app/lab_client.py job/1
python3 /app/lab_client.py status
```

Every field is required (no defaults):

- **`kind`:** `"scf"`, or `"relax"` (BFGS relaxation of all atomic positions at fixed cell, starting from the given geometry; scalar-relativistic only).
- **`relativistic`:** `"scalar"` (scalar-relativistic pseudopotentials) or `"full"` (fully relativistic pseudopotentials, noncollinear spinors with spin–orbit coupling).
- **Geometry:** `lattice_A`, the hexagonal lattice parameter a of the unstrained cell (a1 = a(1, 0), a2 = a(−1/2, √3/2); x is a zigzag direction, y the armchair direction); `strain`, a symmetric 2×2 in-plane strain tensor ε in these axes (cell vectors and ideal sites are mapped r → (1 + ε) r; components within ±0.05); `sulfur_height_A`, the height h of each S plane above or below the Mo plane; `sulfur_shift_A`, an in-plane Cartesian displacement of both S atoms from their strained ideal site; `cell_height_A`, the cell period normal to the layer (15, 20, 25 or 30 Å).
- **Numerics:** `ecutwfc_Ry` (35, 45, 55 or 70; density cutoff 8×), `kmesh` (n×n×1 Γ-centred Monkhorst–Pack mesh for the self-consistent density, n from 6 to 18) and `band_kpoints` (up to 60 points [k1, k2] in crystal coordinates of the strained reciprocal lattice; may be empty).
- **Results:** total energy, stress tensor (GPa, pw.x sign convention, averaged over the whole cell including the vacuum), forces, final positions, S height and shift, the highest occupied and lowest unoccupied levels on the scf mesh, the planar-averaged electrostatic potential along z, all band energies at the requested k points (raw pw.x eigenvalues, eV) and the number of occupied bands; for relaxations, the ionic steps.
- **Cost:** proportional to the number of symmetry-distinct k points of the scf mesh (lower-symmetry cells need more), times 3 for a relaxation, plus a quarter of a k point for each band k point, times 4 for spin–orbit and a factor for the cutoff; minimum 1 unit, maximum 7 units per request. `task` gives the exact rule. The budget is 48 units. Failed jobs stay charged.

## Conventions

- **Protocol:** PBE with the lab's PSlibrary pseudopotentials, converged with respect to the cutoff, k-point sampling and cell height the lab exposes. Geometries (lattice parameter, atomic positions, transverse strain) are those of scalar-relativistic calculations with relaxed atomic positions; band energies are those of fully relativistic calculations (spin–orbit coupling) at those geometries.
- `a0_A`: equilibrium lattice parameter of the unstrained monolayer (zero in-plane stress).
- Strains are relative to a0: a biaxial strain ε means a = a0(1 + ε) in both directions. Report strains in percent.
- `direct_gap_strain_range_percent`: [lower, upper], the biaxial strains that bound the range in which both the valence-band maximum and the conduction-band minimum lie at K.
- `edge_leaving_K_at_lower_end` and `edge_leaving_K_at_upper_end`: which band edge (`"valence"` or `"conduction"`) leaves K at that end of the range, and the crystal coordinates [k1, k2] of the point it moves to (any symmetry-equivalent point is accepted).
- `deformation_potentials_eV`: dE/dε at ε = 0 for biaxial strain, in eV per unit strain: `valence_K` (highest valence state at K) and `conduction_K` (lowest conduction state at K), both on an absolute energy scale referred to the vacuum level; `gap_K`, the direct gap at K.
- `uniaxial_zigzag_compression`: the monolayer with ε_xx = −4 % (relative to a0) along the zigzag direction x, the transverse direction y free (zero transverse stress) and atomic positions relaxed. `transverse_strain_percent` is its ε_yy; `gap_eV` its band gap; `vbm_k_crystal` and `cbm_k_crystal` the crystal coordinates of the band edges in the strained reciprocal lattice (any symmetry-equivalent point is accepted).

## Deliverables

Save `/app/submission.json`:

```json
{"a0_A": 0.0,
 "direct_gap_strain_range_percent": [0.0, 0.0],
 "edge_leaving_K_at_lower_end": {"band": "valence", "k_crystal": [0.0, 0.0]},
 "edge_leaving_K_at_upper_end": {"band": "conduction", "k_crystal": [0.0, 0.0]},
 "deformation_potentials_eV": {"valence_K": 0.0, "conduction_K": 0.0, "gap_K": 0.0},
 "uniaxial_zigzag_compression": {"transverse_strain_percent": 0.0, "gap_eV": 0.0,
                                 "vbm_k_crystal": [0.0, 0.0], "cbm_k_crystal": [0.0, 0.0]},
 "evidence_ids": [1, 2]}
```

Values are placeholders. `evidence_ids` lists the calculations your numbers rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how the numerical settings and the lattice parameter were chosen, and the evidence that each graded quantity is converged;
- how the geometries were obtained at each strain, and what changes if the atoms are not relaxed;
- which band edges compete at each end of the direct-gap range, their orbital character, and why they respond differently to strain; what spin–orbit coupling changes;
- how energies were put on an absolute scale;
- how the uniaxial-stress state was found and how its band edges were located;
- how the protocol compares with experiment, and which conclusions these calculations do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks every field of the submission; physical reasoning is reviewed separately.
