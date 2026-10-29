# Germanium with spin–orbit coupling: strain-split valence bands and the s state at Γ (PBE)

Compressive biaxial strain is what turns germanium into a heavy-hole-dominated channel material: (001) in-plane compression splits the fourfold Γ8 valence-band top into heavy- and light-hole levels, and the light hole mixes with the spin–orbit split-off band. Using new Quantum ESPRESSO calculations with spin–orbit coupling, determine for PBE germanium the spin–orbit splitting and the position of the s-like antibonding level at Γ at the protocol's own lattice parameter, how 1 % biaxial in-plane compression with a relaxed out-of-plane lattice parameter splits and shifts the levels at Γ, the shear deformation potential b, and the compression at which the s-like level and the split-off level become degenerate.

## Calculation service

Every request runs a new self-consistent `pw.x` calculation for the two-atom primitive cell of diamond Ge with conventional lattice parameters a (in plane, x and y) and c (out of plane, z), optionally followed by a band calculation at k points you choose, with optional symmetry labels of the bands from `bands.x`. One calculation runs at a time and takes from about half a minute to ten minutes; `evaluate` returns within 20 s with a job id while it runs, and you poll `job/<id>` until it is `complete` or `failed` (never resubmit a running job).

```sh
python3 /app/lab_client.py task          # full request schema, allowed values, cost rule and budget
python3 /app/lab_client.py evaluate '{"kind": "pw", "a_par_A": 5.70, "c_A": 5.70, "spin_orbit": true,
    "ecutwfc_Ry": 50, "kmesh": 6, "kshift": 0, "occupations": "smearing", "degauss_Ry": 0.01,
    "band_kpoints_inverse_A": [[0, 0, 0], [0, 0, 0.05]], "symmetry_labels": true}'
python3 /app/lab_client.py job/1
python3 /app/lab_client.py status
```

Every field is required (no defaults):

- **Geometry:** `a_par_A` and `c_A`, the conventional lattice parameters in Å (c = a is the cubic crystal; the atoms sit at the ideal diamond sites of the strained cell).
- **Physics:** `spin_orbit` (true: noncollinear spinors with spin–orbit coupling, fully relativistic SG15 pseudopotential; false: scalar-relativistic run with the scalar-relativistic SG15 file).
- **Numerics:** `ecutwfc_Ry` (40, 50, 60 or 80; norm-conserving, density cutoff 4×), `kmesh` (n×n×n Monkhorst–Pack, n = 4–12) and `kshift` (0 contains Γ, 1 is offset by half a step), `occupations` (`"fixed"` or `"smearing"`, Marzari–Vanderbilt) with `degauss_Ry` (null when fixed).
- **Bands:** `band_kpoints_inverse_A`, up to 30 Cartesian k points in Å⁻¹ (including 2π), or `[]`; `symmetry_labels` (bands.x irreducible-representation labels at those k points).
- **Results:** total energy, stress tensor and pressure (GPa, pw.x sign convention: the negative of the Cauchy stress), Fermi energy or highest occupied and lowest unoccupied levels, all band energies at the mesh point nearest Γ, the k-point and symmetry counts, and, when requested, all band energies (and labels) at your band k points. Band energies are on each run's own absolute scale.
- **Cost:** proportional to the number of symmetry-distinct k points of the SCF mesh (lower-symmetry cells need more), with factors for spin–orbit coupling and the cutoff, plus the band k points; about one unit per minute of computing. `task` gives the exact rule. The budget is 60 units. Failed jobs stay charged.

## Conventions

- **Protocol:** PBE with spin–orbit coupling and the fully relativistic file `Ge_ONCV_PBE_FR-1.0.upf`, converged with respect to the cutoff, k-point sampling and occupations the lab exposes; static lattice.
- `a0_A`: the cubic lattice parameter at zero pressure.
- **Strain** is measured from a0: ε∥ = a/a0 − 1 (in plane, ε_xx = ε_yy) and ε⊥ = c/a0 − 1 (ε_zz). A **relaxed** biaxial state has zero stress along z at the given ε∥.
- **Levels** are energies at Γ (k = 0) in eV; each Kramers pair counts as one level. In the cubic crystal the p-like valence states form the fourfold Γ8+ (j = 3/2) and the twofold split-off Γ7+ (j = 1/2) level; `delta0_at_a0_eV` = E(Γ8+) − E(Γ7+) at a0.
- **s-like level:** the twofold odd-parity level derived from the antibonding s state (the Γ2′ state, Γ7− with spin–orbit coupling, that forms the conduction-band minimum at Γ in an ordinary diamond-structure semiconductor). `s_level_at_a0_eV` = E(s) − E(Γ8+) at a0, signed.
- **Strained levels:** under (001) strain Γ8+ splits into the heavy-hole level HH (m_j = ±3/2 character) and the light-hole level LH (m_j = ±1/2 character); SO is the split-off level. At ε∥ = −1 % relaxed: `c_over_a_at_minus1pct` = c/a, `hh_minus_lh_at_minus1pct_eV` = E(HH) − E(LH) and `so_minus_hh_at_minus1pct_eV` = E(SO) − E(HH), signed.
- `b_shear_deformation_potential_eV`: defined by E(HH) − E(LH) = b (ε_xx + ε_yy − 2ε_zz) to first order in the strain (the small-strain limit), in eV.
- `s_so_crossing_strain_percent`: the in-plane strain ε∥, in percent (negative for compression), along the relaxed biaxial path at which E(s) = E(SO).

## Deliverables

Save `/app/submission.json`:

```json
{"a0_A": 0.0, "delta0_at_a0_eV": 0.0, "s_level_at_a0_eV": 0.0,
 "c_over_a_at_minus1pct": 0.0, "hh_minus_lh_at_minus1pct_eV": 0.0, "so_minus_hh_at_minus1pct_eV": 0.0,
 "b_shear_deformation_potential_eV": 0.0, "s_so_crossing_strain_percent": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders. `evidence_ids` lists the calculations your numbers rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how a0 and the numerical settings were chosen, and the evidence that each graded quantity is converged;
- how each level at Γ was identified (not just its energy), at a0 and under strain, and the electronic structure of the protocol near the Fermi level;
- how the out-of-plane relaxation was obtained, and how b was extracted (state the model and its range of validity);
- how the crossing strain was located;
- how the protocol compares with experiment, and which conclusions these calculations do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks every field of the submission; physical reasoning is reviewed separately.
