# GeTe: the ferroelectric double well, its strain coupling and its collapse under pressure (PBE)

Below about 700 K germanium telluride is a rhombohedral (R3m) ferroelectric: relative to the rocksalt structure, Ge and Te shift against each other along a cube diagonal and the cell shears along the same axis. Using new Quantum ESPRESSO calculations, determine how deep the ferroelectric double well of GeTe is in PBE, how much of that depth comes from the coupling of the polar displacement to the cell (shear and volume) rather than from the displacement alone, at what hydrostatic pressure the rocksalt structure becomes stable against the polar displacement, and how spin–orbit coupling changes the well depth.

## Calculation service

Every request runs one new `pw.x` SCF calculation for the two-atom rhombohedral primitive cell of GeTe (one formula unit): Ge at (0, 0, 0) and Te at (x, x, x) in crystal coordinates; α = 60° with x = 0.5 is rocksalt. One calculation runs at a time and takes from about twenty seconds to about ten minutes; `evaluate` returns within 20 s with a job id while it runs, and you poll `job/<id>` until it is `complete` or `failed` (never resubmit a running job).

```sh
python3 /app/lab_client.py task          # full request schema, allowed values, cost rule and budget
python3 /app/lab_client.py evaluate '{"kind": "scf", "pseudopotential_set": "sssp", "a_rh_A": 4.25, "alpha_deg": 60.0,
    "x_Te": 0.5, "ecutwfc_Ry": 50, "kmesh": 8, "kshift": 0, "occupations": "fixed", "degauss_Ry": null}'
python3 /app/lab_client.py job/1
python3 /app/lab_client.py status
```

Every field is required (no defaults):

- **Geometry:** `a_rh_A` (rhombohedral lattice parameter, the edge of the primitive cell, 3.9–4.6 Å; the conventional cubic parameter of rocksalt is √2 a_rh), `alpha_deg` (54–62°), `x_Te` (0.42–0.58).
- **Numerics:** `pseudopotential_set` (`sssp`: the SSSP 1.3 precision pair, GBRV Ge with 3d in valence and a PSlibrary-derived Te; `pslibrary`: PSlibrary 1.0.0 Ge (3d frozen in the core) and Te, scalar-relativistic; `pslibrary_fullrel`: their fully relativistic counterparts, which switch spin–orbit coupling on), `ecutwfc_Ry` (40, 50, 60 or 75; density cutoff 8×), `kmesh` (n×n×n Monkhorst–Pack, n = 6–24) and `kshift` (0 contains Γ, 1 is offset by half a step), `occupations` (`"fixed"` or `"smearing"`, Marzari–Vanderbilt) with `degauss_Ry` (null when fixed).
- **Results:** total energy (Ry per cell), stress tensor and pressure (GPa, pw.x sign convention: the negative of the Cauchy stress, so that the pressure is positive under compression), forces on both atoms (Ry/bohr), the band edges found on the k mesh with insulator filling (valence top, conduction bottom, indirect and smallest direct gap and where they occur), the Fermi energy and smearing term when smearing is used, the number of k points and symmetry operations, and the wall time.
- **Cost:** proportional to the number of symmetry-distinct k points of the request (polar structures have fewer symmetries than rocksalt and need more), times factors for the pseudopotential set (spin–orbit is expensive) and the cutoff, with a minimum of 1 unit and a maximum of 8 per request; `task` gives the exact rule. The budget is 70 units. Failed jobs stay charged.

## Conventions

- **Protocol:** every graded number except the spin–orbit change refers to PBE with the `sssp` set, static lattice, converged with respect to the cutoff, the k-point sampling and the occupations the lab exposes (zero-smearing limit). Energies are per formula unit (GeTe), in meV.
- `cubic_a0_A`: conventional lattice parameter of rocksalt GeTe at zero pressure.
- `r3m_a_rh_A`, `r3m_alpha_deg`, `r3m_x_Te`: the fully relaxed R3m structure at zero pressure (x is the Te coordinate with Ge at the origin, x < 0.5).
- `well_depth_meV_per_fu`: E(rocksalt at zero pressure) − E(R3m fully relaxed at zero pressure).
- `clamped_well_depth_meV_per_fu`: the depth of the double well along x alone, with the cell fixed at the zero-pressure rocksalt cell (α = 60°, a_rh = cubic_a0_A/√2): E(x = 0.5) − min over x of E.
- `soft_mode_critical_pressure_GPa`: the hydrostatic pressure above which rocksalt GeTe is stable against the zone-centre polar displacement, i.e. where the harmonic force constant ∂²E/∂x² of the rocksalt structure at x = 0.5 changes sign.
- `soc_well_depth_change_meV_per_fu`: the change of E(rocksalt) − E(R3m) caused by spin–orbit coupling alone, evaluated at fixed structures, namely the protocol's zero-pressure rocksalt cell and its relaxed R3m structure.

## Deliverables

Save `/app/submission.json`:

```json
{"cubic_a0_A": 0.0, "r3m_a_rh_A": 0.0, "r3m_alpha_deg": 0.0, "r3m_x_Te": 0.0,
 "well_depth_meV_per_fu": 0.0, "clamped_well_depth_meV_per_fu": 0.0,
 "soft_mode_critical_pressure_GPa": 0.0, "soc_well_depth_change_meV_per_fu": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders. `evidence_ids` lists the calculations your numbers rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how the structures were found and the evidence that each graded quantity is converged;
- the electronic structure of the two phases near the gap, at zero pressure and under compression, and what it implied for your numerical choices;
- how the well depth divides into the contribution of the polar displacement alone and that of its coupling to the cell, and the symmetry argument for that coupling;
- how the critical pressure was obtained, and whether your data say if the transition is continuous;
- how the spin–orbit change was isolated;
- how the protocol compares with experiment, and which conclusions these calculations do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks every field of the submission; physical reasoning is reviewed separately.
