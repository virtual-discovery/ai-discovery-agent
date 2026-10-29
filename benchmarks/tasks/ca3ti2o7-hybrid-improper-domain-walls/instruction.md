# Ca3Ti2O7: ferroelectric domain walls in a hybrid improper ferroelectric

In the Ruddlesden–Popper titanate Ca3Ti2O7 the polarization P is not a soft mode: it is induced by the combination of an octahedral rotation about the stacking axis (ψ_z) and an octahedral tilt (φ). Using the authors' released machine-learned potentials, determine which 180° domain walls are cheap, how the energy of charged walls behaves, and how the order parameters evolve across the walls.

## Data and calculation service

`/app/catalog.json` gives the bulk polar ground state (Cmc2₁, 48 atoms, relaxed with the reference potential), the high-symmetry I4/mmm parent structure, and the potentials on offer. There are no domain-wall structures: you build them. The lab evaluates any Ca3Ti2O7 structure you submit (whole formula units, at most 1200 atoms; Å; `cell` rows are the lattice vectors):

```sh
python3 /app/lab_client.py evaluate '{"kind": "single_point", "structure": S}'
python3 /app/lab_client.py evaluate '{"kind": "relax", "structure": S, "cell": "fixed", "fmax": 0.02, "max_steps": 600}'
python3 /app/lab_client.py evaluate '{"kind": "relax", "structure": S, "cell": "normal", "normal_axis": 0}'
python3 /app/lab_client.py evaluate '{"kind": "born_charges", "structure": S}'
python3 /app/lab_client.py job/<id>        # poll a running calculation
python3 /app/lab_client.py status
```

with `S = {"symbols": [...], "positions": [[x, y, z], ...], "cell": [[...], [...], [...]]}`.

- `single_point` (1 unit): energy, forces summary and stress.
- `relax` (⌈N_atoms/400⌉ units; ⌈N_atoms/800⌉ with `"potential": "NEP"`): relaxes the atomic positions, with the cell fixed or with the length of one cell vector along the normal of the plane spanned by the other two (`"cell": "normal"`, `normal_axis` = index of that vector). It returns the energy, convergence information and the final structure.
- `born_charges` (1 unit): per-atom charges and Born effective charge tensors of the reference potential.
- `"potential"` may be `"qNEP-w-BECs"` (default, the reference model), `"qNEP"` or `"NEP"`.

The budget is 24 units. Calculations take seconds to a few minutes; poll running jobs. Results print as JSON; redirect them to files.

## Conventions

- **Wall types.** All walls here reverse P (180° walls). A **ψ wall** reverses the rotation ψ_z and keeps the tilt; a **φ wall** reverses the tilt and keeps the rotation.
- **Neutral wall:** the wall plane contains the stacking axis and P (its normal is the in-plane direction perpendicular to P). **Layer wall:** the wall plane is parallel to the perovskite layers. **Charged wall:** the normal is parallel to P, so walls alternate head-to-head and tail-to-tail.
- **Wall energy** E_DW: the excess energy over the bulk ground state per unit wall area, in mJ/m², for the reference potential, with the in-plane lattice of the wall plane equal to that of the bulk and all atoms relaxed. For each wall type, report the lowest-energy wall position, converged with respect to the separation between walls.
- **Charged-wall slope:** the change of the charged ψ-wall energy E_DW with the domain width d (the distance between adjacent walls), dE_DW/dd in mJ m⁻² nm⁻¹, in the long-domain regime.
- **Profiles:** `"ising"` if the order parameter that reverses passes through (nearly) zero at the wall centre; `"rotational"` if it reverses by rotating, keeping at least half of its bulk magnitude everywhere.

## Deliverables

Save `/app/submission.json`:

```json
{"E_neutral_psi_wall_mJ_m2": 0.0,
 "E_neutral_phi_wall_mJ_m2": 0.0,
 "E_layer_psi_wall_mJ_m2": 0.0,
 "charged_psi_wall_slope_mJ_m2_per_nm": 0.0,
 "psi_wall_profile": "ising",
 "phi_wall_profile": "ising",
 "evidence_ids": [1, 2]}
```

Values are placeholders. The profiles refer to the neutral walls. `evidence_ids` lists the calculations your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how the coupling between P, ψ_z and φ determines which domain states you combined, and how you built and checked the walls;
- how the wall energies were normalized and converged, and which wall positions you compared;
- why the charged walls behave as they do, and which potential you trusted for them;
- the order-parameter and polarization profiles across the walls;
- the ranking of the walls, its structural origin, and which conclusions the model and these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the four energies or slopes and the two profiles; physical reasoning is reviewed separately.
