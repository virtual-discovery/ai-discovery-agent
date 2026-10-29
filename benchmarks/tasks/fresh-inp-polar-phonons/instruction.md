# InP: LO-TO splitting, Born charges and dielectric screening under compression (DFPT-PBE)

Zincblende InP has one polar optical phonon at the zone centre. Its LO-TO splitting is set by the Born effective charge Z*, the electronic dielectric constant ε∞ and the cell volume. Using new Quantum ESPRESSO calculations, determine what this computational protocol (PBE, the lab's ultrasoft pseudopotentials, static lattice) predicts for these quantities at its own equilibrium; how they change under hydrostatic compression, and which factor carries the change; and how a [001] uniaxial compression makes the LO frequency depend on the direction of the phonon wavevector.

## Calculation service

Every request runs new pw.x (and ph.x) calculations for the two-atom primitive cell (In at (0,0,0), P at (¼,¼,¼) in crystal coordinates). Requests run one at a time and take from under a minute to about ten minutes; `evaluate` returns a job id while it runs, and you poll it.

```sh
python3 /app/lab_client.py task          # this interface, the cost table and the budget
python3 /app/lab_client.py evaluate '{"kind": "dfpt", "lattice_A": 5.95, "kmesh": 6, "kshift": 1}'
python3 /app/lab_client.py job/3
python3 /app/lab_client.py status
```

Request fields:

- `kind`:
  - `"scf"`: pw.x self-consistent run. Returns the total energy, pressure and stress tensor, and the highest occupied and lowest unoccupied band energies found on the k mesh.
  - `"dielectric"`: SCF plus the ph.x electric-field response at q = 0. Also returns the electronic dielectric tensor and the Born effective charge tensors of In and P as computed (dF/dE; rows are field directions, columns force directions; no sum rule imposed).
  - `"dfpt"`: SCF plus the ph.x electric-field and phonon response at q = 0. Also returns the Γ force-constant matrix (6×6, Ry/bohr², ordered In x,y,z then P x,y,z; no sum rule imposed), the atomic masses, and the analytic Γ frequencies obtained by diagonalising it (no nonanalytic term).
- `lattice_A`: cubic lattice parameter a in Å (5.70–6.10).
- `strain_zz` (optional, default 0): the z components of all three primitive vectors are multiplied by 1 + strain_zz (−0.03 to 0.03); x and y are unchanged and the atoms stay at their ideal crystal coordinates.
- `kmesh`: n for an n×n×n Monkhorst–Pack mesh; `kshift`: 0 for a Γ-centred mesh, 1 for a mesh offset by half a step in each direction.
- `ecutwfc_Ry` (optional): 50 (default) or 65; the charge-density cutoff is 8× the wavefunction cutoff.

Costs (units) depend on the kind, the mesh and whether the cell is strained; `task` lists the full table and the allowed meshes. Cubic cells: `scf` 1 (2 for k10–12); `dielectric` k4 1, k6 2, k8 4, k10 7; `dfpt` k4 2, k6 4. Strained cells: `dielectric` k4 3, k6 5; `dfpt` k4 4, k6 8. A 65 Ry cutoff doubles the cost. The budget is 40 units. Failed jobs stay charged.

## Conventions

- Report everything for this protocol and the static lattice at its own zero-pressure lattice parameter a0, converged with respect to the numerical parameters the lab exposes.
- `born_charge_In`: the Born effective charge of In (a positive number), after imposing charge neutrality.
- `eps_inf`: the electronic (clamped-ion) dielectric constant; `eps_0`: the static dielectric constant including the polar phonon.
- `omega_TO_cm-1`, `omega_LO_cm-1`: zone-centre optical frequencies, the LO for q → 0.
- Volume derivatives are d ln X / d ln V at a0 under hydrostatic (cubic) strain: X = Z* and ε∞, and the TO frequency, reported as the Grüneisen parameter γ_TO = −d ln ω_TO / d ln V.
- `LOTO_splitting_under_compression`: `"grows"` or `"shrinks"`, the sign of the change of ω_LO − ω_TO under small hydrostatic compression from a0.
- Strained cell: the cubic cell at a0 with `strain_zz` = −0.020. Report the LO frequency for q → 0 along z ([001]) and along x ([100]).

## Deliverables

Save `/app/submission.json`:

```json
{"a0_A": 0.0,
 "omega_TO_cm-1": 0.0, "omega_LO_cm-1": 0.0, "born_charge_In": 0.0, "eps_inf": 0.0, "eps_0": 0.0,
 "dlnZ_dlnV": 0.0, "dln_eps_inf_dlnV": 0.0, "gruneisen_TO": 0.0, "LOTO_splitting_under_compression": "grows or shrinks",
 "strained_omega_LO_q_z_cm-1": 0.0, "strained_omega_LO_q_x_cm-1": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders. `evidence_ids` lists the calculations your numbers rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how a0 and the numerical settings were chosen, and the evidence that each graded quantity is converged;
- how the LO frequencies were obtained from the raw outputs;
- whether the LO-TO splitting grows or shrinks under compression, and how much of the change comes from Z*, from ε∞ and from the volume;
- how far the protocol's splitting is from experiment and how much of that comes from the PBE ε∞;
- why the strained crystal has two different LO frequencies;
- which conclusions these calculations do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the submitted values; physical reasoning is reviewed separately.
