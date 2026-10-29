# MnF2: a two-sublattice exchange model from spin spirals, and what spirals cannot see

Build a Heisenberg exchange model for rutile MnF2 from author-released noncollinear DFT spin-spiral calculations, test how far it transfers to propagation directions you cannot buy, compute the classical Néel temperature it implies, and decide whether the altermagnetic exchange can be determined from these data. Here, altermagnetic exchange means bonds of the same length and direction whose exchange constant differs between the two Mn sublattices.

## Data and calculation service

`/app/catalog.json` describes the following (VASP, PBE, constrained moments, generalized Bloch theorem):

- the structure;
- how the spirals are defined;
- the units as documented by the archive;
- the spirals you may purchase (`spirals_available`, keyed by direction and q);
- the spirals whose energies you must predict (`spirals_to_predict`).

A purchase returns the released calculation: the VASP total energy of the cell, the archive's tabulated energy, the constrained moments, the INCAR, the k-mesh, the warnings VASP printed and the number of electronic steps.

```sh
python3 /app/lab_client.py evaluate '{"kind": "spiral", "key": "[214]-07"}'
python3 /app/lab_client.py status
```

Every purchase costs 1 unit; the budget is 24 units. These are retrieved author calculations, not new DFT.

## Conventions

Report exchange constants for H = −Σ_{i≠j} J_ij ê_i·ê_j over the Mn moments: unit vectors, with every ordered pair counted, so each bond appears twice.

- **J1:** Mn–Mn along c, same sublattice, 3.36 Å.
- **J2:** Mn1–Mn2 along the body diagonal, 3.89 Å.
- **J3:** Mn–Mn along a, same sublattice, 4.95 Å.

The Néel temperature is the thermodynamic-limit transition temperature, in K, of the classical model with these unit-vector spins.

Spiral energies are E(q) − E(q = 0) of the 6-atom primitive cell (2 Mn + 4 F), in meV.

## Deliverables

Save `/app/submission.json`:

```json
{"spiral_predictions_meV_per_cell": {"P111-01": 0.0},
 "J1_meV": 0.0, "J2_meV": 0.0, "J3_meV": 0.0,
 "neel_temperature_K": 0.0,
 "altermagnetic_shell": {"distance_A": 0.0, "sublattice_mean_J_meV": 0.0, "sublattice_difference_J_meV": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **Predictions:** give an energy for every id in `spirals_to_predict`.
- **`altermagnetic_shell`:** describe the shortest Mn–Mn shell in which symmetry allows bonds of the same length and direction to carry different exchange on Mn1 and Mn2:
  - `distance_A`: the length of these bonds;
  - `sublattice_mean_J_meV`: the average exchange (J(Mn1, d) + J(Mn2, d))/2;
  - `sublattice_difference_J_meV`: the difference J(Mn1, d) − J(Mn2, d) for d = a[110], or `null` if these data do not determine it.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the model, conventions and energy normalization;
- the spiral geometry and the magnetic ground state;
- how the configurations and the number of exchange shells were chosen and validated;
- how the Néel temperature was obtained, and its uncertainty;
- what these data determine about the altermagnetic exchange, and what calculation would determine the rest;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the predictions, J1–J3, the Néel temperature and the altermagnetic shell; physical reasoning is reviewed separately.
