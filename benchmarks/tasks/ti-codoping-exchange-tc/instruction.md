# Co-doped Bi2Te3: which magnetic dopant pairs order, and at what temperature?

Mozumder et al. (2024) embedded 3d and 4d transition-metal impurities on the Bi sites of the topological insulator Bi2Te3 with the KKR Green-function method and computed the exchange constants J_ij of impurity dimers at 17 separations. From mean-field estimates they concluded that co-doping Cr with Ti, V, Mo or Nb strengthens ferromagnetic order, as needed for the quantum anomalous Hall effect. Using their released exchange table, determine the Curie temperatures of randomly doped Bi2Te3 for the compositions below, and decide which co-dopant raises the Curie temperature of Cr-doped Bi2Te3 the most.

## Data and calculation service

`/app/catalog.json` describes the host structure, the dopants, the columns of the archived table, the compositions and the two purchase types:

- **`pair`** returns the archived rows of one dopant pair in one shell window, verbatim. Each purchase costs 1 unit.
- **`monte_carlo`** runs a classical Monte Carlo simulation of an exchange model you specify, for randomly placed dopants on the Bi sublattice of a periodic Bi2Te3 supercell. It returns per-sample thermal averages of order-parameter moments. The cost depends on the supercell size and the number of random placements: 1–8 units per 16 placements (see the catalog).

```sh
python3 /app/lab_client.py evaluate '{"kind": "pair", "species": ["Cr", "V"], "shells": "0-6"}'
python3 /app/lab_client.py job/3
python3 /app/lab_client.py status
```

The budget is 60 units, far short of the catalog. Pair purchases are retrieved author data, not new calculations. Monte Carlo runs take up to a few minutes; poll the job rather than resubmitting. You may also run your own simulations with the tools in your container.

## Conventions

- **Doping:** dopants substitute Bi on randomly chosen sites of both Bi layers. Concentrations are fractions of all Bi sites: `V5Cr5` means 5% V and 5% Cr, i.e. (Bi0.90V0.05Cr0.05)2Te3.
- **Spin model:** classical unit-vector spins on the magnetic dopants, H = −(1/2) Σ_{i≠j} J_ij ê_i·ê_j (each pair counted once), J > 0 ferromagnetic. Every pair of magnetic dopants at the separation of an archived shell couples with the archived J of that shell and species pair, whichever species sits on which site. There are no couplings beyond the archived shells.
- **Magnetic dopants:** a dopant carries a spin if its single-impurity moment is at least 10⁻³ μB. Other dopants are non-magnetic substitutions.
- **Curie temperature T_C:** the temperature below which the disorder-averaged model has spontaneous long-range order with a nonzero net magnetization (ferromagnetic or ferrimagnetic), in the thermodynamic limit, in K.
  - Report 0 if there is no such order above 2 K.
  - Report `null` if the released data do not determine T_C.

## Deliverables

Save `/app/submission.json`:

```json
{"curie_temperature_K": {"V5Cr5": 0.0, "Cr5Mo5": 0.0, "Cr5Mn5": 0.0, "V2Cr2": 0.0, "Cr5Co5": 0.0},
 "best_codopant_for_Cr": "X",
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`curie_temperature_K`:** T_C for every composition id in the catalog.
- **`best_codopant_for_Cr`:** the element X among `codopant_candidates_with_Cr` for which Cr 5% + X 5% has the highest T_C.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the model, conventions and how the archived rows were turned into couplings;
- how T_C was obtained, and how finite size and disorder were handled;
- how and why the Monte Carlo results differ from mean-field estimates;
- what decides which co-dopant helps, with the relevant couplings;
- any problems found in the released table, and how they were handled;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the Curie temperatures and the best co-dopant; physical reasoning is reviewed separately.
