# Fluorite superionic conductors: transition, activation energies and cell-size effects

Determine the superionic (order–disorder) transition temperatures and the activation energies of anion diffusion in three fluorite-structure conductors (CaF₂, PbF₂ and UO₂) from released classical molecular-dynamics runs. Establish how the size of the periodic simulation cell changes these quantities and the diffusion coefficients themselves, and contrast the fluorites with α-AgI.

## Data and calculation service

`/app/catalog.json` lists released NVT molecular-dynamics runs (LAMMPS, empirical force fields, fixed lattice parameter, 0.8 ns production per run) by material, supercell (n×n×n conventional cells) and target temperature. A purchase returns, for one run:

- its LAMMPS input deck;
- its released mean-square-displacement table (ANALISI output, averaged over 4 blocks): one row per lag, with the MSD of each atom type and of each type's centre of mass;
- a summary of its thermodynamic log.

```sh
python3 /app/lab_client.py evaluate '{"kind": "md_run", "material": "CaF2", "supercell": 3, "temperature_K": 1400}'
python3 /app/lab_client.py status
```

The cost of a run grows with its number of atoms (field `cost` in the catalog). The budget is 140 units. The runs listed under `predict` are released but withheld; you must predict them. These are retrieved author calculations, not new simulations.

## Conventions

- **D:** tracer self-diffusion coefficient of the mobile ions (F in CaF₂ and PbF₂, O in UO₂, Ag in AgI), in Å²/ps, from MSD = 6Dt in the diffusive regime, measured relative to the immobile sublattice.
- **T_c:** the temperature of the maximum of the isochoric heat capacity c_V(T) of the 4×4×4 (768-atom) cell, in K.
- **Activation energies:** E_a from D = D₀ exp(−E_a/k_BT), in eV, fitted to 768-atom runs at least 100 K below T_c (`ordered`) or at least 100 K above T_c (`superionic`).

## Deliverables

Save `/app/submission.json`:

```json
{"transition_temperature_K": {"CaF2": 0.0, "PbF2": 0.0, "UO2": 0.0},
 "activation_energy_eV": {"CaF2": {"ordered": 0.0, "superionic": 0.0},
                          "PbF2": {"ordered": 0.0, "superionic": 0.0}},
 "diffusion_predictions_A2_per_ps": {"CaF2-4x4x4-1250K": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **Predictions:** give D for every id in the catalog's `predict` list.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how D was extracted from the released tables;
- how T_c was located, and its uncertainty;
- how the cell size changes D, T_c and E_a, and how you decided which cells are large enough;
- which explanation of the size effects the data support, and which they rule out;
- why α-AgI behaves differently;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks T_c, the activation energies and the predictions; physical reasoning is reviewed separately.
