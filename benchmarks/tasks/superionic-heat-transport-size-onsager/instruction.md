# Fluorite superionics: lattice thermal conductivity from the crystal to the superionic state

Determine the lattice thermal conductivity κ of UO₂, CaF₂ and PbF₂ from released classical molecular-dynamics runs, from the crystalline regime through the superionic transition, where the anion sublattice becomes mobile. Establish how κ depends on temperature and on the size of the periodic simulation cell, whether it keeps falling like that of a phonon gas or levels off once the anions diffuse, and which of the released simulations can be trusted for this.

## Data and calculation service

`/app/catalog.json` lists released LAMMPS runs (empirical force fields, fixed volume, 0.8 ns of production per run) by material, ensemble (NVE or NVT), supercell (n×n×n conventional cells) and nominal temperature. A purchase returns, for one run, its input deck, a thermodynamic summary, and the time series of the heat flux and of the summed anion velocities that the deck writes to its thermo log (format and decoding are described in the catalog):

```sh
python3 /app/lab_client.py evaluate '{"kind": "md_run", "run": "PbF2-NVE-3x3x3-700K"}' > /tmp/run.json
python3 /app/lab_client.py status
```

A run costs 1 (2×2×2), 2 (3×3×3) or 4 (4×4×4) units; the budget is 48 units. The states listed under `predict` are released but withheld; you must predict them. These are retrieved author simulations, not new ones.

## Conventions

- κ in W m⁻¹ K⁻¹ is the Green–Kubo thermal conductivity of a state as these simulations give it: this force field, this cell, this fixed volume, isotropic average over the three Cartesian directions.
- A state is (material, supercell, nominal temperature). For a run, use the run's own data.

## Deliverables

Save `/app/submission.json`:

```json
{"kappa_predictions_W_mK": {"UO2-4x4x4-1700K": 0.0},
 "kappa_runs_W_mK": {"UO2-NVE-4x4x4-3500K": 0.0},
 "thermostat_biases_kappa": false,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`kappa_predictions_W_mK`:** κ for every state in the catalog's `predict` list.
- **`kappa_runs_W_mK`:** κ for every run in the catalog's `report_runs` list.
- **`thermostat_biases_kappa`:** `true` if the NVT runs give a κ that differs from that of the corresponding NVE runs beyond statistical error, otherwise `false`.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how κ and its statistical uncertainty were obtained from the series, and how you validated the estimator in each regime;
- how κ depends on temperature, and what carries heat once the anions diffuse;
- how the cell size and the ensemble change κ, and how you decided which runs to trust;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the predictions, the run conductivities and the thermostat verdict; physical reasoning is reviewed separately.
