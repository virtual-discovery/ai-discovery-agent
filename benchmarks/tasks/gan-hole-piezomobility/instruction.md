# Wurtzite GaN: how strain reshapes the hole mobility tensor

The low hole mobility of gallium nitride limits p-channel devices. Determine how each strain component changes the phonon-limited hole mobility tensor of wurtzite GaN, which strain maximizes the in-plane or the c-axis hole mobility, and why.

## Data and calculation service

`/app/catalog.json` describes 37 author-released first-principles calculations: unstrained GaN and 36 states strained along a single component, labelled `compress_X%_eK` or `expand_X%_eK` (X = 1, 2, 3; K = 1 to 6). The catalog gives each state's cell and atomic positions, and the methods. For any state you may purchase:

- the GW quasiparticle band energies on a 10×10×10 k-grid;
- the EPW hole mobility tensors at 100–500 K for a hole density of 10¹⁸ cm⁻³, either in the self-energy relaxation-time approximation (SERTA) or from the iterative Boltzmann equation (IBTE). An IBTE purchase also returns its SERTA starting point and the solver log.

```sh
python3 /app/lab_client.py evaluate '{"kind": "bands", "state": "unstrained"}'
python3 /app/lab_client.py evaluate '{"kind": "mobility", "state": "compress_1%_e1", "solver": "serta"}'
python3 /app/lab_client.py evaluate '{"kind": "mobility", "state": "compress_1%_e1", "solver": "ibte"}'
python3 /app/lab_client.py status
```

Costs are 1 unit for bands, 1 unit for a SERTA mobility and 3 units for an IBTE mobility. The budget is 50 units. These are retrieved author calculations, not new ones. The mobilities of the states listed in `to_predict` are not served; their bands are. Results print as JSON; redirect large ones to a file.

## Conventions

- **Mobility:** the phonon-limited hole drift mobility tensor in cm²/Vs from the IBTE at zero magnetic field and 300 K. Use the Cartesian frame of the catalog cells (x along a1, z along c). Off-diagonal elements are the symmetrized (μ_ij + μ_ji)/2.
- **Strain:** the symmetric small strain ε of a state relative to the unstrained cell, as the Voigt vector (e1, …, e6) = (ε_xx, ε_yy, ε_zz, 2ε_yz, 2ε_xz, 2ε_xy). Components are fractions (0.01 = 1%), negative for compression.
- **The twelve 2% states** are the archived states labelled `compress_2%_*` and `expand_2%_*`.
- **In-plane directional mobility:** the largest eigenvalue of the xy block of the tensor, that is, the largest mobility along any in-plane direction.

## Deliverables

Save `/app/submission.json`:

```json
{"mobility_predictions_cm2_per_Vs": {"compress_3%_e3": {"xx": 0.0, "yy": 0.0, "zz": 0.0, "yz": 0.0, "xz": 0.0, "xy": 0.0}},
 "best_mu_zz_2pct": {"strain_voigt": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "mu_zz_cm2_per_Vs": 0.0},
 "best_inplane_2pct": {"strain_voigt": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "mu_inplane_max_cm2_per_Vs": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`mobility_predictions_cm2_per_Vs`:** the full 300 K IBTE tensor of every state in `to_predict`.
- **`best_mu_zz_2pct`:** the strain of the 2% state with the largest μ_zz, and that μ_zz.
- **`best_inplane_2pct`:** the strain of the 2% state with the largest in-plane directional mobility, and that mobility.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the character of the top valence bands and how each strain component changes their order;
- how the predictions were obtained, and the symmetry relations you relied on and how you checked them;
- over what strain range the response is linear;
- a design recommendation, with the size of the effect and its temperature dependence;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the predictions and the two 2% design answers; physical reasoning is reviewed separately.
