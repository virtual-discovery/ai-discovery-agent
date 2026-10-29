# When does the ferromagnetic reference give the right exchange? NiO, FeCl2, FePS3 and FeP

Tellez-Mora et al. (2024) determine magnetic ground states self-consistently: compute exchange parameters with the magnetic force theorem in a reference state (first the ferromagnet), test that reference in linear spin-wave theory, update the reference and repeat until the reference is stable with its own exchange. Using their released TB2J output, establish for each material whether the exchange computed in the ferromagnetic (FM) reference already describes the converged magnetic state, and quantify how much the exchange depends on the reference: rock-salt NiO, layered FeCl2 and FePS3 (insulators; FePS3 with and without a Hubbard U) and the metallic helimagnet FeP. For FeP, also test whether the converged spiral is consistent with its own exchange.

## Data and calculation service

`/app/catalog.json` describes ten released TB2J calculations: for each material (and for FePS3, each of the DFT and DFT+U workflows) an FM reference and the authors' updated antiferromagnetic (AFM) reference. It also gives the exchange convention and what can be purchased. Purchases return archived author output; they are not new calculations.

```sh
python3 /app/lab_client.py evaluate '{"kind": "exchange", "calculation": "FeP/FM", "window": "0-5"}'
python3 /app/lab_client.py status
```

An `exchange` purchase returns every isotropic pair of one calculation whose archived distance lies in one window, together with that calculation's TB2J header, cell, atoms, charges and moments. A window costs 1–3 units; the catalog lists the windows and their costs. The budget is 68 units, well short of the whole catalog. Results are large: redirect them to files.

## Conventions

- **Exchange:** the archived TB2J convention, E = −Σ_{i≠j} J_ij ê_i·ê_j over ordered pairs, with unit vectors along the local moments, J > 0 ferromagnetic, in meV.
- **Converged calculation:** for each material (and FePS3 workflow), the calculation whose own reference is stable in linear spin-wave theory, so that the loop stops there.
- **Stabilization energy:** ΔE = E(FM) − E(converged configuration) per magnetic ion (Ni or Fe), for classical unit vectors, with both energies evaluated with one exchange set. The converged configuration is the moment arrangement of the converged calculation, in its cell. For the FM-reference exchange, carry its couplings onto that cell by lattice translation and ignore the small structural differences between the calculations.
- **k_min:** for an exchange set on the four-Fe cell of FeP, the wave vector at which the lowest eigenvalue of the linear spin-wave bosonic matrix around the ferromagnetic configuration is smallest (the k_min of the paper). Give it in reciprocal-lattice units of the cell of the `FeP/FM` calculation. Vectors related by sign changes of components or by integer shifts are treated as equivalent.

## Deliverables

Save `/app/submission.json`:

```json
{"fm_reference_kmin_FeP_rlu": [0.0, 0.0, 0.0],
 "converged_exchange_kmin_FeP_rlu": [0.0, 0.0, 0.0],
 "stabilization_energy_meV_per_ion": {
   "NiO": {"fm_reference": 0.0, "converged_reference": 0.0},
   "FeCl2": {"fm_reference": 0.0, "converged_reference": 0.0},
   "FePS3_DFT": {"fm_reference": 0.0, "converged_reference": 0.0},
   "FePS3_DFT+U": {"fm_reference": 0.0, "converged_reference": 0.0},
   "FeP": {"fm_reference": 0.0, "converged_reference": 0.0}},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`fm_reference_kmin_FeP_rlu`:** k_min of the FeP FM-reference exchange, to four decimals.
- **`converged_exchange_kmin_FeP_rlu`:** k_min of the exchange of FeP's converged calculation, mapped onto the four-Fe cell, with bonds that are equivalent under that calculation's supercell lattice translations treated as one coupling.
- **`stabilization_energy_meV_per_ion`:** ΔE of each material's converged configuration from its FM-reference exchange (`fm_reference`) and from the converged calculation's own exchange (`converged_reference`). `FePS3_DFT` is the U = 0 workflow and `FePS3_DFT+U` the U = 3 eV workflow. Give `null` for any value that the released data do not determine.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the conventions, and how the different cells, supercells and moment configurations were handled;
- how each calculation was tested and how the converged one was identified;
- how far the exchange lists must extend, and how that was established within the budget;
- the reference dependence of each material, with numbers, and what distinguishes the cases;
- whether FeP's converged spiral is consistent with its own exchange;
- any problems found in the released files, and how they were handled;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the two FeP k_min vectors and the ten stabilization entries; physical reasoning is reviewed separately.
