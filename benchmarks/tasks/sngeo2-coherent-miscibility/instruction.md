# Rutile Sn1−xGexO2: equilibrium solubility or coherent metastability?

Epitaxial rutile Sn1−xGexO2 films have been grown single-phase over a wide composition range at about 800–1000 K. Using the authors' released DFT data for this alloy, decide whether these films can be equilibrium solid solutions, and which mechanism explains them: Ge–Sn ordering, a small mixing enthalpy, kinetic trapping inside an incoherent miscibility gap, coherency strain between the decomposing phases, substrate strain, or entropy. Predict the incoherent and coherent phase boundaries that support your answer.

## Data and calculation service

`/app/catalog.json` describes the released archive (VASP; PBE and r2SCAN). It stages, free of charge:

- the relaxed end members;
- the authors' 80 small-cell Ge/Sn orderings with their PBE energies.

The other data are purchased:

```sh
python3 /app/lab_client.py evaluate '{"kind": "sqs", "functional": "PBE", "x_Ge": 0.375, "index": 1}'
python3 /app/lab_client.py evaluate '{"kind": "strain_scan", "functional": "PBE", "id": "100-GeO2-03-07"}'
python3 /app/lab_client.py evaluate '{"kind": "ordering_r2SCAN", "id": "O17"}'
python3 /app/lab_client.py evaluate '{"kind": "hse_setting_run", "name": "Alloy_PBE_Values/Sn_25_2"}'
python3 /app/lab_client.py status
```

- **`sqs`:** a relaxed 72-atom special quasirandom structure with its archived energy.
- **`strain_scan`:** one of the archive's coherency-strain scans of GeO2 or SnO2.
- **`ordering_r2SCAN`:** the r2SCAN result for one of the orderings.
- **`hse_setting_run`:** a separate PBE relaxation of an SQS or end member, run with the hybrid-functional settings.

The catalog explains each kind and the id formats. PBE purchases and `hse_setting_run` cost 1 unit; r2SCAN purchases cost 2. The budget is 58 units. These are retrieved author calculations, not new DFT. Results print as JSON; redirect them to files.

## Conventions

- **Composition:** x is the Ge fraction. Energies are per cation (formula unit), and temperatures are in K.
- **Mixing enthalpy:** of the disordered alloy, represented by the SQS, relative to relaxed rutile SnO2 and GeO2.
- **Coherency strain energy Δe_cs(x):** Cahn's approximation, as used by Doak and Wolverton. It is the composition-weighted energy of the two end members, each held at the in-plane lattice of the homogeneous alloy of composition x on the habit plane and relaxed in every other degree of freedom, relative to the relaxed end members.
- **Habit plane:** the released plane with the smallest Δe_cs. Plane labels are the archive's (`001`, `100`, `101`).
- **Coherent free energy:** Δg_coh(x, T) = Δh_dis(x) − Δe_cs(x) − T Δs_ideal(x), with the ideal configurational entropy of the cation sublattice.
  - The coherent spinodal is where ∂²Δg_coh/∂x² = 0, and its critical temperature is the maximum of that curve.
  - The coherent binodal is the common tangent of Δg_coh at fixed T.
- **Functionals:** quantities are PBE unless labelled otherwise. A quantity labelled r2SCAN uses r2SCAN data only.
- **Incoherent critical temperature:** the top of the miscibility gap of the unconstrained (incoherent) alloy, from any method you can justify.

## Deliverables

Save `/app/submission.json`:

```json
{"mixing_enthalpy_x050_meV_per_cation": 0.0,
 "coherency_strain_energy_x050_meV_per_cation": {"001": 0.0, "100": 0.0, "101": 0.0},
 "habit_plane": "<label>",
 "incoherent_Tc_K": 0.0,
 "coherent_spinodal_Tc_K": {"PBE": 0.0, "r2SCAN": 0.0},
 "coherent_boundaries_700K_PBE": {"binodal": [0.0, 0.0], "spinodal": [0.0, 0.0]},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`coherency_strain_energy_x050_meV_per_cation`:** Δe_cs at x = 0.5 for each released dataset.
- **`coherent_boundaries_700K_PBE`:** the x values of the two coherent binodal compositions and of the two coherent spinodal compositions at 700 K.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- whether the alloy tends to order or to cluster, and the evidence;
- which released data your numbers rest on, and why you trust them;
- the incoherent phase diagram and the method behind it, including the effect of short-range order;
- the coherency calculation, the habit plane, and how sensitive the result is to the functional;
- which explanation of the homogeneous films the numbers support, and which they rule out;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the mixing enthalpy, the strain energies and habit plane, the incoherent and coherent critical temperatures, and the 700 K coherent boundaries. Physical reasoning is reviewed separately.
