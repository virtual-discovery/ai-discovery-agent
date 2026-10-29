# Hydrogen evolution on Ag-Au-Pd-Pt alloys: adsorption energetics, coverage and the work function

Hydrogen-evolution activity is usually rationalised through the hydrogen adsorption energy and the Sabatier principle. On a multimetallic surface, however, each composition exposes a distribution of adsorption sites, and the bare-surface work function has been proposed as an additional activity descriptor.

Using author-released measurements and calculations on Ag-Au-Pd-Pt alloys, determine:

- what controls the measured activity;
- whether the work function carries information beyond the adsorption energetics;
- the current densities of a library you cannot measure.

## Data and calculation service

`/app/catalog.json` describes three thin-film materials libraries (Lib-1, Lib-2, Lib-3) screened for acidic hydrogen evolution by scanning electrochemical cell microscopy (SECCM).

- **Measurement areas.** Each library has about 340 measurement areas (MDI), each with an EDX composition, a stage position and a wafer row.
- **Descriptors.** `/app/descriptors.json` gives, for every measurement area of every library, surrogate-model descriptors computed at its composition. They are free; the catalog and the file itself describe them.
- **DFT.** The catalog lists the released DFT compositions (RPBE, fcc(111) slabs) and the released hydrogen-coverage calculations.

Purchases return released data. They are not new measurements or calculations.

```sh
python3 /app/lab_client.py evaluate '{"kind": "seccm_row", "library": "Lib-1", "row": 9, "potential_mV": -200}'   # 1 unit
python3 /app/lab_client.py evaluate '{"kind": "lib3_spot", "mdi": 125}'                                          # 2 units
python3 /app/lab_client.py evaluate '{"kind": "coverage_path", "id": "pt100/variant-00001"}'                     # 1 unit
python3 /app/lab_client.py evaluate '{"kind": "dft_composition", "id": "pt100"}'                                  # 1 unit
python3 /app/lab_client.py status
```

What each purchase returns:

- **`seccm_row`:** the current densities of every measurement area in one wafer row of Lib-1 or Lib-2, at one potential.
- **`lib3_spot`:** the current densities of one Lib-3 measurement area from `calibration_spots_available`, at all four potentials. At most 5 distinct Lib-3 spots can be bought.
- **`coverage_path`:** one released hydrogen-coverage calculation.
- **`dft_composition`:** the site-resolved DFT adsorption energies and the work functions of one DFT composition.

The budget is 60 units. Results print as JSON; redirect them to files.

## Conventions

- **Potentials and currents.** Potentials are vs RHE. Current densities j are in A/cm², and cathodic current is negative. Predictions are log10(|j| / (A/cm²)).
- **Adsorption-energy scale.** ΔE_H is the hydrogen adsorption energy of one site relative to ½E(H₂), from RPBE total energies without zero-point or entropy terms. The dilute value refers to 1/9 ML, which is the scale of the released DFT data and descriptors.
- **`volcano_centre_eV`.** The dilute ΔE_H of a site at which that site's contribution to the hydrogen-evolution rate is largest, as determined from the measured activity.
- **`coverage_shift_eV`.** For the equimolar AgAuPdPt surface, treat the H–H interaction in mean field, linear in the coverage θ. Report how much higher a given site's adsorption energy is at θ = 1 ML than its own dilute (θ = 1/9 ML) value. Average over the released equimolar decorations.
- **`alpha_wf`.** The dimensionless coefficient α_Φ in a rate factor exp[α_Φ(Φ − Φ_ref)/k_BT]. Here Φ is the bare-surface work function of a measurement area (eV), Φ_ref = 5.499 eV (pure Pt) and T = 298.15 K. Take α_Φ from your best description of the measured currents that also accounts for the hydrogen adsorption energetics. A positive value means that a higher work function gives more current.

## Deliverables

Save `/app/submission.json`:

```json
{"lib3_log10_abs_current": {"1": [0.0, 0.0, 0.0, 0.0]},
 "volcano_centre_eV": 0.0,
 "coverage_shift_eV": 0.0,
 "alpha_wf": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`lib3_log10_abs_current`:** for every MDI in Lib-3 `spots_to_predict` (as a string), the predicted log10|j| at −100, −200, −300 and −400 mV, in that order.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- what controls the activity, with numbers, and which competing explanations the data rule out;
- how hydrogen coverage and the distribution of adsorption sites enter your description;
- how differences between libraries were handled, and how the Lib-3 predictions were made and validated;
- how you tested whether the work function adds information, and what it would mean physically;
- which compositions are most active and why;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the Lib-3 predictions, the volcano centre, the coverage shift and α_Φ; physical reasoning is reviewed separately.
