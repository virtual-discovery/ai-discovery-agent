# Oxygen vacancies in ferroelectric Hf1−xZrxO2: what sets their formation energy, and where do they sit?

Using author-released DFT calculations, determine what controls the formation energy of neutral oxygen vacancies in polar orthorhombic (Pca2₁) Hf1−xZrxO2: the cations around the vacated site, the kind of site, or the overall composition. Predict the formation energy at sites you cannot buy, find which sites hold the equilibrium vacancy population as the Zr content changes, and check which conclusions survive a change of exchange-correlation functional.

## Data and calculation service

`/app/catalog.json` describes:

- **Supercells:** three 96-atom supercells at x = 0.25, 0.50 and 0.75 (Hf24Zr8O64, Hf16Zr16O64, Hf8Zr24O64), one random Hf/Zr decoration each, as CIF files in `/app/structures/`. Oxygen sites are identified by their 0-based atom index in the CIF.
- **Released table:** the authors' site table for the x = 0.50 cell, `/app/x050_site_table.csv`, with its columns described in the catalog.
- **Purchasable calculations:** retrieved author results, not new DFT.

```sh
python3 /app/lab_client.py evaluate '{"kind": "vacancy", "composition": "x0.25", "site": 40}'       # VASP-PBE vacancy cell, 1 unit
python3 /app/lab_client.py evaluate '{"kind": "pristine", "composition": "x0.25"}'                  # VASP-PBE pristine cell, 1 unit
python3 /app/lab_client.py evaluate '{"kind": "castep_vacancy", "functional": "LDA", "site": 57}'   # CASTEP vacancy cell, 1 unit
python3 /app/lab_client.py evaluate '{"kind": "castep_333_vacancy", "site": 112}'                   # CASTEP-PBE 3x3x3 cell, 2 units
python3 /app/lab_client.py status
```

The budget is 36 units. The sites listed in `sites_to_predict` cannot be purchased.

## Conventions

- **Composition:** x is the Zr fraction.
- **Formation energy:** for a neutral vacancy, E_f = E(supercell with the vacancy) − E(pristine supercell of the same composition) + μ_O. Use μ_O = −10.09 eV (the O-poor limit), on the absolute scale of the VASP total energies. No further corrections apply.
- **First-shell cations:** the Hf and Zr atoms within 2.3 Å of an oxygen in the staged (unrelaxed) structures, periodic images included.
- **Equilibrium population:** dilute, non-interacting neutral vacancies at T = 700 K (k_B = 8.617333 × 10⁻⁵ eV/K), with O-poor E_f and configurational statistics only. The crystal is macroscopic, and its cation sites are occupied independently at random (Zr with probability x). The vacancy site fraction is the number of vacancies per O site: the average of exp(−E_f/k_BT) over all O sites of that crystal.
- **LDA description:** describe the CASTEP-LDA vacancy enthalpies by E = E₀ + Δ·s + b·n_Zr. Here s = 1 for vacancy type O3 and 0 otherwise, and n_Zr is the number of Zr neighbours, both taken from the authors' labels of each LDA vacancy. Δ and b are the least-squares values over all 64 LDA vacancies.

## Deliverables

Save `/app/submission.json`:

```json
{"predicted_formation_energies_eV": {"x0.25-O34": 0.0},
 "vacancy_Zr_neighbour_distribution": {"x0.25": [0.0, 0.0, 0.0, 0.0, 0.0], "x0.75": [0.0, 0.0, 0.0, 0.0, 0.0]},
 "vacancy_site_fraction": {"x0.25": 0.0, "x0.75": 0.0},
 "lda_common_slope_model": {"offset_O3_eV": 0.0, "slope_per_Zr_eV": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`predicted_formation_energies_eV`:** the O-poor E_f for every id in `sites_to_predict`.
- **`vacancy_Zr_neighbour_distribution`:** for each composition, the fractions of equilibrium vacancies whose vacated oxygen had 0, 1, 2, 3 or 4 first-shell Zr neighbours. The five fractions sum to 1.
- **`vacancy_site_fraction`:** the equilibrium vacancy site fraction at each composition.
- **`lda_common_slope_model`:** your estimates of Δ (`offset_O3_eV`) and b (`slope_per_Zr_eV`), in eV.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- what controls E_f, with numbers, and which candidate explanations the data rule out;
- how the purchases were chosen and how the predictions were validated;
- how the populations were obtained;
- which conclusions are robust to the functional and to the supercell size;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the predictions, the populations and the LDA description; physical reasoning is reviewed separately.
