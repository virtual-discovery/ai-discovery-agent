# n-type Mg₂Si, Mg₂Ge and Mg₂Sn: the mid-infrared band and the conduction states at X

n-type Mg₂Si and Mg₂Ge absorb in a band near 0.40 and 0.56 eV; n-type Mg₂Sn absorbs from about 0.17 eV. Using the author-released calculations of one study, and digitized spectra, establish which transition produces this absorption, how its energy depends on the level of theory and on temperature, and how its strength depends on temperature.

## Data and calculation service

`/app/catalog.json` describes the archive (Zenodo record 23319151, H. Sugawara, 2026): the codes, the conventions of its files, the request formats and every run you may purchase. `/app/experiment.json` holds digitized absorption spectra (Koenig et al. 1961 for Mg₂Si; Lott and Lynch 1966 for Mg₂Ge; interband edges) and is free.

Purchases return the archived results; nothing is recomputed:

- **`bands`:** all-electron LDA or QSGW bands along L–Γ–X–W–Γ;
- **`x_scan`:** fine LDA or QSGW scans from X along X–Γ or X–W;
- **`qsgw_box`:** QSGW energies of the two lowest conduction bands on a box of wave vectors around X;
- **`momentum_matrix`:** pseudopotential PBE momentum matrix elements between the two lowest conduction bands of Mg₂Si on the box grid;
- **`tote`:** rows of the ecalj TOTE files of the G0W0 and QSGW runs;
- **`questaal`:** eigenvalues at Γ and X from the Questaal runs on Mg₂Si, with and without ladder diagrams;
- **`pseudopotential_g0w0`:** BerkeleyGW quasiparticle energies;
- **`electron_phonon`:** electron–phonon renormalizations of the band edges.

```sh
python3 /app/lab_client.py evaluate '{"kind": "x_scan", "compound": "Mg2Si", "level": "QSGW", "direction": "X-W"}'
python3 /app/lab_client.py evaluate '{"kind": "qsgw_box", "compound": "Mg2Ge", "box": "coarse"}' > box.json
python3 /app/lab_client.py evaluate '{"kind": "momentum_matrix", "points": [[0, 0, 0.02], [0.02, 0, 0]]}'
python3 /app/lab_client.py evaluate '{"kind": "tote", "run": "g0w0_pbe/Mg2Ge"}'
python3 /app/lab_client.py evaluate '{"kind": "electron_phonon", "compound": "Mg2Si", "method": "ZG", "supercell": 3, "temperature_K": 300}'
python3 /app/lab_client.py status
```

A `qsgw_box` request costs 5 units and every other request 1 unit (a `momentum_matrix` request takes up to 40 points). The budget is 32 units. Results print as JSON; redirect large ones to a file.

## Conventions

- **X1 and X3** are the symmetry labels of the two lowest conduction states at X. In Mg₂Si, at every level of theory archived here, X1 is the lower of the two. A splitting is E(X3) − E(X1) in eV, with its sign.
- **Levels of theory.** All are all-electron ecalj calculations at the room-temperature lattice constants of the catalog, without spin–orbit coupling.
  - **LDA.**
  - **G0W0:** the one-shot G0W0 from the LDA with production settings. Its quasiparticle energies solve the linearized equation E = ε + Z[Σ(ε) − V_xc], where ε is the LDA eigenvalue and Z the renormalization factor.
  - **QSGW:** the production quasiparticle self-consistent GW calculation at its final, unscaled self-energy.
- **Ladder correction:** the change of X3 − X1 of Mg₂Si, in meV, when ladder diagrams are included in the screened interaction W of a QSGW calculation, all else equal, with the production ladder window.
- **Infrared/Raman ratio.** Electronic Raman scattering by the same transition has an intensity proportional to the electron density n. The infrared absorption integrated over the band is proportional to n⟨|p|²⟩. Here p is the momentum matrix element between the lower and the upper conduction state at the electron's wave vector, taken at its lowest non-vanishing order in the wave vector measured from X. The average ⟨…⟩ runs over the conduction electrons in thermal equilibrium at density n, all of which are in the X valleys. Report Q = (IR/Raman at 300 K)/(IR/Raman at 77 K) at fixed n, on the QSGW band structure.
- **Predicted splitting:** the X3 − X1 splitting that an optical measurement on the compound at the given temperature should find. Start from the QSGW splitting and apply the archived calculations of how the splitting changes with the temperature of the measurement. Treat these changes as additive, and take the lattice constant below 100 K to equal its 100 K value. `T0` means T → 0.

## Deliverables

Save `/app/submission.json`:

```json
{"x3_minus_x1_eV": {"Mg2Si": {"LDA": 0.0, "G0W0": 0.0, "QSGW": 0.0},
                    "Mg2Ge": {"LDA": 0.0, "G0W0": 0.0, "QSGW": 0.0},
                    "Mg2Sn": {"LDA": 0.0, "G0W0": 0.0, "QSGW": 0.0}},
 "vertex_correction_Mg2Si_meV": 0.0,
 "ir_raman_ratio_300K_over_77K": {"Mg2Si_1e16": 0.0, "Mg2Si_3e18": 0.0, "Mg2Sn_1e16": 0.0},
 "predicted_splitting_eV": {"Mg2Si_T0": 0.0, "Mg2Si_300K": 0.0, "Mg2Ge_300K": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`x3_minus_x1_eV`:** the signed splitting at each level of theory for each compound.
- **`vertex_correction_Mg2Si_meV`:** the ladder correction.
- **`ir_raman_ratio_300K_over_77K`:** Q for Mg₂Si at n = 1 × 10¹⁶ and 3 × 10¹⁸ cm⁻³, and for Mg₂Sn at 1 × 10¹⁶ cm⁻³.
- **`predicted_splitting_eV`:** the predicted splitting of Mg₂Si at T → 0 and at 300 K, and of Mg₂Ge at 300 K.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- which transition produces the absorption, its selection rule, and the evidence against other explanations;
- how X1 and X3 were identified in each compound and at each level of theory;
- how the levels of theory compare with each other and with the measured bands, and what the ladder correction implies;
- how Q was obtained, and why it differs from the ratio of the temperatures;
- how the predicted splittings were assembled;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the splittings, the ladder correction, the infrared/Raman ratios and the predicted splittings; physical reasoning is reviewed separately.
