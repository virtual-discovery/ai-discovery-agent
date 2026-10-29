# 1D wires: thermal contraction, anharmonic heat capacity and the limits of stability

Exfoliable one-dimensional wires (CuC2, TaSe3, AuSe2) are candidate interconnects. Using the authors' archived stochastic self-consistent harmonic approximation (SSCHA) and DFT data, determine whether these wires contract on heating and by how much, how strongly anharmonicity changes their heat capacity relative to the harmonic (Dulong–Petit) picture, and up to which temperature the TaSe3 wire remains dynamically stable.

## Data and calculation service

`/app/catalog.json` describes the archive: the SSCHA runs of each wire (temperature, configurations per population, and which files exist), the zero-temperature DFPT phonons, the DFT input/output files, the single strained calculation, the file formats and units, and the request formats. Purchases return archived files or verbatim excerpts; nothing is recomputed:

```sh
python3 /app/lab_client.py evaluate '{"kind": "sscha_log", "wire": "TaSe3", "temperature_K": 300, "configurations": 10000}'
python3 /app/lab_client.py evaluate '{"kind": "sscha_matrices", "wire": "CuC2", "temperature_K": 300, "configurations": 10000, "matrix": "auxiliary"}' > aux.json
python3 /app/lab_client.py evaluate '{"kind": "dfpt", "wire": "TaSe3", "asr": "1d"}'
python3 /app/lab_client.py evaluate '{"kind": "dft_file", "wire": "TaSe3", "file": "scf.out"}'
python3 /app/lab_client.py evaluate '{"kind": "strained_calculation"}'
python3 /app/lab_client.py status
```

Every purchase costs 1 unit; the budget is 36 units. Results print as JSON; redirect large ones to a file. Some files of the archive are withheld by this task, as marked in the catalog.

## Conventions

- All SSCHA runs use the same fixed cell (the exfoliated DFT geometry, vacuum included) at every temperature. Lengths and volumes refer to that cell.
- **Axial stiffness** k = ∂²E/∂ε² per unit cell (eV), for an infinitesimal axial strain ε of the zero-temperature harmonic DFT (DFPT) wire with all internal coordinates relaxed.
- **Linear thermal expansion** α_L = (1/L)(∂L/∂T) at zero axial force, obtained from the fixed-length SSCHA runs with 10,000 configurations, with the length response of the wire given by k. Report it in units of 10⁻⁶ K⁻¹ (negative means contraction).
- **Heat capacity** C_V = T(∂S/∂T) at fixed length, where S is the SSCHA vibrational entropy of the 1 × 1 × 6 simulation supercell (10,000-configuration runs). Report the ratio C_V / ((3N − n₀) k_B), with N the number of atoms in the supercell and n₀ the number of its zero-frequency rigid-body modes, so that a classical harmonic crystal gives 1.
- **Instability temperature:** the lowest temperature at which the TaSe3 wire is dynamically unstable in the SSCHA sense, judged on the commensurate q grid of the simulation supercell. Between archived temperatures, interpolate the signed squared frequency (sign(ω) ω²) linearly in T.

## Deliverables

Save `/app/submission.json`:

```json
{"TaSe3": {"axial_stiffness_eV": 0.0, "alpha_L_300K_1e-6_per_K": 0.0, "cv_ratio_300K": 0.0,
           "instability_temperature_K": 0.0},
 "CuC2": {"cv_ratio_300K": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`TaSe3`:** k; α_L at 300 K; the heat-capacity ratio at 300 K; the instability temperature (K).
- **`CuC2`:** the heat-capacity ratio at 300 K.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the mechanism and sign of the thermal expansion, and how you established the sign from the archived data;
- how k was obtained, and why you trusted that route over others available in the archive;
- how the heat capacity follows from SSCHA quantities, and the size and sign of its deviation from the harmonic value for each wire;
- how you judged stability, and which features of the archived phonons you did not treat as instabilities, and why;
- how stochastic sampling and convergence of the archived runs limit your numbers;
- which conclusions these data do not support (for example, about the quasi-harmonic approximation and about AuSe2).

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks k, α_L, the two heat-capacity ratios and the instability temperature; physical reasoning is reviewed separately.
