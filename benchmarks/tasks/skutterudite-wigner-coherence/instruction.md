# Yb-filled skutterudites: when does wave-like heat tunnelling rival particle-like propagation?

Skutterudite cages can host loosely bound filler atoms ("rattlers") that lower the lattice thermal conductivity. Using the authors' archived first-principles Wigner-transport data for CoSb3, FeSb3 and IrSb3 and their Yb-filled counterparts, determine whether Yb filling drives these crystals from particle-like (Boltzmann) to wave-like (Wigner) heat conduction, and at what temperature the coherence conductivity becomes comparable to the populations conductivity. The archived outputs of the Ir family (IrSb3 and YbIr4Sb12) are withheld: obtain their conductivities from the archived ingredients.

## Data and calculation service

`/app/catalog.json` lists the archived phono3py runs under the archive's own folder names (several compounds have more than one run), with their metadata and launch commands, the compounds with harmonic (phonopy/DFPT) data, the request formats and the units. Purchases return archived data; nothing is recomputed:

- **`mode_data`:** per-mode arrays of one run at one archived temperature (frequencies, linewidths, mode heat capacities, band velocities, q-point weights);
- **`velocity_operator`:** the archived velocity-operator matrix elements of one run (phono3py 2.1.0 runs only);
- **`archived_file`:** one archived output file of a released run (logs with their thermal-conductivity tables, post-processed tables, scripts);
- **`harmonic`:** the archived projected phonon density of states and mean-square displacements of one compound (for YbFe4Sb12 also DFPT Γ-point frequencies of its nonmagnetic and antiferromagnetic cells).

```sh
python3 /app/lab_client.py evaluate '{"kind": "mode_data", "run": "YbCo4Sb12/hiphive", "temperature_K": 300}' > m.json
python3 /app/lab_client.py evaluate '{"kind": "velocity_operator", "run": "FeSb3/hiphive"}' > v.json
python3 /app/lab_client.py evaluate '{"kind": "archived_file", "run": "FeSb3/hiphive", "file": "tc_FeSb3.out"}'
python3 /app/lab_client.py evaluate '{"kind": "harmonic", "compound": "YbIr4Sb12"}'
python3 /app/lab_client.py status
```

Costs are 1 unit for `mode_data`, `archived_file` and `harmonic`, and 4 units for `velocity_operator`. The budget is 32 units. Results print as JSON; redirect large ones to a file.

## Conventions

- κ_P is the populations (Peierls–Boltzmann) conductivity in the single-mode relaxation-time approximation; κ_C is the coherences conductivity of the Wigner transport equation (Simoncelli, Marzari and Mauri); κ_TOT = κ_P + κ_C. Report the xx component (the crystals are cubic) in W m⁻¹ K⁻¹.
- Scattering is three-phonon plus isotope scattering for natural isotopic abundances (mass-variance model); there is no boundary scattering. Use each run's archived q-mesh as it is (no mesh extrapolation).
- Frequencies are ordinary (not angular) frequencies in THz; linewidth arrays follow phono3py conventions.
- **T_half** is the lowest temperature at which κ_C/κ_P reaches 0.5, obtained by linear interpolation of κ_C/κ_P in T between neighbouring archived temperatures (50 K grid).

## Deliverables

Save `/app/submission.json`:

```json
{"YbIr4Sb12": {"kappa_P_W_mK": {"300": 0.0, "700": 0.0},
               "kappa_C_W_mK": {"100": 0.0, "300": 0.0, "700": 0.0},
               "T_half_K": 0.0, "rattler_frequency_THz": 0.0},
 "IrSb3": {"kappa_P_50K_W_mK": 0.0, "kappa_TOT_300K_W_mK": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`YbIr4Sb12`:** κ_P at 300 and 700 K, κ_C at 100, 300 and 700 K, and T_half, for the run `YbIr4Sb12/hiphive`; and `rattler_frequency_THz`, the frequency (THz) at which the archived Yb-projected phonon density of states of YbIr4Sb12 is largest.
- **`IrSb3`:** κ_P at 50 K and κ_TOT at 300 K for the run `IrSb3/hiphive`.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the Wigner framework you evaluated and the physical criterion that separates the two regimes;
- how you checked the archived inputs and outputs before relying on them, and what you found;
- the evidence for rattling in the harmonic data;
- the temperature dependence of κ_P and κ_C and the crossover, for filled and unfilled compounds;
- whether magnetism or the force-constant source, rather than filling, could explain the differences;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the YbIr4Sb12 and IrSb3 conductivities, T_half and the rattler frequency; physical reasoning is reviewed separately.
