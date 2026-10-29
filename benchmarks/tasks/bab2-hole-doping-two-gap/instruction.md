# Hole-doped BaB2: stability threshold, Eliashberg Tc and two gaps

BaB2 in the AlB2 structure has been predicted to become an MgB2-like superconductor once it is hole doped. Using the authors' archived DFPT, EPW and Migdal–Eliashberg calculations, determine the smallest hole concentration at which the lattice is dynamically stable, the isotropic and anisotropic Eliashberg Tc and the gap structure at and above that threshold, and whether replacing Ba by Cs reproduces the effect of homogeneous hole doping.

## Data and calculation service

`/app/catalog.json` lists the archived runs under the archive's own folder names: `sensitivity/0` to `sensitivity/6` (primitive cell at several hole dopings) and `331`, `331_cs` (3×3×1 supercells). It also quotes the archive README, and describes the methods, the request formats and the solver. Purchases return archived files, or lines from them, together with each run's input files:

- **`dfpt_phonons`:** harmonic frequencies at the 21 irreducible DFPT q-points of a primitive-cell run;
- **`a2f`:** the archived EPW Eliashberg spectral function file;
- **`iso_log`:** the archived isotropic Eliashberg input and log;
- **`aniso_gap`:** the archived anisotropic gap distribution at one archived temperature, on the imaginary axis or continued to the real axis;
- **`eliashberg_iso`:** a new numpy solution of the isotropic Migdal–Eliashberg equations on a run's archived α²F, at up to 8 temperatures;
- **`supercell_phonons`:** the archived supercell frequencies.

```sh
python3 /app/lab_client.py evaluate '{"kind": "dfpt_phonons", "run": "sensitivity/0"}'
python3 /app/lab_client.py evaluate '{"kind": "aniso_gap", "run": "sensitivity/0", "temperature_K": 15, "axis": "imaginary"}'
python3 /app/lab_client.py evaluate '{"kind": "eliashberg_iso", "run": "sensitivity/0", "mu_star": 0.13, "cutoff_meV": 500, "temperatures_K": [15, 20, 25]}'
python3 /app/lab_client.py evaluate '{"kind": "supercell_phonons", "run": "331"}'
python3 /app/lab_client.py status
```

Costs are 2 units for `dfpt_phonons`, `aniso_gap` and `supercell_phonons`, and 1 unit for `a2f`, `iso_log` and `eliashberg_iso`. The budget is 34 units. The DFPT results of `sensitivity/3` and `sensitivity/5` are not served. Results print as JSON; redirect large ones to a file.

## Conventions

- **Hole doping x:** holes per formula unit (one BaB2 per primitive cell), added as a homogeneous background charge.
- **Dynamically stable:** all harmonic frequencies are real, apart from the three acoustic modes at Γ. The **critical doping** x_c is the smallest x at which the primitive cell is dynamically stable, estimated from the archived DFPT results.
- **Isotropic Tc:** the temperature at which the gap Δ(iω₀) of the isotropic Migdal–Eliashberg equations vanishes, for the lab's solver on the archived α²F of that doping (default `a2f_column`), with μ* = 0.10 and `cutoff_meV` = 500.
- **Anisotropic Tc:** the temperature at which the gap of the authors' archived anisotropic Eliashberg solution for x = 0.065 closes.
- **Gap ratio:** in the archived imaginary-axis distribution of Δ_nk(iω₀) for x = 0.065 at 25 K, the gaps form two clusters. The ratio is the weight-averaged gap of the large-gap cluster divided by that of the small-gap cluster.
- **Soft mode:** the lowest harmonic frequency, in cm⁻¹ (negative for imaginary), at the wave vector where the undoped primitive cell is unstable. Report it (i) for a homogeneous hole doping of 1/9 per formula unit, and (ii) for the `331_cs` supercell, in which one Ba in nine is replaced by Cs (also 1/9 hole per formula unit), at the supercell wave vector onto which that primitive wave vector folds.

## Deliverables

Save `/app/submission.json`:

```json
{"critical_hole_doping": 0.0,
 "iso_tc_K": {"0.065": 0.0, "0.10": 0.0},
 "aniso_tc_K": 0.0,
 "gap_ratio_25K": 0.0,
 "soft_mode_cm1": {"homogeneous_one_ninth": 0.0, "cs_substituted": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`critical_hole_doping`:** x_c.
- **`iso_tc_K`:** the isotropic Tc at x = 0.065 and x = 0.10.
- **`aniso_tc_K`:** the anisotropic Tc at x = 0.065.
- **`gap_ratio_25K`:** the gap ratio at x = 0.065 and 25 K.
- **`soft_mode_cm1`:** the soft-mode frequency for homogeneous doping and for Cs substitution.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- why the lattice is unstable at low doping, which mode is involved and how x_c was located;
- what the electron–phonon coupling and Tc values below x_c mean;
- the pairing mechanism and the gap structure, and why the anisotropic and isotropic Tc differ;
- how each Tc was obtained and checked;
- whether Cs substitution can reach the stable phase, and the evidence;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks x_c, the Tc values, the gap ratio and the soft-mode frequencies; physical reasoning is reviewed separately.
