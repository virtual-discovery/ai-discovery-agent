# Superconducting hydrides: how strong is the Coulomb repulsion, and what Tc follows?

First-principles Tc predictions for hydrides usually assume a Coulomb pseudopotential μ* ≈ 0.10–0.16. From the authors' released calculations of the static RPA-screened Coulomb interaction W(ε,ε′) and the Eliashberg Tc of Mg₂IrH₆ and H₃S at several pressures, work out how large the Coulomb repulsion really is, why it differs between compounds and pressures, and what Tc follows at pressures where the screened interaction is not available.

## Data and calculation service

`/app/catalog.json` lists the archived calculations (system, pressure), what can be bought for each, the request formats and the solver conventions. You can purchase:

- **`inputs`:** the Eliashberg spectral function α²F(ω), the electronic DOS file, the cell and the target pressure;
- **`coulomb_table`, `w_kernel`:** the Mantel tables of μ versus electronic smearing, and the screened interaction W(ε,ε′);
- **`archived_run`:** the summary and settings of one archived IsoME Eliashberg run;
- **`eliashberg`:** a new solution of the linearized isotropic Migdal–Eliashberg equations on archived inputs, with a fixed μ* (constant or energy-dependent DOS) or with the archived W(ε,ε′).

```sh
python3 /app/lab_client.py evaluate '{"kind": "inputs", "system": "Mg2IrH6", "pressure_GPa": 20}'
python3 /app/lab_client.py evaluate '{"kind": "eliashberg", "system": "Mg2IrH6", "pressure_GPa": 20, "solver": "vdos", "mu_star": [0.1, 0.16], "omega_c_meV": 7000}'
python3 /app/lab_client.py evaluate '{"kind": "eliashberg", "system": "H3S", "pressure_GPa": 200, "solver": "vdos_W"}'
python3 /app/lab_client.py status
```

`eliashberg` costs 2 units; every other request costs 1 unit. The budget is 38 units. A running calculation returns a job id; poll it with `python3 /app/lab_client.py job/<id>`. Results print as JSON; redirect large ones to a file.

At Mg₂IrH₆ 15 and 25 GPa and H₃S 180 and 190 GPa, only `inputs` and the fixed-μ* solvers are available.

## Conventions

- **Full-W Tc:** the Tc of the isotropic full-bandwidth Eliashberg equations with the energy-dependent DOS and the static screened interaction W(ε,ε′), with the archive's settings (the α²F column and DOS file served by `inputs`, Matsubara cutoff 7 eV). The lab's `vdos_W` solver is this calculation. The graded values are the authors' archived IsoME results.
- **μ\*:** in `cdos` and `vdos`, a constant Coulomb pseudopotential applied to all Matsubara frequencies up to the cutoff `omega_c_meV`.

## Deliverables

Save `/app/submission.json`:

```json
{"tc_full_W_K": {"Mg2IrH6_15GPa": 0.0, "Mg2IrH6_25GPa": 0.0, "H3S_180GPa": 0.0, "H3S_190GPa": 0.0},
 "mu_star_cdos_1eV_Mg2IrH6_20GPa": 0.0,
 "N_EF_H3S_180GPa_per_fu_per_spin": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`tc_full_W_K`:** the full-W Tc at the four pressures where W is not available to you.
- **`mu_star_cdos_1eV_Mg2IrH6_20GPa`:** the μ\* that, in the constant-DOS solver with a Matsubara cutoff of 1 eV, gives the full-W Tc of Mg₂IrH₆ at 20 GPa.
- **`N_EF_H3S_180GPa_per_fu_per_spin`:** the density of states at the Fermi level of H₃S at 180 GPa from its DOS file, in states per eV per spin per H₃S formula unit.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how large μ is in each compound and what controls it across compounds and pressures;
- the difference between μ and μ\*, and why conventional μ\* values fail here;
- how the Tc at the withheld pressures were obtained and checked;
- which approximations and data limit these conclusions.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the four Tc, μ\* and N(E_F); physical reasoning is reviewed separately.
