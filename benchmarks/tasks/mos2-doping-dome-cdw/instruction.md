# Electron-doped MoS₂ monolayer: the superconducting dome and the 1×1 → 2×2 CDW transition

Gated or intercalated MoS₂ becomes superconducting, and its Tc first rises and then falls with electron doping. From the authors' released first-principles data for the 1×1 H, 2×2 charge-density-wave (CDW) and 1T′ phases, explain why Tc rises, where the 1×1 H phase stops being the ground state, and how high Tc gets before the decline. Treat phonon softening, nonadiabatic phonon renormalization, spin–orbit coupling and the structural transition consistently.

## Data and calculation service

`/app/catalog.json` describes the released data, the simulation cells, the request formats and the methods. You can purchase:

- **`spectrum`, `dos`, `moments_table`:** Eliashberg spectral functions α²F(ω) (adiabatic A, nonadiabatic NA, fully nonadiabatic FA), electronic densities of states, and the authors' table of λ and ω_log;
- **`pes`, `soft_mode`, `mo_displacement`:** potential-energy surfaces between the 1×1 H and 2×2 CDW structures, soft-mode frequencies and relaxed Mo displacements;
- **`archived_tc_table`, `a1g_raman`, `experiment_dome`:** the authors' Tc tables, A1g frequencies and the measured dome;
- **`eliashberg`, `eliashberg_custom`:** a linearized isotropic Migdal–Eliashberg solver run on an archived spectrum or on a spectrum you supply. It reproduces the authors' archived Eliashberg results.

```sh
python3 /app/lab_client.py evaluate '{"kind": "spectrum", "phase": "1x1H", "label": "0.14", "flavor": "FA"}'
python3 /app/lab_client.py evaluate '{"kind": "eliashberg", "phase": "1x1H", "label": "0.14", "flavor": "FA", "dos": "full", "mu_star": 0.13}'
python3 /app/lab_client.py evaluate '{"kind": "pes", "label": "pes_074"}'
python3 /app/lab_client.py status
```

Eliashberg runs cost 2 units; every other request costs 1 unit. The budget is 34 units. A running calculation returns a job id; poll it with `python3 /app/lab_client.py job/<id>`. Results print as JSON; redirect large ones to a file.

## Conventions

- **Doping x:** extra electrons per MoS₂ formula unit, i.e. per 1×1 cell.
- **Tc:** the critical temperature of the linearized isotropic Migdal–Eliashberg equations as solved by the lab, with μ* = 0.13 and the fully nonadiabatic (FA) spectra. For the 1×1 H phase keep the energy dependence of the electronic DOS (`"dos": "full"`). For phases without an archived DOS use the constant-DOS limit (`"dos": "constant"`).

## Deliverables

Save `/app/submission.json`:

```json
{"tc_1x1_K": {"0.15": 0.0, "0.17": 0.0},
 "critical_doping_e_per_cell": 0.0,
 "tc_1x1_soc_at_critical_doping_K": 0.0,
 "tc_cdw_high_K": {"0.35": 0.0, "0.40": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`tc_1x1_K`:** Tc of the 1×1 H phase without spin–orbit coupling at x = 0.15 and 0.17. Their spectra and DOS are withheld.
- **`critical_doping_e_per_cell`:** the doping x_c above which the 1×1 H phase is no longer the ground state.
- **`tc_1x1_soc_at_critical_doping_K`:** the highest Tc that the 1×1 H phase reaches while it is the ground state, with spin–orbit coupling included.
- **`tc_cdw_high_K`:** Tc of the high-doping 2×2 CDW phase at x = 0.35 and 0.40.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- why Tc rises with doping in the 1×1 H phase, and how much nonadiabatic effects and spin–orbit coupling change it;
- how the withheld Tc values were estimated and checked;
- how x_c was determined and how robust it is;
- how the Tc of the high-doping CDW phase was obtained;
- the shape of the dome across the phases, compared with experiment;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the five Tc values and x_c; physical reasoning is reviewed separately.
