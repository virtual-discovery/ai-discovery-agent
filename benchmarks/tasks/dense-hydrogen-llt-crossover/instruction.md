# Dense hydrogen: liquid–liquid transition or supercritical crossover?

Compressed fluid hydrogen turns from a molecular into an atomic (dissociated) liquid near 150–200 GPa and 1000–2000 K. Small, short first-principles simulations reported pressure plateaus and density jumps, read as a first-order liquid–liquid transition ending at a critical point. Using released simulations driven by one machine-learning potential (MACE trained on PBE), decide where the transition is first order and where it is a continuous crossover, locate the crossover (Widom) line, establish what the low-density state is where the transition is first order, and quantify the effect of nuclear quantum effects.

## Data and calculation service

`/app/catalog.json` lists the released runs. They are retrieved author results, not new simulations:

- **Classical NVT molecular dynamics** with N = 128, 256, 512 and 2048 atoms, grouped by isotherm (900–1800 K) and density (Wigner–Seitz radius rs).
  - `thermo` returns a run's thermodynamic time series, its input deck and job history.
  - `structure` returns observables computed from frames sampled out of the run's released trajectory. It is available at 900–1100 K.
- **Path-integral MD** (`pimd`), with N = 512 and several independent replicas per density, archived under a temperature label.

```sh
python3 /app/lab_client.py evaluate '{"kind": "thermo", "N": 2048, "T_K": 1000, "rs": "1.43"}'
python3 /app/lab_client.py evaluate '{"kind": "structure", "N": 512, "T_K": 900, "rs": "1.425"}'
python3 /app/lab_client.py evaluate '{"kind": "pimd", "T_label": "1000K", "rs": "1.48"}'
python3 /app/lab_client.py status
```

Costs:

- `thermo`: 1 unit for N ≤ 256 and 2 units for N = 512 or 2048;
- `structure`: 2, 2, 3 and 4 units for N = 128, 256, 512 and 2048;
- `pimd`: 1 unit per density.

The budget is 135 units. Responses are large JSON objects; redirect them to files.

## Conventions

- rs is in bohr, with (4π/3)(rs a₀)³ = V/N. Report pressures in GPa.
- The specific heat per particle is c_v/k_B = (⟨E²⟩ − ⟨E⟩²)/(N k_B² T²), with E the total energy of the classical NVT run.
- **Widom line:** at a given temperature, the rs of the maximum of c_v along the N = 2048 isotherm.
- **Crossover pressure** of an isotherm: the pressure at its maximal isothermal compressibility, i.e. where P(rs) is flattest.
- **NQE shift** at temperature T: the crossover pressure of the path-integral isotherm minus that of the classical N = 2048 isotherm.
- **Classification:** the nature of the molecular-to-atomic transition at each temperature in the thermodynamic limit, either `first_order` or `crossover`.

## Deliverables

Save `/app/submission.json`:

```json
{"classification": {"900": "crossover", "950": "crossover", "1000": "crossover", "1100": "crossover"},
 "low_density_phase_900K": "liquid",
 "widom_rs_N2048": {"1100": 0.0, "1200": 0.0, "1400": 0.0},
 "pimd_crossover_rs_1000K": 0.0,
 "nqe_crossover_pressure_shift_GPa": {"900": 0.0, "1000": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`low_density_phase_900K`:** whether the molecular (low-density) side of the 900 K transition is a `crystal` or a `liquid`.
- **`pimd_crossover_rs_1000K`:** the rs at the crossover of the path-integral isotherm at 1000 K.
- **`nqe_crossover_pressure_shift_GPa`:** the NQE shift at 900 K and at 1000 K. Give `null` for a temperature where the released data cannot determine it.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how the time series were treated and c_v obtained;
- how first-order behaviour was distinguished from a crossover, and how cell size enters;
- what the low-density state at 900 K is, and why;
- why small, short simulations suggested a first-order liquid–liquid transition;
- how nuclear quantum effects change the picture;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the classification, the 900 K phase, the Widom line, the PIMD crossover and the NQE shifts; physical reasoning is reviewed separately.
