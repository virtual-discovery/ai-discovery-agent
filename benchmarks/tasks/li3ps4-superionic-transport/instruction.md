# Li3PS4: the superionic transition, ionic conductivity and the paddle-wheel question

Lithium ortho-thiophosphate (Li3PS4) is a solid electrolyte whose high-temperature phases conduct Li⁺ fast. Machine-learning molecular dynamics (ML-MD) with potentials trained on three DFT references (PBEsol, r2SCAN and PBE0) has been released for runs started from the ideal α, β and γ structures between 200 and 1200 K. From these data, locate the superionic transition for each reference and quantify the ionic conductivity, its ion–ion correlations and its activation energies. Then decide whether a paddle-wheel coupling of PS4 rotations to Li hops is what makes the material superionic.

## Data and calculation service

`/app/catalog.json` describes the simulations (768-atom cell with 288 Li, 96 P and 384 S, NpT at zero pressure), what can be purchased and the cost of each request. Purchases return the authors' archived analyses of the released runs; they are not new simulations.

```sh
python3 /app/lab_client.py evaluate '{"kind": "conductivity", "functional": "PBE0", "initial_phase": "gamma", "T_K": 700, "estimator": "green_kubo"}'
python3 /app/lab_client.py evaluate '{"kind": "sublattice_dynamics", "functional": "r2SCAN", "T_K": 650}'
python3 /app/lab_client.py evaluate '{"kind": "heat_capacity", "functional": "PBEsol", "T_low_K": 550, "T_high_K": 575}'
python3 /app/lab_client.py status
```

The other request kinds are `li_msd`, `current_breakdown` (with a `frame`), `li_pair_correlation` and `order_parameters`; the catalog documents them.

Costs:

- a Nernst–Einstein conductivity: 1 unit;
- a Green–Kubo conductivity: 2 units;
- `li_msd`: 2 units;
- every other request: 1 unit.

The budget is 60 units. Results print as JSON; redirect large ones to a file.

## Conventions

- σ is the ionic (electrical) conductivity in S/m. The Haven ratio is H_R = σ_NE/σ, where σ_NE is the Nernst–Einstein conductivity.
- Activation energies E_a are defined by σT ∝ exp(−E_a/k_BT).
- The transition temperature of a reference is the temperature at which its simulations started from γ transform from the γ structure into the superionic state. When it falls between two simulated temperatures, report their midpoint.
- Hop rates are per Li ion; flip rates are per PS4 tetrahedron, where a flip is a reorientation of the tetrahedron.

## Deliverables

Save `/app/submission.json`:

```json
{"transition_temperature_K": {"PBE0": 0.0, "r2SCAN": 0.0},
 "sigma_predictions_S_per_m": {"PBEsol/gamma/625": 0.0, "PBE0/gamma/750": 0.0, "r2SCAN/gamma/800": 0.0},
 "haven_ratio_PBE0_gamma_750K": 0.0,
 "activation_energy_eV_PBE0_gamma": {"below_transition": 0.0, "superionic_to_1000K": 0.0},
 "li_s_cross_fraction": 0.0,
 "log10_hop_to_flip_rate_ratio_PBE0_800K": 0.0,
 "paddle_wheel_supported": null,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`sigma_predictions_S_per_m`:** σ of the runs listed in `to_predict` (functional/initial phase/temperature). These runs are not served.
- **`haven_ratio_PBE0_gamma_750K`:** H_R of the ML-PBE0 run started from γ at 750 K.
- **`activation_energy_eV_PBE0_gamma`:** E_a of σ for the ML-PBE0 runs started from γ, below the transition and from the transition up to 1000 K.
- **`li_s_cross_fraction`:** the fraction of σ that arises from dynamical Li–S coupling, from the released current breakdown.
- **`log10_hop_to_flip_rate_ratio_PBE0_800K`:** log10 of the Li hop rate over the PS4 flip rate in the ML-PBE0 runs started from γ at 800 K. An order-of-magnitude estimate is enough.
- **`paddle_wheel_supported`:** `true` if PS4 rotations directly assist individual Li hops in the superionic phase, otherwise `false`.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how the transition was located for each reference;
- which conductivity estimator answers which question, and how the predictions and activation energies were obtained;
- the role of ion–ion correlations;
- the mechanism test and its verdict;
- how the DFT reference changes the results;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the transition temperatures, predictions, Haven ratio, activation energies, Li–S fraction, rate ratio and verdict; physical reasoning is reviewed separately.
