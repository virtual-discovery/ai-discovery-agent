# Gated HgTe topological-insulator films: where does the gate put the carriers?

Thick, tensile-strained HgTe films on CdTe are three-dimensional topological insulators with a narrow bulk gap. A top gate fills and empties the top and bottom topological surface states, bulk conduction subbands and bulk valence-band hole pockets. Determine how the gate redistributes carriers among these subbands in a 70 nm film, from its raw Shubnikov–de Haas (SdH) data and from the authors' released self-consistent 8-band k·p (kdotpy) calculations. Then judge how well the self-consistent model reproduces the measured densities.

## Data and calculation service

`/app/catalog.json` describes three top-gated Hall bars (26, 70 and 107 nm HgTe) and what can be purchased for each:

- **`hall_scan`**: the low-field Hall slope and the zero-field longitudinal signal at every gate voltage.
- **`transport`**: a block of the raw lock-in record (VG, B ≤ 2 T, VI, Vxx, Vxy) for all gate voltages within ±0.025 V of a catalog centre.
- **`calculation`**: one released self-consistent kdotpy run at a total carrier density `n_nm2`. It returns the archived `.kpout` energies, the potential profile and the band-resolved integrated density of states.
- **`pocket_calculation`**: one released kdotpy run in the same self-consistent potential, restricted to an in-plane momentum wedge (`kphi_deg`).

The exact kdotpy command lines are in the catalog.

```sh
python3 /app/lab_client.py evaluate '{"kind": "hall_scan", "sample": "107nm"}'
python3 /app/lab_client.py evaluate '{"kind": "transport", "sample": "26nm", "vg_center": 0.1}'
python3 /app/lab_client.py evaluate '{"kind": "calculation", "sample": "26nm", "n_nm2": 0.002}'
python3 /app/lab_client.py evaluate '{"kind": "pocket_calculation", "sample": "26nm", "n_nm2": 0.0, "kphi_deg": [105, 165]}'
python3 /app/lab_client.py status
```

Every purchase costs 1 unit, and the budget is 25 units. These are retrieved measurements and archived calculations, not new ones.

## Conventions

- Densities are in 10^11 cm^-2, holes counted positive. Each density is that of a single spin-resolved Fermi pocket, which is what one SdH frequency measures.
- The top surface is the gate-side HgTe interface; the bottom surface is the substrate side.
- **Gate calibration.** n_tot(Vg) is the straight line fitted to the single-carrier low-field Hall density n_H = 1/(e dR_xy/dB) over 0.05 V < Vg ≤ 0.29 V, using |B| ≤ 0.2 T. It is positive for electrons. Model values are the released 70 nm runs evaluated at n_tot(Vg), interpolated linearly in density between runs.

## Deliverables

Save `/app/submission.json`:

```json
{"experiment_1e11_cm2": {"top_surface_at_+0.20V": 0.0, "bottom_surface_at_+0.05V": 0.0,
                         "first_bulk_subband_at_+0.30V": 0.0, "top_surface_at_-0.45V": 0.0,
                         "bottom_surface_at_-0.45V": 0.0, "hole_pocket_at_-0.70V": 0.0},
 "model_1e11_cm2": {"bottom_surface_at_+0.20V": 0.0, "top_surface_at_-0.45V": 0.0,
                    "hole_pocket_at_-0.60V": 0.0},
 "hole_onset_gate_V": 0.0,
 "majority_pocket_multiplicity": 0,
 "evidence_ids": [1, 2]}
```

Values are placeholders. All quantities refer to the 70 nm sample.

- **`experiment_1e11_cm2`**: densities measured by SdH at the stated gate voltages.
  - `first_bulk_subband` is the most populated bulk (non-surface) electron subband branch.
  - `hole_pocket` is one pocket of the most populated bulk valence-band hole-pocket type.
- **`model_1e11_cm2`**: the same kinds of density in the released self-consistent calculations, under the gate calibration above.
- **`hole_onset_gate_V`**: the gate voltage at which the measured majority hole-pocket density extrapolates to zero.
- **`majority_pocket_multiplicity`**: the number of symmetry-equivalent copies of that majority pocket in the in-plane Brillouin zone.
- **`evidence_ids`**: the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how the SdH frequencies were converted and assigned to subbands, and how other spectral features were recognised;
- how the top and bottom surfaces were told apart and why their gate dependences differ;
- how the hole-pocket densities and their multiplicity were obtained from the calculations and the data;
- which electrostatic and band-structure ingredients the model needs, and where it fails to match the data;
- what these data cannot determine.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the densities, the onset voltage and the multiplicity. Physical reasoning is reviewed separately.
