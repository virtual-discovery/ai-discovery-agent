# Mn5Si3: which spin model describes the noncollinear AFM1 spin waves?

Mn5Si3 orders antiferromagnetically below about 100 K (AFM2) and enters a second phase, AFM1, below about 60 K. Two Mn Wyckoff sites, Mn1 and Mn2, may carry moments. Using author-archived first-principles exchange tables and inelastic neutron data, decide which exchange model and classical ground state reproduce the AFM1 spin waves at 10 K, and which Mn sublattice carries ordered moments and how those moments are coupled. Then predict measurements that the lab withholds.

## Data and calculation service

`/app/catalog.json` describes the orthorhombic cell (a × √3a × c), all atomic sites and the items for sale. All purchases are retrieved author-archived data, not new measurements or calculations:

- **`exchange_set`:** a KKR exchange table (bond list with Cartesian bond vectors) or the authors' Mn1 couplings.
- **`dispersion_table`:** archived 10 K branch positions along (Qh,1,0) or (0,Qk,1).
- **`constant_energy_cut`:** a raw 10 K constant-energy scan along (Qh,1,0) or (0,Qk,1). Costs 2 units.
- **`energy_scan`:** a raw energy scan (gap region at 10 K, (0,3,0) at several temperatures, (1,2,0) at 50 K in field).
- **`gap_table`:** fitted gap-mode energies near (1,2,0) at 10 K (AFM1) and 80 K (AFM2).
- **`q_scan`:** raw constant-energy scans along (Qh,2,0) across temperature and field.

```sh
python3 /app/lab_client.py evaluate '{"kind": "exchange_set", "name": "noncollinear_r4.5"}'
python3 /app/lab_client.py evaluate '{"kind": "dispersion_table", "direction": "h"}'
python3 /app/lab_client.py evaluate '{"kind": "constant_energy_cut", "direction": "k", "energy_meV": 7}'
python3 /app/lab_client.py evaluate '{"kind": "energy_scan", "name": "gap_10K_dh0.00"}'
python3 /app/lab_client.py status
```

Each purchase costs 1 unit unless stated otherwise. The budget is 25 units. Nothing measured along (0,3,Ql) is sold, and neither is the 5–20 meV energy scan at (2,1,0).

## Conventions

- Momenta Q = (H, K, L) are in reciprocal-lattice units of the orthorhombic cell. All predictions refer to T = 10 K and zero field.
- **Branch positions:** for each energy in `to_predict.l_dispersion`, the two positions Ql < 1 and Ql > 1 that a constant-energy scan along (0,3,Ql) would give when fitted with two Gaussians of common width on a constant background. This is the method behind the dispersion tables.
- **Zone-centre excitation:** the energy in meV of the strongest magnetic excitation between 5 and 20 meV in an energy scan at Q = (2,1,0) (centre of a Gaussian fitted above the background).
- **Canting:** θ = arcsin |m̂(Mn2a)·m̂(Mn1)|, averaged over the Mn2a sites, in degrees, for the zero-field classical AFM1 ground state of your model. It is the tilt of the Mn2a moments out of the plane perpendicular to the Mn1 moments. Report 0 if your model has no Mn1 moment.

## Deliverables

Save `/app/submission.json`:

```json
{"l_branch_positions_rlu": {"2": [0.0, 0.0], "20": [0.0, 0.0]},
 "zone_centre_excitation_meV": 0.0,
 "mn1_carries_moment": true,
 "mn2_canting_deg": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`l_branch_positions_rlu`:** one `[Ql below 1, Ql above 1]` pair for every energy in `to_predict.l_dispersion.energies_meV`, keyed by the energy as written there (for example `"2"`, `"11"`).
- **`mn1_carries_moment`:** whether the Mn1 sites carry an ordered moment in AFM1.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the spin Hamiltonian, its normalisation and the spin length, and how they were established;
- how the ground state was found and what it looks like;
- how spin-wave energies and neutron intensities were computed, and how the measurement (resolution, peak extraction) was modelled;
- which exchange table and which magnetic sublattices you adopted, and why;
- how well the model describes the data in each direction, and what the competing models predict;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the branch positions, the zone-centre excitation, the Mn1 answer and the canting; physical reasoning is reviewed separately.
