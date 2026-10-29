# CeIn3: nearest-neighbour local moments or long-range Kondo-lattice exchange?

CeIn3 is a cubic heavy-fermion antiferromagnet (T_N = 10 K, ordering wave vector R = (1/2,1/2,1/2)) whose Ce³⁺ ions carry a Γ7 crystal-field doublet. Decide from the released inelastic neutron scattering data whether its spin-wave spectrum is captured by a nearest-neighbour Heisenberg model of local Γ7 moments, or whether it requires exchange of longer range, such as the oscillatory RKKY interaction mediated by the conduction electrons. Then predict data that the lab withholds.

## Data and calculation service

`/app/catalog.json` describes the archived measurements (CNCS, T = 1.8 K), how they were reduced and binned, the coordinate conventions, the purchasable items and the withheld data you must predict. The data are the imaginary part of the dynamic susceptibility χ″(Q,E) from two settings:

- **12-meV data:** maps along R–Γ–X–M–Γ, sold as single constant-Q energy spectra.
- **High-resolution data** (E_i = 3.315 meV): maps along lines through R, sold as single constant-q spectra along R–Γ and R–X, plus constant-energy cuts along R–Γ and R–X at 0.6 and 1.4 meV.

```sh
python3 /app/lab_client.py evaluate '{"kind": "he_spectrum", "path": "XM", "index": 10}'
python3 /app/lab_client.py evaluate '{"kind": "hr_spectrum", "path": "RX", "index": 20}'
python3 /app/lab_client.py evaluate '{"kind": "hr_cut", "path": "RG", "energy_meV": 1.4}'
python3 /app/lab_client.py evaluate '{"kind": "band_structure"}'
python3 /app/lab_client.py status
```

Costs:

- a constant-energy cut: 2 units;
- every other request: 1 unit.

The budget is 60 units. These are retrieved author measurements, not new experiments.

## Conventions

- Momenta are in reciprocal-lattice units (r.l.u., 2π/a); the catalog defines the path coordinate s and the signed distance q from R.
- Spin model: H = Σ_⟨ij⟩ J_ij S_i·S_j, each bond counted once, J > 0 antiferromagnetic, with S the effective spin of the Ce³⁺ ground doublet.
- **W:** the magnon bandwidth, the maximum magnon energy over the Brillouin zone, in meV.
- **v:** the magnon velocity at R, |∇_q ω| as q → R, in meV per r.l.u.
- **η:** the dimensionless ratio v/(2πW).
- **exchange_range:** `nearest_neighbour` if J1 alone describes the spectrum; `short_range` if couplings up to third neighbours (|r| ≤ √3 a) suffice; `long_range` if couplings beyond third neighbours are required; `no_spin_waves` if the response is an incoherent continuum rather than spin waves.

## Deliverables

Save `/app/submission.json`:

```json
{"rm_map_prediction": [[0.0]],
 "rm_cut_1p4meV_prediction": [0.0],
 "integrated_chi_prediction": [0.0],
 "bandwidth_meV": 0.0,
 "velocity_bandwidth_ratio": 0.0,
 "J1_nearest_neighbour_meV": 0.0,
 "exchange_range": "",
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`rm_map_prediction`:** χ″ on the grid of `to_predict.rm_map`, as a list of rows (one per q row, one value per energy column).
- **`rm_cut_1p4meV_prediction`:** χ″ at the points of `to_predict.rm_cut_1p4meV`.
- **`integrated_chi_prediction`:** the energy-integrated χ″ at the points of `to_predict.integrated_chi`.
- **Scale of the predictions:** they are compared with the withheld data after the best non-negative overall scale. For the R–M map and cut, the best constant offset is also applied.
- **`J1_nearest_neighbour_meV`:** the exchange constant of the best nearest-neighbour (J1-only) model of the 12-meV data.
- **`bandwidth_meV`** and **`velocity_bandwidth_ratio`:** W and η of the model you conclude describes CeIn3.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the model, its conventions and how it was compared with the data;
- which measurements decide between the competing descriptions, and why;
- how the exchange range was chosen and how overfitting was avoided;
- how broadening from resolution, binning and intrinsic damping were treated;
- the physical interpretation in terms of local moments, Kondo screening and RKKY exchange;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the three predictions, W, η, J1 and the exchange range; physical reasoning is reviewed separately.
