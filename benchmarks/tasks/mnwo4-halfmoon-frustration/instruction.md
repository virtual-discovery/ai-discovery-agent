# MnWO4: which competing exchanges make the paramagnetic half-moons, and is AF1 really the ground state?

Polarised neutron scattering on multiferroic MnWO4 shows half-moon-shaped diffuse scattering that bridges the magnetic Bragg positions and survives well above T_N. Natsch, Tošić, Soh and Spaldin built a classical spin-dynamics model of MnWO4 (twelve Heisenberg exchange paths plus a uniaxial single-ion anisotropy) and archived its simulated outputs. Treat that model as the system under study: its parameters are not given. From its archived outputs, determine which exchange interactions produce the half-moons and why they persist above T_N, and decide whether the same Hamiltonian makes the collinear AF1 state the ground state and how close the multiferroic spiral comes. Then predict archived outputs that the lab withholds.

## Data and calculation service

`/app/catalog.json` describes the structure, the twelve exchange paths (bond lists, without values), the easy axis, the AF1 configuration, the simulation protocols, the purchasable outputs and the withheld outputs you must predict. All purchases are the authors' archived UppASD outputs; they are not new calculations.

- **`spectrum`:** the 1.5 K dynamical structure factor (AF1 state) at one of 101 points of the archived q-path. Costs 1 unit.
- **`map`:** the static structure factor at 20, 30 or 60 K in the plane q = (−l/2, k, l). Costs 4 units.
- **`reduced_model_map`:** the 20 K static structure factor of reduced Hamiltonian 5b or 5c (see the catalog). Costs 4 units.
- **`zone_slice`:** one constant-k slice of the 11 K static structure factor on the full q grid of a 36×24×36 or a 48×24×24 supercell. Costs 1 unit.

```sh
python3 /app/lab_client.py evaluate '{"kind": "spectrum", "index": 50}'
python3 /app/lab_client.py evaluate '{"kind": "map", "temperature_K": 20}' > map20.json
python3 /app/lab_client.py evaluate '{"kind": "reduced_model_map", "model": "5c"}' > m5c.json
python3 /app/lab_client.py evaluate '{"kind": "zone_slice", "supercell": "36x24x36", "k_index": 23}'
python3 /app/lab_client.py status
```

The budget is 33 units. Map results are large: redirect them to files.

## Conventions

- **Hamiltonian:** H = −Σ_{i≠j} J_ij ê_i·ê_j + D Σ_i (ê_i·n)², with unit vectors ê, every ordered pair counted (each bond twice), n the catalog's easy axis, and J and D in meV. J < 0 is antiferromagnetic. J1…J12 label the catalog's exchange paths.
- **Propagation vector:** the wave vector q, in reciprocal-lattice units of the catalog cell, at which the largest eigenvalue of the 2×2 exchange matrix J_ss′(q) = Σ_R J e^{2πi q·(R + r_s′ − r_s)} is maximal (Luttinger–Tisza) for the infinite crystal. Vectors related by sign or by reciprocal-lattice vectors are equivalent.
- **Spiral versus AF1:** the classical energy per Mn of the flat (circular) spiral at that propagation vector, with its spin plane containing n and the optimal relative phase of the two Mn sites, minus the classical energy per Mn of the catalog's AF1 state. Both energies include the anisotropy term.
- **Map predictions:** each withheld map is compared with its archived block averages of |S| = √(Sxx² + Syy² + Szz²) by the Pearson correlation, so overall scale and offset do not matter.

## Deliverables

Save `/app/submission.json`:

```json
{"exchange_meV": {"J1": 0.0, "J2": 0.0, "J3": 0.0, "J4": 0.0, "J5": 0.0, "J6": 0.0,
                  "J7": 0.0, "J8": 0.0, "J9": 0.0, "J10": 0.0, "J11": 0.0, "J12": 0.0},
 "anisotropy_D_meV": 0.0,
 "propagation_vector_rlu": [0.0, 0.0, 0.0],
 "spiral_minus_af1_meV_per_Mn": 0.0,
 "map_14K_prediction": [[0.0]],
 "reduced_model_predictions": {"5d": [[0.0]], "5e": [[0.0]]},
 "half_moon_origin": "",
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`exchange_meV`** and **`anisotropy_D_meV`:** your best estimate of the model's parameters.
- **`map_14K_prediction`:** the 27 × 29 block values of the withheld 14 K map defined in the catalog (`to_predict.map_14K`).
- **`reduced_model_predictions`:** the 32 × 32 block values of the withheld 20 K maps of reduced models 5d and 5e (`to_predict.reduced_models`).
- **`half_moon_origin`:** one of `spin_liquid_local_constraint`, `quasi_1d_zigzag_chain`, `competing_3d_exchange`, `single_ion_anisotropy`.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the model, its conventions and how the archived outputs were modelled;
- how the parameters were determined from the purchases and how they were validated;
- how the paramagnetic structure factor was computed, and how its temperature dependence was validated;
- which interactions produce the half-moons, and why they persist above T_N;
- the ground state, the spiral–AF1 energy balance and its link to multiferroicity;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the propagation vector, the three map predictions and the half-moon origin; the exchange constants, D, the spiral–AF1 energy difference and the physical reasoning are reviewed separately.
