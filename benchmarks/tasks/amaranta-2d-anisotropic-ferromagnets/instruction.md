# Cr-based 2D magnets: where does the magnetic anisotropy come from?

Six Cr-based monolayers often discussed as two-dimensional ferromagnets, CrCl3, CrBr3, CrI3, CrGeTe3, CrSiSe3 and CrSBr, were treated with spin–orbit-coupled first-principles energy mapping by the AMaRaNTA workflow. From energies of spin configurations, determine each monolayer's classical magnetic ground state, the direction of its ordered moments and the anisotropy energy that pins them, and how much of that anisotropy comes from the single-ion term rather than from the exchange interaction. Determine the anisotropic nearest-neighbour exchange of CrI3 in the frame of its octahedra and the ordering temperature of the CrI3 model, the antisymmetric exchange and the spin-wave gap of CrSBr, and predict energies of configurations you cannot buy.

## Data and calculation service

`/app/catalog.json` gives, for each compound, the released monolayer structure (cell, all atoms, the magnetic sites) and the archive's moment label. It also defines the supercells, the site order, the request format and the configurations whose energies you must predict (`to_predict`).

The service returns total energies (eV) of classical spin configurations (unit vectors, one per magnetic site) in a diagonal n1 × n2 supercell of a compound's catalog cell, with at most 36 magnetic sites and 1–6 configurations per request:

```sh
python3 /app/lab_client.py evaluate '{"kind": "energies", "compound": "CrBr3", "supercell": [1, 1], "spins": [[[0,0,1], [0,0,1]]]}'
python3 /app/lab_client.py status
```

The energies come from the authors' spin model, fitted to their noncollinear LSDA+U+SOC calculations, with its parameters withheld; deterministic noise of DFT-like size is added (details in the catalog). They are not new DFT calculations. Each configuration costs ceil(n_sites/12) units; the budget is 245 units.

## Conventions

- **Spins and frame:** classical unit vectors e_i in the Cartesian frame of the catalog cell (x along a1, z normal to the layer). E denotes the service's energy function of the spins (without its arbitrary constant). A bond's coupling matrix J_ij is defined by its contribution e_i·J_ij·e_j to E, each pair counted once.
- **`ground_states`:** one of `"ferromagnet"`, `"Neel antiferromagnet"` (all nearest neighbours antiparallel), `"zigzag antiferromagnet"`, `"stripe antiferromagnet"`, `"noncollinear"`.
- **`easy_axis_polar_deg`:** the angle in degrees between the ordered moments of the classical ground state and the layer normal (0: perpendicular to the layer; 90: in the plane).
- **`anisotropy_meV_per_Cr`:** take the ground-state pattern (every moment ±n, relative signs fixed) and rotate it rigidly. `total` is [E(n = x) − E(n = z)] per Cr. `single_ion` is the part of `total` that comes from the on-site terms of E, the energy of each Cr moment as a function of its own direction.
- **`CrI3_bond_exchange_meV`:** for a nearest-neighbour Cr–Cr bond of CrI3, write the symmetric part of its coupling matrix as J·1 + K e3e3ᵀ + Γ(e1e2ᵀ + e2e1ᵀ) + Γ′(e1e3ᵀ + e3e1ᵀ + e2e3ᵀ + e3e2ᵀ). Here e1, e2, e3 are the ideal cubic axes of the CrI6 octahedra: unit vectors at arccos(1/√3) from +z, with azimuths 120° apart, each exactly perpendicular to one of the three bonds of a Cr. Of the two such sets, take the one nearest to the Cr→I directions toward the three I atoms above that Cr. e3 is the axis perpendicular to the bond. Report J, K, Gamma and Gamma_prime.
- **`CrI3_Tc_K`:** the ordering temperature, in K, of the classical model E (unit-vector spins) in the thermodynamic limit; k_B = 0.086173 meV/K.
- **`CrSBr_dm_vector_meV`:** the vector D of the bond from magnetic site 0 of the catalog cell to the same site translated by a1. It is defined by the antisymmetric part of that bond's coupling matrix: e_i·[(J − Jᵀ)/2]·e_j = D·(e_i × e_j), with i the site at the origin of the bond.
- **`CrSBr_spin_wave_gap_meV`:** the lowest zone-centre magnon energy of linear spin-wave theory about the classical ground state. Use spins of length S = 3/2 with the quantum Hamiltonian H = E({S_i/S}), so that classical energies are unchanged. Treat single-ion terms in the same (classical) way as exchange, as SpinW does. ħ = 1; energies in meV.
- **Predictions:** for each `to_predict` id, (E(configuration) − E(all spins along +z)) divided by the number of magnetic sites, in meV, for the same supercell.

## Deliverables

Save `/app/submission.json`:

```json
{"ground_states": {"CrCl3": "ferromagnet"},
 "easy_axis_polar_deg": {"CrCl3": 0.0},
 "anisotropy_meV_per_Cr": {"CrCl3": {"total": 0.0, "single_ion": 0.0}},
 "CrI3_bond_exchange_meV": {"J": 0.0, "K": 0.0, "Gamma": 0.0, "Gamma_prime": 0.0},
 "CrI3_Tc_K": 0.0,
 "CrSBr_dm_vector_meV": [0.0, 0.0, 0.0],
 "CrSBr_spin_wave_gap_meV": 0.0,
 "energy_predictions_meV_per_site": {"CrI3-P1": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`ground_states`, `easy_axis_polar_deg`, `anisotropy_meV_per_Cr`:** one entry for each of the six compounds.
- **`energy_predictions_meV_per_site`:** a value for every id in `to_predict`.
- **`evidence_ids`:** the calculations your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the model you inferred for each compound, its conventions and how it was validated;
- how the configurations and supercells were designed, and why;
- the ground states, with the criteria and numbers that decide them;
- for each compound, how the anisotropy splits into single-ion and exchange parts, how the split was measured, and the trend with the ligand;
- the CrI3 bond exchange and the frame used;
- how the ordering temperature was obtained, and what role the anisotropy plays in it;
- the antisymmetric exchange of CrSBr, and the anisotropy that sets its gap;
- which conclusions about the real materials these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the ground states, easy axes, anisotropy energies, the CrI3 bond exchange (J, K, Γ) and ordering temperature, the CrSBr DM vector and gap, and the predictions; physical reasoning is reviewed separately.
