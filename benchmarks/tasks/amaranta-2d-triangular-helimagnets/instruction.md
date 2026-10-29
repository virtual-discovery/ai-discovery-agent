# Triangular dihalide monolayers: ferromagnets, 120° antiferromagnets or spirals?

Eight transition-metal dihalide monolayers whose magnetic ions form (nearly) triangular lattices, NiCl2, NiBr2, VCl2, VBr2, VI2, MnCl2, MnBr2 and MnI2, were treated with spin–orbit-coupled first-principles energy mapping by the AMaRaNTA workflow. Determine from energies of spin configurations which magnetic order each monolayer adopts (ferromagnet, 120° antiferromagnet, collinear stripe or incommensurate spiral) and with which wave vector, establish what decides the ordered-moment direction in VI2, predict energies of configurations you cannot buy, and state what these data imply about ordering at finite temperature.

## Data and calculation service

`/app/catalog.json` gives, for each compound, the released monolayer structure (cell, all atoms, the magnetic sites) and the archive's moment label. It also defines the supercells, the site order, the request format and the configurations whose energies you must predict (`to_predict`).

The service returns total energies (eV) of classical spin configurations (unit vectors, one per magnetic site) in a diagonal n1 × n2 supercell of a compound's catalog cell, with at most 36 magnetic sites and 1–6 configurations per request:

```sh
python3 /app/lab_client.py evaluate '{"kind": "energies", "compound": "VCl2", "supercell": [3, 3], "spins": [[[0,0,1], [0,0,1], [0,0,1], [0,0,1], [0,0,1], [0,0,1], [0,0,1], [0,0,1], [0,0,1]]]}'
python3 /app/lab_client.py status
```

The energies come from the authors' spin model, fitted to their noncollinear LSDA+U+SOC calculations, with its parameters withheld; deterministic noise of DFT-like size is added (details in the catalog). They are not new DFT calculations. Each configuration costs ceil(n_sites/12) units; the budget is 240 units.

## Conventions

- **Spins:** classical unit vectors in the Cartesian frame of the catalog cell (x along a1, z normal to the layer).
- **q_rlu:** the propagation vector of the classical ground state of the energy function the service implements, in reciprocal-lattice units of the catalog cell (b_i·a_j = 2π δ_ij; for VBr2, VI2, MnCl2, MnBr2 and MnI2 this is the rectangular cell with two magnetic sites). Any vector equivalent by a reciprocal-lattice translation, by sign or by a symmetry of the lattice is accepted.
- **type:** one of `"ferromagnet"`, `"120-degree"`, `"collinear stripe"`, `"incommensurate spiral"`.
- **Predictions:** for each `to_predict` id, (E(configuration) − E(all spins along +z)) divided by the number of magnetic sites, in meV, for the same supercell.
- **VI2_moment_direction:** the Cartesian direction of the ordered moments in the VI2 ground state (sign irrelevant).

## Deliverables

Save `/app/submission.json`:

```json
{"ground_states": {"NiCl2": {"type": "ferromagnet", "q_rlu": [0.0, 0.0]},
                   "NiBr2": {"type": "ferromagnet", "q_rlu": [0.0, 0.0]}},
 "energy_predictions_meV_per_site": {"NiCl2-P1": 0.0},
 "VI2_moment_direction": [0.0, 0.0, 1.0],
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`ground_states`:** one entry for each of the eight compounds.
- **`energy_predictions_meV_per_site`:** a value for every id in `to_predict`.
- **`evidence_ids`:** the calculations your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the model you inferred for each compound, its conventions and how it was validated;
- how the configurations and supercells were designed, and why;
- the criteria that decide each ground state, with numbers;
- what sets the VI2 moment direction, and how it was established within the noise;
- what the data imply for magnetic order at finite temperature;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the eight ground states, the predictions and the VI2 moment direction; physical reasoning is reviewed separately.
