# InMnTi2: census of the Weyl nodes of a magnetic full-Heusler compound

The inverse full-Heusler compound InMnTi2 (space group 216, no inversion centre) has been proposed as a magnetic Weyl semimetal once Hubbard corrections with first-principles U are included. Using the authors' archived calculations, establish how many Weyl nodes the magnetic state has near the Fermi level, how symmetry organises them, where they are in energy and momentum, how they compare with plain PBE, and what anomalous Hall response the magnetic state carries.

## Data and calculation service

The archive (Materials Cloud record hkprp-xw234) holds two electronic structures of the same cell, under the archive's folder names: `standard` (PBE with spin–orbit coupling) and `hubbard` (DFT+U with spin–orbit coupling and linear-response U). Inputs and selected outputs are staged unchanged in `/app/archive` (scf runs, vc-relax, non-self-consistent inputs, Wannier90 inputs and logs, the Kubo input, README). `/app/catalog.json` lists the staged files, the purchasable runs, the request formats and the conventions of the service.

Purchases return archived results or new evaluations of the two archived Wannier90 Hamiltonians (`InMnTi2_hr.dat` of each model). Wannier bands are numbered from 1 within each model.

- **`wannier_bands`** (1 unit): band energies of a Wannier model at up to 500 k-points.
- **`gap_scan`** (1 unit for mesh ≤ 24, 2 units for mesh 25–36): local minima of E_{n+1} − E_n of a Wannier model on a uniform mesh.
- **`find_nodes`** (1 unit): minimisation of E_{n+1} − E_n from up to 8 seeds.
- **`berry_flux`** (1 unit): flux of the Berry curvature of Wannier bands 1..n through up to 8 spheres.
- **`dft_bands`** (3 units): archived band energies of one non-self-consistent DFT run: `hubbard/nscf`, `hubbard/nscfXYZ`, `standard/nscf` or `standard/nscfXYZ`.
- **`kubo`** (2 units): archived Wannier90 Kubo output files of one folder `hubbard/wannier/kubo/<folder>`, with folders `10`, `30`, `58`, `120`, `210` and `300`.
- **`wannier_tools`** (2 units): archived WannierTools outputs: `hubbard/find`, `hubbard/weyl`, `standard/find` or `standard/weyl`.

```sh
python3 /app/lab_client.py evaluate '{"kind": "wannier_bands", "model": "hubbard", "k_frac": [[0, 0, 0], [0.5, 0, 0.5]]}'
python3 /app/lab_client.py evaluate '{"kind": "gap_scan", "model": "standard", "lower_band": 10, "mesh": 20}'
python3 /app/lab_client.py evaluate '{"kind": "find_nodes", "model": "hubbard", "lower_band": 10, "seeds_cart_invA": [[0.1, 0.2, 0.3]]}'
python3 /app/lab_client.py evaluate '{"kind": "berry_flux", "model": "hubbard", "occupied_bands": 10, "spheres": [{"center_cart_invA": [0.1, 0.2, 0.3], "radius_invA": 0.02, "n_theta": 16, "n_phi": 32}]}'
python3 /app/lab_client.py evaluate '{"kind": "dft_bands", "run": "standard/nscf"}'
python3 /app/lab_client.py evaluate '{"kind": "kubo", "run": "30", "components": ["S_xx", "A_xy"]}'
python3 /app/lab_client.py evaluate '{"kind": "wannier_tools", "run": "standard/find"}'
python3 /app/lab_client.py status
```

The budget is 28 units. A request that takes longer than 20 s returns a job id; poll `python3 /app/lab_client.py job/<id>`. Results print as JSON; redirect large ones to a file.

## Conventions

- **k-points:** Cartesian components in Å⁻¹ including 2π, along the axes of `CELL_PARAMETERS`; `k_frac` are coordinates in the reciprocal basis of that cell.
- **Bands:** DFT bands are counted from the bottom of each pw.x spectrum. The cell has 52 valence electrons, and the census concerns the point crossings between DFT bands 52 and 53 anywhere in the Brillouin zone.
- **Energies:** relative to the Fermi energy of the self-consistent run of the same calculation, in meV.
- **Chirality:** χ of a node is the outward flux of Ω = ∇ × A, A = i⟨u|∇u⟩, of bands 1–52 through a small sphere enclosing only that node, divided by 2π.
- **Hall conductivity:** intrinsic, at the Fermi level, ω → 0, in S/cm, with the sign convention of the archived Wannier90 Kubo outputs.

## Deliverables

Save `/app/submission.json`:

```json
{"node_count": {"dftu": 0, "pbe": 0},
 "orbit_sizes": [0],
 "node_energy_meV": 0.0,
 "pair_separation_invA": 0.0,
 "kz_offset_plus_invA": 0.0,
 "pbe_node_energy_meV": 0.0,
 "ahc_S_per_cm": {"xy": 0.0, "yz": 0.0, "zx": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`node_count`:** the number of Weyl nodes (not pairs) between bands 52 and 53, for the magnetic DFT+U state (`dftu`) and for the PBE calculation (`pbe`).
- **`orbit_sizes`:** the sizes of the sets of DFT+U nodes that are related to each other by the symmetry of the archived magnetic state, in any order.
- **`node_energy_meV`:** the energy of the DFT+U nodes of the family that contains the nodes closest to the k_z = 0 plane, for the DFT+U electronic structure itself.
- **`pair_separation_invA`:** for that family, the distance in Å⁻¹ between a node and the nearest node of opposite chirality.
- **`kz_offset_plus_invA`:** the k_z coordinate, in Å⁻¹, of the χ = +1 nodes of that family.
- **`pbe_node_energy_meV`:** the energy of the PBE nodes.
- **`ahc_S_per_cm`:** σ_yz, σ_zx and σ_xy of InMnTi2 in its archived magnetic state, converged.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- which archived calculation describes the magnetic state, and what the PBE calculation describes;
- how the Weyl nodes were found, verified and counted, and how the symmetry of the magnetic state organises them;
- how reliable each source of energies and momenta is, and which you used;
- the Weyl type, the expected interband optical onset and the role of trivial pockets;
- how the Hall conductivities were obtained and checked;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the node counts, the orbit sizes, the node energy, separation and k_z offset, the PBE node energy and the Hall conductivities; physical reasoning is reviewed separately.
