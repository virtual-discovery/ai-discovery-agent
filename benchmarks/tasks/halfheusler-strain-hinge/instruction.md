# Half-Heusler LiSbZn and LiBiZn: hinge states and the archived strain models

Cubic half-Heusler LiSbZn (normal band order) and LiBiZn (Γ6–Γ8 band inversion) have been proposed as a platform for hinge states, possibly of higher-order topological origin and tunable by uniaxial strain. Using the authors' archived Wannier Hamiltonians, establish the bulk electronic structure of both compounds and what the archived models labelled "strained" contain. Then determine whether an infinite {111}-faceted prism of LiSbZn hosts in-gap hinge states, where they sit, how they disperse and converge, and whether a local perturbation can remove them.

## Data and calculation service

The archive (Zenodo record 16927281) holds four FPLO Wannier Hamiltonians (8 Wannier functions: Zn-4s and pnictogen-p spinors; Li carries none). They are served under these names:

| Name | Archived file |
|---|---|
| `LiSbZn_orig` | `wannier/+hamdata_lisbzn_eq_orig_setup.txt` |
| `LiSbZn_strained_setup` | `wannier/+hamdata_lisbzn_eq_strained_setup.txt` |
| `LiBiZn_orig` | `wannier/+hamdata_libizn_eq_orig_setup.txt` |
| `LiBiZn_strained_setup` | `wannier/+hamdata_libizn_eq_strained_setup.txt` |

`/app/catalog.json` lists each model's lattice vectors, Wannier functions and centres, the request formats and the conventions of the service. `/app/archive` holds the record description, the authors' code (`code/strained_hH.py`) and the header block of each Wannier file. The hopping blocks are evaluated by the service only. Every purchase is a new evaluation of an archived Hamiltonian:

- **`bulk_bands`** (1 unit): band energies and Zn-4s weights at up to 300 k-points.
- **`band_edges`** (1 unit for mesh ≤ 24, 2 for 25–48): the highest values of band n and the lowest of band n+1 on a uniform mesh.
- **`fermi_level`** (1 unit for mesh ≤ 32, 2 for 33–64): the energy that holds a given number of electrons per cell, with the band fillings.
- **`slab`** (1 unit): states of a slab that is finite along one basis vector, with their weight near its two faces.
- **`prism`** (1, 2 or 3 units for size N ≤ 12, 13–20 or 21–28; up to 8, 4 or 2 k-points): states of an infinite prism of N × N cells, with their weight in each corner, facet and bulk region. An onsite potential can be added to a region.

```sh
python3 /app/lab_client.py evaluate '{"kind": "bulk_bands", "model": "LiSbZn_orig", "k_frac": [[0, 0, 0], [0.5, 0, 0]]}'
python3 /app/lab_client.py evaluate '{"kind": "band_edges", "model": "LiBiZn_orig", "band": 4, "mesh": 20}'
python3 /app/lab_client.py evaluate '{"kind": "fermi_level", "model": "LiBiZn_orig", "electrons_per_cell": 2, "mesh": 24}'
python3 /app/lab_client.py evaluate '{"kind": "slab", "model": "LiSbZn_orig", "basis": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "layers": 20, "k_phases": [[0, 0]], "window_eV": [-0.5, 0.5]}'
python3 /app/lab_client.py evaluate '{"kind": "prism", "model": "LiSbZn_orig", "basis": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "size": 10, "k_phases": [0, 1.5], "window_eV": [-0.2, 0.6], "potential": [{"region": "c00", "depth": 1, "V_eV": 0.1, "species": "Zn"}]}'
python3 /app/lab_client.py status
```

The budget is 30 units. A request that takes longer than 20 s returns a job id; poll `python3 /app/lab_client.py job/<id>`. Results print as JSON; redirect large ones to a file.

## Conventions

- **Energies:** meV. For LiSbZn, energies are relative to the bulk valence-band maximum (VBM) of the same model, for the electron count of the neutral compound.
- **Strain:** δ = (c₀ − c)/c₀ in percent. Here c is the model's lattice parameter along its z axis (the tetragonal axis of a strained cell), and c₀ is the cubic lattice parameter of the `orig` model of the same compound.
- **Fundamental gap:** the minimum over the Brillouin zone of the lowest unoccupied band minus the maximum of the highest occupied band, for the neutral electron count. It is negative when the two bands overlap in energy.
- **LiBiZn Fermi level:** the charge-neutral Fermi level of `LiBiZn_orig`, relative to the fourfold Γ8 level at Γ.
- **Prism P_N:** the cubic frame is that of the `orig` files, in which the Zn site lies at (a/4)(1,1,1) from a pnictogen (X) site.
  - P_N is periodic along (a/2)[011], with a cross-section spanned by N·(a/2)[101] and N·(a/2)[110].
  - It contains every Zn and X site whose reduced coordinates along these two vectors, measured from an X site, lie in [0, N). In the `orig` lattice basis this is the identity basis matrix.
  - Report every hinge quantity in the large-N limit, that is, converged in N.
- **Hinge state:** the in-gap state at k = 0 that is localized at one hinge (an edge where two facets of the prism meet). Name its hinge by the outermost atomic species of the two facets: `Zn/Zn`, `X/X` or `Zn/X`.
- **Hinge band:** the hinge state and its continuation over the one-dimensional Brillouin zone, k ∈ [−π, π] (Bloch phase per period), including both spin-split branches.
- **Hinge potential:** a potential V added to the onsite energies of all eight orbitals of the 2 × 2 cells at the hinge that carries the state. In the service this is a corner region with depth 2 and species `all`.

## Deliverables

Save `/app/submission.json`:

```json
{"strain_percent": {"LiSbZn_orig": 0.0, "LiSbZn_strained_setup": 0.0, "LiBiZn_orig": 0.0, "LiBiZn_strained_setup": 0.0},
 "band_inverted": {"LiSbZn": false, "LiBiZn": false},
 "fundamental_gap_meV": {"LiSbZn_orig": 0.0, "LiSbZn_strained_setup": 0.0, "LiBiZn_orig": 0.0},
 "libizn_fermi_level_meV": 0.0,
 "hinge": {"termination": "Zn/X",
           "energy_k0_meV": {"LiSbZn_orig": 0.0, "LiSbZn_strained_setup": 0.0},
           "band_window_meV": [0.0, 0.0],
           "zn_s_weight_k0": 0.0,
           "removal_potential_eV": null},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`strain_percent`:** δ of each archived model.
- **`band_inverted`:** whether the Zn-4s Γ6 doublet lies below the Γ8 quartet, for each compound.
- **`fundamental_gap_meV`:** the fundamental gap of the three models named.
- **`libizn_fermi_level_meV`:** the LiBiZn Fermi level defined above.
- **`hinge.termination`:** the hinge that carries the in-gap hinge state of P_N of `LiSbZn_orig`.
- **`hinge.energy_k0_meV`:** the energy of the k = 0 hinge state of P_N, relative to the VBM of its own model. Give it for `LiSbZn_orig` and for the same physical prism built from `LiSbZn_strained_setup`.
- **`hinge.band_window_meV`:** the lowest and highest energies of the `LiSbZn_orig` hinge band, relative to the VBM.
- **`hinge.zn_s_weight_k0`:** the probability of the `LiSbZn_orig` k = 0 hinge state on the Zn-4s orbitals.
- **`hinge.removal_potential_eV`:** the hinge potential V at which the `LiSbZn_orig` k = 0 hinge level reaches the bulk VBM. Use `null` if no V between −1.5 and 1.5 eV brings it down to the VBM.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- what each archived model contains (geometry, frame, lattice basis, energy reference), and what the "strained" files are;
- the bulk band order, band inversion, gaps and insulating or semimetallic character, and how they were established;
- how the hinge state was found, separated from facet and bulk states, assigned to a hinge and converged;
- how the hinge band disperses relative to the projected bulk bands;
- whether the hinge states are of higher-order topological origin, and which calculation decides it;
- why hinge states are more visible in one compound than in the other;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the strains, the band inversions, the gaps, the LiBiZn Fermi level and the hinge quantities; physical reasoning is reviewed separately.
