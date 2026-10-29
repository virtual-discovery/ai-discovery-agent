# Plasma-facing materials: hydrogen uptake, surface binding and what the archive supports

A divertor armour must retain little tritium and resist sputtering. From the authors' released first-principles calculations for 21 refractory candidates (W and its alternatives, including carbides, nitrides and borides), determine how strongly selected hosts dissolve hydrogen and where it sits, how strongly their surfaces hold their atoms, and whether the surface binding energy simply follows the bulk cohesive energy. Then judge which candidates combine low hydrogen uptake with high sputtering resistance, and which conclusions the archive cannot support.

## Data and calculation service

`/app/catalog.json` describes the archive (Quantum ESPRESSO, PBEsol): its README quoted verbatim, the units, and every run you may purchase, listed per material. Material keys are the archive's folder labels; for example `C-h-194` is graphite and `C-c-227` is diamond. For each material the catalog lists:

- **Hydrogen interstitials:** the pristine bulk supercell, and the relaxations of that supercell with one H atom started at different interstitial positions.
- **Surface binding:** the relaxed bulk unit cell, the isolated atoms, the relaxed slabs of every facet and termination, the slab supercells of the facet the authors selected, and unrelaxed slab supercells with one surface atom removed.

The H₂ molecule and the isolated H atom are listed separately.

A purchase returns one released run: its input settings, cell and positions, run header, energies, convergence messages and final positions.

```sh
python3 /app/lab_client.py evaluate '{"kind": "calculation", "id": "H-IFE/H2"}'
python3 /app/lab_client.py evaluate '{"kind": "calculation", "id": "SBE/WC-h-187/unitcell"}'
python3 /app/lab_client.py status
```

Every purchase costs 1 unit; the budget is 61 units. These are retrieved author calculations, not new ones. Results print as JSON; redirect large ones to a file.

## Conventions

Report energies in eV, with 1 Ry = 13.605693122994 eV.

- **Hydrogen formation energy:** E_f = E(host supercell + H) − E(host supercell) − ½E(H₂), using the relaxed archived runs. No zero-point or thermal terms. The preferred site is the archived relaxation with the lowest E_f.
- **Distinct sites:** the formation energies of the symmetry-distinct relaxed hydrogen configurations the archive contains, one value per distinct configuration, in ascending order.
- **Preferred-site geometry:** in the relaxed structure of the preferred site, the chemical symbol of the host atom nearest to H and its distance in Å, periodic images included.
- **Surface binding energy (SBE):** E_sb = E(isolated atom) + E(slab supercell with that surface atom removed, unrelaxed) − E(intact slab supercell). Use the archived surface-vacancy runs. For each element, report the smallest E_sb among the archived sites of that element on the termination with the lowest surface energy, within the facet the authors selected.
- **Cohesive energy:** [Σ E(isolated atoms) − E(bulk crystal)] / N, in eV per atom, for the relaxed bulk crystal.

## Deliverables

Save `/app/submission.json`:

```json
{"h_formation_energy_eV": {"WC-h-187": 0.0},
 "distinct_site_energies_eV": {"WC-h-187": [0.0, 0.0]},
 "preferred_site": {"WC-h-187": {"nearest_host_element": "X", "nearest_host_distance_A": 0.0}},
 "surface_binding_energy_eV": {"WC-h-187": {"W": 0.0, "C": 0.0}},
 "cohesive_energy_eV_per_atom": {"WC-h-187": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **Hydrogen fields:** give every material in `graded_materials.hydrogen`.
- **`surface_binding_energy_eV`:** give every material and element in `graded_materials.surface_binding`.
- **`cohesive_energy_eV_per_atom`:** give every material in `graded_materials.cohesive_energy`.
- **Undetermined values:** any surface-binding or cohesive entry may be `null` if you conclude that these data do not determine it.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the hydrogen reservoir and how the formation energies were assembled;
- how the sites were identified and compared, and how robust the site order is;
- how each surface binding and cohesive energy was assembled;
- whether the surface binding energy tracks the cohesive energy;
- which candidates combine endothermic hydrogen uptake with high surface binding;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the hydrogen formation energies, distinct sites, preferred-site geometry, surface binding energies and cohesive energies; physical reasoning is reviewed separately.
