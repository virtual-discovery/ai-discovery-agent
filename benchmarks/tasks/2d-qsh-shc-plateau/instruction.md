# Two-dimensional quantum spin Hall candidates: spin Hall plateaus and their spin direction

A high-throughput study of non-magnetic monolayers (Zhou, Poncé and Charlier, npj 2D Mater. Appl. 2025) released fully relativistic Wannier Hamiltonians of 216 insulators. For eight of these insulators, establish which are quantum spin Hall insulators (Z2 = 1), what in-gap plateau of the spin Hall conductivity each shows, and why the plateaus are or are not quantized. For two of them, Bi2_P-3m1 and Hg4O2_P1, also determine the direction of the spin carried by the in-gap spin Hall current.

## Data and calculation service

The eight materials are named by their archive folders in `Wannier_Hamiltonians/insulators/`. For each of them, the archived inputs and outputs are staged in `/app/archive/<material>/`: the self-consistent fully relativistic Quantum ESPRESSO run (`aiida.in`, and `aiida.out.excerpt`: the output without its final band listing) and the Wannier90 input and log (`aiida.win`, and `cut_aiida.wout.excerpt`: the log without its iteration lines). `/app/catalog.json` lists the materials, the request formats and the service conventions.

Purchases evaluate the archived Wannier Hamiltonians (`aiida_hr.dat` of each material) or return archived results:

- **`bands`** (1 unit): Wannier band energies, optionally with spin expectation values, at up to 1000 k-points.
- **`band_edges`** (1 unit): the maximum of band n and the minimum of band n + 1 on a uniform mesh.
- **`wcc`** (1–3 units, by resolution): hybrid Wannier charge centres of bands 1..n from Wilson loops.
- **`shc`** (1 unit per 40 000 k-points, plus 1 unit per extra spin component): the Kubo spin Hall conductivity σ^γ_xy(μ) of a Wannier model, on a chosen mesh (optionally with adaptive refinement), broadening, velocity-operator convention and list of chemical potentials.
- **`ribbon`** (2 units): the spectrum of a finite-width ribbon with edge weights.
- **`archived_shc`** (1 unit): one archived Wannier90 spin Hall Fermi scan of the study's `Top_materials/` folder, with its Wannier90 input.

```sh
python3 /app/lab_client.py evaluate '{"kind": "band_edges", "material": "Bi2_P-3m1", "band": 10, "mesh": 120}'
python3 /app/lab_client.py evaluate '{"kind": "wcc", "material": "Bi2_P-3m1", "bands": 10, "loop_axis": 1, "k": [0, 0.1, 0.2, 0.3, 0.4, 0.5], "n_loop": 60}'
python3 /app/lab_client.py evaluate '{"kind": "shc", "material": "Bi2_P-3m1", "spin": ["z"], "mesh": [60, 60], "eta_eV": 0.005, "position": "<lattice or wannier_centres>", "mu_eV": [-3.5]}'
python3 /app/lab_client.py evaluate '{"kind": "archived_shc", "folder": "Top_materials/Bi2/SHC-xyz"}'
python3 /app/lab_client.py status
```

The budget is 40 units. A request that takes longer than 20 s returns a job id; poll `python3 /app/lab_client.py job/<id>` (a dense spin Hall calculation can take several minutes). Results print as JSON; redirect large ones to a file.

## Conventions

- **k-points:** reduced coordinates of the reciprocal basis of `CELL_PARAMETERS` in each `aiida.in`; k_z = 0.
- **Energies:** eV, on the energy scale of the archived calculation of each material.
- **Spin Hall conductivity:** σ^γ_xy is the response of the spin current flowing along x, with spin along γ, to an electric field along y. x, y and z are the Cartesian axes of the archived cell. Signs follow the service output.
- **Plateau:** the intrinsic spin Hall conductivity at T = 0, in the clean limit (broadening η → 0), with the chemical potential inside the bulk gap of the archived Wannier model, converged in k. Report it per sheet in units of (ħ/e)(e²/h) = e/2π, in which an insulator that conserves S_z and has one pair of helical edge states has |σ^z_xy| = 1.
- **Spin direction:** the unit vector along (σ^x_xy, σ^y_xy, σ^z_xy) of the plateau, given by its polar angle θ from +z (0° to 180°) and its azimuth φ from +x towards +y (−180° to 180°).
- All graded values refer to the archived Wannier Hamiltonians of `Wannier_Hamiltonians/insulators/` as evaluated by the service.

## Deliverables

Save `/app/submission.json`:

```json
{"z2": {"Bi2_P-3m1": 0, "Hf2Br2_P-3m1": 0, "CuLi2As_P-6m2": 0, "Ti2I2N2_Pmm2": 0,
        "Li2Tl2_P3m1": 0, "Bi2In2_P4_nmm": 0, "Hg4O2_P1": 0, "KSnAs_P3m1": 0},
 "plateau_sigma_z": {"Bi2_P-3m1": 0.0, "Hf2Br2_P-3m1": 0.0, "CuLi2As_P-6m2": 0.0, "Ti2I2N2_Pmm2": 0.0,
                     "Li2Tl2_P3m1": 0.0, "Bi2In2_P4_nmm": 0.0, "Hg4O2_P1": 0.0, "KSnAs_P3m1": 0.0},
 "spin_direction": {"Bi2_P-3m1": {"theta_deg": 0.0, "phi_deg": 0.0},
                    "Hg4O2_P1": {"theta_deg": 0.0, "phi_deg": 0.0}},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`z2`:** the Z2 invariant (0 or 1) of the occupied bands of each archived Wannier model.
- **`plateau_sigma_z`:** the plateau of σ^z_xy of each material, in (ħ/e)(e²/h).
- **`spin_direction`:** θ and φ of the plateau spin direction for Bi2_P-3m1 and Hg4O2_P1.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how the topological invariant of each material was established and checked;
- where the chemical potential was placed, and why;
- how the spin Hall conductivities were obtained, converted and converged;
- why each plateau is or is not quantized, and what sets the spin direction;
- which conclusions these data and models do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the invariants, the plateaus and the spin directions; physical reasoning is reviewed separately.
