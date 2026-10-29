# Polar phonons in one dimension: the LO frequency of BN nanotubes near the zone centre

In a bulk polar crystal the macroscopic electric field of a longitudinal optical (LO) phonon gives a finite LO-TO splitting at Γ. Determine how this dielectric shift behaves near Γ in truly one-dimensional materials, how periodic images in a standard 3D-periodic supercell calculation change it, and what LO dispersion the armchair BN nanotubes (5,5) and (6,6) have at small wavevector along the axis when they are isolated. Answer from the authors' archived DFPT data (N. Rivano, N. Marzari, T. Sohier, npj Comput. Mater. 9, 194 (2023)): a BN atomic chain, BN nanotubes (4,4), (5,5) and (6,6), and a GaAs nanowire, each computed with standard 3D-periodic Quantum ESPRESSO and with the authors' version with a 1D Coulomb cutoff, plus bulk GaAs. For BN(5,5) and BN(6,6), the explicit 1D calculations below 0.15 |ΓZ| and the archived interpolated dispersions are withheld.

## Data and calculation service

`/app/catalog.json` quotes the archive README, lists every folder with its input files, and lists what you may buy:

- **`dfpt` (1 unit):** one archived ph.x dynamical-matrix file, as text (the Γ files also carry the dielectric tensor of the simulation cell and the Born effective charges).
- **`interpolated` (2 units):** an archived q2r/matdyn interpolated dispersion along Γ–Z (frequencies and the axial share of each mode).
- **`archived_file` (1 unit):** an archived pw.x, ph.x, q2r.x or matdyn.x output, as text.

```sh
python3 /app/lab_client.py evaluate '{"kind": "dfpt", "file": "BN_nanotubes/44/1D_code/dyn1"}' > d.json
python3 /app/lab_client.py evaluate '{"kind": "interpolated", "folder": "BN_atomic_chain/3D_code"}' > i.json
python3 /app/lab_client.py evaluate '{"kind": "archived_file", "file": "BN_nanotubes/44/3D_code/ph.out"}' > f.json
python3 /app/lab_client.py status
```

The budget is 30 units. These are retrieved author calculations; nothing is recomputed.

## Conventions

- **Wavevector:** q along the tube axis, given as a fraction of |ΓZ| = π/L, where L is the axial period of the tube.
- **LO branch:** at each q, the phonon mode with the largest longitudinal mode effective charge |Σ_a Σ_β Z*_a,zβ e_aβ / √M_a|, where e is the normalized eigenvector of the mass-scaled dynamical matrix and Z* are the Born charges (archived E-U convention, [field][displacement]) of the Γ run of the same tube and boundary conditions.
- **Boundary conditions:** "1D" means the authors' 1D-cutoff calculations (`1D_code` folders), i.e. the isolated tube; "3D" means standard 3D-periodic Quantum ESPRESSO on the archived supercell (`3D_code` folders).
- Frequencies in cm⁻¹.

## Deliverables

Save `/app/submission.json`:

```json
{"BN55_LO_1D_cm-1": {"0.01": 0.0, "0.025": 0.0, "0.05": 0.0, "0.075": 0.0, "0.125": 0.0},
 "BN66_LO_1D_cm-1": {"0.01": 0.0, "0.025": 0.0, "0.05": 0.0, "0.075": 0.0, "0.1": 0.0, "0.125": 0.0},
 "BN55_LO_3D_q0_cm-1": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`BN55_LO_1D_cm-1`, `BN66_LO_1D_cm-1`:** the LO frequency of the isolated tube (1D boundary conditions) at each listed q/|ΓZ|.
- **`BN55_LO_3D_q0_cm-1`:** the LO frequency of BN(5,5) in the 3D-periodic description of its archived supercell, in the limit q → 0 along the axis.
- **`evidence_ids`:** the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how the dielectric shift behaves near Γ in 3D, 2D and 1D, and the evidence for the form the 1D data follow;
- how and why the 3D-periodic results differ from the 1D ones, and what sets the 3D q → 0 value;
- how the BN(5,5) and BN(6,6) frequencies were obtained and validated, and how any model parameters were fixed;
- the uncertainty of the predictions and where your method stops being valid;
- the consequences for infrared and Raman spectroscopy of single tubes and of arrays;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the eleven 1D LO frequencies and the 3D q → 0 LO frequency; physical reasoning is reviewed separately.
