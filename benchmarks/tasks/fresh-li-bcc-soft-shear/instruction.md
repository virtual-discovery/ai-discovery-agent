# bcc Li: the soft tetragonal shear, elastic anisotropy and compression (PBE)

Body-centred-cubic lithium is famously soft against the tetragonal (Bain) shear that connects bcc to fcc, and it transforms martensitically at low temperature. Using new Quantum ESPRESSO calculations, determine what this computational protocol (PBE, the lab's ultrasoft pseudopotential, static lattice at T = 0) predicts at its own zero-pressure lattice parameter a0: the bulk modulus, the tetragonal shear modulus C′ = (C11 − C12)/2, the trigonal shear modulus C44, the Zener anisotropy A = C44/C′, and whether the crystal is mechanically stable. Then determine how the shear stiffnesses respond to compression, in the cubic cell with a = 3.30 Å.

## Calculation service

Every request runs a new pw.x self-consistent calculation for the one-atom bcc primitive cell. Before strain, the primitive vectors are (a/2)(−1,1,1), (a/2)(1,−1,1) and (a/2)(1,1,−1). Requests run one at a time and take from a few seconds to about ten minutes; `evaluate` returns a job id while it runs, and you poll it.

```sh
python3 /app/lab_client.py task          # this interface, the full cost table and the budget
python3 /app/lab_client.py evaluate '{"lattice_A": 3.45, "strain_voigt": [0, 0, 0.01, 0, 0, 0], "kmesh": 32, "smearing_Ry": 0.01}'
python3 /app/lab_client.py job/3
python3 /app/lab_client.py status
```

Request fields:

- `lattice_A`: cubic lattice parameter a in Å (3.20–3.60).
- `strain_voigt` (optional, default all zero): a homogeneous strain [e1, e2, e3, e4, e5, e6] in Voigt notation with engineering shear strains, i.e. ε_xx = e1, ε_yy = e2, ε_zz = e3, ε_yz = e4/2, ε_xz = e5/2, ε_xy = e6/2. Every lattice vector r becomes (1 + ε) r. Each component must satisfy |e| ≤ 0.04.
- `kmesh`: n for a Γ-centred n×n×n Monkhorst–Pack mesh of the (strained) primitive cell.
- `smearing_Ry`: Marzari–Vanderbilt cold-smearing width, 0.005, 0.01 or 0.02 Ry.
- `ecutwfc_Ry` (optional): 40 (default) or 60; the charge-density cutoff is 8× the wavefunction cutoff.

Each result gives the free energy F = E − TS (the variational functional, in Ry and eV per atom) and the internal energy E, the Cauchy stress tensor and pressure in GPa (Quantum ESPRESSO sign convention: positive = compressive), the cell vectors and volume, the Fermi energy, the number of irreducible k-points and symmetry operations, and the SCF iteration count.

Costs (units, about one per minute of computing) depend on the k mesh and on the point symmetry of the strained lattice, which sets the number of irreducible k-points. Cubic cells: k16–32 1, k40–48 2, k56 3. Tetragonal (e.g. e1 = e2 ≠ e3, no shear): k16–24 1, k32 2, k40 4, k48 6, k56 9. Trigonal (e1 = e2 = e3, e4 = e5 = e6): k16–24 1, k32 3, k40 4, k48 7. Lower symmetries cost more and allow only coarser meshes; `task` lists the full table. A 60 Ry cutoff doubles the cost. The budget is 90 units. Failed jobs stay charged.

## Conventions

- Report everything for this protocol and the static lattice, converged with respect to every numerical parameter the lab exposes (k mesh, smearing width, cutoff).
- Elastic constants are in GPa, in Voigt notation with engineering shear strains, defined at zero pressure by the energy density ΔF/V = ½ Σ C_IJ e_I e_J for small strains e of the cubic cell at a0.
- `mechanically_stable_at_a0`: whether the cubic crystal at a0 satisfies the Born stability criteria in this protocol.
- `compressed_*`: the cubic cell with a = 3.30 Å. Report its pressure and its stress–strain coefficients B′ = (B11 − B12)/2 and B44, defined by σ_I(e) − σ_I(0) = Σ_J B_IJ e_J to first order in a small strain e applied to that cell (σ is the Cauchy stress with tension positive, in Voigt order xx, yy, zz, yz, xz, xy). These are the coefficients that govern long-wavelength stability and sound propagation in the compressed crystal; at zero pressure they equal C′ and C44.

## Deliverables

Save `/app/submission.json`:

```json
{"a0_A": 0.0, "bulk_modulus_GPa": 0.0,
 "Cprime_GPa": 0.0, "C44_GPa": 0.0, "zener_anisotropy": 0.0, "mechanically_stable_at_a0": true,
 "compressed_pressure_GPa": 0.0, "compressed_Bprime_GPa": 0.0, "compressed_B44_GPa": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders. `evidence_ids` lists the calculations your numbers rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how a0 and the numerical settings were chosen, and the evidence that each graded quantity (not only the total energy) is converged;
- which strains, amplitudes and fits were used for each modulus, and why;
- how the pressure in the compressed cell was handled;
- whether bcc Li is mechanically stable in this protocol, and why C′ is so much smaller than C44;
- how the results compare with experiment, and which conclusions these calculations do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the submitted values; physical reasoning is reviewed separately.
