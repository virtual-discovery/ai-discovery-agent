# 135 kagome metals: what makes the Cr compounds strongly correlated?

Cr-based kagome metals CsCr₃Sb₅ are far more strongly correlated than their Ti and V analogues. Using the authors' released Wannier Hamiltonians and a slave-spin mean-field service, decide whether this comes from the d-shell filling (Hund's-metal physics), from a reduced kinetic energy of the d electrons, or from both, and predict how correlated the not-yet-synthesized CsCr₃Bi₅ is.

## Data and calculation service

`/app/catalog.json` describes the staged files and the service.

- `/app/archive/` holds the files of Zenodo record 17951035 for six models, one folder each (`TiSb`, `VSb`, `CrSb`, `TiBi`, `CrBi`, `CrBi_9.2`): Wannier90 inputs and outputs, DFT band structures and `Fermi_level.dat`, plus the record's example DFT input and output (`DFT/`).
- `/app/models/<model>.npz` holds each archived Wannier Hamiltonian (`*_hr.dat`), converted losslessly to numpy arrays.

The six models are CsTi₃Sb₅, CsV₃Sb₅, CsCr₃Sb₅, CsTi₃Bi₅ and two models of CsCr₃Bi₅: one at the lattice constants of CsCr₃Sb₅ and one at those of CsTi₃Bi₅. The archive does not say which folder is which.

The service runs a slave-spin mean-field (SSMF) calculation on one model at a chosen electron count, U and J/U, and returns the quasiparticle weight Z of every d Wannier function:

```sh
python3 /app/lab_client.py evaluate '{"kind": "ssmf", "model": "VSb", "electrons": 30, "U_eV": 4.0, "J_over_U": 0.18}'
python3 /app/lab_client.py status
```

Each evaluation costs 1 unit; the budget is 20 units. An evaluation takes from under a minute to several minutes; if it returns a job id, poll `python3 /app/lab_client.py job/<id>`.

## Conventions

- **Electron count:** electrons per unit cell in the 33-orbital Wannier model, both spins. The **pristine** count of a model is the number of electrons that its charge-neutral DFT calculation places in the 33-orbital manifold (an integer).
- **Protocol:** U = 6 eV and J = 0.18 U unless a field says otherwise; always at the pristine count of the model whose Hamiltonian is used, unless a field says otherwise.
- **Z̄ and Z_min:** the arithmetic mean and the minimum of Z over the five d orbitals of a transition-metal site (the three sites are equivalent).
- **Orbital names:** the d orbitals in the site-local frame of each Wannier function as defined by the Wannier90 projections of the archive (local z along c, local x as given in the projections), one of `dxy`, `dyz`, `dz2`, `dxz`, `dx2-y2`.
- **Kinetic-energy scale σ of an orbital:** the square root of the second centred moment of that orbital's projected density of states of the full 33-orbital model at U = 0.

## Deliverables

Save `/app/submission.json`:

```json
{"electron_counts": {"TiSb": 0, "VSb": 0, "CrSb": 0, "TiBi": 0, "CrBi": 0, "CrBi_9.2": 0},
 "CsCr3Sb5_U6": {"Z_mean": 0.0, "Z_min": 0.0},
 "most_renormalized_orbital": "dxy",
 "filling_share": 0.0,
 "hund_ratio": 0.0,
 "CsCr3Bi5_U4p5": {"Sb_lattice_Z_mean": 0.0, "TiBi_lattice_Z_mean": 0.0},
 "sigma_eV": {"CsCr3Bi5_Sb_lattice": 0.0, "CsTi3Sb5": 0.0},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`electron_counts`:** the pristine electron count of each model.
- **`CsCr3Sb5_U6`:** Z̄ and Z_min of pristine CsCr₃Sb₅.
- **`most_renormalized_orbital`:** the orbital with the smallest Z in pristine CsCr₃Sb₅.
- **`filling_share`:** the share of the change in ln Z̄ from pristine CsTi₃Sb₅ to pristine CsCr₃Sb₅ that is due to the electron count rather than to the Hamiltonian. Exchange the two ingredients (electron count and Hamiltonian) one at a time, and average the electron-count contribution over the two possible orders. Divide by the total change.
- **`hund_ratio`:** Z̄ of pristine CsCr₃Sb₅ at U = 6 eV and J = 0, divided by Z̄ at U = 6 eV and J = 0.18 U.
- **`CsCr3Bi5_U4p5`:** Z̄ of pristine CsCr₃Bi₅ at U = 4.5 eV (J = 0.18 U), for the model at the CsCr₃Sb₅ lattice constants and for the model at the CsTi₃Bi₅ lattice constants.
- **`sigma_eV`:** σ of the orbital named in `most_renormalized_orbital`, in CsCr₃Bi₅ at the CsCr₃Sb₅ lattice constants and in CsTi₃Sb₅.
- **`evidence_ids`:** the evaluations your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how the pristine electron counts were established, and how they compare with nominal valence counts;
- which archive folder is which CsCr₃Bi₅ model, and how you know;
- how the orbitals were identified in the local frames;
- how filling and kinetic energy were separated, and what the decomposition shows;
- what the J dependence says about Hund's versus Mott physics;
- why Bi substitution changes the correlations, separating chemistry from lattice expansion, and what happens to CsCr₃Bi₅ at U = 6 eV;
- which conclusions these calculations do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks every field above. Physical reasoning is reviewed separately.
