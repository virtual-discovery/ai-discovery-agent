# θ-TaN under pressure: which scattering channels control the lattice thermal conductivity?

Semimetallic θ-TaN (WC-type, two atoms per cell) has an unusually large lattice thermal conductivity, and its pressure dependence is not monotonic. Using the authors' archived first-principles force constants and phonon–electron scattering rates, determine which scattering channels (three-phonon, four-phonon, isotope, phonon–electron) control the conductivity and its pressure dependence, and predict the conductivity where parts of the archive are withheld.

## Data and calculation service

`/app/catalog.json` lists the six archived pressures (0–100 GPa), what the archive contains, what is withheld (the fourth-order force constants at 40 and 80 GPa, the phonon–electron table at 80 GPa), the 16×16×16 irreducible q-point list, the request formats and the units. The lab is a phonon Boltzmann-transport service built on the archived force constants (ShengBTE/FourPhonon conventions):

- **`archived_file`** (1 unit): one archived text file of one pressure, verbatim (POSCAR, CONTROL, README.md, and BTE.sr_phel_metal where released).
- **`harmonic`** (1 unit): phonon frequencies, group velocities and density of states on the irreducible points of an n×n×n grid.
- **`rates`** (2 units): mode-resolved scattering rates on the irreducible points: three-phonon (computed for your mesh, broadening and temperature), isotope, and the precomputed four-phonon rates where they are available.
- **`kappa`** (1 unit; 2 units for meshes above 16): the conductivity tensor for your choice of pressure, mesh, temperature, Gaussian broadening factor (`scalebroad`), scattering channels and solver (`RTA` or `iterative`). On the 16×16×16 grid you may add your own diagonal scattering rates (`extra_rates_ps`, 1/τ in ps⁻¹, catalog irreducible order, 6 branches sorted by frequency).

```sh
python3 /app/lab_client.py evaluate '{"kind": "archived_file", "pressure_GPa": 20, "file": "CONTROL"}'
python3 /app/lab_client.py evaluate '{"kind": "rates", "pressure_GPa": 20, "mesh": 16, "T_K": 300, "scalebroad": 0.1}' > r.json
python3 /app/lab_client.py evaluate '{"kind": "kappa", "pressure_GPa": 20, "mesh": 12, "T_K": 300, "scalebroad": 1.0, "channels": ["3ph", "isotope"], "solver": "RTA"}'
python3 /app/lab_client.py job/7
python3 /app/lab_client.py status
```

A `kappa` or `rates` call takes up to a few minutes; if it returns a running job, poll `job/<id>`. The budget is 30 units. Results print as JSON; redirect large ones to a file.

## Conventions

- Conductivities are lattice (phonon) thermal conductivities in W m⁻¹ K⁻¹; xx lies in the hexagonal basal plane and zz along c.
- "All channels" means three-phonon, four-phonon, isotope (natural abundance) and phonon–electron scattering, evaluated with the authors' archived calculation settings, solved as the authors did: the converged iterative solution for three-phonon scattering, with four-phonon, isotope and phonon–electron scattering as diagonal (relaxation-time) terms.
- The archived phonon–electron rates are used unchanged at every temperature.
- The phonon–electron share of the thermal resistivity is 1 − κ_xx / κ_xx(all channels except phonon–electron), at the same pressure and temperature.

## Deliverables

Save `/app/submission.json`:

```json
{"kappa_xx_80GPa_300K_W_mK": 0.0,
 "kappa_zz_80GPa_300K_W_mK": 0.0,
 "kappa_xx_40GPa_300K_W_mK": 0.0,
 "kappa_xx_60GPa_600K_W_mK": 0.0,
 "phel_resistivity_fraction_80GPa_300K": 0.0,
 "evidence_ids": [1, 2]}
```

Values are placeholders. All conductivities are for all channels. `evidence_ids` lists the lab jobs your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- the role of each scattering channel and why the pressure dependence is non-monotonic;
- how you established the calculation settings and checked their adequacy;
- how you read and checked the archived phonon–electron tables before relying on them;
- how you obtained the withheld ingredients at 40 GPa, 80 GPa and 600 K, and how uncertain they are;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the five values; physical reasoning is reviewed separately.
