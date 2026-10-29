# Wurtzite AlN and (Al,B)N: how does the polarization reverse?

Wurtzite AlN becomes switchable when alloyed, and B substitution lowers the coercive field. From author-archived DFT calculations of AlN and B-substituted AlN supercells, determine how the polarization reverses: collectively, through the nonpolar layered-hexagonal (h-BN-like) structure, or by a non-collective process. Report the switching barriers, the spontaneous polarizations, and the polarization and structure along the (Al,B)N switching path, and explain what they imply for the coercive field.

## Data and calculation service

`/app/catalog.json` describes the archive (Quantum ESPRESSO 6.7, PBE). For each system it lists the authors' nudged-elastic-band runs under their folder names, the Berry-phase calculations under their folder names, the vc-relax runs, and the authors' Born effective charges. A purchase returns one archived result:

- `neb_run` (2 units): the settings and final iteration of one neb.x run, with the image energy and error table, recent activation-energy lines, the convergence message if any, and the end of the log;
- `neb_image` (1 unit): the final structure of one image of one run;
- `berry_phase` (1 unit): one Berry-phase calculation (`lberry`, `gdir = 3`), with its settings, scf total energy, structure and the polarization section of the pw.x output;
- `born_charges` (1 unit) and `vc_relax` (1 unit).

```sh
python3 /app/lab_client.py evaluate '{"kind": "neb_run", "system": "AlBN", "run": "NEB"}'
python3 /app/lab_client.py evaluate '{"kind": "berry_phase", "system": "AlN", "folder": "00_image"}'
python3 /app/lab_client.py status
```

The budget is 28 units. These are retrieved author calculations, not new ones. Results print as JSON; redirect them to files.

## Conventions

- **Formula unit (f.u.):** one cation (Al or B) plus one N.
- **Initial state:** the first image of the NEB runs of a system (all runs of a system share their end states). The final state is the same structure with the polarity reversed.
- **Barrier:** the highest energy along a path minus the energy of the initial state, in meV per f.u. of the supercell. `barrier_meV_per_fu` is the barrier along the minimum-energy path of each system. `AlBN_collective_barrier_meV_per_fu` is the barrier of (Al,B)N for homogeneous switching, in which every cation passes through its planar, layered-hexagonal site at the same time.
- **Polarization:** the component along the third cell vector (+z), in C/m².
  - `spontaneous_polarization_C_m2` is the polarization of the initial state relative to the nonpolar layered-hexagonal reference structure of the same composition and cell.
  - `AlBN_saddle_polarization_change_C_m2` is P(saddle) − P(initial) along the (Al,B)N minimum-energy path. The saddle is the highest-energy image.
- **Intermediate:** the lowest-energy local minimum strictly between the end states of the (Al,B)N minimum-energy path. Its energy is measured from the initial state, in meV/f.u.
  - A cation has inverted when the sign of its z offset from the plane of its three basal (in-plane) N neighbours is opposite to that in the initial state.
  - `inverted_cations` is the number of inverted cations at the intermediate.
- **`AlBN_switching_mode`:** `"collective"` if all cations invert together through the layered-hexagonal state along the minimum-energy path, otherwise `"non-collective"`.

## Deliverables

Save `/app/submission.json`:

```json
{"barrier_meV_per_fu": {"AlN": 0.0, "AlBN": 0.0},
 "AlBN_collective_barrier_meV_per_fu": 0.0,
 "spontaneous_polarization_C_m2": {"AlN": 0.0, "AlBN": 0.0},
 "AlBN_saddle_polarization_change_C_m2": 0.0,
 "AlBN_intermediate": {"energy_meV_per_fu": 0.0, "inverted_cations": 0},
 "AlBN_switching_mode": "collective",
 "evidence_ids": [1, 2]}
```

Values are placeholders. `evidence_ids` lists the purchases your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- which runs you trusted for each barrier, and why;
- how the polarizations were obtained from the Berry-phase output, including the choice of branch and reference structure;
- the switching mechanism in each system and the structural evidence for it;
- whether B lowers the barrier through the mechanism or through homogeneous effects such as the change of lattice;
- what the barriers imply for the coercive field, and which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the barriers, polarizations, the intermediate and the switching mode; physical reasoning is reviewed separately.
