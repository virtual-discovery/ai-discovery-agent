# Germanium under pressure: diamond, β-tin and the metastable decompression phases

Under compression cubic-diamond (CD) Ge transforms to β-tin (β-Sn). On decompression it usually does not return to diamond: depending on the decompression rate it ends up in metastable ST12 or BC8/R8. From author-released PBE calculations and the authors' machine-learned potential, establish the equilibrium transition pressures, where the metastable phases become favoured over β-Sn, how the metastable phases are ordered, and what thermodynamics alone can and cannot explain about the decompression products.

## Data and calculation service

`/app/catalog.json` lists the released VASP (PBE) calculations of the training database of a DeePMD potential for Ge (Fantasia et al., 2024) for five phases: CD, β-Sn (`bSn`), ST12, BC8 and R8. Each record is one archived calculation, either a single point or a relaxation. The catalog gives its phase, number of atoms, number of ionic steps and the volume per atom of its last step. Two request kinds are available:

- **`record`** (1 unit): the archived calculation, as parsed from its OUTCAR. It returns the settings (cutoff, smearing, k mesh, `PSTRESS`, `ISIF`, `IBRION`, `NSW`, …), whether the relaxation reached the required accuracy, and every ionic step. For each step you get the energy (σ→0) and free energy in eV per cell, the cell, the volume, the stress (kB, VASP order XX YY ZZ XY YZ ZX, positive = compressive), the Cartesian positions and the forces.
- **`dp`** (1 unit for up to 20 structures, at most 64 atoms each): the authors' released DeePMD potential evaluated on structures you supply. It runs in one of three modes:
  - `single`: energy, stress and forces;
  - `relax_ions`: relaxes the ions at a fixed cell;
  - `relax_cell`: relaxes the cell and ions at the given `pressure_GPa`.

```sh
python3 /app/lab_client.py evaluate '{"kind": "record", "id": "CD-098"}'
python3 /app/lab_client.py evaluate '{"kind": "dp", "structures": [{"cell": [[...],[...],[...]], "positions": [[...], ...], "mode": "relax_ions"}]}'
python3 /app/lab_client.py job/<id>        # poll a running request
python3 /app/lab_client.py status
```

The budget is 34 units. Records are retrieved author calculations, not new DFT. Potential evaluations take seconds to a few minutes. Results print as JSON; redirect them to files.

## Conventions

- All graded quantities refer to the PBE energy landscape of the released DFT calculations, at zero temperature and without zero-point energy. Every phase is taken fully relaxed at each pressure, meaning both its cell and its internal coordinates.
- **Transition pressure A→B:** the pressure at which the enthalpies per atom of A and B are equal, in GPa.
  - `CD_to_bSn` is the CD/β-Sn equilibrium.
  - `bSn_to_ST12` and `bSn_to_BC8` are the pressures below which β-Sn has a higher enthalpy than ST12 or BC8, respectively.
- **Zero-pressure energy differences:** differences between the minimum energies per atom of the fully relaxed phases, in meV/atom.
- **`metastable_order_at_0GPa`:** ST12, BC8 and R8 ordered from lowest to highest energy at 0 GPa.

## Deliverables

Save `/app/submission.json`:

```json
{"transition_pressures_GPa": {"CD_to_bSn": 0.0, "bSn_to_ST12": 0.0, "bSn_to_BC8": 0.0},
 "zero_pressure_energy_differences_meV_per_atom": {"ST12_minus_BC8": 0.0, "R8_minus_BC8": 0.0},
 "metastable_order_at_0GPa": ["ST12", "BC8", "R8"],
 "evidence_ids": [1, 2]}
```

Values are placeholders. `evidence_ids` lists the requests your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- which archived calculations you used and which you rejected, and why;
- how the energies, volumes and pressures were obtained and fitted;
- how far the released potential can be trusted here, and for what you used it;
- the sequence of phases on compression and on decompression, and which parts of the observed decompression behaviour are thermodynamic and which must be kinetic;
- how the PBE pressures relate to experiment, and which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the three transition pressures, the two energy differences and the order; physical reasoning is reviewed separately.
