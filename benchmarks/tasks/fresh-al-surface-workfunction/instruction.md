# Aluminium surfaces: surface energies, work functions and the Smoluchowski ordering

Determine the PBE surface energies and work functions of Al(111) and Al(100) from slab calculations that are converged with respect to slab thickness, and decide whether aluminium follows the Smoluchowski ordering φ(111) > φ(100) of a close-packed face or the reverse. Quantify the outward or inward relaxation of the Al(111) surface, and predict two slabs thicker than the calculation service can run.

## Calculation service

The lab runs new Quantum ESPRESSO 6.7 calculations (pw.x; pp.x for the potentials) on two cores: PBE, the PAW pseudopotential `Al.pbe-n-kjpaw_psl.1.0.0.UPF` (SSSP 1.3 precision), ecutwfc 30 Ry and ecutrho 240 Ry, and an fcc lattice constant a = 4.04 Å for every cell. Each request is written into a pw.x input, which the result echoes; a parameter you do not set takes the default listed below, or pw.x's own default.

```sh
python3 /app/lab_client.py evaluate '{"kind": "slab", "face": "111", "layers": 7, "kpoints": 18}'
python3 /app/lab_client.py evaluate '{"kind": "bulk", "cell": "primitive", "kpoints": 24}'
python3 /app/lab_client.py evaluate '{"kind": "bulk", "cell": "oriented", "face": "100", "kpoints": 18, "kz": 18}'
python3 /app/lab_client.py job/3
python3 /app/lab_client.py status
```

**Slab** (`"kind": "slab"`): a symmetric slab of `layers` bulk-terminated planes of the 1×1 surface cell, centred in the cell along z.

- `face`: `"111"` or `"100"`.
- `layers`: from 3 up to a maximum that depends on `kpoints` (see `limits` in `GET /task`, printed by `python3 /app/lab_client.py task`).
- `kpoints`: n of the in-plane Γ-centred n×n×1 mesh, one of 12, 18 or 24.
- `vacuum_A`: the gap between the outermost planes across the vacuum, 8–20 Å (default 15).
- `smearing_Ry`: Marzari–Vanderbilt smearing width, 0.01, 0.02 or 0.04 Ry (default 0.02).
- `nbnd`: number of Kohn–Sham bands (optional; pw.x default if omitted).
- `relax`: `"none"` (default) or `"outer_two_both_sides"`, which relaxes the z positions of the two outermost planes on each side with BFGS (the slab stays symmetric; 18×18 mesh, 5–10 layers).

A slab result gives:

- the free energy E − TS of the cell, in eV;
- the Fermi energy, the number of bands and the minimum over k of the highest computed band, in eV;
- the atomic z positions, interlayer spacings and forces (final values after a relaxation);
- planar averages on the FFT grid along z of the electrostatic potential (ionic plus Hartree) and of the total local potential (which also contains exchange–correlation), in eV;
- SCF iterations, the k-point count and the wall time.

**Bulk** (`"kind": "bulk"`): either `"cell": "primitive"`, the one-atom fcc cell with an n×n×n Γ-centred mesh (`kpoints` 8–40), or `"cell": "oriented"`, the 1×1 surface cell of `face` stacked with 3 (111) or 2 (100) planes into a periodic bulk cell with an n×n×kz mesh (`kpoints` 12, 18 or 24; `kz` 8–40). It returns the free energy per cell and per atom and the Fermi energy. Both kinds accept `smearing_Ry`.

**Costs** (units): a slab SCF costs ceil(layers × f / 4) units, at least 1, with f = 0.5, 1 and 1.6 for `kpoints` 12, 18 and 24; a relaxation costs three times its SCF; a bulk cell costs 1. For example, a 7-layer slab on the 18×18 mesh costs 2 units and a 12-layer one 3. The budget is 58 units. Calculations take from seconds to about nine minutes; while one runs, poll its `job/<id>` and do not start another.

## Conventions

- **Reference protocol.** Graded slab quantities refer to the in-plane 18×18 mesh, Marzari–Vanderbilt smearing 0.02 Ry and the settings above, in the limit of thick slabs and large vacuum.
- **Surface energy σ:** energy per unit area of one surface of the bulk-terminated (unrelaxed) slab, in J/m².
- **Work function φ:** vacuum level minus Fermi level of the bulk-terminated slab, in eV.
- **Relaxation Δd12:** the relaxed spacing between the first and second planes of Al(111) relative to the bulk spacing, in percent and positive for expansion, for the thick-slab limit of the symmetric `"outer_two_both_sides"` relaxation.
- **Thick slabs:** 18-layer Al(111) and 18-layer Al(100) slabs, bulk-terminated, computed with the reference protocol and a 15 Å vacuum. For each, predict its work function and its free energy F18 exactly as a slab result would report it, given as the offset F18 + 18 × 537.46 eV (an arbitrary anchor of −537.46 eV per atom, to keep the number short).

## Deliverables

Save `/app/submission.json`:

```json
{"surface_energy_J_per_m2": {"111": 0.0, "100": 0.0},
 "work_function_eV": {"111": 0.0, "100": 0.0},
 "work_function_ordering": "phi111>phi100",
 "relaxation_d12_percent_111": 0.0,
 "thick_slab_predictions": {"111_N18": {"free_energy_offset_eV": 0.0, "work_function_eV": 0.0},
                            "100_N18": {"free_energy_offset_eV": 0.0, "work_function_eV": 0.0}},
 "evidence_ids": [1, 2]}
```

Values are placeholders.

- **`work_function_ordering`:** `"phi111>phi100"`, `"phi100>phi111"`, or `"unresolved"` if your converged values cannot distinguish them.
- **`evidence_ids`:** the calculations your conclusions rest on.

Save `/app/report.md` (at least 400 characters) explaining:

- how σ and φ were obtained from the calculations;
- how convergence with slab thickness was established, and what limits it;
- which numerical parameters you checked and what the checks showed;
- the relaxation result, the work-function ordering and its significance;
- which conclusions these data do not support.

Then submit with `python3 /app/lab_client.py submit "$(cat /app/submission.json)"` and poll `python3 /app/lab_client.py submission` until grading completes. Submission is terminal.

The numerical grade checks the surface energies, work functions, ordering, relaxation and thick-slab predictions; physical reasoning is reviewed separately.
