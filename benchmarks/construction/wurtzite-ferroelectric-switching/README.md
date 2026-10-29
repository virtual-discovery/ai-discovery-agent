# Construction of wurtzite-ferroelectric-switching

This task is a released start. It is built from the Calderon et al. (Science 2023) DFT archive, Zenodo record 7916565 (CC BY 4.0).

## Rebuild

1. `./fetch.sh` downloads `AlBN-system-Science.zip` (16.0 MB, md5 814d04e76cfb237221c961154c5e8716) into `~/.cache/ai-theorist/wurtzite-ferroelectric-switching/` and unpacks it there (94 MB).
2. `python3 -B extract.py` writes `pool.json`, which holds the parsed neb.x runs, final image structures, Berry-phase outputs, vc-relax results and Born charges. It asserts that the polarization runs are the refined_NEB2 images and that the raw P values equal the authors' table.
3. `python3 -B make_task.py` writes `catalog.json`, the lab `data.json`, `private/hidden*.json` and `solution/ablation_values.json`. It prints the reference, ablation and robustness values.
4. `.venv/bin/python -m benchmarks.tools.check_task wurtzite-ferroelectric-switching --cpuset 6,7,22,23` writes `check.json`.

## What the archive contains (findings that became decision points)

- **AlN.** All three NEB stages are converged and collective, with the same barrier: 8.368 eV per 32-atom cell, or 523 meV/f.u.
- **(Al,B)N first-stage NEB.** It is converged (path_thr 0.05 eV/Å) but collective: every cation is planar at its saddle. Its barrier is 7.354 eV, or 460 meV/f.u., only 12% below AlN. Images 1–5 and 20–25 are padding at the end-state energy.
- **(Al,B)N refined_NEB.** It never converged. The job was killed with "Disk quota exceeded" after 86 iterations, and the image errors are up to 0.55 eV/Å.
- **(Al,B)N refined_NEB2.** A 75-image restart that converged (max error 0.048 eV/Å). Its barrier is 3.198 eV, or 200 meV/f.u. The path is non-collective:
  - columns of cations invert in sequence;
  - there are four main humps, with local minima at 0-based images 16, 40, 49 and 58;
  - the lowest intermediate, image 40, has 8 of 16 cations inverted;
  - at the saddle (image 28 0-based, image 29 in neb.x numbering) the B column is exactly planar. The B column inverts second, not first.
- **Berry-phase runs.** They exist only for the refined_NEB2 images, in 0-based folders. pw.x prints P modulo 0.4723 C/m² (AlN) and 0.4814 C/m² (AlBN), the e·c/Ω quantum of the supercell.
  - AlN steps reach 0.119 C/m² per image. Nearest-branch unwrapping fails at stride 2 on AlN and at stride 6 on (Al,B)N.
  - A Born-charge estimate (the archived Z33, ground state) underestimates ΔP by 10–15%. It picks the wrong branch end to end but works over segments of ≤ 0.7 C/m².
- The nscf atom order differs from the NEB order for (Al,B)N (B is listed last). The lab serves each file's own order.

## Deviations from the catalog design

- **No fresh pw.x mode.** No zinc-blende calculation is archived, and building one would need pseudopotentials and Berry-phase runs outside the archive. The zinc-blende convention is therefore a rubric point. Its numerical analogue is ablation dp4: raw end-state Berry values of about 0.1 C/m², which resemble zinc-blende-referenced literature values.
- **The barrier ratio AlBN/AlN (0.38) is not graded separately.** It follows from the two graded barriers. The categorical mode is graded, together with two quantitative descriptors of the mechanism: the intermediate's energy and its inverted-cation count.
- **Added graded quantities.** These are the collective (Al,B)N barrier (which separates mechanism from homogeneous effects) and the polarization change at the (Al,B)N saddle (which tests branch handling, folder mapping and saddle choice).
- **The coercive field is a rubric point only.** Homogeneous estimates depend on convention (factor 2, volume per f.u. or per cell).
- **Tolerances.** Polarization ±0.015 C/m² instead of the design's ±0.05: the data give an exact value, and the nearest wrong path, the half-switched reference, is 0.027 away. Barriers ±8% instead of ±12%: the nearest wrong path, the AlN collective barrier reused for (Al,B)N, is 13.8% away.

## Calibration

| field | reference | robustness | nearest ablation | tolerance |
|---|---|---|---|---|
| barrier AlN (meV/f.u.) | 522.98 | 523.06 (first-stage run) | 261.5 (per atom) | ±8% |
| barrier AlBN (meV/f.u.) | 199.84 | 199.87 (relaxed endpoint), 199.07 (unconverged run) | 153.3 (first hump) | ±8% |
| collective barrier AlBN | 459.64 | 459.64 (QE spline) | 522.98 (AlN value) | ±8% |
| Ps AlN (C/m²) | −1.3197 | −1.3196 (midpoint reference) | −1.0835 (end-to-end Born charges) | ±0.015 |
| Ps AlBN (C/m²) | −1.3623 | same for any stride ≤ 5 or segmented Born charges | −1.3891 (half-switched reference) | ±0.015 |
| saddle ΔP AlBN (C/m²) | 0.9939 | 0.9919 (parabolic saddle) | 1.0401 / 0.9508 (off-by-one folder) | ±0.025 |
| intermediate energy (meV/f.u.) | 100.94 | 102.9 (unconverged run) | 134.7 (first minimum) | ±10% |
| inverted cations | 8 | — | 4 (first minimum) | exact |
| switching mode | non-collective | — | collective | exact |

## Oracle

The oracle makes 17 evaluations for a charge of 21 of the 28 units (75%), and runs in about 2 s in the gate container. Its route uses no literature values:

- 4 neb_run logs: AlN refined_NEB2, plus AlBN NEB, refined_NEB and refined_NEB2;
- 1 neb_image (the AlBN NEB saddle, to confirm the collective structure);
- the Born charges;
- AlN Berry folders 00, 06 and 12;
- AlBN Berry folders 00, 10, 20, 28, 40, 50, 60 and 74.

The branches come from segment Born-charge estimates. Each estimate is asserted to lie within 0.3 quanta of the chosen branch.

- **Gate (2026-10-26):** passed. All 13 ablations fail, the 4 robustness variants pass, and rounding to 4 significant digits passes.
- **Harbor oracle:** reward 1.0 (campaign `build-wurtzite-ferroelectric-switching`).

## Open risks

- **Memorization.** The qualitative mechanism (non-collective switching in (Al,B)N) is published, so the binary mode can be recalled. The numbers can't: the paper plots them, but the per-f.u. and branch-specific values here are protocol-specific. A recalled AlN Ps of about −1.35 C/m² falls outside ±0.015. A recalled value between −1.305 and −1.335 would pass that one field, but the (Al,B)N fields still need the data.
- **Folder mapping.** Off by one changes the saddle ΔP by 0.046 against a 0.025 tolerance. The mapping is discoverable from the scf energies served with each Berry result.
- **Budget.** At 28 units, buying all six NEB logs (12 units) leaves 16 units for structures and Berry runs. That is enough with segmented Born-charge or displacement estimates, but not for dense unwrapping of the (Al,B)N path.
