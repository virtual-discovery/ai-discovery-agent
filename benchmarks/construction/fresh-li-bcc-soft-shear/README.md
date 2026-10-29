# Construction of fresh-li-bcc-soft-shear

Fresh-start the benchmark task (see `benchmarks/tools/BUILDING.md`): the lab runs new Quantum ESPRESSO 6.7 pw.x
calculations for (strained) one-atom bcc Li that the agent requests; the hidden truth comes from reference
calculations made through the same `task_lab.run` code path.

## Rebuild

1. `python3 -B fetch_pseudos.py --install` downloads the SSSP 1.3.0 PBE-precision archive (Materials Cloud
   rcyfm-68h65) into `~/.cache/ai-theorist/fresh-li-bcc-soft-shear/`, checks the archive md5s, extracts
   `li_pbe_v1.4.uspp.F.UPF` (GBRV 1.4, GPL v3), checks it against the SSSP metadata md5 and the pinned sha256
   (`pseudo_checksums.json`) and copies it and the archive `LICENSE.txt` (as `SSSP_LICENSE.txt`) into `environment/lab/`.
2. `python3 -B reference_jobs.py` writes `jobs_ref.json` (141 lab requests in six batches); `./run_batch.sh jobs_ref.json`
   runs every one through `environment/lab/task_lab.run` in the runtime image (`--cpuset-cpus 12,13,28,29`,
   `mpirun -np 2`), skipping finished ones; results land in `~/.cache/.../runs/<name>/result.json`. About 13 h on
   2 cores (the k56/k64 convergence runs dominate). `run()` does not validate, so the reference runs can use
   meshes denser than the lab offers (tetragonal k64, trigonal k56/k64) with the same input writer and parser.
   The parsed results are kept in `reference_results.json` so steps 3-4 can be rerun without QE.
3. `python3 -B make_task.py` writes `private/hidden.json`, `private/hidden_values.json`, `solution/ablation_values.json`
   and `calibration.json`, and prints the calibration. It uses the oracle's own analysis functions (`solution/solve.py`).
4. `python3 -B write_design.py` writes `design.json` from `design_notes.json` and `calibration.json`.
5. `.venv/bin/python -m benchmarks.tools.check_task fresh-li-bcc-soft-shear --cpuset 12,13,28,29` writes `check.json`.

## Lab

- `environment/lab/qe_driver.py`: one-atom bcc cell with primitive vectors (a/2)(-1,1,1), (a/2)(1,-1,1), (a/2)(1,1,-1),
  deformed as r' = (1 + eps) r from a Voigt strain with engineering shears; MV smearing; Gamma-centred n^3 mesh;
  ecutrho = 8 ecutwfc; conv_thr 1e-12; `mpirun -np 2 pw.x -npool 2`. Returns F (= QE "total energy") and E,
  the stress tensor from the Ry/bohr^3 block (QE sign, positive = compressive; 1e-6 GPa resolution), pressure,
  Fermi energy, irreducible k-points, symmetry operations, SCF iterations, cell, wall time. Scratch is deleted.
- `task_lab.validate`: field and range checks (a 3.20-3.60 A, |e| <= 0.04, smearing 0.005/0.01/0.02 Ry, ecut 40/60);
  the charge comes from the point symmetry of the strained lattice (counted in `qe_driver.symmetry_operations` with
  integer matrices that preserve the metric; agrees with pw.x for every class) and the mesh, from the cost table in
  `public.json`. Charges track measured wall time (about 1 unit per minute on 2 cores, minimum 1):

  | class (ops) | k16 | k24 | k32 | k40 | k48 | k56 |
  |---|---|---|---|---|---|---|
  | cubic (48) measured min | 0.1 | 0.3 | 0.7 | 1.3 | 2.2 | 3.5 |
  | tetragonal (16) | 0.3 | 0.8 | 1.7 | 3.1 | 5.7 | 8.0 |
  | trigonal (12) | 0.3 | 0.8 | 1.9 | 3.7 | 6.6 | (10.5, not offered) |

  Orthorhombic (8; e.g. a single e6 shear), monoclinic (4) and triclinic (2) costs are scaled from the counted
  irreducible k-points (14725 at orthorhombic k48 = 10 min, not offered). The slowest offered request is a
  tetragonal k56 (8.0 min measured; 9 units). ecutwfc 60 roughly doubles the time (measured 2.0-2.1x) and is offered
  where the 40 Ry cost is at most 4.

## Physics found during construction (a0 = 3.4348 A unless stated)

C' from symmetric +-1 % volume-conserving tetragonal strains (stress route) and C44 from +-1 % trigonal shear:

| mesh | 0.005 Ry | 0.01 Ry | 0.02 Ry | (0.04 Ry, design survey at 3.43 A) |
|---|---|---|---|---|
| k16 | C' -2.06, C44 6.16 | (3.43 A: 0.40) | C' 1.53, C44 10.82 | 0.96 |
| k24 | 0.78 | (3.43 A: 0.77) | 0.77, C44 11.61 | 0.94 |
| k32 | (3.43 A: 0.80) | 0.65, C44 11.55 | 0.84, C44 11.33 | 0.94 |
| k40 | (3.43 A: 0.74) | 0.87, C44 11.15 | (3.43 A: 0.85) | 0.94 |
| k48 | 0.930, C44 11.73 | 0.945, C44 11.42 | 0.851 | 0.94 |
| k56 | 0.974, C44 11.95 | 0.938, C44 11.53 | (3.43 A: 0.848) | |
| k64 | (3.43 A: 0.955); C44 11.53 | (3.43 A: 0.936); C44 11.50 | | |

- Total energies and pressures converge on the coarsest meshes (pressure at a0 within 0.01 GPa at k16/0.02 Ry;
  a0 from the k16/0.02 EOS 3.4350 A vs 3.4348 A), so they say nothing about C'.
- One-sided moduli at k48/0.01 Ry: C' = 1.27 (+1 %) and 0.63 (-1 %), i.e. the third-order term equals the
  harmonic one at 1 % strain. Symmetric differences are flat in amplitude: 0.959 (0.5 %), 0.945 (1 %), 0.945 (2 %),
  0.980 (3 %); trigonal 2 %: 11.44.
- 0.02 Ry is k-converged by k32 (0.842, 0.849, 0.851 at k32/40/48) but 10 % low; 0.04 Ry (design survey, not
  offered) is converged by k24 at 0.94, close to the converged value by accident, so it was dropped from the lab.
- Compressed cell a = 3.30 A (P = 2.042 GPa, converged at every mesh): B' = 0.809/0.784/0.787 (k48/56/64, 0.01 Ry),
  0.870/0.837 (k48/56, 0.005 Ry); B44 = 12.47/12.67/12.64 (0.01) and 12.86/13.26/12.80 (0.005). The tetragonal
  shear softens under compression while C44 and B stiffen. Energy curvature of the linear (d, d, -2d) strain gives
  1.90 GPa and of the trigonal strain 13.49 = B44 + P/2; the exactly volume-conserving tetragonal energy curvature
  gives 0.87 (no P dV term).
- Experimental lattice (3.477 A, 78 K): P = -0.48 GPa, C' 0.945 (insensitive), C44 11.01.
- Uniaxial e3 = +-1 %: C11 = 15.22, C12 = 13.28 (stress), C' = 0.971; C11 from energy 15.18, so 3(C11 - B)/4 = 0.946.
  The design's "C11 - C12 noise amplification" wrong path does not materialise at converged settings (see deviations).
- Cutoff: from 40 to 60 Ry the pressure at a0 changes by 0.002 GPa and C' (k40, 0.01 Ry) by 1e-5 GPa (0.87148 vs 0.87148); 40 Ry is converged.

## Reference values (private/hidden.json)

a0 = 3.4348 A (Birch-Murnaghan, k48/0.01 Ry; zero pressure 3.4347), B = 13.92 GPa; C' = 0.947 GPa (mean of k48/k56 x
0.01/0.005 Ry: 0.945, 0.938, 0.930, 0.974); C44 = 11.61 GPa (k48/56/64 x 0.01/0.005: 11.42, 11.73, 11.53, 11.95, 11.50,
11.53); A = 12.26; stable. Compressed a = 3.30 A: P = 2.042 GPa; B' = 0.818 GPa (k48/56/64 at 0.01, k48/56 at 0.005);
B44 = 12.78 GPa (k48/56/64 x both widths). Derived: C11 = 15.18, C12 = 13.29, C12 - C44 = +1.7 GPa.

## Decision points and calibration

| dp | fork | wrong path (ablation) | consequence |
|---|---|---|---|
| dp1 | convergence of the shear response vs total energy | k16/0.02 Ry, where E and P are converged | C' 1.53, B' 1.34, C44 10.82 |
| dp2 | small smearing on coarse meshes | k16 or k24 at 0.005 Ry | C' -2.06 (unstable) or 0.78; B' -0.26 |
| dp3 | largest smearing, k-converged | 0.02 Ry at k48 | C' 0.851, B' 0.674 |
| dp4 | one-sided vs symmetric strain on the Bain path | single +1 % or -1 % cell | C' 1.27 / 0.63, B' 1.11 / 0.51 |
| dp5 | reference volume | experimental 3.477 A | a0 +0.042, C44 11.01 (C' silent) |
| dp6 | finite pressure | energy curvature of linear strains; Lagrangian coefficients | B' 1.90 / 2.86, B44 13.49 / 14.83 |
| dp7 | engineering shear, three trigonal components | factor 1/4 or 3 | C44 2.90 / 34.8 |
| dp8 | transfer to the compressed cell | reuse a0 moduli; scale with B | B' 0.947 / 1.38, B44 11.61 / 16.9 |

| field | reference | tolerance | max robustness error | nearest failing ablation (error) |
|---|---|---|---|---|
| a0_A | 3.4348 | 0.005 | 0.0001 | experimental lattice (0.042) |
| bulk_modulus_GPa | 13.92 | 3 % | 0.02 % | none (no wrong path found) |
| Cprime_GPa | 0.947 | 0.065 | 0.033 (3 % strains) | 0.02 Ry smearing (0.096) |
| C44_GPa | 11.61 | 0.42 | 0.344 (k56/0.005) | one-sided +1 % (0.50) |
| zener_anisotropy | 12.26 | 6 % | 3.4 % | 0.02 Ry smearing (8.7 %) |
| compressed_pressure_GPa | 2.042 | 0.05 | 0.0004 | none |
| compressed_Bprime_GPa | 0.818 | 0.09 | 0.053 (k48/0.005) | reuse a0 C' (0.129) |
| compressed_B44_GPa | 12.78 | 0.56 | 0.477 (k56/0.005) | one-sided +1 % (0.65) |

Robustness variations (all pass): oracle-like k48 at both widths; k56 at 0.01 and 0.005; k64 at 0.01; k48 at 0.005;
volume-conserving energy route; 0.5 %, 2 % and 3 % symmetric strains; a0 from zero pressure; 60 Ry cutoff correction;
uniaxial-strain stress route; C' = 3(C11 - B)/4 with C11 from the uniaxial energy curvature. Full table:
`calibration.json` and `design.json` -> `calibration`.

## Build notes and deviations from the design

- **Smearing list 0.005/0.01/0.02 Ry (design: 0.005/0.01/0.02 and a survey at 0.04).** At 0.04 Ry C' is k-converged
  by k24 at 0.94 GPa, within the tolerance of the converged value by a coincidental cancellation, which would
  have made the expensive dense small-smearing ladder unnecessary for no physical reason, so 0.04 is not offered.
- **k meshes up to 56 (design: 16-48).** Cubic and tetragonal k56 are offered (3 and 9 units), trigonal only up to
  k48 (k56 is 10.5 min). The reference uses k48-k64 (k64 and trigonal k56 only in construction).
- **Graded quantities.** The design's hidden-SCF energy prediction at 0.9 V0 (+-0.02 meV/atom) was replaced by
  grading the stress-strain coefficients B', B44 and the pressure of an explicitly specified cell (a = 3.30 A,
  V/V0 = 0.887, P = 2.04 GPa). Specifying a fixed cell keeps the target unique and does not reveal a0; defining
  B' and B44 by the stress-strain relation in the prompt makes the finite-pressure route a genuine decision point
  (energy curvatures of volume-changing strains give B + P/2). `cauchy_violation_sign` was dropped as a graded flag
  (C12 - C44 = +1.7 GPa is robust and follows from graded numbers); it is a rubric item.
- **dp "C' from C11 - C12" dropped as an ablation.** With converged symmetric strains, C' from a uniaxial strain
  (stress route 0.971, or C11 from energy minus B, 0.946) is accurate, so these are robustness variations. The
  difficulty of that route is cost (tetragonal-cost cells) and the same anharmonicity.
- **Strain amplitude.** Symmetric differences are amplitude-independent between 0.5 and 3 % (the 4th-order term is
  small), so large amplitudes are robustness variations, not an ablation; the anharmonic trap is the asymmetry of
  one-sided strains (dp4).
- **Relaxed vs clamped internal coordinates** does not arise for a one-atom Bravais lattice and is not a decision point.
- **Experimental-lattice path (dp5)** is silent in C' (C' is flat near a0) and fails through a0 and C44.

## Open risks

- Small-smearing trigonal values scatter by +-0.3 GPa (C44 11.42-11.95 over k48-k56); the C44/B44 tolerances
  (0.42/0.56 GPa) are set by that spread, so they do not discriminate C44-only paths finely.
- The converged PBE C' of Li is protocol-specific; agents may anchor on the experimental 1.2 GPa or on k16/0.02-like
  values (1.5), which is the intended trap, but a lucky coarse level could land inside the C' window by accident
  (none of the levels tested does: k24/0.02 0.77, k32/0.01 0.65, k40/0.01 0.87).
- Tetragonal k56 requests take 8 min on 2 dedicated cores; on an oversubscribed host they can exceed 10 min (the
  driver timeout is 50 min).

## Results

- **Gate** (`check_task --cpuset 12,13,28,29`, 2026-10-25T17:59Z): passed, no issues; 8 decision points, 13 ablations
  all failing, 13 robustness variations all passing, 4-significant-digit rounding passes, grading deterministic;
  oracle 18 evaluations, 79 of 90 units, 0 failed; 4052 s container time (67 min of pw.x).
- **Harbor oracle** (`runs/build-fresh-li-bcc-soft-shear/oracle-r1-fresh-li-bcc-soft-shear`): reward 1.0, all 10
  numerical checks true; lab pinned by the coordinator's daemon to cores 2,18,3,19; 91 min trial wall
  (21:35-23:06Z). Oracle errors: C' -0.009, C44 +0.12, B' -0.009, B44 -0.31 GPa, A +2.1 %.
- **Disk:** peak about 70 MB in `~/.cache/ai-theorist/fresh-li-bcc-soft-shear/` (63 MB archive, deleted after the
  gate); 5 MB of QE inputs/outputs from the 141 reference runs kept there; task package 0.54 MB.
