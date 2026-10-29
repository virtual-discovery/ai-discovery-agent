# Source and transformations

This is a **fresh-start** task: every number the agent sees comes from a new Quantum ESPRESSO 6.7 `pw.x` calculation that the lab runs on request. No archived calculation results are served or used. The only external inputs are two pseudopotential files.

The physics question is paired, at the level of the scientific question only, with the released r2SCAN lattice-dynamics study of the group-IV elements: Haxhijaj, Riemelmoser, Pasquarello, *Lattice dynamics and structural phase stability of group IV elemental solids with the r2SCAN functional*, Phys. Rev. B (2026), doi:10.1103/znkf-3g37, Materials Cloud Archive record z7tes-s0y64 (DOI 10.24435/materialscloud:dp-1f, CC BY 4.0). No file of that record is used: the task grades protocol-specific PBE values computed here.

## Pseudopotentials

| file | source | valence | md5 | sha256 | license |
|---|---|---|---|---|---|
| `ge_pbe_v1.4.uspp.F.UPF` | GBRV 1.4 (K. F. Garrity, J. W. Bennett, K. M. Rabe, D. Vanderbilt, Comput. Mater. Sci. 81, 446 (2014)); the SSSP 1.3.0 PBE precision choice for Ge (recommended cutoffs 45/360 Ry) | 3d10 4s2 4p2 (Z = 14) | 9c9eaa91e581c3f09632fb3098b2c6b2 (matches the SSSP metadata) | 1b95046ab0f3228e0404719ca66d0309af5a62313be0b22171b963566c90430b | GNU GPL v3 (GBRV family, per the SSSP archive `LICENSE.txt`) |
| `Ge.pbe-n-rrkjus_psl.1.0.0.UPF` | PSlibrary 1.0.0 (A. Dal Corso, Comput. Mater. Sci. 95, 337 (2014)), from the Quantum ESPRESSO pseudopotential library, https://pseudopotentials.quantum-espresso.org/upf_files/Ge.pbe-n-rrkjus_psl.1.0.0.UPF | 4s2 4p2, 3d frozen in the core with nonlinear core correction (Z = 4) | a23a7f3f7456476e534e54484b27d8e8 | b756d20fe98992cb3afef5fcde78515720a89336d412c13b0ceddf794d54cfcd | GNU GPL v2 or later |

- **SSSP archive:** G. Prandini, A. Marrazzo, I. E. Castelli, N. Mounet, N. Marzari, *Precision and efficiency in solid-state pseudopotential calculations*, npj Comput. Mater. 4, 72 (2018); Materials Cloud Archive record rcyfm-68h65, DOI 10.24435/materialscloud:f3-ym, https://archive.materialscloud.org/records/rcyfm-68h65. Files fetched: `SSSP_1.3.0_PBE_precision.tar.gz` (62,963,841 bytes, md5 fde94756886f32ada7bf597547557eb5), `SSSP_1.3.0_PBE_precision.json` (22,485 bytes, md5 1692c5c9ce89e1c7c783f8f0eee0cbfa), `LICENSE.txt` (3,673 bytes, md5 c0f9df05fb4cd4921f4caa25606815ef); MD5s as published by the record API. The record is CC BY 4.0 except the pseudopotentials, which keep their family licenses; `LICENSE.txt` is vendored as `environment/lab/SSSP_LICENSE.txt`.
- Both UPF files are copied unchanged, with their embedded generation information, into `environment/lab/` (the lab image) only. The agent image contains neither Quantum ESPRESSO nor the pseudopotentials.
- **Refetch:** `benchmarks/construction/fresh-ge-internal-strain/fetch_pseudos.py` downloads the files, checks every MD5/SHA256 above and installs them.

## Calculations

- `environment/lab/task_lab.py` writes every `pw.x` input from the agent's structured request (kind, pseudopotential, lattice parameter, strain tensor, sublattice displacement, cutoff, k mesh and shift, occupations, smearing width, relaxation threshold; ranges in `public.json`), runs `mpirun -np 2 pw.x -npool 2` with one thread per rank, and parses the output (energy, stress, forces, band information, relaxation steps). Fixed choices: PBE from the pseudopotential, ecutrho = 8 ecutwfc, conv_thr 1e-10 Ry, mixing_beta 0.5, four empty bands, Marzari–Vanderbilt smearing when smearing is requested, BFGS with etot_conv_thr 1e-5 Ry for relaxations, Ge mass 72.63 amu. `pw.in` and `pw.out` of every job are kept as evidence; wavefunctions are deleted.
- The charge of a request is computed before it runs from its canonical size: the number of symmetry-distinct k points pw.x will use (the lab reproduces pw.x's `kpoint_grid` + `irreducible_BZ` count from the strain and displacement symmetry; checked against pw.x on every construction run), times a pseudopotential/cutoff factor and ×3 for relaxations.
- **Hidden references** were produced with the same `task_lab.run` code path, in the runtime image `ai-theorist-benchmark-runtime:20261007` on 2 cores (`benchmarks/construction/fresh-ge-internal-strain/qe_batch.py`, requests in `reference_plan.py`, raw results in `reference_runs.jsonl`). `make_task.py` applies the estimators of `solution/solve.py` to them; the derived graded values are stored only in `environment/lab/private/`.

## What is hidden

The graded values (see `design.json` → `calibration`): the zero-pressure lattice parameter, C11, C12, the clamped- and relaxed-ion C44, the Kleinman parameter, the zone-centre optical frequency at a0 and its Grüneisen parameter, all for PBE with the 3d-in-valence pseudopotential, converged with respect to the lab's numerical parameters. None of them is returned by the lab as such: each is a fit, a strain or displacement derivative, or a combination of several calculations.

## Scope

New static-lattice PBE calculations (no spin–orbit coupling, no zero-point or thermal effects) with one pseudopotential and one code version: protocol-specific values, not experimental or literature numbers.
