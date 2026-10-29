# Source and transformations

This is a **fresh-start** task: every number the agent sees comes from a new Quantum ESPRESSO 6.7 `pw.x` calculation that the lab runs on request. No archived calculation result is served or used. The only external inputs are six pseudopotential files.

The physics question is paired, at the level of the scientific question only, with the released GeTe study of Dangic, Fahy and Savic: *Gaussian approximation potentials (GAP) for germanium telluride*, Materials Cloud Archive record x2tdq-3bt13 (DOI 10.24435/materialscloud:pd-g9, CC BY 4.0), supplement to npj Comput. Mater. (2021), doi:10.1038/s41524-021-00523-7, arXiv:2012.08414 (ferroelectric transition and thermal conductivity of GeTe). No file of that record is used: the task grades protocol-specific PBE values computed here.

## Pseudopotentials

| file | set | source | valence | md5 | sha256 | license |
|---|---|---|---|---|---|---|
| `ge_pbe_v1.4.uspp.F.UPF` | `sssp` | GBRV 1.4 (Garrity, Bennett, Rabe, Vanderbilt, Comput. Mater. Sci. 81, 446 (2014)); SSSP 1.3.0 PBE precision choice for Ge (45/360 Ry) | 3d10 4s2 4p2 (Z = 14) | 9c9eaa91e581c3f09632fb3098b2c6b2 | 1b95046ab0f3228e0404719ca66d0309af5a62313be0b22171b963566c90430b | GNU GPL v3 |
| `Te.us.z_6.ld1.psl.v1.0.0-low.upf` | `sssp` | PSlibrary 1.0.0 Te ultrasoft regenerated with ld1.x by the SSSP authors; SSSP 1.3.0 PBE precision choice for Te (50/400 Ry) | 5s2 5p4 (Z = 6) | 44eeed7063cd65092dfe4811eda022ea | 66e19d27182ab8f64e558f3d998fa3da3ca161ff830dcfeeac2cbe9cbd6b2f7c | GNU GPL v2 or later (PSlibrary) |
| `Ge.pbe-n-rrkjus_psl.1.0.0.UPF` | `pslibrary` | PSlibrary 1.0.0 (A. Dal Corso, Comput. Mater. Sci. 95, 337 (2014)), QE pseudopotential library | 4s2 4p2, 3d frozen, NLCC (Z = 4) | a23a7f3f7456476e534e54484b27d8e8 | b756d20fe98992cb3afef5fcde78515720a89336d412c13b0ceddf794d54cfcd | GNU GPL v2 or later |
| `Te.pbe-n-rrkjus_psl.1.0.0.UPF` | `pslibrary` | PSlibrary 1.0.0, QE pseudopotential library | 5s2 5p4 (Z = 6) | 203a1e86ef55f0cd9de2762a27c58f54 | 1a21c8b1d595a3f0165eb1a6ad3711d78fa512d8f583215e898cd5a1716791ba | GNU GPL v2 or later |
| `Ge.rel-pbe-n-rrkjus_psl.1.0.0.UPF` | `pslibrary_fullrel` | PSlibrary 1.0.0 fully relativistic counterpart | as above | e5dcc699b05aed9f3c83ff0819d69beb | 5bcfc408bbf02c1a614b2805239f9bbcd99c42e27f85feaa7c072c87a8136b38 | GNU GPL v2 or later |
| `Te.rel-pbe-n-rrkjus_psl.1.0.0.UPF` | `pslibrary_fullrel` | PSlibrary 1.0.0 fully relativistic counterpart | as above | 7fd30b6833abd190471e4dcc2367351b | 937b711a653d4b8faffe7e46326a98c1e4ad6ccb17fb4e37e1288b39f6cd4a37 | GNU GPL v2 or later |

- **SSSP archive:** G. Prandini, A. Marrazzo, I. E. Castelli, N. Mounet, N. Marzari, *Precision and efficiency in solid-state pseudopotential calculations*, npj Comput. Mater. 4, 72 (2018); Materials Cloud Archive record rcyfm-68h65, DOI 10.24435/materialscloud:f3-ym. Files fetched: `SSSP_1.3.0_PBE_precision.tar.gz` (62,963,841 bytes, md5 fde94756886f32ada7bf597547557eb5), `SSSP_1.3.0_PBE_precision.json` (md5 1692c5c9ce89e1c7c783f8f0eee0cbfa), `LICENSE.txt` (md5 c0f9df05fb4cd4921f4caa25606815ef, vendored as `environment/lab/SSSP_LICENSE.txt`). The record is CC BY 4.0 except the pseudopotentials, which keep their family licenses. The SSSP UPF MD5s match the SSSP metadata.
- **PSlibrary files:** https://pseudopotentials.quantum-espresso.org/upf_files/<name>. The scalar and fully relativistic files of each element were generated with the same parameters (same author, date and cutoffs in their headers). The SSSP Te file has the same total pseudo-energy as `Te.pbe-n-rrkjus_psl.1.0.0.UPF` (it is a regeneration of it), so the `sssp` and `pslibrary` sets differ essentially through Ge (3d in valence versus frozen).
- All UPF files are copied unchanged into `environment/lab/` (the lab image) only. The agent image contains neither Quantum ESPRESSO nor the pseudopotentials.
- **Refetch:** `benchmarks/construction/fresh-gete-ferroelectric-well/fetch_pseudos.py` downloads the files, checks every checksum (`pseudo_checksums.json`) and installs them.

## Calculations

- `environment/lab/task_lab.py` writes every `pw.x` input from the agent's structured request (pseudopotential set, a_rh, α, x, cutoff, k mesh and shift, occupations, smearing width; ranges in `public.json`), runs `mpirun -np 2 pw.x -npool 2` with one thread per rank and parses the output. Fixed choices: ibrav = 0 with the threefold axis along Cartesian [111] (so rocksalt keeps all 48 operations), ecutrho = 8 ecutwfc, conv_thr 1e-10 Ry, mixing_beta 0.4, six empty bands (twelve with spin–orbit), Marzari–Vanderbilt smearing when requested, masses 72.63 and 127.6 amu. `pw.in` and `pw.out` of every job are kept; the pw.x scratch is deleted, also for failed jobs.
- The charge is computed before the run from the canonical request size: the number of symmetry-distinct k points pw.x will use (reproduced from pw.x's `kpoint_grid` + `irreducible_BZ` logic and checked against pw.x on every construction run), times a pseudopotential-set factor and a cutoff factor.
- **Hidden references** were produced with the same `task_lab.run` code path in the runtime image `ai-theorist-benchmark-runtime:20261007` on 2 cores (`qe_batch.py`, requests in `reference_plan.py`, raw results in `reference_runs.jsonl`). `make_task.py` applies the estimators of `solution/solve.py` to them; the derived graded values are stored only in `environment/lab/private/`.

## What is hidden

The graded values (see `design.json` → `calibration`): the zero-pressure rocksalt lattice parameter, the relaxed R3m structure (a_rh, α, x), the relaxed and clamped well depths, the soft-mode critical pressure and the spin–orbit change of the well depth. None of them is returned by the lab: each is a fit, a stationary point, a sign change or a difference of several calculations.

## Scope

New static-lattice PBE calculations with one code version: protocol-specific values, not experimental or literature numbers.
