# Source and transformations

This is a **fresh-start** task: every number the agent sees comes from a new Quantum ESPRESSO 6.7 calculation (`pw.x`, and `bands.x` for symmetry labels) that the lab runs on request. No archived calculation result is served or used. The only external inputs are two pseudopotential files.

The physics question is paired, at the level of the scientific question only, with the released study M. Bonacci, M. Schüler, G. Colonna, N. Marzari, *Ab initio study of strain and quantum confinement shaping the valence-band structure of Ge quantum wells*, arXiv:2603.18753 (2026); Materials Cloud Archive record vz6hg-tvd13, DOI 10.24435/materialscloud:yg-hd, https://archive.materialscloud.org/records/vz6hg-tvd13 (CC BY 4.0). That study uses HSE with spin–orbit coupling at the experimental lattice parameter. No file of the record is used or staged: the task grades protocol-specific PBE values computed here. (The record ships the same fully relativistic Ge potential; its copy differs from the SG15 file only by one header line holding the file name.)

## Pseudopotentials

| file | source | used for | size (bytes) | md5 | sha256 | license |
|---|---|---|---|---|---|---|
| `Ge_ONCV_PBE_FR-1.0.upf` | SG15 ONCV, fully relativistic: P. Scherpelz, M. Govoni, I. Hamada, G. Galli, J. Chem. Theory Comput. 12, 3523 (2016), built with the SG15 parameters of M. Schlipf, F. Gygi, Comput. Phys. Commun. 196, 36 (2015) and ONCVPSP (D. R. Hamann, PRB 88, 085117 (2013)); http://www.quantum-simulation.org/potentials/sg15_oncv/upf/Ge_ONCV_PBE_FR-1.0.upf | `spin_orbit: true` (the graded protocol) | 167,286 | 046fb2884f640e00dddb7bc679dd92bf | ee520682a4103e45bff9360d87b9fb8f15a2045eb206da5ef412372984fd53af | **CC BY-SA 4.0** (stated in the file's PP_INFO; copyright 2015 The Regents of the University of California, modifications 2016 P. Scherpelz) |
| `Ge_ONCV_PBE-1.0.upf` | SG15 ONCV, scalar-relativistic (Schlipf and Gygi 2015); http://www.quantum-simulation.org/potentials/sg15_oncv/upf/Ge_ONCV_PBE-1.0.upf | `spin_orbit: false` | 115,422 | 39f69ca9d3fedd83738352a2cbce5cbd | 1023970eba91ac075d8a5d85a4817afe98f98affcde8a0dd378a4365105548d7 | **CC BY-SA 4.0** (stated in the file's PP_INFO; copyright 2015 The Regents of the University of California) |

- Both files are norm-conserving with 3d10 4s2 4p2 in valence (Z = 14) and are copied unchanged, with their license notices, into `environment/lab/` (the lab image) only. The agent image contains neither Quantum ESPRESSO nor the pseudopotentials. As CC BY-SA files they are redistributed under the same license; nothing derived from them is distributed beyond calculation results.
- **Why two files.** Quantum ESPRESSO 6.7 does not j-average the fully relativistic file when spin–orbit coupling is switched off: in construction tests that mode gave a pressure of −16.5 GPa near a0 and placed the s-like level 1.3 eV above the p triplet (both with `noncolin = .false.` and with `noncolin = .true., lspinorb = .false.`). The scalar-relativistic runs therefore use the scalar SG15 file of the same family (its p-triplet level agrees with the j-weighted average of the spin–orbit Γ8+/Γ7+ levels to 2 meV).
- **Refetch:** `benchmarks/construction/fresh-ge-soc-strained-valence/fetch_pseudo.py` downloads both files and checks size, MD5, SHA256 and the license line.

## Calculations

- `environment/lab/task_lab.py` writes every input from the agent's structured request (a, c, spin–orbit switch, cutoff, k mesh and shift, occupations, smearing, band k points, symmetry labels; ranges in `public.json`), runs `mpirun -np 2 pw.x -npool 2` with one thread per rank (SCF, then optionally a `bands` run on the converged density), runs `bands.x` with `lsym = .true.` when labels are requested, and parses energy, stress, Fermi level or band edges, the band energies at the mesh point nearest Γ and at the requested k points, and the labels. Fixed choices: PBE, ecutrho = 4 ecutwfc, 36 spinor (18 scalar) bands, conv_thr 1e-10 Ry, mixing_beta 0.5, Marzari–Vanderbilt smearing when smearing is requested. Inputs and outputs of every job are kept as evidence; wavefunctions are deleted.
- The charge of a request is computed before it runs from its canonical size: the number of symmetry-distinct k points pw.x will use for the cubic or tetragonal cell (reproduced by the lab and checked against pw.x on every construction run), times a factor for spin–orbit coupling and cutoff, plus the band k points.
- **Hidden references** were produced with the same `task_lab.run` code path in the runtime image `ai-theorist-benchmark-runtime:20261007` on 2 cores (`benchmarks/construction/fresh-ge-soc-strained-valence/qe_batch.py`; requests in `reference_plan.py`; raw results in `reference_runs.jsonl`). `make_task.py` applies the estimators of `solution/solve.py` to them; the derived graded values are stored only in `environment/lab/private/`.

## What is hidden

The graded values (see `design.json` → `calibration`): a0, Δ0 and the s-like level at a0, c/a, E(HH) − E(LH) and E(SO) − E(HH) at ε∥ = −1 % relaxed, b in the small-strain limit, and the s/split-off crossing strain. None of them is returned by the lab as such: each is a pressure root, a level difference that needs identification by character, a stress interpolation, a strain derivative or an interpolated crossing.

## Scope

New static-lattice PBE calculations with spin–orbit coupling, one pseudopotential family and one code version: protocol-specific values, not experimental or literature numbers.
