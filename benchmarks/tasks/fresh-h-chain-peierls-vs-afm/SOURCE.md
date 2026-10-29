# Source and transformations

This is a **fresh-start** task: every number the agent sees comes from a new Quantum ESPRESSO 6.7 (pw.x) calculation that the lab runs on request. No archived energies are served.

## Pseudopotential (the only external input)

- **Library:** SSSP 1.3.0 PBE precision. G. Prandini, A. Marrazzo, I. E. Castelli, N. Mounet, N. Marzari, *Precision and efficiency in solid-state pseudopotential calculations*, npj Comput. Mater. 4, 72 (2018); Materials Cloud Archive record rcyfm-68h65, DOI 10.24435/materialscloud:f3-ym, https://archive.materialscloud.org/records/rcyfm-68h65.
- **Archive files fetched:** `SSSP_1.3.0_PBE_precision.tar.gz` (62,963,841 bytes, md5 fde94756886f32ada7bf597547557eb5), `SSSP_1.3.0_PBE_precision.json` (22,485 bytes, md5 1692c5c9ce89e1c7c783f8f0eee0cbfa), `LICENSE.txt` (3,673 bytes, md5 c0f9df05fb4cd4921f4caa25606815ef). The MD5s are those published by the record API.
- **File vendored:** `environment/lab/H_ONCV_PBE-1.0.oncvpsp.upf` (SG15 ONCV norm-conserving, M. Schlipf and F. Gygi, Comput. Phys. Commun. 196, 36 (2015)); md5 1790becc920ee074925cf490c71280fe (matches the SSSP metadata), sha256 69d631780c0178c5239ed3dcb4a3e22c13bd13cde9beb97e123d036949da5c09. SSSP recommended cutoffs for H: 80 Ry (wavefunctions), 320 Ry (density).
- **License of the vendored file:** SG15 pseudopotentials are distributed under **CC BY-SA 4.0** (stated per family in the archive's `LICENSE.txt`); the rest of the SSSP record is CC BY 4.0. The UPF is copied unchanged, with its embedded generation notice, into the lab image only. The agent image does not contain it.
- **Refetch:** `benchmarks/construction/fresh-h-chain-peierls-vs-afm/fetch_pseudo.py` downloads the three archive files, checks the MD5s, extracts the UPF and checks its MD5 and SHA256.

## Calculations

- `environment/lab/task_lab.py` writes every pw.x input from the agent's structured request (bond lengths, starting magnetic moments, k-points, smearing, box, cutoff; ranges in `public.json`), runs `mpirun -np 2 pw.x` and parses the output. Fixed choices: PBE, the SG15 H pseudopotential, ecutrho = 4 ecutwfc, conv_thr = 1e-10 Ry per atom, mixing_beta = 0.3, Γ-centred k-mesh along the chain only.
- **Hidden references** were produced with the same `task_lab.run` code path, in the runtime image `ai-theorist-benchmark-runtime:20261007` on 2 cores, by `benchmarks/construction/fresh-h-chain-peierls-vs-afm/run_reference.py` (requests listed in `reference_plan.py`); the raw reference results are kept in that construction folder (`reference_runs.jsonl`). `make_task.py` applies the estimators of `solution/solve.py` to them, and the derived graded values are stored only in `environment/lab/private/`.

## What is hidden

The graded values (see `design.json` → `calibration`): the equilibrium dimerization and the dimerization energy at R = 1.25 and 1.75 Å, the antiferromagnetic energy gain of the uniform chain at R = 1.75 and 2.0 Å, the critical bond alternation above which the magnetic solution disappears at R = 1.75 and 2.0 Å, the converged band gap of the lowest uniform chain at R = 1.25 Å, and the ground-state labels at R = 1.0 and 2.0 Å. None of them is returned by the lab as such: each is an energy difference, a fitted minimum or boundary, or (the gap) a converged-limit value that needs smearing and k-point convergence tests.

## Scope

PBE mean-field treatment of an idealized infinite chain at fixed mean spacing (collinear spins, clamped nuclei, no zero-point motion). The task does not reproduce any published number; it is related to the hydrogen-chain benchmark literature (Motta et al., Phys. Rev. X 7, 031059 (2017); Phys. Rev. X 10, 031058 (2020)) only through the physical question.
