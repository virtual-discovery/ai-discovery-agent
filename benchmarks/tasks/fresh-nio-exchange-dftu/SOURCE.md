# Source and transformations

This is a **fresh-start** task: every number the agent sees comes from a new Quantum ESPRESSO 6.7 (pw.x) calculation that the lab runs on request. No archived energies are served.

## Pseudopotentials (the only external input)

- **Library:** SSSP 1.3.0 PBE precision. G. Prandini, A. Marrazzo, I. E. Castelli, N. Mounet, N. Marzari, *Precision and efficiency in solid-state pseudopotential calculations*, npj Comput. Mater. 4, 72 (2018); Materials Cloud Archive record rcyfm-68h65, DOI 10.24435/materialscloud:f3-ym, https://archive.materialscloud.org/records/rcyfm-68h65.
- **Archive files fetched:** `SSSP_1.3.0_PBE_precision.tar.gz` (62,963,841 bytes, md5 fde94756886f32ada7bf597547557eb5), `SSSP_1.3.0_PBE_precision.json` (22,485 bytes, md5 1692c5c9ce89e1c7c783f8f0eee0cbfa), `LICENSE.txt` (3,673 bytes, md5 c0f9df05fb4cd4921f4caa25606815ef). The MD5s are those published by the record API.
- **Files vendored into `environment/lab/` (lab image only; the agent image never contains them), copied unchanged:**
  - `ni_pbe_v1.4.uspp.F.UPF`: GBRV 1.4 ultrasoft (K. F. Garrity, J. W. Bennett, K. M. Rabe, D. Vanderbilt, Comput. Mater. Sci. 81, 446 (2014)); 612,014 bytes, md5 1ee80287db30b12d2bc1f57a5b5d6409 (matches the SSSP metadata), sha256 f76b86ce60cde3d83dfcc8df79ba05478db573d158289f1b226919442d977d25. **License: GNU GPL v3** (per the archive's `LICENSE.txt`). SSSP cutoffs 50/400 Ry.
  - `O.pbe-n-kjpaw_psl.0.1.UPF`: PSlibrary 0.3.1 PAW (A. Dal Corso, Comput. Mater. Sci. 95, 337 (2014)); 930,566 bytes, md5 0234752ac141de4415c5fc33072bef88 (matches the SSSP metadata), sha256 7c4b6ed541f83d0afdf5c1d3a8c611340f073f3e2b11e95b7963bb2ee26929aa. **License: GNU GPL v2 or later.** SSSP cutoffs 75/600 Ry.
  - `SSSP_LICENSE.txt`: the archive's `LICENSE.txt` (per-family licenses; the rest of the record is CC BY 4.0).
- **Refetch:** `benchmarks/construction/fresh-nio-exchange-dftu/fetch_pseudos.py --install` downloads the three archive files, checks the MD5s, extracts the two UPF files, checks them against the SSSP metadata MD5 and the SHA256 pinned in `pseudo_checksums.json`, and installs them.

## Calculations

- `environment/lab/task_lab.py` writes every pw.x input from the agent's structured request (cell, starting moment per Ni, U, k mesh and shift, cutoff, smearing; schema in `public.json`, no defaults), runs `mpirun -np 2 pw.x -npool 2` with one thread per rank and parses pw.out and the XML. Fixed choices: rocksalt NiO at a = 4.17 Å; PBE; `lda_plus_u_kind = 0` (Dudarev) with `U_projection_type = 'atomic'` on every Ni species (QE 6.7 still uses the `Hubbard_U(i)` namelist input, not the HUBBARD card of QE ≥ 7.1); one Ni species per distinct starting magnetization; ecutrho = 8 ecutwfc; nbnd = 14 per formula unit + 4; conv_thr = 1e-9 Ry per formula unit; mixing_beta = 0.3; at most 200 SCF iterations.
- **Hidden references** were produced with the same `task_lab.run` code path, in the runtime image `ai-theorist-benchmark-runtime:20261007` on 2 cores (cpuset 8,9,24,25), by `benchmarks/construction/fresh-nio-exchange-dftu/run_reference.py` (requests listed in `reference_plan.py`, runner `qe_batch.py`); the raw results are kept in that folder (`reference_runs.jsonl`). `make_task.py` applies the estimators of `solution/solve.py` to them, and the derived graded values are stored only in `environment/lab/private/`.

## What is hidden

J1 and J2 at U = 4 and 8 eV (converged same-cell energy differences), the molecular-field Néel temperatures derived from them, E(FM) − E(AFM-II) of plain PBE in the converged zero-smearing limit, and the U at which the ferromagnetic gap closes. None is returned by the lab as such: each is a combination of several calculations (energy differences of converged k meshes in two different cells, a smearing-corrected metallic energy, an interpolated gap-closing point).

## Context (not used as data)

The physics question mirrors G. dos Santos, L. Binci, G. Menichetti, R. Mahajan, N. Marzari, I. Timrov, *Comparative study of magnetic exchange parameters and magnon dispersions in NiO and MnO from first principles*, Phys. Rev. B (2025), doi:10.1103/gtxm-6vtg, arXiv:2508.12153 (Materials Cloud record kg706-sbx37, DOI 10.24435/materialscloud:jw-h1, CC BY 4.0). That study uses PseudoDojo pseudopotentials, linear-response U and other protocols; nothing from its archive is served or graded. The graded values are protocol-specific to this lab and do not reproduce any published number.
