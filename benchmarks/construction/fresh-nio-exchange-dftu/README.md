# Construction of fresh-nio-exchange-dftu

Fresh-start the benchmark task (see `benchmarks/tools/BUILDING.md`): the lab runs new Quantum ESPRESSO 6.7 PBE+U
`pw.x` SCF calculations of rocksalt NiO requested by the agent; the hidden truth comes from reference calculations
made through the same `task_lab.run` code path.

## Rebuild

All QE work runs in `ai-theorist-benchmark-runtime:20261007` with `docker run --rm --cpuset-cpus 8,9,24,25 --user UID:GID -e OMPI_MCA_btl_vader_single_copy_mechanism=none` and `mpirun -np 2` (done by `qe_batch.py`).

1. `python3 -B fetch_pseudos.py --install` downloads SSSP 1.3.0 PBE precision (record rcyfm-68h65) into
   `~/.cache/ai-theorist/fresh-nio-exchange-dftu/`, checks the record MD5s, extracts `ni_pbe_v1.4.uspp.F.UPF`
   (GBRV, GPL v3) and `O.pbe-n-kjpaw_psl.0.1.UPF` (PSlibrary 0.3.1, GPL v2+), checks them against the SSSP metadata
   MD5 and the SHA256 in `pseudo_checksums.json`, and copies them plus `SSSP_LICENSE.txt` into `environment/lab/`.
2. `python3 -B run_reference.py` runs the 129 requests of `reference_plan.py` through `task_lab.run` and appends
   them to `reference_runs.jsonl` (resumable). Total ≈ 2.8 h on 2 cores. One request fails by design
   (`U0.0_rh_fm_k444`: the metallic PBE FM state does not converge in 200 iterations on that mesh).
3. `python3 -B make_task.py` applies the estimators of `solution/solve.py` to the reference runs and writes
   `private/hidden.json`, `private/hidden_values.json`, `solution/ablation_values.json`, `calibration.json` and
   `design.json` → `calibration`.
4. `.venv/bin/python -B replay.py` (from the repository root) replays `solve.py` against `reference_runs.jsonl`
   with a stand-in lab client (seconds, no QE) and grades the ablations and robustness variants as check_task does.
5. `.venv/bin/python -m benchmarks.tools.check_task fresh-nio-exchange-dftu --cpuset 8,9,24,25` writes `check.json`.

## Lab

- **Protocol (fixed):** a = 4.17 Å; PBE; `lda_plus_u_kind = 0` (Dudarev) with `U_projection_type = 'atomic'` on every
  Ni species (QE 6.7 namelist input; the HUBBARD card arrived in QE 7.1); ecutrho = 8 ecutwfc; nbnd = 14/f.u. + 4;
  conv_thr 1e-9 Ry/f.u.; mixing β 0.3; ≤ 200 SCF iterations; one Ni species per distinct starting magnetization.
- **Request (no defaults):** cell (`primitive_1fu`, `rhombohedral_2fu` = AFM-II cell, `tetragonal_2fu` = AFM-I cell,
  `cubic_4fu`), starting magnetization per Ni (or null), U (0–9 eV, 0.1 steps), MP mesh (per-cell maxima 16/7/10/5),
  shift, ecutwfc (40/50/60/75), smearing (gaussian/mv/fd) and degauss (0.002–0.05 Ry).
- **Returned:** F, −TS and E per cell (Ry), Fermi energy, gap on the mesh, Ni/O moments, total/absolute
  magnetization, Ni 3d occupation traces and eigenvalues per spin, iterations, irreducible k, wall time; pw.in and
  pw.out stay in the job folder.
- **Charge:** ceil(w × n1n2n3/1000 × f), w = 0.4/10/2.5/20 (primitive/rhombohedral/tetragonal/cubic),
  f = 0.72/1/1.31/1.84 for 40/50/60/75 Ry; about one unit per minute on 2 cores (the rhombohedral cell has a
  large, low-symmetry FFT grid and is ~3× slower per k-point than the tetragonal cell). Budget 36.
- **Measured wall times (2 cores):** rhombohedral 4³ 70–130 s, 6³ 120–210 s, 7³ 170–250 s at 50 Ry; tetragonal
  8×8×6 55–65 s, 10×10×8 90–105 s; primitive 12³ 25–40 s, 16³ 50–65 s. Slowest allowed requests measured: rhombohedral 7³ at 75 Ry for the
  metallic PBE FM state 403 s (6.7 min); cubic 5³ at 75 Ry 124 s. The pw.x timeout is 1500 s. np = 4 on the hyperthread siblings was only 9% faster than np = 2, so the lab uses 2 ranks.

## Physics found in construction (50 Ry, Gaussian 0.005 Ry unless stated)

| U (eV) | E_FM − E_AFM-II (meV/Ni) | E_AFM-I − E_FM (meV/Ni) | J1 (meV) | J2 (meV) | FM gap / AFM-II gap (eV) | Ni moment FM / AFM-II |
|---|---|---|---|---|---|---|
| 0 | 264.0 (FM metal, primitive 16³) | 23 (k886, both metallic) | — | — | 0 / 0.89 | 1.24 / 1.33 |
| 2 | 232 | 26.7 | — | — | 0 / 1.71 | 1.63 / 1.49 |
| 4 | 169.98 | 17.48 | −2.186 | 30.52 | 0.77 / 2.47 | 1.678 / 1.579 |
| 6 | 121.5 | 11.0 | −1.37 | 21.6 | 1.58 / 3.0 | 1.722 / 1.651 |
| 8 | 83.74 | 6.34 | −0.792 | 14.75 | 2.18 / 3.12 | 1.763 / 1.712 |

- J2 dominates and roughly halves from U = 4 to 8 eV (×0.48), close to 1/U; J1 is weak, ferromagnetic and falls
  faster (×0.36). At U = 6 eV the protocol gives J1 = −1.37 meV, within 10% of the experimental −1.4 meV, which is
  why U = 6 eV is not graded (recall risk); J2 is never near the experimental 19 meV at U = 4 or 8 eV.
- **Cross-cell artefact (dp2).** The FM energy per f.u. converges at different rates in different cells: at 4³
  the rhombohedral FM energy is 4.95 (U = 4) / 6.2 meV (U = 8) above the converged value, whereas the same-cell
  difference FM − AFM-II is within 1.2% at 4³. Pairing that FM with a converged tetragonal AFM-I gives J1 −1.57
  instead of −2.19 (U = 4) and −0.01 instead of −0.79 meV (U = 8).
- **k convergence.** FM − AFM-II: 3³ +10% (U = 4), +12% (U = 8); 4³ −1.2/+0.9%; 6³ within 0.1%. AFM-I − FM:
  4×4×4 −10/−7%, 6×6×4 −5.3/−3.5%, 6×6×6 −1.1%, 8×8×6 −0.5%, 8×8×8 and 10×10×8 converged.
- **Cutoff** 40 → 50 → 75 Ry changes FM − AFM-II by 0.1–0.2% and AFM-I − FM by ≤ 1% (J1 at U = 8: 0.015 meV):
  not a decision point, only a budget trade-off.
- **Smearing at U = 4 eV (dp7).** Fermi–Dirac 0.02 Ry and Gaussian 0.05 Ry close the 0.77 eV FM gap and lower the
  FM moment by 0.03 μB; J1 +8–10%, J2 +5–6% even with the (E+F)/2 zero-smearing estimate. MV 0.02 keeps the gap
  but shifts F by up to 0.6 mRy per f.u. (J2 −3%). Fermi–Dirac 0.05 Ry is far worse (J2 −17%).
- **Plain PBE (dp8–dp10).** AFM-II is an insulator (gap 0.8–0.9 eV, needs 6³; 4³ is 21 meV/f.u. too high). FM is a
  metal: its moment drifts with mesh and smearing (1.32 at 8³, 1.26 at 12³, 1.24 at 16³), the −TS term is several
  meV/f.u., and Gaussian 0.02 Ry F is 29 meV/f.u. below the zero-smearing estimate. A starting magnetization of 0.1
  lands in a low-moment FM solution (0.69 μB, +46 meV/Ni) in the rhombohedral cell and fails to converge in the
  primitive cell; 0.5 and 1.0 give the same high-moment state. In the rhombohedral 4³ cell the FM SCF does not
  converge at all. At U ≥ 4 eV every start (0.1, 0.5, 1.0; one-site-only start) gives the same high-spin
  t2g⁶e_g² state (occupation eigenvalues ≈ 1 ×5 majority, 3 × ~1 + 2 × ~0.2 minority), so no metastable DFT+U
  occupation state was found there; the design's metastability fork moved to U = 0.
- **FM gap closing (dp11).** FM gaps on a 12³ primitive mesh: 0 at U ≤ 2.0, 0.005 (2.1), 0.047 (2.2), 0.088 (2.3),
  0.170 (2.5), 0.37 (3.0), 0.57 (3.5), 0.771 (4.0), 2.176 (8.0) eV; 16³ agrees within 0.002 eV. Reference
  U_c = 2.087 eV (line through the first two insulating points, bracketed by the last metallic one).
  Shifted meshes (kshift = 1) miss the band extrema: shifted 12³ gives 0.077 (U = 2.0), 0.158 (2.2; 0.110 at
  shifted 16³) and 0.276 eV (2.5), hence U_c = 1.81 eV (dp12).
- **Non-spin-polarized NiO** at U = 6 eV is 1.8 eV/f.u. above AFM-II (metallic) — the lab allows it, nothing uses it.
- The cubic 4-f.u. cell gives FM, AFM-I and the 3:1 ferrimagnet (E0 + 3J2 per Ni), an alternative J1 route
  (exploration at U = 6 eV, 4³: AFM-I − FM = 10.1 meV/Ni against 11.0 converged; its absolute energies per f.u.
  at 4³ are 16 meV above the converged ones, another instance of the cross-cell trap).

## Deviations from the candidate design

- **Graded U = 4 and 8 eV instead of 6 eV**, because J1(U = 6) of this protocol (−1.37 meV) coincides with the
  famous experimental value (recall risk).
- **No hidden SCF at a withheld U (4.3 eV).** U is purchasable in 0.1 eV steps and interpolation between
  purchasable U values is accurate to < 1%, so a withheld-U target would not test the J(U) law. The U dependence is
  graded through J at two U values plus the FM gap-closing U (the regime boundary below which the FM reference is
  metallic).
- **Fixed lattice parameter** (4.17 Å) instead of an agent-chosen one: J depends steeply on volume and an
  agent-chosen a would make the targets non-unique.
- **"heisenberg_valid_at_U0: false" flag replaced** by two quantitative U = 0 targets (converged
  E_FM − E_AFM-II with the metallic FM, and U_c); the rigidity argument is a rubric item.
- **Metastability fork moved to U = 0** (no metastable occupation state exists at U ≥ 4 eV for these starts).
- **Cutoff fork dropped** as a graded decision point (≤ 1% effect); AFM-III not offered (not needed).
- **T_N formula:** k T_N = S(S+1)/3 × 6 J2 for AFM-II (the design's (2/3) S(S+1) × 6 J2 counts pairs twice).
- **θ_CW sign not graded** (trivially negative); it appears only as the Curie–Weiss ablation of T_N.

## Calibration (full table: `calibration.json`, `design.json` → `calibration`)

| field | reference | tolerance | largest robustness error | nearest failing ablation |
|---|---|---|---|---|
| J1 U = 4 / 8 (meV) | −2.1858 / −0.7919 | 3% | 0.6% / 2.0% (40 Ry) | 6×6×4 tetragonal −5.3% / −3.5%; cross-cell −28% / −98% |
| J2 U = 4 / 8 (meV) | 30.516 / 14.750 | 2% | 1.4% / 0.9% (4³ rhombohedral at 75 Ry) | cross-cell −3.2% / −4.5%; FM+AFM-II only −7.2% / −5.4%; FD 0.02 +5.2%; 3³ +9% / +12% |
| T_N U = 4 / 8 (K) | 1416.5 / 684.66 | 2% | as J2 | Curie–Weiss −14% / −11%; S² −50%; plus all J2 paths |
| PBE E_FM − E_AFM-II (meV/Ni) | 264.05 | 2% | 0.5% (12³, Gaussian 0.002 / MV 0.005) | 8³ primitive FM +3.9%; coarse meshes −4.0%; FD 0.02 +9.3%; Gaussian 0.02 F −11%; small start +16% |
| FM gap-closing U (eV) | 2.087 | ±0.1 | 0.011 | shifted 12³ meshes −0.28; linear from U = 4, 8: −0.28; 2/4 eV bracket +0.91 |

Robustness variants: oracle meshes (rhombohedral 6³, tetragonal 8×8×6), rhombohedral 4³ + tetragonal 8×8×8,
the converged primitive-cell FM as a common reference, 40 and 75 Ry cutoffs, PBE FM with Gaussian 0.002 (12³, 16³)
and MV 0.005, PBE AFM-II at 6³, and three gap-closing estimators. Every ablation fails at least one check; the
tightest margins are J1 at U = 8 for the 6×6×4 ablation (3.5% against 3%; the same ablation fails J1 at U = 4 by
5.3%) and J2 at U = 4 for the 4³-at-75-Ry robustness variant (1.4% against 2%).

## Oracle

`solution/solve.py`: rhombohedral FM/AFM-II at 4³ and 6³ and tetragonal FM/AFM-I at 8×8×6 and 10×10×8 (U = 4);
rhombohedral 6³ and tetragonal 8×8×6 pairs at U = 8; PBE AFM-II at 6³, PBE FM in the primitive cell at 12³ and 16³
plus a starting-moment check; FM gaps at U = 2.0, 2.2, 2.5 (primitive 12³). 19 evaluations, 32 of 36 units (89%),
≈ 27 min of SCF time. The offline replay reproduces the references to ≤ 0.2% and 0.003 eV.

## Gate and Harbor

- **Gate** (`check_task --cpuset 8,9,24,25`, final full run 2026-10-25T22:02Z on the finished package): passed, no
  issues. 12 decision points, 18 ablations all failing, 12 robustness variants and 4-significant-digit rounding
  passing, grading deterministic, 8 hidden values, task 1.55 MB (the two UPF files). Live oracle: 19 evaluations,
  32 of 36 units, 0 failed, 1513 s container time; its errors against the references are ≤ 0.2% and 0.0024 eV
  (identical requests reproduce the construction energies exactly).
- **Harbor oracle** (`runs/build-fresh-nio-exchange-dftu/oracle-r1-fresh-nio-exchange-dftu`, 21:36–22:19Z):
  reward 1.0, all 9 numerical checks true, charge 32; the coordinator's daemon pinned the lab to cores 4,20,5,21;
  43 min wall. (An earlier Harbor run, 17:11–17:55Z, also scored 1.0; it was repeated after the last package
  changes: the dp12 ablation, rubric wording, and `task_lab.run` now deleting the pw.x scratch of failed jobs too.)
- **Disk:** peak ≈ 160 MB (62 MB SSSP archive in `~/.cache/ai-theorist/fresh-nio-exchange-dftu/`, deleted after
  the gate; ≈ 95 MB of exploration and reference job folders under `runs/`, of which 6 MB of reference pw.in/pw.out
  are kept in `runs/build-fresh-nio-exchange-dftu/reference_jobs/`).

## Open risks

- **Recall.** NiO exchange is famous (experiment J1 ≈ −1.4, J2 ≈ 19 meV; T_N = 523 K). The graded values are at
  U = 4 and 8 eV, where this protocol gives J2 = 30.5 and 14.8 meV and J1 = −2.19 and −0.79 meV, so a recalled value
  fails every J check.
- **Dev-easy risk.** The energy mapping, pair counting and S(S+1) are textbook (current models handle them). The
  discriminating forks are the discovered ones: the cross-cell FM reference (dp2), the slow tetragonal k
  convergence at U = 4 eV (dp4), metal-style smearing closing the 0.77 eV FM gap (dp7), the metallic plain-PBE FM
  (dp8–dp10: zero-smearing limit, dense meshes, metastable low-moment start) and shifted meshes hiding the FM gap
  closing (dp12). A careful agent that converges every absolute energy per formula unit passes dp2 legitimately.
- **Tight margins.** J1 at U = 8 eV for the 6×6×4 ablation (3.5% against 3%; it also fails U = 4 by 5.3%); J2 at
  U = 4 eV for the 4³-rhombohedral robustness variant (1.2–1.4% against 2%). Tolerances of 2–3% are far above the
  run-to-run noise (identical requests reproduce energies to 1e-8 Ry) and the cutoff spread.
- **Timing.** The slowest allowed request (rhombohedral 7³ at 75 Ry, metallic PBE FM) took 403 s; metallic U = 0
  SCFs need up to 33 iterations, and some do not converge (rhombohedral 4³ FM, primitive FM with a 0.1 start), which
  is part of the task (failed jobs stay charged).
- The gap-closing U is defined on the lab's meshes as the U where the FM gap of the Γ-centred-converged band
  structure vanishes; 12³ and 16³ agree within 0.004 eV, so the definition is unique to well inside the tolerance.
