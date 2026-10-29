# Construction of tmd-moire-chern-twist

Released-start task built from X.-W. Zhang, C. Wang, X. Liu, Y. Fan, T. Cao and D. Xiao, *Polarization-driven band topology evolution in twisted MoTe2 and WSe2*, Nat. Commun. 15, 4223 (2024) (arXiv:2311.12776). It uses Zenodo record 10888769 (CC BY 4.0; one file, `source-data.tar.gz`, 20.8 MB). Provenance and checksums are in the task's `SOURCE.md` and in `members.md5`.

## Rebuild

1. `python3 -B fetch.py [DEST]` downloads `source-data.tar.gz` (md5 038b59d5…) and extracts the 76 members listed in `members.md5` (fig1, fig2, fig3). It never extracts `ML/` (which holds a VASP POTCAR). The default DEST is `~/.cache/ai-theorist/tmd-moire-chern-twist`.
2. `bash run_in_image.sh reference.py hartree` and `bash run_in_image.sh reference.py continuum` write `reference/hartree.json` and `reference/continuum.json` (log: `reference/continuum.log`; about 8 min on 4 CPUs).
3. `bash run_in_image.sh make_task.py` writes the staged archive, the catalog, the lab data, the private targets, the ablation values and `calibration.json` (log: `make_task.log`).
4. `python3 -B write_design.py`, then `.venv/bin/python -m benchmarks.tools.check_task tmd-moire-chern-twist --cpuset 4,5,20,21` (writes `check.json`). Then `python3 -B record_oracle.py` (writes `oracle_run.json`), and rerun `write_design.py` to record the oracle length.

The other files:

- `continuum.py`: the numpy plane-wave continuum model of Eq. (1) and the Fukui–Hatsugai–Suzuki Chern evaluator, with a periodic gauge across the zone boundary.
- `c3_labels.json`: the hand transcription of the Fig. 2 C3z labels.

## What the task is

The question: do the Chern numbers of the top K-valley spin-up moiré valence bands of tWSe2 and tMoTe2 change with twist angle, and is the change explained by a fixed continuum model or by the angle-dependent moiré potential?

The lab serves only archived results:

| Request | Cost | Returns |
|---|---|---|
| `dft_bands` | 1 | SIESTA moiré bands along γ–m–κ–γ (8 calculations) |
| `c3` | 1 | spin and C3z eigenvalue of the six highest states at γ, κ or κ′ (the Fig. 2 labels; K↑ and K′↓ entries in a fixed but arbitrary order within each doublet) |
| `hartree_map` | 3 | tWSe2 top/bottom-surface Hartree maps (3.15/1.70/1.54/1.25°) |
| `hartree_linecut` | 1 | archived `pot_diag-top.dat` |
| `hartree_fourier` | 1 | archived `pot_top_bot_G.dat`: magnitudes only |
| `structure` | 1 | relaxed tWSe2 POSCAR |

The budget is 27 units; buying everything costs 56. The agent computes the continuum model in its own container (numpy). Staged files:

- both fig1 parameter maps, as set P (= fig1b, DFT-fitted) and set Q (= fig1a, local stacking), with the m* and a₀ of their source papers in the catalog;
- the 3.15° tWSe2 and 3.89° tMoTe2 bands;
- the 3.15°/1.25° relaxed structures;
- the 3.15° Hartree maps, line cut and scripts.

Graded values:

- the DFT Chern pairs (C1, C2) of all 8 calculations (exact);
- the continuum Chern pairs of P and Q at 3.89° (exact);
- which set matches the DFT bands;
- set P with fixed parameters at 2.14° (exact);
- Δv_H at XM at the 4 angles (±1.0 meV);
- the first-star part of Δv_H at XM at 1.70° and 1.25° (±max(0.25 meV, 5 %)).

## What the archive shows (and the task is built on)

- **Chern numbers exist only in the figure.** The `c3_labels.json` transcription reproduces the printed C_K↑ of every Fig. 2 panel through exp(2πiC/3) = −ξ_γ ξ_κ ξ_κ′. Results:
  - tWSe2: (+1,+1) at 3.15° and 1.70°, (0,0) at 1.54°, (−1,−1) at 1.25°. The paper text mentions only "(+1,+1) to (0,0)", so (−1,−1) is figure-only.
  - tMoTe2: (+1,−1) at 3.89° and 3.15°, (+1,+1) at 2.14°, (0,0) at 1.25°.
- **DFT bands are K↑/K′↓ doublets.** The two members are degenerate at γ to 10⁻³ meV and within 0.15 meV at κ. Energies are in eV relative to the VBM; path indices are 0 (γ), 20 (m), 40 (κ), 60 (γ). The archived SIESTA `.fdf` has `#SPIN spin-orbit` commented out, although the bands are clearly spin-split doublets. It is staged as found; no target depends on it.
- **fig1 maps.**
  - The grids include the periodic endpoint (100² and 10²). The first-star harmonics are (1,0), (0,1) and (1,1), so the grid basis is 60° apart. The authors' `plot-sky.py` plots them with 120° vectors.
  - Δ_T contains (0,0), (0,−1) and (−1,−1) with coefficient w. Together with right-handed axes, the stated forms fix b1 = g1 and b2 = g3.
  - Fitted values: set Q has V = 8.0 meV, ψ = −89.6°, w = −8.5 meV (the local-stacking set); set P has V = 20.8 meV, ψ = 107.7°, w = −23.8 meV (the 3.89° DFT fit).
  - Conjugating Δ_T (as `plot-sky.py` reads it) changes nothing, because it is gauge-equivalent. It was therefore not used as a fork.
- **Continuum results (FHS, 61/91/169 plane waves per layer, 12²–30² meshes, all identical).**
  - Q gives (−1,−1) and P gives (+1,−1) at 3.89°. This reproduces the paper's statement: K↑ top band −1 with local stacking and +1 with the DFT parameters.
  - RMS distance of the top two bands to the 3.89° SIESTA bands: P 1.6 meV, Q 7.6 meV.
  - P with fixed parameters at 2.14° stays at (+1,−1) (minimum direct gap 0.86 meV between bands 1 and 2), whereas DFT gives (+1,+1).
  - At 1.25° the fixed P model has a gap below 0.01 meV, so it is not graded.
- **Hartree potential.**
  - In all four relaxed structures, fractional (1/3,1/3) has top-layer Se above bottom-layer W (XM) and (2/3,2/3) has top W above bottom Se (MX).
  - Δv_H(XM) is −18.58 meV (3.15°), +8.92 (1.70°), +11.73 (1.54°) and +16.08 (1.25°): it changes sign between 3.15° and 1.70°.
  - The first-star coefficient goes from +4.13i (3.15°) through −0.07i (1.70°, nearly zero) to −0.98i (1.54°) and −2.37i (1.25°). The first-star part at XM is −21.5, +0.34, +5.04 and +12.25 meV.
  - At 1.70° the drop is carried by the second and third stars ((2,0)-type amplitude 1.63 meV).
  - The archived Fourier file holds magnitudes only. The 3.15° bottom map is archived as `pos_2d_bot_sm.dat`, a typo kept as found.

## Deviations from the design (candidates.json / DETAILS.md section 25)

- **Targets.**
  - The Chern pairs of all 8 calculations are graded, not 4.
  - The angle bracket and the inversion momenta are rubric items, not graded fields: they follow directly from the purchased C3 data.
  - The γ-point E1–E2 splitting was dropped, because the lab serves the band energies.
  - In their place the task grades the continuum Chern pairs (P, Q at 3.89°; P at 2.14°), the matching set, Δv_H(XM) at 4 angles and the first-star part at XM (1.70°, 1.25°).
- **The hole/electron sign fork (design dp3) was replaced.** A sign-unambiguous quantity (Δv_H = top − bottom as archived) is graded. Its sign forks are now the XM/MX location (dp8), first star vs full potential (dp9) and the magnitude-only Fourier file (dp10).
- **Interlayer tunnelling fork (design dp7) dropped.** The archive has no angle-dependent Δ_T. The fixed-model transfer fork (dp7 here) tests the same competing explanation (twist-independent continuum model).
- **No continuum-scoring request.** The agent has numpy and the staged DFT bands, so it can score models itself.
- **No wavefunction (atom_proj) request.** Not needed for any graded value; it remains a rubric topic only through the paper's argument.
- **The pot_top_bot_G.dat of 3.15° is not staged.** A coincidental substring matched a graded value at 4 significant digits; the file is sold by `hartree_fourier`.

## Decision points (10)

| dp | Fork | Kind |
|---|---|---|
| dp1 | K↑ vs K′↓ entries of each doublet | sector; signs flip |
| dp2 | spinful indicator minus sign | convention; mod-3 shift |
| dp3 | κ′ eigenvalue bought or assumed conjugate | discriminating purchase |
| dp4 | doublet pairing in the energy-ordered state list | band indexing |
| dp5 | continuum map basis (60° vs plot script's 120°, handedness) and layer assignment | evidence in the data vs scripts |
| dp6 | which parameter set matches DFT (literature prior vs bands) | model choice |
| dp7 | transfer of the 3.89°-fitted model to 2.14° vs buying DFT symmetry data | regime / budget |
| dp8 | XM vs MX location from the structures | position convention |
| dp9 | first-star vs full potential at small angles | model class |
| dp10 | magnitude-only Fourier file vs signed maps | budget / evidence |

Cascades:

- dp5 propagates into both continuum pairs.
- dp8 propagates into all six Hartree values.
- dp3/dp4 and dp7 compete with the Hartree maps for budget: maps cost 3 units each, and the oracle needs 17 C3 queries plus 7 units of Hartree data.

## Calibration

The full table is in `design.json` → `calibration`; `calibration.json` lists every ablation value.

| Field | Reference | Defensible variants (pass) | Nearest wrong paths (fail) | Tolerance |
|---|---|---|---|---|
| DFT Chern pairs | see above | — | spin-down sector (all signs flipped), no minus sign, κ′ = conj(κ), state index as band | exact |
| continuum 3.89° | P (+1,−1), Q (−1,−1) | 4/5/7 shells, 12²/18²/24² | plot-script basis Q (−1,0), left-handed Q (−1,+1), swapped layers P (−1,+1) Q (+1,+1) | exact |
| matching set | P | — | Q | exact |
| P at 2.14° | (+1,−1) | 5/7 shells, 18²/24²/30² | (+1,+1) (= DFT) | exact |
| Δv_H(XM) 3.15/1.70/1.54/1.25 | −18.577 / 8.916 / 11.735 / 16.083 | nearest grid, bilinear, archived line cut (≤ 0.035 meV away) | MX values 18.00 / −9.15 / −12.07 / −16.46; first-star only −21.5 / 0.34 / 5.04 / 12.25; no sign flip | ±1.0 meV |
| first star at XM 1.70/1.25 | 0.337 / 12.250 | FFT keeping the duplicated endpoint 0.305 / 12.195 | MX −0.33 / −12.31; magnitudes with 3.15° sign −0.30 / −12.23; full drop 8.92 / 16.08 | ±max(0.25, 5 %) |

## Oracle

The oracle uses no literature value beyond the catalog's m* and a₀:

- continuum Chern numbers of P and Q from the staged maps (91 plane waves per layer, 18² mesh), and P at 2.14°;
- the set choice from the RMS distance to the staged 3.89° bands;
- C3 at γ and κ for all 8 calculations plus κ′ once (17 units);
- Hartree maps at 1.70° and 1.25° (6 units) and the 1.54° line cut (1 unit);
- XM located from the staged CONTCARs.

That totals 20 evaluations and 24 of 27 units (89 %). Its Hartree errors are at most 0.003 meV. Lab job time is 2.3 s, and the container ran for 51 s, most of it the continuum Chern calculations. The estimated agent length is 45–75 tool calls.

- **Gate (2026-10-26, cpuset 4,5,20,21):** passed. All 15 ablations fail, the 4 robustness variants pass, rounding to 4 significant digits passes and grading is deterministic. The task has 10 decision points and 8 rubric items, is 15.2 MB and has 6 hidden values. After `design.json` was rewritten with the oracle length, the gate was re-graded with `--reuse-state` (passed).
- **Harbor oracle (`runs/build-tmd-moire-chern-twist/oracle-r1-tmd-moire-chern-twist`):** reward 1.0, 1 min 4 s.

## Open risks

- **Figure transcription.** The C3 labels and the Chern numbers they imply come from a figure, not an archived file. Every panel was cross-checked against the printed C_K↑. The labels for band 3 are not graded.
- **Memorisation.** The paper text states the tWSe2 (+1,+1)→(0,0) and tMoTe2 (+1,−1)→(+1,+1)→(0,0) sequences. The tWSe2 1.25° pair (−1,−1) is figure-only, and the other graded values (continuum, Hartree) are protocol-specific. A model that recalls the text and assumes "(0,0) at small angle" fails tWSe2 1.25°.
- **Mechanical indicator.** Once the indicator is applied correctly, the DFT Chern part is mechanical; its difficulty rests on dp1–dp4 and the budget. Development models may solve it; the continuum basis (dp5) and the Hartree forks (dp8–dp10) are expected to discriminate more.
- **Fixed-model 2.14° gap is small** (0.86 meV direct gap between bands 1 and 2), but the result is stable over cutoffs and meshes. A careless agent with a very coarse mesh still gets (+1,−1), because FHS is exact once the gap is open.
- **Thin margin.** The first-star part at 1.70° (0.34 meV) has a 0.25 meV tolerance. The estimator spread is 0.03 meV, and the nearest wrong path is 0.64 meV away.

## Disk

- **Peak:** about 110 MB of downloads, extraction and scratch, plus the 1.15 GB gate lab image layer.
- **After the build:** raw downloads deleted (refetch with `fetch.py`, `members.md5`); gate image removed; no Harbor `tmd-moire-chern-twist__*` images remained.
