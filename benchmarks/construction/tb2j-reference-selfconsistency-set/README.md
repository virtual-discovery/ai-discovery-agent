# Construction of tb2j-reference-selfconsistency-set

This is a released-start task for the benchmark (see `benchmarks/tools/BUILDING.md`). The source is Tellez-Mora, He, Bousquet, Wirtz and Romero, npj Comput. Mater. 10, 20 (2024), with data in Materials Cloud record rbrcp-0tc34 (CC BY 4.0). Provenance and checksums are in the task's `SOURCE.md`.

## Rebuild

1. `python3 benchmarks/construction/tb2j-reference-selfconsistency-set/fetch.py` downloads `NiO.tar.gz`, `FeCl2.tar.gz`, `FePS3.tar.gz` and `FeP.tar.gz` (12.1 MB) to `~/.cache/ai-theorist/tb2j-reference-selfconsistency-set/`. It verifies md5 and sha256 and extracts the 18 members used, checking each member's sha256.
2. Inside the runtime image (`ai-theorist-benchmark-runtime:20261007`, `--cpuset-cpus <yours>`, `--user $(id -u):$(id -g)`, `HOME` set so that the cache path resolves), run `python3 -B benchmarks/construction/tb2j-reference-selfconsistency-set/make_task.py`. It takes about 1 minute and writes:
   - the catalog and lab `data.json` (2.4 MB);
   - `private/hidden.json` and `hidden_values.json`;
   - `solution/ablation_values.json`;
   - the budget in both `public.json` files;
   - `calibration.json` in this folder.

   It asserts that the paper's FM-reference results and the authors' archived h(k) files are reproduced (see Validation), and that the FePS3 DFT+U FM list duplicates the U = 0 list.
3. `.venv/bin/python -m benchmarks.tools.check_task tb2j-reference-selfconsistency-set --cpuset 6,7,22,23` writes `check.json`.

`tb2jmodel.py` holds the parser, the isotropic model in each calculation's own cell, the noncollinear Toth–Lake bosonic matrix, the classical relaxation, the minimum-image remapping, the folding and unfolding between the FM cell and the AFM supercells, the k_min search and the stabilization energies. `design.json` was written from `calibration.json`. `solution/solve.py` is a stand-alone re-implementation and does not import this folder.

## Archive files used (shared source group)

From record rbrcp-0tc34 only `NiO.tar.gz`, `FeCl2.tar.gz`, `FePS3.tar.gz` and `FeP.tar.gz` are used. The served data are the ten `exchange.out` files (isotropic pairs, header, cell, atoms). Four `dispersion_*` files are used for validation only. The sibling task `tb2j-groundstate-cuo-multiferroic` uses only `CuO.tar.gz`. MnF2, CrSBr and Mn5Si3 are not used. MnF2 also appears in `heis-configdesign-mnf2-altermagnet`, which draws on a different archive.

## What the data show (found during construction)

- **The FM-reference results of the paper are reproduced.**
  - FePS3 (U = 0): Z lies 2.00 meV below Γ.
  - FeP: k_min = (0, 0.0306, 0.2179), against the paper's (0.00, 0.03, 0.21).
  - The archived NiO couplings are the PBE+U rows of Table 1.
- **FeP needs very long-range exchange.**
  - The FM-reference k_min jumps from (0, 0.35, 0) at 10 Å, and (0, 0.33, 0.09) at 12 Å, to (0, 0.021, 0.224) at 14 Å. It converges only beyond about 16 Å.
  - The FM-reference stabilization energy is still 15 % too high with the 0–18 Å list; it converges with 0–23 Å.
  - In the insulators everything that matters lies within 11–14 Å. FeCl2's stacking energy needs the interlayer couplings at 6.4–10 Å, which are 0.003–0.03 meV.
- **The FeCl2 and FeP AFM supercell lists are aliased.** They were computed with R = 0 along the long axis only, so bonds that cross the cell boundary are listed at the far image:
  - FeP's 3.81 Å bond (−4.17 meV) appears at 25.69 Å;
  - FeCl2's 7.34 Å interlayer bonds appear at 30.5 Å.

  Consequences:
  - The total K(0), and so every classical energy, is correct only if every listed window is kept.
  - Any k-dependence along the long axis is lost. The authors' archived FeP AFM h(k) is exactly flat along Γ–Z (to 1e-8 meV), and the as-listed model reproduces that file to 0.1 meV.
  - Folding the spiral list onto the four-Fe cell at the listed vectors gives a spurious k_min = (0, 1/2, 0). With minimum-image vectors it gives (0.002, 0.027, 0.1885).
- **The FePS3 DFT+U FM file is the U = 0 FM calculation.** It has the same cell, atoms (moment 3.377 μB) and 1196 couplings. Only the TB2J version, date, working directory and pair order differ. The DFT+U AFM run has 3.642 μB and a 1 % larger cell. Using the file as a U = 3 eV FM reference gives a spurious 44 % "reference effect" in an insulator.
- **The FeP spiral moments are not a stationary point** of their own exchange (torques 0.3–4.3 meV).
  - The bosonic matrix at the archived moments has a −0.99 meV eigenvalue at Γ, also visible (−1.10 meV) in the authors' file.
  - Relaxing the directions in the 1×1×5 cell (a unique minimum: all 40 random starts reach it) removes it with as-listed couplings. With minimum-image couplings a −0.077 meV instability remains near k_z = 0.1 of the supercell, consistent with the converged exchange's preference for q_c ≈ 0.19. This lies inside the paper's 0.5 meV ε window.
  - Relaxation changes the stabilization energies by only 0.2 % (converged) and 2 % (FM reference), so it is a rubric point rather than a graded fork.
- **Reference dependence.** The ratio of converged to FM-reference stabilization energy is:
  - NiO 0.988;
  - FeCl2 0.970;
  - FePS3 (U = 0) 1.113, with the AFM structure zigzag-distorted (J1' bond 3.47 → 3.64 Å), so lattice and reference effects are confounded;
  - FeP 7.61.

  In FeP the net coupling sum changes from +15.7 to −48.5 meV per Fe, and the short antiferromagnetic bonds strengthen by factors of 2–5.
- **FeP phase arrangement.** The FM-reference Luttinger–Tisza helix of FeP has a different intra-cell phase arrangement from the converged double helix. Its own optimum is stabilized by 29.1 meV/Fe, but it gives only 13.4 meV/Fe on the converged arrangement.

## Deviations from the design (candidates.json / DETAILS.md section 22)

- **Categorical ground states not graded.** They are the experimental structures, the AFM metadata reveal them, and they are memorizable. Reference dependence is graded instead through the stabilization energies from both exchange sets for every material.
- **"Fractional change of dominant J within 20 %" became stabilization energies within 4 %.** Per-bond comparisons need a bond mapping across differently relaxed cells, and 20 % would separate neither NiO nor FeCl2 from their wrong paths.
- **"LSWT magnon energies at hidden q" dropped.**
  - For FeP, the spiral is not stationary and its list is aliased along the helix axis. Magnon energies along c then depend on how the images are assigned (0 against 5.8 meV at Z), and the minimum-image model is marginally unstable. They are not uniquely defined.
  - For the insulators, magnons would repeat the CuO sibling's targets.
- **FeP helix target split in two:** the FM-reference k_min and the k_min of the converged (spiral-reference) exchange, which is the self-consistency test.
- **New decision points from the data:** aliased supercell lists (dp4, dp5), the duplicated DFT+U file (dp6) and observable-dependent range convergence (dp3). The design's "budget allocation" fork is realized through the binding budget (68 of 160 units) rather than as a separate graded item.

## Decision points (8)

| dp | Fork | Kind |
|---|---|---|
| dp1 | FeP exchange range for k_min | convergence under a binding budget |
| dp2 | commensurate approximant or paper vector quoted as k_min | regime judgement; anti-memorization |
| dp3 | range that converges k_min reused for energies | cascade |
| dp4 | long windows of the aliased AFM supercell lists skipped | evidence in raw output |
| dp5 | spiral list folded at the archived (far-image) vectors | aliasing artefact; cross-check |
| dp6 | FePS3 DFT+U FM file used as a U = 3 eV reference | metadata inconsistency |
| dp7 | FM reference's own helix energy, or reference independence assumed | competing explanations |
| dp8 | each bond counted once | convention |

## Calibration (full table in `design.json` → `calibration`)

| Field | Reference | Defensible spread | Nearest wrong path | Tolerance |
|---|---|---|---|---|
| FeP FM k_min | (0, 0.0306, 0.2179) | ≤ 0.0028 (18 and 23 Å lists) | paper (0, 0.03, 0.21) 0.0080; 14 Å list 0.0118; (0, 0, 0.2) 0.0354 | 0.005 r.l.u. |
| FeP converged k_min | (0.0021, 0.0268, 0.1885) | ≤ 0.0058 (sum over five, dropping half-box pairs, 0–30 Å) | as listed, sum over five 0.0186; (0, 0, 0.2) 0.0293; FM-reference value 0.0297; as listed 0.509 | 0.012 r.l.u. |
| NiO ΔE fm / conv (meV/Ni) | 139.17 / 137.42 | −1.5 % / −0.4 % (0–5 Å) | bond counted once −50 % | 4 % |
| FeCl2 ΔE fm / conv (meV/Fe) | 0.4932 / 0.4785 | −0.2 % (0–11 Å) | conv without 30–40 Å −16 %; fm 0–8 Å −33 % | 4 % |
| FePS3 U = 0 fm / conv | 48.48 / 53.96 | −2.4 % / −1.6 % (0–14 Å) | conv 0–11 Å −14 %; fm 0–11 Å −3.8 % (passes alone) | 4 % |
| FePS3 U = 3 fm / conv | null / 27.04 | conv 0–11 Å −0.3 % | duplicate file 48.48 | null / 4 % |
| FeP fm / conv (meV/Fe) | 13.39 / 101.91 | +1.9 % / +0.2 % (relaxed spiral) | fm 0–18 Å +15 %; fm own LT helix +117 %; conv without 23–30 Å −5.4 % | 4 % |

The tightest margins are:

- FeP FM k_min: 0.003 r.l.u. to the paper's rounded vector;
- FeP converged energy: 1.4 points to the skipped 23–30 Å window;
- FeP converged k_min: 0.0066 to the as-listed sum-over-five folding.

## Oracle and length

- **Oracle:** 45 evaluations, 60 of 68 units (88 %). The catalog is worth 160 units, so the budget covers 42.5 % of it. It buys:
  - NiO FM/AFM 0–8 Å (2 + 2 units);
  - FeCl2 FM 0–14 Å (4 units);
  - FePS3 U = 0 FM/AFM 0–18 Å (6 + 6 units);
  - FePS3 DFT+U AFM 0–11 Å (3 units);
  - FePS3 DFT+U FM 0–5 Å (1 unit, the provenance check);
  - FeCl2 AFM 0–40 Å (13 units, after detecting R = 0 along c);
  - FeP FM 0–30 Å (10 units, with a 14/18/23/30 Å k_min convergence test);
  - FeP AFM 0–40 Å (13 units).
- **Graded errors:** both k_min vectors exact; energies within 0.75 %.
- **Runtime:** about 30 s in the gate container; the Harbor oracle took 74 s, with reward 1.0.
- **Expected agent length:** 50–80 tool calls.

## Open risks

- **Dev-easy risk.** An agent that buys every window of the two expensive supercell lists avoids dp4 by brute force. The budget allows this only if it economizes elsewhere. dp5 and dp6 still require noticing the aliasing and the duplicate.
- **Lucky passes on the converged FeP k_min.** Some wrong paths pass by chance: folding as listed from a 0–23 Å purchase gives 0.182, within tolerance. The k_min of the converged exchange is a delicate long-range balance, and partial lists scatter between 0.18 and 0.22.
- **Null answer.** The FePS3 DFT+U FM entry must be null. The evidence is strong (identical cell, moments and 1196 couplings to 1e-4 meV), and the instruction explicitly allows null.
- **Model-derived truth.** The graded numbers are computed from the archived pair lists, not archived themselves. The implementation is pinned by the paper's FePS3 and FeP FM results and by the authors' archived h(k) for FePS3, FeCl2 and FeP.

## Disk

- **Peak:** about 76 MB: 66 MB in the cache (six downloaded archives with a full unpack, including MnF2 and CrSBr, which were inspected and not used) and 10 MB of paper and SI PDFs in the scratch folder.
- **After the gate:** the raw downloads were deleted. Kept: the refetch script with checksums, `calibration.json` and `check.json`.
- **Task package:** 2.4 MB.

## Gate

- **check_task (2026-10-19):** passed. 8 decision points, all 14 ablations fail, all 7 robustness variants pass, rounding to 4 significant digits passes, 13 hidden values with no leaks.
- **Harbor oracle:** reward 1.0, in `runs/build-tb2j-reference-selfconsistency-set/oracle-r1-tb2j-reference-selfconsistency-set`.
