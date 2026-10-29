# Construction of tb2j-groundstate-cuo-multiferroic

This is a released-start task for the benchmark (see `benchmarks/tools/BUILDING.md`). The source is Tellez-Mora, He, Bousquet, Wirtz and Romero, npj Comput. Mater. 10, 20 (2024), with data in Materials Cloud record rbrcp-0tc34 (CC BY 4.0). Provenance and checksums are in the task's `SOURCE.md`.

## Rebuild

1. `python3 benchmarks/construction/tb2j-groundstate-cuo-multiferroic/fetch.py` downloads `CuO.tar.gz` (6.5 MB, md5 b57a770e…) to `~/.cache/ai-theorist/tb2j-groundstate-cuo-multiferroic/`. It verifies md5 and sha256 and extracts the 11 members used, checking each member's sha256.
2. Inside the runtime image (`ai-theorist-benchmark-runtime:20261007`, `--cpuset-cpus <yours>`, `HOME` set so that the cache path resolves), run `python3 -B benchmarks/construction/tb2j-groundstate-cuo-multiferroic/make_task.py`. It takes about 1 minute and writes:
   - the catalog and lab `data.json` (3.1 MB);
   - `private/hidden.json` and `hidden_values.json`;
   - `solution/ablation_values.json`;
   - the budget in both `public.json` files;
   - `calibration.json` in this folder.

   It asserts that the authors' archived h(k) spectra of the FM/ligand and AFM/ligand calculations are reproduced from the archived pair lists, and that the k_min values equal the paper's SI values.
3. `.venv/bin/python -m benchmarks.tools.check_task tb2j-groundstate-cuo-multiferroic --cpuset 8,9,24,25` writes `check.json`.

`cuomodel.py` holds the parser, the native (as-archived) and symmetrised models, Luttinger–Tisza, the authors' bosonic matrix and the paraunitary magnon solver. `design.json` was written from `calibration.json`. `solution/solve.py` is a stand-alone re-implementation; it does not import this folder.

## Archive files used (shared source group)

From record rbrcp-0tc34 only `CuO.tar.gz` is used, and from it:

- `README`;
- the four `exchange.out`;
- the four `dispersion_*` files;
- `ligand/AFM/CuO_AFM.xsf`, for the full-precision AFM cell in the path check.

The `TB2J.yaml` files duplicate `exchange.out`; they were cross-checked to 5e-5 meV during exploration and are not used. The sibling task `tb2j-reference-selfconsistency-set` uses `FeP.tar.gz`, `FeCl2.tar.gz`, `FePS3.tar.gz` and `NiO.tar.gz` from the same record. This task touches none of them, and its targets (CuO k_min, stability verdicts, AF1 magnons, chain coupling) do not overlap with that design's FeP helix, magnon or reference-change targets.

## What the data show (found during construction)

- **Both FM references are unstable.** The ligand-corrected k_min is (0.4911, 0, −0.4736) r.l.u. in the staged axes, which is the paper's SI (0.491, 0.0, 0.474) with the opposite c sign. The uncorrected k_min is (0, 0.5042, 0); the SI gives (0, 0.5, 0). Both need couplings beyond about 8 Å. At a 6 Å cutoff the vectors move by 0.010 and 0.006 r.l.u.; at 4–5 Å the ligand-corrected minimum moves to another point of the zone.
- **The AFM/ligand calculation is self-consistent.** It uses a 2×1×2 cell whose origin is shifted by c/2 relative to the FM cell, with cell lengths differing by up to 0.4 %. Its moments are AF1, a single Fourier component at (1/2, 0, 1/2) with the two plaquette orientations in antiphase. With its own pair list, h(k) is positive semi-definite over the zone (Goldstone zero only), and AF1 is the top eigenvector of its K(0).
- **The AFM/no_ligand calculation is broken.** It uses a {2a, b, a+c} cell, not the 1×2×1 cell the SI describes. Its structure has two inequivalent Cu sites, with different Cu–O bond sets (1.992/2.021 Å and 2.004/2.028 Å, against 2.000/2.030 Å in the FM cell) and two moment sizes, 0.686 and 0.660 μB. Nominally equivalent bonds take −27.8 and +27.9 meV, and its lowest h eigenvalue at its own reference is −5.2 meV. The SI says this branch never converged. Its archived dispersion file cannot be reproduced from its pair list by any collinear configuration, and the FM/no_ligand file is off by 0.4 %. They are served as archived, and nothing graded depends on them.
- **The AFM-reference pair list carries reference dependence.** Bonds equivalent by crystal symmetry differ with the relative spin orientation: antiparallel pairs are more antiferromagnetic, by 1.07 meV at 2.95 Å and 0.80 meV at 3.04 Å. Symmetrising the list does three things:
  - it lowers the AF1 stabilisation from 44.33 to 40.18 meV/Cu;
  - it makes an incommensurate spiral at (0.502, 0, 0.516) the Luttinger–Tisza optimum, only 0.04 meV/Cu below AF1 (and close to the experimental AF2 vector, a coincidence at the level of the method's resolution);
  - it lowers the magnon energies by about 14 %.

  This is dp5.
- **The archived "dispersion" files are h(k) eigenvalues, not magnon energies.** The AFM file was evaluated at the FM path's Cartesian k-vectors in the AFM frame, which is rotated by about 2°. Its path is therefore slightly off the high-symmetry points.

## Deviations from the design (candidates.json / DETAILS.md section 21)

- **Anisotropic/DMI fork dropped.** Only isotropic exchange is served. The archive README declares the anisotropic constants invalid, and they come from collinear or single-axis runs; TB2J's own header warns about this for one file. Offering them would either be an instructed trap or would require serving data the authors disown.
- **Classical Monte Carlo T_N not graded.** The self-consistent exchange depends on the reference: it is defined pair by pair in the AF1 supercell. A finite-temperature simulation must extend it to disordered states, and that choice moves T_N by more than 20 %. It is a rubric item instead: mean field overestimates T_N for this quasi-1D magnet, and T_N is not supported by these data.
- **"Magnon energies at hidden q" became zone-averaged quantities.** The graded quantities are the mean magnon energy and the band maximum. Single-q energies depend on whether the 16-site or the 8-site magnetic cell is used, through band folding; these two quantities do not.
- **"Dominant coupling within 15 %" became within 4 %.** The FM-reference value is only 9.7 % off, so 15 % would not test the reference fork. The defensible spread over single bonds is ±0.3 %.
- **Ground-state target split into three parts:**
  - k_min of both FM references (graded vectors);
  - stability verdicts for all four calculations;
  - the converged propagation vector.

  The Luttinger–Tisza constraint fork (design dp5) does not bite: both plaquette sublattices carry equal weight in every optimum.
- **Spin-wave specification.** The instruction states that the spin waves use the archived pairwise exchange in the calculation's own cell, "as the authors' workflow does". Without this the magnon targets would be ambiguous between the native and symmetrised models. The data, and the authors' archived h(k) at Γ, show which is meant. The symmetrising temptation remains (dp5) but is partly instructed.

## Decision points (7)

| dp | Fork | Kind |
|---|---|---|
| dp1 | FM-reference exchange as final | self-consistency loop; reference dependence |
| dp2 | uncorrected workflow accepted without testing | metadata inconsistency found in the purchase |
| dp3 | exchange range (long windows cost 2–3 units) | convergence under a binding budget |
| dp4 | commensurate approximant, paper SI or AF2 quoted as k_min | regime judgement; anti-memorisation (the SI's c sign differs from the staged axes) |
| dp5 | symmetrising the AFM-reference list | evidence in the data; cascade into q and magnons |
| dp6 | h(k) eigenvalues as magnon energies | confounder in the archived files |
| dp7 | S = 1 / pair counting | convention |

## Calibration (full table in `design.json` → `calibration`)

| Field | Reference | Defensible spread | Nearest wrong path | Tolerance |
|---|---|---|---|---|
| k_min, ligand-corrected | (0.4911, 0, −0.4736) | ≤ 0.0023 (8–17 Å cutoffs) | 6 Å cutoff 0.0102; paper SI verbatim 0.0179; AF2 0.0176; commensurate 0.0278 | 0.006 r.l.u. |
| k_min, uncorrected | (0, 0.5042, 0) | ≤ 0.0007 | commensurate 0.0042; 6 Å cutoff 0.0063 | 0.003 r.l.u. |
| converged q | (1/2, 0, 1/2) | exact | symmetrised 0.0165; FM reference 0.0278 | 0.006 r.l.u. |
| stability verdicts | F, T, F, F | robust to any cutoff ≥ 4 Å | AFM/no_ligand judged converged | exact |
| chain J (S = 1/2, meV) | 129.74 | 129.39–130.11 | FM reference 117.11 (−9.7 %); uncorrected 170.08; pair counting 64.87 / 32.44 | ±4 % |
| mean magnon (meV) | 149.12 | 148.89–149.90 (6 Å: 145.69) | symmetrised 128.43; FM reference 132.43; h eigenvalues 177.32; S = 1 74.56 | ±5 % |
| band maximum (meV) | 190.49 | 188.36–191.01 (6 Å: 184.50) | FM reference 178.56 (−6.3 %); symmetrised 174.54; h eigenvalues 327.48; S = 1 95.24 | ±4 % |

The tightest margins are:

- k_min uncorrected: 0.0012 r.l.u. to the commensurate (0, 1/2, 0);
- band maximum: 2.3 percentage points to the FM-reference path;
- k_min ligand-corrected: 0.004 r.l.u. to the 6 Å cutoff.

## Oracle and length

- **Oracle:** 20 evaluations, 21 of 24 units (87.5 %). The catalogue is worth 68 units, so the budget covers 35 % of it. It buys:
  - FM/ligand 0–17 Å (8 units), with a 14 Å versus 17 Å convergence check;
  - FM/no_ligand 0–12 Å (5 units);
  - AFM/ligand 0–14 Å (6 units);
  - AFM/no_ligand 0–4 Å (1 unit);
  - the AFM/ligand author dispersion (1 unit), to check its own h(0) against the archive.
- **Graded errors:** k_min 1.4e-5 and 1.8e-4 r.l.u.; q exact; chain 0; magnons −0.013 % and −0.010 %.
- **Runtime:** 26 s in the gate container; the Harbor oracle took 65 s, with reward 1.0.
- **Expected agent length:** 45–75 tool calls.

## Open risks

- **Dev-easy risk.**
  - The pipeline (Luttinger–Tisza, a 16-site paraunitary LSWT) is standard for strong models.
  - The converged q (AF1) is the well-known experimental vector.
  - dp5 is partly instructed by the uniqueness clause.
  - The difficulty rests on doing every step correctly under the budget: exchange range, testing all four calculations, native supercell LSWT with S = 1/2 and the factor conventions.
- **Tight uncorrected-k_min tolerance (0.003).** It is needed to reject the commensurate (0, 1/2, 0), the paper's rounded value. An agent that rounds its own minimiser to two decimals fails; the instruction asks for three.
- **Model-derived truth.** The graded numbers are computed from the archived pair lists, not archived themselves. The authors' archived h(k) for FM/ligand and AFM/ligand is reproduced to 1e-5–1e-4 eV, which pins the implementation, including the native use of the AFM pair list.
- **Uncorrected files are inconsistent.** The AFM/no_ligand data contradict the paper's description, and the uncorrected dispersion files do not match their pair lists. They are served as archived. Only the stability verdict uses AFM/no_ligand, and it is robust to the exchange cutoff (lowest eigenvalue −5 to −9 meV).

## Disk

- **Peak:** about 50 MB: 39 MB in the cache (`CuO.tar.gz` and a full unpack including the `TB2J.yaml` files) and 10 MB of paper and SI PDFs in the scratch folder.
- **After the gate:** the raw downloads were deleted. Kept: the refetch script with checksums, `calibration.json` and `check.json`.
- **Task package:** 3.1 MB.

## Gate

- **check_task (2026-10-19):** passed. 7 decision points, all 12 ablations fail, all 7 robustness variants pass, rounding to 4 significant digits passes, 9 hidden values with no leaks.
- **Harbor oracle:** reward 1.0, in `runs/build-tb2j-groundstate-cuo-multiferroic/oracle-r1-tb2j-groundstate-cuo-multiferroic`.
