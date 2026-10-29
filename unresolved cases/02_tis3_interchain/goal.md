# TIS3-INTERCHAIN-1 — Is the 1.28 eV low-temperature absorption peak of TiS3 for E ⊥ chains an interchain exciton?

**Source finding.** Boldyrev et al., JETP Letters 120, 565–572 (2024), DOI 10.1134/S002136402460335X (shortlist entry 03): polarized IR/optical
transmission of TiS3 whiskers, room temperature to liquid helium. Below ≈ 150 K an absorption peak at 1.28 eV appears, predominantly for light
polarized perpendicular to the Ti–S chains (chains ∥ b). The authors propose an interchain exciton and note that earlier monolayer
calculations show a different polarization trend. No first-principles calculation in the paper.

**Thread.** Standalone campaign; pipeline parent = RuCl3 round 8 (chain driver the round-8 chain (stage list + templates in ../SHARED_RUNBOOK.md); real-space pipeline
../env/run_realspace.sh + pick_states.py, plotxct rules in SHARED_RUNBOOK.md §5; GPU recipe ../env/setup_gpu_vast.sh).
Sibling card: ZRSE3-STRAIN-1 (same P2_1/m 8-atom structure type; one chain driver serves both). Literature anchors: Molina-Mendoza et al.,
Adv. Electron. Mater. 1, 1500126 (2015): bulk electronic gap 1.2 eV (STS), optical gap 1.07 eV, E_b ≈ 0.13 eV, bulk RPA/BSE in the paper;
monolayer GW-BSE (e.g. Chem. Phys. 2022, DOI 10.1016/j.chemphys.2022.111678 family): strongly anisotropic, lowest bright exciton E ∥ b,
E ∥ a exciton higher with a different binding energy. The per-polarization edge energies and the 1.28 eV peak width of Boldyrev et al. are
behind a paywall on this machine: **S0 freezes them from the paper's figures** (PI has access) before any production run.

**Status vocabulary (terminal, exactly one):** `INTERCHAIN_CONFIRMED` | `INTRACHAIN_ASSIGNED` | `MIXED_CHARACTER` | `PEAK_NOT_REPRODUCED` |
`LITCHECK_FAIL` | `PROCESS_BLOCKED`.

## 1. Why this campaign exists

The observable is exactly what round 8 produces: polarization-resolved BSE absorption (momentum operator b1/b2/b3), the identity of the
leading states behind each peak (ratio-of-ratios, state overlaps), and real-space electron densities with the hole fixed (plotxct on the
Born–von Kármán box). The hypothesis "interchain" is a statement about the electron density on chains other than the hole's chain, which the
maps measure directly. Round-8 lessons that shape the card: (i) gate on the peak's energy relative to an anchored edge exciton after ONE
global shift, with width-scaled tolerance, not on absolute GW-BSE energies; (ii) report E_QP and E_b separately so a gap error and
under-binding cannot compensate; (iii) a negative control makes every character measure non-vacuous (the E ∥ b edge exciton must come out
intrachain by the same measure); (iv) a dimensionality control (monolayer chain) is required before any "theory vs experiment" statement,
because the earlier monolayer trend is what the authors flag.

## 2. Read this and nothing more

- This file, `acceptance_draft.json` (this folder) (frozen before S2 as `acceptance.json` + sha256, with the S0 anchors filled in).
- Structure: ICSD TiS3, P2_1/m, 8 atoms (a ≈ 4.96, b ≈ 3.40, c ≈ 8.78 Å, β ≈ 97.3°; chains ∥ b, layers in ab, stacked along c). Monolayer =
  one ab layer, vacuum ≥ 15 Å, BerkeleyGW slab truncation.
- Pseudopotentials: PseudoDojo v0.4 NC PBE scalar-relativistic Ti (12 e, semicore) and S (6 e), standard.
- Round-8 driver, realspace pipeline, GPU recipe, Vast skill. The three papers above, anchors only.

## 3. Definitions

- **Channels.** b1 = E ∥ b (chain axis), b2 = E ∥ a (in-plane ⊥ chains), b3 = E ∥ c* (recorded only). The experiment's "⊥ chains" on a
  whisker face is b2 unless S0 shows the measured face contains c*; S0 records which.
- **Edge exciton X_b.** Leading bright stick of the b1 channel; anchor = the measured E ∥ b edge/exciton energy frozen at S0 (fallback: optical
  gap 1.07 eV, Molina-Mendoza).
- **Global shift s.** ONE shift, |s| ≤ 0.15 eV, that puts X_b on its anchor. Applied to every channel.
- **Candidate peak P_a.** In the b2 channel after the shift: the strongest local maximum of ε2 (gaussian 0.03 eV, matching a low-T line) in
  [1.28 − 0.10, 1.28 + 0.10] eV; tolerance 0.10 = width class of a 1.28 eV peak measured at 4–150 K (replace by half the measured FWHM if S0
  gives it, floor 0.06).
- **Selectivity R.** ε2_b2(E_P) / ε2_b1(E_P) at the peak energy, same broadening.
- **Identity.** Leading state(s) carrying ≥ 50 % of P_a's stick weight within ±0.03 eV; their eigenvector overlap with X_b's leading state; their
  binding E_b = E_QP,min(k of dominant transition) − E.
- **Interchain weight w_ic.** From plotxct with the hole on a Ti site (and, as a check, on the bonding S site) of chain 0: w_ic = electron
  density integrated over all chains of the same layer other than chain 0, divided by the in-layer total; w_il = density in other layers.
  Box = Born–von Kármán cell of the fine grid, hole translated to the box centre, OMP_NUM_THREADS = 1, log-scale renders, containment
  (face density / max) ≤ 1e−5 for a bound state. For a weakly bound bulk exciton (E_b ≈ 0.13 eV, radius a few nm) the fine grid must be
  ≥ 8×12×4 (box ≈ 40×41×35 Å); the containment check decides, not the assumption.
- **Chain assignment of a grid point:** nearest Ti–S chain axis in the ab plane (chains at fixed fractional a positions in P2_1/m; two chains
  per cell). The assignment script and a figure of it are deliverables.

Parameter count: zero fitted parameters. Covariates: coarse/fine k ladder (two levels), band window of the kernel (two levels).

## 4. Frozen discriminators

| gate | statement | PASS | FAIL |
|---|---|---|---|
| G0 LITCHECK | bulk chain vs the Li = Molina-Mendoza et al. 2015 bulk G0W0+BSE (Yambo, PPA, 10×10×10, 30v×70c): E_QP 1.15, E_opt 1.05, E_b 0.10 eV, and the presence of their second feature at ~1.2–1.3 eV | all three within 0.10 eV and the second feature present | else `LITCHECK_FAIL`, stop and audit (we were the outlier once) |
| G1 edge anchor | bulk: X_b after the shift; E_QP (G0W0 direct) vs 1.2 eV (STS); E_b(X_b) vs 0.13 eV | \|s\| ≤ 0.15; E_QP within ±0.25; E_b within ±0.07 (both components separately) | E_b off by more → record `BINDING_OFF`, proceed; \|s\| > 0.15 → `PROCESS_BLOCKED` (chain not anchored) |
| G2 peak exists | P_a in the b2 channel in [1.18, 1.38] after the shift, prominence ≥ 0.20 of the b2 maximum in [0.9, 1.8] | yes | no local maximum → `PEAK_NOT_REPRODUCED` |
| G2b selectivity | R at E_P | R ≥ 3 ("predominantly perpendicular") | 1 < R < 3 → recorded as weak selectivity (non-blocking); R ≤ 1 → `PEAK_NOT_REPRODUCED` |
| G3 identity | leading state(s) of P_a | bound (E_b ≥ 0.03 eV) and distinct from X_b (overlap ≤ 0.1, k/band composition reported) | unbound (resonance) or same state as X_b → `INTRACHAIN_ASSIGNED` is not available; verdict `MIXED_CHARACTER` with the composition |
| G4 interchain character | w_ic of P_a's leading state(s), hole on Ti and on S (both must agree in class) | w_ic ≥ 0.50 → `INTERCHAIN_CONFIRMED`; w_ic ≤ 0.25 → `INTRACHAIN_ASSIGNED`; between → `MIXED_CHARACTER` | — |
| G4c control | w_ic of X_b by the same measure | ≤ 0.35 (the measure discriminates) | > 0.35 → the measure is non-discriminating → `PROCESS_BLOCKED`, report both |
| G5 dimensionality (control only, not a LITCHECK) | monolayer chain, same settings: energy of the lowest b2 exciton relative to X_b (Δ_mono) vs the bulk Δ_bulk = E_P − E(X_b) | the sign/size difference is stated with numbers; if \|Δ_mono − Δ_bulk\| ≥ 0.10 eV the "dimensionality explains the different monolayer trend" sentence is allowed | otherwise that sentence is forbidden |
| G6 convergence | k ladder level 1 → 2: E_P − E(X_b), R, w_ic | ≤ 0.03 eV, ≤ 25 % relative, ≤ 0.10 | else one refinement, then `PROCESS_BLOCKED` |
| G7 window | kernel band window ×1.5 | E_P − E(X_b) ≤ 0.03 eV change | else report as a range (round-8 window lesson) |

Non-vacuity: G2 can fail (no b2 peak in the window), G4 returns three classes, G4c can fail, G0 can fail. Recorded before reading numbers.

## 5. Verdict and stopping rule

| outcome | status |
|---|---|
| G0 fails | `LITCHECK_FAIL` |
| G1 shift cap fails or G4c fails | `PROCESS_BLOCKED` |
| G2 or G2b(R ≤ 1) fails | `PEAK_NOT_REPRODUCED` (report what the b2 channel does show in 0.9–1.8 eV) |
| G2, G3 pass, G4 class | `INTERCHAIN_CONFIRMED` / `INTRACHAIN_ASSIGNED` / `MIXED_CHARACTER` |

Stop at the verdict. The temperature onset (peak visible only below ≈ 150 K) is NOT predicted by this chain (no electron–phonon): the record
says so in one sentence and does not speculate beyond linewidth.

## 6. Out of scope

Temperature dependence, phonons and exciton–phonon coupling, defects/surface states, evGW, SOC (Ti/S light; declared off), the far-IR phonon
findings of the same paper, few-layer thickness series, substrate.

## 7. Deliverables — exactly these 9

1. `acceptance.json` + sha256 with the S0 anchors (edge energies per polarization, 1.28 eV peak FWHM, measured face/channel mapping).
2. `results/spectra.json`: b1/b2/b3 ε2 and stick lists, bulk and monolayer, both k levels, both windows, shift s.
3. `results/states.json`: leading states of X_b and P_a (energy, strength, E_b, k/band composition, overlap).
4. `results/realspace.json`: w_ic, w_il, containment, box, hole site, for P_a and X_b (control), both hole sites.
5. Fig. 1 eye test: b1 and b2 spectra vs the measured transmission/absorption traces at the same broadening, 1.28 eV marked.
6. Fig. 2 real-space maps (log scale, BvK box) for P_a and X_b, hole on Ti; chain assignment overlay.
7. Fig. 3 bulk vs monolayer polarization trend (Δ_bulk vs Δ_mono) with the published monolayer values.
8. `RESULT.md` (verdict box, gate table, captions, provenance), `make_figures.py`, `MANIFEST.sha256`.
9. `vast_run.json` per instance + destroy verification.

## 8. Plan

| stage | what | where | est. |
|---|---|---|---|
| S0 | anchors from the paper (PI/user supplies the figure readings); cell + pseudos + cutoff ladder; bulk PBE-D3 relax (internal coordinates, experimental lattice); monolayer cell; which crystal face was measured → channel mapping | local | 0.5 day |
| S1 | DFT bands bulk + monolayer; KS gaps; check the chain geometry assignment script on the DFT charge density | local | 2 h |
| S2 | freeze the card; `tis3_chain.py` = ZrSe3 driver + slab truncation branch + b1/b2/b3 absorption + eigenvectors output; frozen scorer `interchain_score.py` (spectra gates) + `realspace_score.py` (w_ic) | local | 0.5–1 day |
| S3 | pilot bulk chain, level-1 k (coarse 4×6×2, ≈ 300 bands, 10 Ry) + fine 8×12×4 absorption; time it; choose local vs Vast | local | 2–4 h |
| S4 | bulk level-1 chain = G0 LITCHECK against the Molina-Mendoza numbers (stop if it fails); then the monolayer control chain (G5) | local/Vast | 3–6 GPU-h |
| S5 | bulk production level 1: G1–G3; eigenvectors kept WITH their WFN (gauge trap); plotxct on the BvK box for P_a and X_b, two hole sites (G4, G4c) | Vast (maps run remotely; local disk 21 GB) | 4–8 GPU-h |
| S6 | level-2 k (6×8×3 coarse, 12×16×6 fine) and window ×1.5 on bulk (G6, G7) | Vast | 4–8 GPU-h |
| S7 | scorer, figures, RESULT.md, MANIFEST; destroy in the same pipeline as pull + md5 | local | 0.5 day |

Compute budget: ≈ 15–25 GPU-h ⇒ ≈ $45–75 at a verified ≤ $3/h H100; together with ZRSE3-STRAIN-1 this exceeds the $74.87 credit
recorded 2026-10-12 — run ZrSe3 first (cheaper, no maps) and the TiS3 pilot locally; top-up is a user decision. Toolchain: rebuild from source
on the instance unless the CPU model matches the tarball's build host (Illegal-instruction trap); wait for `make` to finish before
collecting bin/ (setup.sh race fix is in gpu_b/setup.sh).

## 9. Honesty constraints

- All gates scored from files by the frozen scorers; coverage stated as "N of M planned stages ran".
- A bound-state claim requires the containment check to pass on the box used; otherwise w_ic is reported as "box-limited" and G4 is not scored.
- E_QP and E_b are always reported separately next to the shifted peak energy; the shift is printed on every figure.
- The monolayer control (G0, G5) runs BEFORE the bulk verdict is written; a verdict without it is `PROCESS_BLOCKED`.
- Nothing is tuned to put a peak at 1.28 eV; the only freedom is the single shift fixed by X_b.
