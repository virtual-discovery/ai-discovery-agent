# ZRSE3-STRAIN-1 — Which term drives the 136 meV/% exciton redshift of ZrSe3 under biaxial compression: E_QP(ε) or E_b(ε)?

**Source finding.** Li et al., Nano Research 19, 94908438 (2026), DOI 10.26599/NR.2026.94908438 (shortlist entry 06): multilayer ZrSe3 on
polycarbonate, cryogenic cooling contracts the substrate biaxially, the A exciton (differential reflectance) redshifts with a gauge factor
≈ 136 meV/%. No first-principles calculation in the paper. The authors' own proposal: compute E_opt(ε) = E_QP(ε) − E_b(ε) and say which
term carries the shift.

**Thread.** Standalone campaign; pipeline parent = RuCl3 round 8 (RuCl3-study/reproduce_li/round_8, chain driver
the round-8 chain (stage list + templates in ../SHARED_RUNBOOK.md), GPU recipe ../env/setup_gpu_vast.sh). Predecessor literature: Li et al., Adv. Mater. 2022
(DOI 10.1002/adma.202103571) — bulk GW (Yambo, 6×6×6, plasmon pole): direct gap 1.29 eV at Γ, indirect 0.66 eV; +1 % uniaxial tension
shifts the direct gap +90 meV (b, chain axis) and +45 meV (a); measured uniaxial gauge 60–95 (b) and 0–15 (a) meV/%; bulk E_b ≈ 49 meV
(cited). Published A exciton ≈ 1.75 eV (RT), A_2s ≈ 1.89 eV at low T (shortlist entry 05).

**Status vocabulary (terminal, exactly one):** `QP_DRIVEN` | `MIXED` | `BINDING_DRIVEN` | `SLOPE_OFF` | `LITCHECK_FAIL` | `PROCESS_BLOCKED`.

## 1. Why this campaign exists

The simplest prior is additivity of the published uniaxial QP shifts: −(90 + 45) = −135 meV/% for biaxial compression, i.e. the whole
redshift is a quasiparticle-gap effect and E_b (≈ 0.05 eV) cannot move by 136 meV per %. That prior has never been computed, and the
published GW gap (1.29 eV) sits 0.46 eV below the measured exciton (1.75 eV), so the published numbers are not a trustworthy baseline
(round-8 LITCHECK lesson: audit the literature chain before building on it). The campaign computes both derivatives explicitly, on one
consistent chain, with the derivatives (not the absolute energies) as the gated observables. Round-8 lesson applied: absolute GW-BSE
energies carry ±0.2–0.3 eV compensation between gap error and under-binding; spacings and derivatives converge far earlier than absolutes,
so the gates are on derivatives, and the two components are reported separately so compensation cannot hide (offset1 lesson).

## 2. Read this and nothing more

- This file, `acceptance_draft.json` (this folder) (to be frozen as `acceptance.json` + sha256 before S2).
- Structure: ICSD/Materials Project ZrSe3, P2_1/m, 8 atoms (a ≈ 5.41, b ≈ 3.75, c ≈ 9.44 Å, β ≈ 97.5°). Record the source and cell used.
- Pseudopotentials: PseudoDojo v0.4 NC PBE scalar-relativistic Zr, Se (standard); fully-relativistic copies only for the SOC covariate.
- Round-8 driver + GPU recipe (paths above); SHARED_RUNBOOK.md §6 (cloud GPU rules).
- The two papers above (for anchors only; no fitted parameters enter the chain — rule the no-fitted-parameters rule (SHARED_RUNBOOK §1) applies by analogy).

## 3. Definitions

- **Strain protocol.** Biaxial in-plane strain ε applied equally to a and b (isotropic substrate contraction); c, β and all internal
  coordinates relaxed at fixed a, b (free stacking direction, as for a flake on a substrate). PBE-D3 relaxation. Points:
  ε ∈ {−1.0, −0.5, 0, +0.5} %. Two uniaxial control points: +1 % along b only and +1 % along a only (same relaxation rule).
- **E_QP(ε)** = G0W0 direct gap at the k-point of the A exciton's dominant transition (Γ expected; record the actual k). Scalar-relativistic,
  plasmon-pole or static-remainder as in round 8 (declare one, keep it fixed across ε).
- **E_A(ε)** = energy of the first bright exciton in the E ∥ b channel (momentum operator, b1 = b axis) of the BSE on the same chain;
  stick position of the leading state, not an envelope maximum (round-8 lesson: envelope vs member trap).
- **E_b(ε)** = E_QP(ε) − E_A(ε).
- **Slopes** s_opt, s_QP, s_b = linear-fit derivatives dE/dε over the four biaxial points (meV/%), with ε > 0 tensile. The measured
  "redshift under compression" therefore means g_opt = s_opt > 0; the published uniaxial QP shifts (+90, +45 meV/%) are in this convention.
- **QP share** f_QP = s_QP / s_opt; binding share f_b = −s_b / s_opt (so f_QP + f_b = 1).

Parameter count: zero fitted parameters. Covariates declared: k-grid ladder (two levels), SOC on/off at ε = 0 and −1 %.

## 4. Frozen discriminators (numeric, fixed before S2)

| gate | statement | PASS | FAIL |
|---|---|---|---|
| G0 LITCHECK | our scalar G0W0 at ε = 0 vs the published Yambo GW: direct-gap difference recorded; uniaxial QP shifts for +1 % b and +1 % a | both uniaxial QP shifts within ±30 meV of +90 / +45 | either misses by > 30 meV → `LITCHECK_FAIL` (stop, audit the chain; no production) |
| G1 anchor (recorded, not gated) | E_A(0) vs 1.75 eV; E_b(0) vs ≈ 0.05 eV (cited) and vs the A_2s–A_1s spacing 0.1–0.15 eV | — | if \|E_A(0) − 1.75\| > 0.30 eV write `ABS_OFF` in the record and proceed |
| G2 gauge factor | g_opt from the four biaxial points | 109 ≤ g_opt ≤ 218 meV/% (0.8×–1.6× measured; lower bound tight because partial strain transfer can only make the measured value smaller than the intrinsic one) | outside → `SLOPE_OFF` |
| G2b linearity | residual of the linear fit over the four points | max residual ≤ 10 meV | else report the quadratic term and gate G2 on the local slope at ε = −0.5 % |
| G3 decomposition | f_QP from the same four points | f_QP ≥ 0.70 → `QP_DRIVEN`; 0.30 < f_QP < 0.70 → `MIXED`; f_QP ≤ 0.30 → `BINDING_DRIVEN` | — (one of the three is always returned if G2 passes) |
| G4 convergence | coarse-k ladder (level 1 → level 2): change of s_QP and of s_b | each ≤ 15 meV/%; and f_QP class unchanged | else refine once; still changing → `PROCESS_BLOCKED` with the two values reported |
| G5 SOC covariate | fully-relativistic chain at ε = 0 and −1 % | change of the local slope ≤ 20 meV/% → scalar chain stands | larger → SOC chain becomes primary for G2/G3 (rerun the four points) |
| G6 uniaxial exciton check | g_b, g_a for the exciton (not the gap) from the two uniaxial points | g_b within 60–95, g_a within 0–15 meV/% (the measured ranges, 11 samples) | miss → reported next to the verdict, does not block (strain-transfer scatter) |

Non-vacuity check (honesty constraint 10): G3 can return each of its three classes; G2 can fail (e.g. if E_b is strongly strain-dependent in
the opposite sense); G0 can fail (a 30 meV miss on a 90 meV shift is realistic). Record which outcomes were possible before reading the numbers.

## 5. Verdict and stopping rule

| outcome | status |
|---|---|
| G0 fails | `LITCHECK_FAIL` — stop; the deliverable is the audit |
| G2 fails | `SLOPE_OFF` — stop; report g_opt, f_QP, G5, G6 |
| G2 passes, G4 passes | the G3 class (`QP_DRIVEN` / `MIXED` / `BINDING_DRIVEN`) |
| G4 fails after one refinement | `PROCESS_BLOCKED` |

Stop at the verdict. No further strain points, no evGW, no thickness series without a new card.

## 6. Out of scope

Thickness dependence (few-layer vs bulk, strain-transfer efficiency), temperature dependence of the exciton (electron–phonon), substrate
optics, the B exciton and the A_2s fine structure (entry 05 — a separate card that reuses this ground state), evGW, exciton–phonon
linewidths. Shear or uniaxial-c strain.

## 7. Deliverables — exactly these 8

1. `acceptance.json` + sha256 (frozen before S2), with the pre-registered expectation f_QP ≈ 1.0 ± 0.15, g_opt ≈ 135 meV/% written in.
2. `results/strain_table.json`: per point ε → a, b, c, β, relaxed coordinates, E_QP, E_A, E_b, oscillator strength, SOC flag, k-level.
3. `results/slopes.json`: g_opt, s_QP, s_b, f_QP, fit residuals, G4 differences, G5/G6 numbers, gate verdicts.
4. Fig. 1 E_QP, E_A, E_b vs ε with the linear fits and the measured 136 meV/% line (eye test).
5. Fig. 2 ε2 (E ∥ b, E ∥ a) at the four biaxial points, same broadening, stacked, A exciton marked.
6. Fig. 3 uniaxial controls vs the published QP shifts and measured gauge ranges.
7. `RESULT.md` (verdict box, gate table, captions, provenance of every number, compute log with $/h verified), `make_figures.py`, `MANIFEST.sha256`.
8. `vast_run.json` per instance (offer, dph_total read from `--raw`, hours, total $) and the destroy verification line.

## 8. Plan (stages; each stage writes a state file and is resumable)

| stage | what | where | est. |
|---|---|---|---|
| S0 | cell + pseudos + cutoff convergence (ecutwfc ladder to 5 meV on the KS gap); PBE-D3 relax of the reference; local BGW flavor check (HDF5 needed for eigenvectors/plotxct? not needed here) | local, 16 ranks | 1–2 h |
| S1 | six strained relaxations (4 biaxial + 2 uniaxial), KS gaps, bands along Γ–X–Y | local | 2–3 h |
| S2 | freeze `acceptance.json`; write `zrse3_chain.py` from `umid_chain.py` with the SOC/U stages removed (SCF → NSCF → pw2bgw → epsilon → sigma G0W0 → kernel → absorption b1/b2/b3, momentum operator); scorer `strain_score.py` frozen with the card | local | 0.5 day |
| S3 | pilot at ε = 0 on level-1 k (e.g. 4×6×2 coarse, ≈ 300 bands, ε cutoff 10 Ry): time one chain; if ≤ 2 h on the local 32 cores, run everything locally; else rent 1× H100 (verify dph_total ≤ $3/h; MAX_HOURS = est × 1.5 + 1) | local or Vast | 2 h |
| S4 | the seven chains at level 1 (G0 → stop if LITCHECK fails) | local/Vast | 7 × 1–2 h |
| S5 | level-2 k ladder on ε = 0, −1.0 % (G4); SOC covariate at ε = 0, −1.0 % (G5) | same | 4 × 2–4 h |
| S6 | scorer, figures, RESULT.md, MANIFEST; destroy the instance in the same pipeline as the pull + md5 check | local | 0.5 day |

Compute budget: ≈ 15–25 GPU-h ⇒ ≈ $45–75 at ≤ $3/h, or ≈ 2–3 days on the local node. Vast credit was $74.87 at the last check
(2026-10-12) — read `show user` first; abort the Vast route if dph_total × est_hours > 60 % of the credit. Local disk is 21 GB free:
keep WFN/epsmat/bsemat on the instance, pull only eqp tables, absorption files, eigenvalues and logs; record checksums of what was left.

## 9. Honesty constraints

- Every gate is scored by the frozen scorer from files on disk; "N of M planned stages ran" is stated explicitly (lesson:
  EXECUTION_PROTOCOL.md §4).
- Report f_QP with its G4 uncertainty, never as a bare number; report the scalar and SOC values side by side.
- The additivity prior (135 meV/%) is written into the card before any number is computed; agreement with it is not evidence by itself —
  the decomposition is.
- No parameter is tuned to hit 136 meV/%; the strain protocol (what is relaxed) is fixed in §3 and not varied.
- The record states what the measured number cannot constrain (strain transfer, thickness) and does not claim to explain it.
