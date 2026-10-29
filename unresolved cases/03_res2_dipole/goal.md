# RES2-DIPOLE-1 — Does AB stacking give ReS2 excitons a permanent out-of-plane dipole that AA stacking forbids, and is it the 13 meV Stark shift?

**Source finding.** Paul et al., Nano Letters 26, 8482–8487 (2026), DOI 10.1021/acs.nanolett.6c01130 (shortlist entry 08): back-gated
few-layer (≈ 6–8 layers) ReS2 FETs; AB-stacked samples show a ≈ 13 meV exciton Stark shift (linear + quadratic contributions) in the
off-state, AA-stacked samples a negligible shift. The paper interprets this as a stacking-dependent exciton dipole, with electrostatic
estimates and fits; no first-principles calculation (main text not verified — provisional).

**The Li (closest published model).** Zhou et al., ACS Photonics 8, 2021 (DOI 10.1021/acsphotonics.0c01225): G0W0-BSE optical absorption
of AA- and AB-stacked few-layer ReS2 along and perpendicular to the Re chain (PL Exc. I/II/III: AB 1.48/1.52/1.57 eV, AA 1.46/1.50/1.56 eV).
Supporting experiment: van der Laan et al., ACS Photonics 10, 2023 (bulk AA vs AB excitons; AB broader, weaker, interpreted as lower
binding). The Li already gives the exciton ladder per stacking and its polarization. It does not give the exciton's charge distribution
across layers, which is what the Stark shift measures.

**Thread.** Standalone campaign; pipeline parent = the shared runbook (../SHARED_RUNBOOK.md); real-space pipeline ../env/run_realspace.sh.
Sibling cards: 01 (ZrSe3), 02 (TiS3). This is the only card with a genuine physics hypothesis added on top of the Li.

**Status vocabulary (terminal, exactly one):** `DIPOLE_CONFIRMED` | `DIPOLE_TOO_SMALL` | `NO_DIPOLE_BY_SYMMETRY` | `LITCHECK_FAIL` |
`PROCESS_BLOCKED`.

## 1. Why this campaign exists

Hypothesis (one rung, one gate): in AA stacking the bilayer/slab has a symmetry operation exchanging the layers (inversion or a mirror
through the interlayer plane), so every exciton has ⟨z_e − z_h⟩ = 0 exactly and only a quadratic (polarizability) Stark shift is possible.
AB stacking (adjacent layers shifted in-plane by the Zhou/van der Laan definition) removes that operation, so an exciton can have a
permanent out-of-plane dipole p_z = e⟨z_e − z_h⟩ and a LINEAR Stark shift ΔE = −p_z F_in. The linear coefficient needs no finite-field
calculation: ⟨z_e − z_h⟩ follows from the BSE eigenvectors (layer-resolved weights of the electron and hole), i.e. from the round-8
real-space machinery. The quadratic part is an optional rung. Round-8 lessons applied: a symmetry-protected null (AA) makes the AB number
non-vacuous; the measured 13 meV is converted with the device field, which is an input from the paper, never a fitted quantity; the exciton
identity (which of Exc. I/II carries the dipole, and its polarization) is reported with the dipole, since the experiment resolves the
exciton by polarization.

## 2. Read this and nothing more

- This file, `acceptance_draft.json` (freeze + sha256 before S2, after the S0 inputs are filled).
- Structures: ReS2 distorted-1T (P-1), 12 atoms per layer; AA and AB stackings as defined in Zhou 2021 / van der Laan 2023 (interlayer
  shift); S0 reproduces both from those papers' structural data (or Qiao et al. 2016, DOI 10.1039/C6NR04212A, for the stacking geometries)
  and records which definition is used. Bilayer and 4-layer slabs, vacuum ≥ 15 Å, slab truncation. PBE-D3 relaxation of positions at
  fixed in-plane lattice (experimental a, b, γ).
- Pseudopotentials: PseudoDojo v0.4 NC PBE scalar-relativistic Re, S (standard). SOC: Re is heavy; the Li (Zhou 2021) states its SOC
  choice — S0 records it; a fully-relativistic bilayer AB chain at one k-level is the declared covariate (G6).
- S0 inputs from the paper (user/PI supplies if paywalled): gate dielectric and thickness, gate-voltage range → external field F_ext;
  sample thickness; the exciton (I or II) whose shift is 13 meV; the linear and quadratic parts if separated.
- ../SHARED_RUNBOOK.md, ../EXECUTION_PROTOCOL.md.

## 3. Definitions

- **Layer weights.** For exciton S with BSE amplitude A^S_{vck}: electron layer weight w_e(l) = Σ_{vck} |A|² ρ_{ck}(l), hole weight
  w_h(l) = Σ |A|² ρ_{vk}(l), where ρ_{nk}(l) is the fraction of |ψ_nk|² in layer l (z-slabs between the interlayer mid-planes). Both from
  eigenvectors.h5 + the WFN (same gauge: keep them together). Cross-check on one state with plotxct (hole at the Re site of layer 1,
  electron density integrated per layer).
- **Dipole length.** d_S = Σ_l z_l [w_e(l) − w_h(l)] (Å), with z_l the layer-centre heights; p_z = e d_S.
- **Internal field.** F_in = F_ext / ε_⊥, with ε_⊥ the slab's out-of-plane dielectric constant computed at the DFT level (finite sawtooth
  field, S1) — not taken from the literature. F_ext from the S0 device inputs.
- **Predicted linear shift.** ΔE_lin = e d_S F_in at the paper's maximal field.
- **Exciton identity.** Exc. I / II per stacking identified as in the Li (energy order and polarization along / across the Re chain).
- **Symmetry test.** The AA slab's space group must contain an operation mapping layer 1 ↔ layer 2 (checked with spglib on the relaxed
  structure, tolerance 1e-3 Å); the AB slab must not. If both or neither do, see G1.

Parameter count: zero fitted parameters. Covariates: k ladder (two levels), vacuum/truncation (two levels), SOC (one AB bilayer chain),
thickness (bilayer vs 4-layer AB).

## 4. Frozen discriminators

| gate | statement | PASS | FAIL |
|---|---|---|---|
| G0 LITCHECK | bilayer AA and AB G0W0-BSE vs Zhou 2021: Exc. I and II energies, their polarization, and the AB−AA sign of Exc. I | energies within 0.15 eV, polarization and sign reproduced | else `LITCHECK_FAIL` (stop, audit) |
| G1 symmetry | spglib on the relaxed AA and AB slabs | AA has a layer-exchange operation, AB does not → hypothesis testable | AB also has one → `NO_DIPOLE_BY_SYMMETRY` (the hypothesis is dead at the structural level; report which operation survives); AA lacks one → `PROCESS_BLOCKED` (wrong AA structure) |
| G2 null | \|d_S\| for Exc. I and II in AA | ≤ 0.05 Å (numerical zero) | larger → `PROCESS_BLOCKED` (gauge/layer-assignment defect; do not proceed to G3) |
| G3 dipole exists | \|d_S\| for the measured exciton in AB (bilayer) | ≥ 0.20 Å | smaller → `DIPOLE_TOO_SMALL` |
| G4 magnitude | ΔE_lin at the paper's maximal field vs the linear part of the 13 meV (or the full 13 meV if not separated) | within a factor 2 → `DIPOLE_CONFIRMED` | outside → `DIPOLE_TOO_SMALL` (or `_TOO_LARGE`, reported with d_S and F_in) |
| G5 convergence | k ladder and vacuum ladder: change of d_S | ≤ 0.10 Å each, G2 class unchanged | else one refinement, then `PROCESS_BLOCKED` |
| G6 SOC covariate | fully-relativistic AB bilayer, one k-level | change of d_S ≤ 0.10 Å → scalar chain stands | larger → SOC chain becomes primary for G3/G4 |
| G7 thickness (recorded, not gated) | 4-layer AB: d_S of the lowest two bright excitons, interior vs surface | — | reported next to the verdict; the 6–8-layer sample is not reproduced |
| G8 quadratic part (optional) | two-field sawtooth DFT → GW-BSE at ±F for AA: polarizability α_S from the energy shifts | compared with the measured quadratic part | optional rung, only if G0–G5 are done and budget remains |

Non-vacuity: G0, G1, G2, G3, G4 can each fail; G2 is a true null (zero by symmetry) and G3 is the positive; the possible outcomes are
recorded before any eigenvector is read.

## 5. Verdict and stopping rule

| outcome | status |
|---|---|
| G0 fails | `LITCHECK_FAIL` |
| G1 AB symmetric | `NO_DIPOLE_BY_SYMMETRY` |
| G2 fails or G5 fails after one refinement | `PROCESS_BLOCKED` |
| G3 fails, or G4 outside factor 2 | `DIPOLE_TOO_SMALL` / `DIPOLE_TOO_LARGE` with numbers |
| G3 and G4 pass | `DIPOLE_CONFIRMED` |

Stop at the verdict. No finite-field BSE beyond G8, no moiré/twist, no 6–8-layer slabs without a new card.

## 6. Out of scope

The 6–8-layer sample itself; defect- or trap-mediated fields; exciton–phonon effects; the broadening difference AA vs AB (van der Laan);
moiré or twisted stackings; interlayer charge-transfer excitons beyond those appearing in the BSE window; the ReS2 excitonic Raman (entry 01).

## 7. Deliverables — exactly these 9

1. `acceptance.json` + sha256 (S0 inputs filled: field, thickness, exciton label, stacking definition, SOC choice of the Li).
2. `results/structures.json`: AA/AB bilayer and 4-layer AB cells, relaxed positions, spglib symmetry report (G1).
3. `results/spectra.json`: along/across-chain ε2 and sticks per stacking and thickness, both k-levels, shift vs the Li.
4. `results/dipoles.json`: w_e(l), w_h(l), d_S for Exc. I/II per stacking/thickness/level/SOC; ε_⊥; F_in; ΔE_lin; gate verdicts.
5. Fig. 1 eye test: AA and AB bilayer spectra along/across the chain vs the Li's spectra and the measured exciton energies.
6. Fig. 2 layer-resolved electron and hole weights for Exc. I/II in AA and AB (bar charts) + one plotxct side view per stacking.
7. Fig. 3 predicted ΔE_lin vs field for AB (line) with the measured 13 meV point and the AA null.
8. `RESULT.md` (verdict box, gate table, captions, provenance, "N of M stages ran"), `make_figures.py`, `MANIFEST.sha256`.
9. `vast_run.json` per instance + destroy verification line.

## 8. Plan

| stage | what | where | est. |
|---|---|---|---|
| S0 | structures (both stackings, from the Li's definition), pseudos, cutoff ladder, device inputs from the paper, Li's SOC choice; PBE-D3 relax at fixed in-plane cell | local CPU | 0.5–1 day |
| S1 | DFT: bands AA/AB bilayer, spglib (G1), ε_⊥ by sawtooth field (two fields, linear fit) | local | 0.5 day |
| S2 | freeze card; `res2_chain.py` from the shared stage list with slab truncation; scorer `dipole_score.py` (reads dipoles.json only) + `layer_weights.py` (eigenvectors.h5 + WFN → w_e, w_h, d_S) validated on a synthetic two-layer test (hand-built amplitudes) | local | 1 day |
| S3 | pilot AA bilayer, level 1 (coarse 6×6×1, ~600 bands, 10 Ry, truncation); time it → GPU plan | 1–2× H100 | 4–8 GPU-h |
| S4 | AB bilayer level 1; G0 on both; G2 (AA null) then G3/G4 (AB) | GPU | 4–8 GPU-h |
| S5 | level-2 k and vacuum ladder on AB (G5); SOC AB bilayer (G6) | GPU | 10–20 GPU-h |
| S6 | 4-layer AB (G7), 48 atoms: one level only | GPU | 15–30 GPU-h |
| S7 | optional G8 | GPU | 10 GPU-h |
| S8 | scorer, figures, RESULT.md, MANIFEST; destroy in the pull+md5 pipeline | local | 0.5 day |

Compute budget: ≈ 40–80 GPU-h ⇒ ≈ $120–250 at a verified ≤ $3/h; the largest of the three cards. Run after 01 and 02; S0–S2 and the
G1 symmetry test cost nothing and may be done first (a `NO_DIPOLE_BY_SYMMETRY` outcome would end the card before any GW).

## 9. Honesty constraints

- d_S is reported with its G5 uncertainty and next to the AA null from the same code path; a nonzero AA value is a bug, not a result.
- F_in uses the computed ε_⊥ and the paper's field; both numbers and their sources are printed with ΔE_lin. No quantity is adjusted to
  land on 13 meV.
- The exciton carrying the dipole is named with its polarization and energy; if the measured exciton cannot be identified, G4 is reported
  for both Exc. I and II and the verdict is `PROCESS_BLOCKED` on identity, not `DIPOLE_CONFIRMED`.
- Thickness (G7) is recorded as a limitation in the first paragraph of RESULT.md.
- Coverage is stated as "N of M planned stages ran"; a reassembled record is audited stage by stage against the files on disk.
