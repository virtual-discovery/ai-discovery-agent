# Unresolved cases — experimental findings with no first-principles explanation, packaged for an executing agent

Three campaigns selected on 2026-10-28 from an 8-entry shortlist of 2024–2026 experimental papers ("First-Principles Research Opportunities",
22 Oct 2026). Selection rule, taken from how the α-RuCl3 round-8 campaign succeeded: **start from the closest published first-principles model
("the Li"), reproduce it to a declared tolerance, then add one hypothesis per rung, each with a gate that can fail.** The per-entry Li,
additions and difficulty for all eight entries are in [LI_ANALYSIS.md](LI_ANALYSIS.md).

| # | case | source | the Li | addition | difficulty | est. compute |
|---|---|---|---|---|---|---|
| 01 | [ZrSe3 strain-driven exciton redshift](01_zrse3_strain/goal.md) — which of E_QP(ε), E_b(ε) carries 136 meV/% | Li et al., Nano Res. 2026 | Li et al., Adv. Mater. 2022 bulk GW-BSE (uniaxial +90/+45 meV/%) | biaxial strain axis + E_b(ε); derivative gates | LOW | 15–25 GPU-h |
| 02 | [TiS3 1.28 eV peak for E ⊥ chains](02_tis3_interchain/goal.md) — interchain exciton or not | Boldyrev et al., JETP Lett. 2024 | Molina-Mendoza et al. 2015 bulk GW-BSE (E_QP 1.15, E_opt 1.05, E_b 0.10 eV, second feature ~1.2–1.3 eV) | polarization channels, state identity, real-space interchain weight; monolayer as control | LOW–MOD | 15–25 GPU-h |
| 03 | [ReS2 stacking-dependent exciton dipole](03_res2_dipole/goal.md) — AB 13 meV Stark shift, AA none | Paul et al., Nano Lett. 2026 | Zhou et al., ACS Photonics 2021 G0W0-BSE AA/AB few-layer | hypothesis: layer-asymmetry exciton dipole d_S = ⟨z_e − z_h⟩, AA null by symmetry | MOD | 40–80 GPU-h |

Run order 01 → 02 → 03. Case 03's S0–S2 and its symmetry gate cost nothing and may be done at any time. Successor (not packaged):
ZrSe3 A_2s/B fine structure (Alcazar Ruano 2025) reuses case 01's ground state.

## Folder layout

```
unresolved cases/
  README.md                 this page
  LI_ANALYSIS.md            closest model, additions, difficulty for all 8 shortlist entries; the three-pick decision
  EXECUTION_PROTOCOL.md     how an agent picks up a card, runs it, reports, audits coverage, hands back
  SHARED_RUNBOOK.md         the stripped round-8 chain: rules, software, stage list, trap list, real-space rules, GPU rules, record format
  templates/                QE and BerkeleyGW input templates (scalar, no U/SOC) with placeholders and comments
  env/                      GPU build script + arch file, fit_scissors.py, real-space pipeline (run_realspace.sh, pick_states.py)
  0N_<slug>/PI_GOAL.md      the PI goal statement (verbatim, authoritative)
  0N_<slug>/goal.md         the card (house style: status vocabulary, frozen discriminators, stopping rule, deliverables, honesty constraints)
  0N_<slug>/acceptance_draft.json   freeze + sha256 before stage S2 (fill TO_FILL fields first)
  0N_<slug>/results/        scorer inputs and the record (RESULT.md, MANIFEST.sha256, figures) — written by the executing agent
```

## What the three cards share

- Pipeline: QE 7.4.1 + BerkeleyGW 4.0 (G0W0, RPA screening, BSE with the momentum operator), scalar-relativistic, no U; SOC only as a
  declared covariate. Same 8-atom P2_1/m structure type for 01 and 02 (one chain driver), slabs with truncation for 03.
- Gates (from round 8): LITCHECK against the Li before production; one global shift with width-scaled tolerances where an absolute energy
  is compared; components reported separately (no compensation); k/window/vacuum ladders as explicit gates; a negative control that makes
  every character measure non-vacuous (edge exciton for 02, AA stacking for 03, uniaxial published shifts for 01).
- Terminal status vocabulary per card, "stop at the verdict", and "N of M planned stages ran" in every record.

## Provenance

Written 2026-10-28 from the RuCl3 round-8 record (RuCl3-study/reproduce_li/round_8, repro package of 2026-10-21) and the shortlist PDF.
Inputs still needed from the papers (paywalled on the authoring machine) are marked `TO_FILL` in each `acceptance_draft.json`: case 02 needs
the per-polarization edge energies, the 1.28 eV linewidth and the measured face from Boldyrev et al.; case 03 needs the device field,
thickness and the exciton label from Paul et al.
