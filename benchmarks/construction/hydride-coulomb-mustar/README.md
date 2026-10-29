# Construction of hydride-coulomb-mustar

Released-start task for the benchmark. Source: Bozier and Pickard, *Enhanced Coulomb pseudopotentials in high-pressure hydrides from first principles* (manuscript, 2026), Zenodo record 23235290 (DOI 10.5281/zenodo.23235290, CC BY 4.0), `mantel_data_files_v1.tar.gz`. Provenance, checksums and transformations are in the task's `SOURCE.md`.

## Rebuild

All compute runs in `ai-theorist-benchmark-runtime:20261007` through `run_in_image.sh` (cpuset `CPUSET`, default 4,5,20,21; cache mounted at `/data`).

1. `bash fetch.sh [DEST]` streams the 892 MB tarball once (nothing large is stored), checks its md5 and writes and verifies the 459 members in `members.md5` (about 100 MB). Default DEST: `~/.cache/ai-theorist/hydride-coulomb-mustar`.
2. `bash run_in_image.sh extract.py /data` writes `/data/pool.json.gz` (13 MB, cache only).
3. `bash run_in_image.sh make_task.py` (about 15 min) writes the catalog, `lab/data.json` (5 MB), `lab/eliashberg.py` (copied from here), `private/hidden.json`, `hidden_values.json`, `solution/ablation_values.json`, the budget in both `public.json`, and `calibration.json` here. It asserts the solver validation listed below. Solver results are cached in `/data/tc_cache.json`.
4. From the repository root: `.venv/bin/python -m benchmarks.tools.check_task hydride-coulomb-mustar --cpuset <set>` (writes `check.json`).
5. `python3 -B record_oracle.py` (writes `oracle_run.json`), then `python3 -B write_design.py` (writes `design.json`).

## Solver

`eliashberg.py` re-implements the IsoME imaginary-axis equations in numpy, linearized at Tc: constant DOS with a fixed μ\* (`cdos`); full bandwidth with N(ε), energy shift χ and chemical-potential update (`vdos`); and the same with the static W(ε,ε′) (`vdos_W`, φ = φ_ph(iω_n) + φ_c(ε)). Conventions follow IsoME: α²F resampled on a uniform grid, energy grid 1/5/50 meV inside ±100/±500/±5000 meV, Matsubara cutoff 7 eV, DOS file divided by 2 for spin. IsoME.jl (MIT) was read for conventions; no code was copied.

Validation against the archive (all asserted in `make_task.py`, details in `calibration.json`):

- 15 archived full-W runs (Mg₂IrH₆ 0–30 GPa, H₃S 140–240 GPa): within +1.6 % (Mg₂IrH₆ 0–10 GPa +1.2 to +1.6 %; all others −0.2 to +0.8 %). At the withheld pressures that have W (Mg₂IrH₆ 15/25, H₃S 190) the deviation is +0.09 %, −0.20 % and +0.25 %.
- Archived fixed-μ\* runs (vDOS and cDOS, μ\* 0.1–0.2): within 0.5 %.
- H₃S 180 GPa has no W in the archive. With W transferred from 170 or 190 GPa the solver gives 257.7 / 257.3 K against the archived 256.5 K.

The low-pressure Mg₂IrH₆ offset traces to μ: IsoME's printed μ is 0.8–3 % below N(E_F)·W(E_F,E_F) interpolated from `W_ee_raw.npy` (the IsoME version that was used evidently interpolated W differently from the current release).

## What the archive shows (and the task is built on)

- **μ tracks N(E_F).** W(E_F,E_F) is nearly constant with pressure (Mg₂IrH₆ 0.365–0.377 eV; H₃S 0.49–0.51 eV), while μ = N(E_F)W(E_F,E_F) falls from 0.62 to 0.52 in Mg₂IrH₆ (0–30 GPa) and rises from 0.25 to 0.31 in H₃S (140–240 GPa).
- **In the full-bandwidth vDOS solver with a 7 eV cutoff, the bare μ is almost the right fixed μ\*.** μ\* = κμ with κ = 0.981 (Mg₂IrH₆, 20 GPa) and 0.919 (H₃S, 200 GPa) reproduces the full-W Tc. In a constant-DOS calculation with a conventional 1 eV cutoff, μ\* = 0.5875 is needed for Mg₂IrH₆ at 20 GPa; with a 7 eV cutoff, no μ\* reaches the full-W Tc (the constant-DOS Coulomb term saturates).
- **σ matters in `mu_results.dat`.** The smeared N(E_F) at σ = 0.3 eV gives μ ≈ 0.33 for Mg₂IrH₆ (dense-grid "mu rescaled": 0.55); at σ = 0.01 eV W(0,0) blows up (μ rescaled 0.64–0.75).
- **H₃S Tc(P) peaks at 180 GPa** (archived 193.5, 228.5, 256.5, 244.5, 232.5 K at 140–200 GPa). λ is 3.8 and 8.8 at 140 and 160 GPa (harmonic phonon instability; the paper masks these pressures), 3.6 at 180 and 3.0 at 190.
- **Archive inconsistencies:**
  - Folder names contradict Info.txt: Mg₂IrH₆ `vdos_muStar0.1_1` used μ\* = 0.2 (all pressures); H₃S `vdos_muStar0.2_1` used μ\* = 0.1 (140, 160, 180, 220, 240 GPa); CaH₆ `vdos_W_1`–`_4` are constant-DOS fixed-μ\* runs; LaH₁₀ `vdos_W_1` is a vDOS fixed-μ\* run.
  - Different α²F smearing columns across systems and runs (CaH₆ 10; LaH₁₀ 10 for cDOS and 20 for vDOS; others 4). YSbH₆ `vdos_W` used the PBEsol α²F and ω_c = 2 eV (its README).
  - H₃S 200 GPa: `W_ee_raw.npy` (500 points, σ = 0.5 eV) is not the W of its Mantel log (200 points, σ = 0.2 eV) or of its IsoME run (μ 0.311 against 0.295 from the staged W); the solver still lands at +0.8 %.
  - The H₃S cell is the conventional bcc cell (H₆S₂); DOS files are per cell, spin-summed.
  - Info.txt paths of H₃S 140–190 GPa name run-time folders `H3S-40GPa` … `H3S-90GPa`.
  - H₃S 170 GPa has W and a dense DOS but no α²F and no N_F in `mu_results.dat`; H₃S 180 GPa has α²F and IsoME runs but no W.

## Deviations from the design (candidates.json / DETAILS.md section 3)

- **Graded fields.**
  - Tc (full W) at Mg₂IrH₆ 15/25 GPa and H₃S 180/190 GPa is kept, at ±3.5 % instead of ±8 %. ±8 % would accept interpolation across the H₃S peak at 190 GPa (−5.3 %) and Allen–Dynes. The design's H₃S 170/220 GPa were replaced: 170 has no α²F (no Tc), and 220 is trivially interpolable (210/230 give 222.5 against 221.5). 180/190 sit on the non-monotonic peak, and 180 has no W in the archive at all.
  - The equivalent μ\* is kept but defined protocol-specifically: the constant-DOS μ\* at a 1 eV Matsubara cutoff matching the full-W Tc of Mg₂IrH₆ at 20 GPa (0.5875, ±0.025 instead of ±0.03). The vDOS equivalent (0.537) is almost μ itself and would make the field trivial; at a 7 eV cutoff no constant-DOS μ\* exists.
  - "Tc reduction relative to μ\* = 0.16 at 20 GPa" was dropped: it is 42 K, within ±5 K of the abstract's memorable "40 K".
  - Added: N(E_F) of H₃S at 180 GPa per formula unit per spin (tests the two-formula-unit cell).
- **Lab.** There is no W transfer between pressures (the agent must transfer μ, not W), no Morel–Anderson/Allen–Dynes request (the agent can compute them), and no full W matrix for CaH₆/LaH₁₀ (diagonal only). The `eliashberg` request takes up to 6 μ\* values per call.
- **Decision points.** The design's dp1–dp6 are kept (as dp1, dp2+dp3, dp4, dp5, dp7, and the units fork as dp8). The design's dp7 (low-pressure H₃S anharmonicity) is a rubric item, not an ablation: it changes no graded number because no graded pressure is below 170 GPa. Added: dp6 (H₃S Tc peak, interpolation fails) and dp9 (which solver defines the equivalent μ\*).

## Decision points (9)

| dp | Fork | Kind |
|---|---|---|
| dp1 | Bare μ vs μ\* at a cutoff (μ itself or Morel–Anderson rescaling as the cDOS μ\*) | convention / regime |
| dp2 | Smeared vs dense-grid N(E_F) in μ | data-quality, found in the table |
| dp3 | Which σ row of `mu_results.dat` (σ → 0 diverges) | convergence judgement |
| dp4 | Constant vs energy-dependent DOS | model class |
| dp5 | Allen–Dynes with the archived `muc_AD` vs Eliashberg | model class, archive metadata |
| dp6 | Interpolating Tc(P) across the H₃S peak vs computing from each α²F | regime, discovered only by computing |
| dp7 | Conventional μ\* or the 0 GPa μ held fixed vs μ tracking N(E_F) | transfer |
| dp8 | Two formula units per H₃S cell | metadata (scf input) |
| dp9 | Which solver defines the equivalent μ\* | protocol |

dp2, dp3 and dp7 cascade into both Mg₂IrH₆ Tc; dp6 into both H₃S Tc. The route to the withheld Tc must be validated against full-W runs at released pressures, which costs budget.

## Calibration

Full table in `design.json` → `calibration` and `calibration.json`.

| Field | Reference | Oracle | Defensible variants (pass) | Nearest wrong paths (fail) | Tolerance |
|---|---|---|---|---|---|
| Tc Mg₂IrH₆ 15 GPa | 88.45 K | +0.5 % | bare μ (no κ) −0.3 %, nearest-pressure W +1.0 % | 0 GPa μ −3.9 %, σ = 0.01 −9.6 %, smeared +22 %, AD +32 %, cDOS +45 %, μ\* 0.16 +47 % | ±3.5 % |
| Tc Mg₂IrH₆ 25 GPa | 91.65 K | −0.5 % | −1.3 %, −0.2 % | 0 GPa μ −6.9 %, σ = 0.01 −15 %, smeared +20 % | ±3.5 % |
| Tc H₃S 180 GPa | 256.5 K | +0.1 % | −1.4 %, +0.1 % | interpolation −10.1 %, cDOS +6.2 %, μ\* 0.16 +7.4 % | ±3.5 % |
| Tc H₃S 190 GPa | 244.5 K | +0.2 % | −1.4 %, +0.3 % | interpolation −5.3 %, AD +6.3 %, cDOS +8.3 %, μ\* 0.16 +9.1 % | ±3.5 % |
| μ\* (cDOS, 1 eV), Mg₂IrH₆ 20 GPa | 0.5875 | 0.5865 | 0.5885 (matched to the lab's full-W Tc) | bare μ 0.5486, vDOS equivalent 0.5367, Morel–Anderson 0.228 | ±0.025 |
| N(E_F) H₃S 180 GPa (per f.u., per spin) | 0.2704 | 0.2704 | nearest grid point | per cell 0.5409 | ±3 % |

Interpolation passes at the two Mg₂IrH₆ pressures (+0.3 %, −0.9 %); the ablation fails through H₃S. Allen–Dynes at H₃S 180 GPa (+2.7 %) passes; that ablation fails through the other three pressures.

## Oracle

18 evaluations, 30 of 38 units (79 %): full-W runs at Mg₂IrH₆ 10/20/30 and H₃S 200 GPa; the W kernel at H₃S 170 GPa; κ calibration at 20 and 200 GPa (vDOS μ\* lists); `inputs` and one vDOS run at each withheld pressure; one archived run (target Tc at 20 GPa); two cDOS μ\* lists for the equivalent μ\*. The catalogue holds about 115 retrieval units plus any number of solver runs.

- **Gate (2026-10-26, cpuset 4,5,20,21):** passed. 9 decision points, all 11 ablations fail, 4 robustness variants and 4-significant-digit rounding pass, 7 rubric items, 6 hidden values, no leaks, package 5.1 MB, oracle 64 s in the gate container.
- **Harbor oracle (`runs/build-hydride-coulomb-mustar/oracle-r1-hydride-coulomb-mustar`):** reward 1.0, runtime 1 min 30 s.

## Open risks

- **Dev-easy risk.** An agent that calibrates a μ route against `vdos_W` runs at released pressures and computes each withheld pressure from its own α²F passes. The discriminating steps are: noticing that the bare μ (not a rescaled μ\*) is the right vDOS input; transferring μ through N(E_F) rather than holding it fixed; not interpolating across the H₃S peak; and finding μ\* by solving rather than by Morel–Anderson.
- **Thin margins.**
  - The 0 GPa-μ ablation fails at Mg₂IrH₆ 15 GPa by 0.4 points (at 25 GPa by 3.4).
  - Bare μ as the cDOS μ\* misses by 0.039 against ±0.025 (margin 0.014).
  - H₃S 190 GPa interpolation fails by 1.8 points; interpolating the lab's own full-W Tc at 160/200 GPa gives −4.7 % (1.2 points).
  - A vDOS route with μ held at the H₃S high-pressure value (0.30) gives −3.4 % at 180 GPa and passes narrowly; it is not listed as an ablation.
- **Reference mix.** The Tc references are archived IsoME values (coarse ±0.5 K for H₃S); the μ\* reference is a validated-solver value. The solver is up to 1.6 % high against IsoME at Mg₂IrH₆ 0–10 GPa (not graded pressures).
- **Not peer reviewed.** The manuscript was released in October 2026.
- **Memorisation.** The abstract's μ = 0.62, μ\* ≈ 0.39 and 40 K (0 GPa) are not graded; recalling 0.39 for the cDOS μ\* fails.

## Disk

- **Peak:** about 1.9 GB in `~/.cache/ai-theorist/hydride-coulomb-mustar`: the full exploratory extraction (all of the tarball except the 3.2 GB of `W_ee.dat`), three `W_ee.dat` files for the equivalence check, the 13 MB pool and the solver cache. The 892 MB tarball itself was only streamed, never stored.
- **After the Harbor run:** the cache and the gate's lab image were deleted. `fetch.sh` with `members.md5` (verified against the extraction before deletion) refetches the 100 MB actually used.
