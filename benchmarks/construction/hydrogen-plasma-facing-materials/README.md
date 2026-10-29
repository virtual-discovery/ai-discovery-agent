# Construction of hydrogen-plasma-facing-materials

Released-start task for the benchmark (see `benchmarks/tools/BUILDING.md`).

- **Source:** Fedrigucci, Marzari and Ricci, Materials Cloud record gpxwh-dx645 (DOI 10.24435/materialscloud:9k-09, CC BY 4.0), accompanying PRX Energy 3, 043002 (2024).
- **Provenance and checksums:** the task's `SOURCE.md`.

## Rebuild

Run compute in `ai-theorist-benchmark-runtime:20261007` (`docker run --cpuset-cpus <set> --user $(id -u):$(id -g) ...`), from this folder unless noted.

1. `python3 fetch.py`
   - Downloads `PFMs_screening_QE_data.tar.xz` (46.9 MB, md5 56aadb3e…), `README.md` and `files_description.md` into `~/.cache/ai-theorist/hydrogen-plasma-facing-materials/`.
   - Checks the md5 values and unpacks the 490 runs (765 MB).
2. `python3 -B extract.py`
   - Writes the cache's `pool.json` (2.5 MB): every run's settings, structure, energies, convergence messages and header.
   - Writes `derived.json` (energies and settings only), kept here.
3. `python3 -B make_task.py`
   - Writes the catalog, the lab data, the private references and `solution/ablation_values.json`.
   - Prints the reference, ablation and robustness values.
4. `cutoff_check/run_cutoff_check.sh` (optional, about 1.5 h on 4 logical CPUs) repeats the construction-time check of isolated-atom cutoffs. Its outcome is in `cutoff_check/results.json`.
5. From the repository root, `.venv/bin/python -m benchmarks.tools.check_task hydrogen-plasma-facing-materials --cpuset <set>` writes `check.json`.

`model.py` holds the shared analysis:

- energies of runs;
- clustering of relaxed H configurations (energy within 5 meV and the eight nearest host distances within 0.05 Å);
- the minimum-image H environment;
- surface energies (the paper's Eq. 11);
- SBEs and cohesive energies.

`solution/solve.py` re-implements the oracle logic stand-alone and reproduces the reference exactly (0.0000 error on every field).

## What the data showed (design inputs)

- **The hydrogen work is internally consistent.** All 21 hosts use 60/480 Ry. Each host and host+H pair shares the cell, k-mesh and pseudopotentials. The H₂ and H references are also at 60/480 Ry; the isolated atom is spin-polarized. ½E(H₂) − E(H) = −2.312 eV per H. The paper's headline Table VI column uses the H atom; the ½H₂ column reproduces our minimum-site values for all 21 hosts.
- **"Inequivalent sites" are random starting points.** The README calls every folder an inequivalent site. The inputs are randomly displaced starting structures, and many relax into the same configuration (energies equal to 10⁻⁷–10⁻⁵ Ry, identical coordination shells):
  - TaN: 6 runs give 1 site.
  - WC: 7 runs give 2 sites.
  - TaC, VB₂ and Mo also collapse.
  - Graphite run 39 hit the BFGS step limit (50) inside the basin of run 17: 0.3 meV higher, same first shell.

  Genuinely different configurations can be close: TiN (Ti₄ site against N–H, 51 meV) and HfB₂ (B-top against B₂Hf₂ bridge, 55 meV).
- **Site competition is not textbook.**
  - hcp Re: the octahedral site beats the tetrahedral one by 0.56 eV.
  - TaC and HfC: H binds to C (C–H 1.15 Å), 0.34 and 0.49 eV below the metal-coordinated void.
  - TiN: the Ti₄ site wins.
- **SBE references have material-specific cutoffs.**
  - The W atom of the W set is at 30/240 Ry, and the W atoms of the WC/W₂C/TaW sets are at 45/360 Ry. They differ by 0.1294 eV, which the construction check reproduces as 0.1289 eV.
  - The paper's Table V W values for WC, W₂C and TaW, and the Ta value for TaN, are reproduced only with the elemental-set atoms (30 Ry W, 45 Ry Ta), so the published compound table mixes cutoffs.
  - The TiB₂ Ti atom is a byte-identical copy of the TiN Ti atom (60/480 Ry; same start time stamp), while all TiB₂ slabs use 35/280 Ry. The construction check gives E_Ti(35) − E_Ti(60) = 0.53 meV, so that mix is numerically harmless.
- **Cohesive energies.** The H-IFE supercells (60/480 Ry) combined with the lower-cutoff SBE atoms shift the cohesive energy by +0.21 eV (W), +0.065 eV (Re), +0.010 eV (Mo) and +0.014 eV (WC). TiB₂ computed fully at 60 Ry (H-IFE supercell, 60 Ry B atom from the BN set) differs from the mixed 35/60 Ry value by 1.4 meV.
- **Terminations and selected facets.**
  - The authors' selected facet is the lowest-surface-energy facet for all 21 materials, and the relaxed-slab and unrelaxed-supercell orderings of terminations agree.
  - WC (100): t2 at 0.251 eV/Å² against t1 at 0.472. TiB₂ (101): t2 at 0.219 against t1 at 0.370.
  - The top atom is the most weakly bound one on these terminations.
- **Data issues outside the graded set.**
  - Several BN, C₃N₄, BC₂N, graphite and WC slab relaxations stopped at the step limit.
  - BN (111)_t2 transforms during relaxation (−4.6 eV); its anomalously low γ selects that termination.
  - The C₃N₄ (100) termination order is undetermined (t3 unconverged and still dropping, 0.08 eV above t7).

  These materials were kept out of the graded SBE set.

## Decision points

| DP | Kind | Wrong path (ablation) |
|---|---|---|
| dp1 | convention (the paper's headline uses the other reservoir) | isolated H atom: all E_f −2.31 eV |
| dp2 | site competition, purchases | textbook site: Re tetrahedral (+0.56 eV), TaC Ta₄ void (+0.34 eV, nearest host Ta) |
| dp3 | README claim against raw outputs; unconverged run | every run a site (WC 7, TaC 4, TiN 3, Mo 2, graphite 3); energy-only 0.1 eV binning merges graphite, HfB₂, TiN |
| dp4 | metadata inconsistency, cascade by reuse | one W atom for all: WC W SBE +0.13, WC cohesive +0.065; reverse: W −0.13 |
| dp5 | discriminating purchase, protocol | other termination (WC W −0.65, C −2.61; TiB₂ B −2.89); paper average (WC W +1.85, TiB₂ B −0.35); site mean |
| dp6 | cascade from H purchases, metadata | cohesive energy with the H-IFE supercell: Re +0.065 |
| dp7 | raw-output reading | starting coordinates: HfB₂ −0.061 Å, Mo +0.047 Å |
| dp8 | competing explanation | cohesive energy as SBE: W −2.09, WC C +1.16, TiB₂ B −0.89 |

The decision points that rest on evidence found in the archive are dp3, dp4, dp6 and dp7. dp2 and dp5 require spending budget on the right runs. dp1 and dp8 are conventions or model choices.

## Build notes and deviations from the design

- **No fresh QE in the lab.**
  - The design wanted a hybrid lab with zero-point energies by finite differences and supercell-size tests. The archived relaxations took 3–12 h on 32–48 cores each, and the runtime image has no pseudopotentials, so a 10-minute evaluation is impossible.
  - The ZPE-corrected E_f, the ZPE site reordering, the Sieverts-law ranking at 1000 K (which needs vibrational entropies and site multiplicities) and the supercell-convergence decision point were dropped. They are rubric items under "what these data cannot support".
- **Nothing withheld.**
  - The design withheld 4–5 materials and asked for predictions; per-host H formation energies cannot be predicted from other hosts to 0.12 eV.
  - Instead every run is purchasable. The graded quantities are derived values that the lab never serves, with a budget of 61 against a catalog of 490 runs.
- **New graded quantities.**
  - The distinct-site energy lists, which test site identity.
  - The preferred-site geometry.
  - Cohesive energies, for the "SBE tracks cohesion" explanation.
  - Per-element minimum SBEs on the lowest-γ termination, instead of the paper's site and termination average.
- **Tolerances.**
  - H E_f is graded at ±0.03 eV, not the design's ±0.12, because the defensible spread is below 1 meV and the nearest wrong configurations are 51–55 meV away.
  - SBE is graded at ±0.06 eV, not ±0.3, because the cutoff trap is 0.13 eV.
- **W excluded from the hydrogen set** (the risk noted in the design). Its SBE and cohesive energy are graded, because the cutoff trap lives there.
- **TiB₂ Ti null.** It was first planned as a required `null` (inconsistent reference). The construction-time QE check showed the cutoff error is 0.5 meV, so the archived combination is sound, and the archive alone cannot show that. The grader accepts the number or `null` for that one entry, and the rubric asks for the assessment.

## Calibration (full table in `design.json` → `calibration`)

| Field | Reference | Robustness max dev | Nearest wrong path | Tolerance |
|---|---|---|---|---|
| E_f (7 hosts) | graphite 1.368, HfB₂ 1.643, Mo 0.509, Re −0.076, TaC 1.774, TiN 2.183, WC 0.543 eV | 0.0001 | second site TiN +0.051, HfB₂ +0.055; textbook +0.34/+0.56; atom −2.31 | ±0.03 eV |
| distinct sites | counts 2, 4, 1, 2, 2, 2, 2 | 0.0002 | counts 3/2/4/3/7 (runs as sites); 1/3/1 (binning) | exact count, ±0.03 eV |
| preferred site | C 1.132, B 1.407, Mo 1.824, Re 1.979, C 1.150, Ti 1.795, W 1.746 Å | 0.0004 | Mo +0.047, HfB₂ −0.061 (starting coords); Re −0.22, TaC Ta | symbol exact, ±0.04 Å |
| SBE | W 11.165; WC W 11.018, C 7.784; TiB₂ B 8.510, Ti 11.489 or null eV | 0.003 (60 Ry C atom) | ±0.129 (W atom reuse); −0.355 (paper average, TiB₂ B) | ±0.06 eV |
| cohesive | Re 8.801, W 9.076, WC 8.943 eV/atom | < 0.0001 | +0.065 (H-IFE supercell Re; W-set atom in WC); −0.129 (WC-set atom in W) | ±0.03 eV |

Rounding every submitted number to 4 significant digits passes.

## Oracle and length

- **Purchases:** 55 of 61 units (90 %).
  - H₂;
  - the pristine supercell and all 25 interstitial runs of the 7 hosts;
  - W: slab supercell, vacancy and atom;
  - WC and TiB₂: both slab supercells, all vacancies of the lower one, and the atoms;
  - unit cells of W, Re and WC, and the Re atom.
- **Runtime:** 22 s in the gate container; the Harbor run took 54 s.
- **Expected agent length:** 45–75 tool calls. The work is parsing about 55 structured runs, clustering configurations, minimum-image geometry, surface energies and bookkeeping of references.

## Results

- **Gate (2026-10-19):** passed. 8 decision points; all 12 ablations fail; the 5 robustness variants and 4-significant-digit rounding pass. 30 hidden values with no leaks; task size 2.6 MB.
- **Harbor oracle:** reward 1.0 (`runs/build-hydrogen-plasma-facing-materials/oracle-r1-hydrogen-plasma-facing-materials/*/verifier/reward.txt`), solver charge 55.
- **Disk:** peak about 1.2 GB in the cache: the tarball, the 735 MB unpacked runs, and about 0.4 MB of cutoff-check inputs plus 408 MB of QE scratch. The raw downloads and scratch were deleted after the Harbor run; `pool.json` (2.5 MB) and `README.md` are kept in the cache for `make_task.py`, and `fetch.py` refetches the rest.

## Open risks

- **Dev-easy risk.** A careful agent that checks the cutoffs of every run it combines, deduplicates relaxed configurations by energy and geometry, and buys every run of each graded host passes. The traps are a 0.13 eV W-atom mismatch, the H-supercell reuse for Re, the README's "inequivalent sites" claim, graphite's unconverged duplicate and Mo's starting coordinates. Strong models may find all of them.
- **Narrow margins.**
  - Mo starting-coordinate distance: +0.047 against a ±0.04 Å tolerance. The HfB₂ entry fails that ablation more clearly (−0.061).
  - Re cohesive with the H-IFE supercell: +0.065 against ±0.03 eV.
- **Paper recall.** The paper tabulates the minimum-site E_f (rounded to 0.01 eV). Recalling them passes the E_f field only; the distinct sites, geometries, minimum SBEs and cohesive energies are not in the paper. The paper's averaged and cutoff-mixed SBEs fail.
- **TiB₂ Ti leniency.** Accepting both the number and `null` makes that entry a check on termination only. The cutoff-mismatch judgement is reviewed through the rubric.
