# Construction of gan-irradiation-defect-levels

Released-start task for the benchmark (see `benchmarks/tools/BUILDING.md`).

- **Source:** Hauck, Jin & Tuttle, Zenodo record 21893631 (DOI 10.5281/zenodo.21893631, CC BY 4.0), supporting J. Appl. Phys. 139, 245702 (2026).
- **Provenance and checksums:** the task's `SOURCE.md` and `members.sha256` here.

## Rebuild

Run every step in `ai-theorist-benchmark-runtime:20261007` (`docker run --cpuset-cpus <set> --user $(id -u):$(id -g) ...`), from this folder unless noted.

1. `python3 fetch.py --pseudo`
   - Downloads `gan_defects.zip` (25.1 MB, md5 85d4a96b…) into `~/.cache/ai-theorist/gan-irradiation-defect-levels/` and unpacks the 67 calculations (115 MB).
   - Fetches the two pslibrary potentials used for the dielectric tensors.
2. `python3 -B extract.py` writes `pool.json` (1.6 MB): settings, energies, structures, eigenvalues, core potentials and site moments of every run.
3. `dfpt/run_dfpt.sh` (mount the potentials at `/pseudo`; about 1 h on 2 shared physical cores) runs pw.x, ph.x and dynmat.x (with `lperm`) for wurtzite GaN. `dfpt/out/` keeps the inputs and outputs.
4. `python3 -B make_dielectric.py` writes `dielectric.json`.
5. `python3 -B make_task.py`
   - Writes the catalog, the lab data (`data.json`, 1.1 MB), the private references, `solution/ablation_values.json` and `calibration.json`.
   - Sets the budget in both `public.json` files.
6. `python3 -B write_design.py` writes the decision points, graded fields and calibration into the task's `design.json`.
7. From the repository root, `.venv/bin/python -m benchmarks.tools.check_task gan-irradiation-defect-levels --cpuset <set>` writes `check.json`.

Code layout:

- `model.py`: site mapping, the anisotropic Ewald point-charge model, the eFNV correction and the level finder.
- `reference.py`: the reference protocol and every variant: band-edge sources, host energy, dielectric tensor, correction scheme, alignment region and sign, site choice, charge states.
- `solution/solve.py`: re-implements the protocol stand-alone from purchased records. It reproduces the reference to 1e-7 eV.

## What the archive shows (design inputs)

- **k-points.** Every defect cell and `perfect`/`perfect/kpt_spec` use one special k-point, (½,½,½). `perfect/kpt_gam` is a Γ-only relaxation (4.708 eV higher in total energy). At (½,½,½) the supercell HOMO/LUMO are 1.2855/7.8185 eV (gap 6.53 eV); at Γ they are 1.9425/5.5080 eV (gap 3.5655 eV). The Fermi energy of `perfect` (1.6165 eV, Gaussian smearing) lies inside the gap.
- **PAW sets.** `perfect/WAV_ref` (Γ, on the `perfect` geometry) and `perfect/kpt_gam/d10` use Ga_d (864 electrons). Their VBM (4.27 eV) is on another eigenvalue scale, and their gap (3.51 eV) looks plausible. All other runs use Ga with 3d in the core (384 electrons).
- **VBM alignment.** The core potentials of `perfect` and `kpt_gam` differ by +0.004 eV (Ga) and −0.105 eV (N). The reference uses the Γ eigenvalue unaligned. All-site alignment (+0.05 eV in every level) is a robustness variant; N-only alignment (+0.105 eV) is recorded in the calibration, but is neither an ablation nor a robustness variant.
- **Sign of the core potentials.** Far from localized charged defects, V_def − V_perf anticorrelates with the point-charge potential (r ≈ −0.85), so they are electron potential energies. Neutral defects shift the far potential by −0.12 (V_Ga), −0.06 (V_N) and −0.18 eV (V_Ga–V_N). The alignment therefore moves every level of a defect by roughly the same amount.
- **Delocalization diagnostics.** Ga_N +5/+6 and V_Ga–V_N +4 have weak potential correlation (r = −0.15 to −0.31), alignment terms of −0.85 to −1.37 eV, and spin spread over many N atoms. N_i(i2) −2/−3 occupy states above the Γ CBM. None of them is a ground state inside the gap.
- **Spin.** NUPDOWN is fixed (0/1). Several even-electron states (V_Ga −1/+1, Ga_N 0/+2, N_i +1, V_Ga–V_N +2/+4) are broken-symmetry singlets.
- **Sites.** Ga_i: i1 and i2 relax to the same structure for q ≤ 0 and q = +3. i2 is lower by 1.035 eV (q = +1) and 0.984 eV (q = +2). N_i, i1 − i2 in eV: −2.26 (−3), +0.10 (−2), +2.94 (−1), +3.70 (0), −0.07 (+1), +1.86 (+2), +1.70 (+3).
- **U values.** V_N: U = −0.067 eV, a robust (3+/+). Ga_i: U = −0.001 eV, effectively zero. The grader therefore merges stability windows narrower than 0.1 eV, and the rubric asks for that judgement.
- **Vertical data.** Only Ga_N −2 starts from another state's relaxed geometry (Ga_N −1). The divacancy runs are single steps on pre-relaxed geometries.
- **Dielectric (builders' DFPT, PBE).** ε∞ = diag(5.962, 5.962, 6.091) and ε0 = diag(10.075, 10.075, 11.052). The point-charge term is 0.189 eV × q². Total |q| = 3 corrections are 1.38–2.01 eV.

## Reference values

- **Gap:** 3.5655 eV.
- **V_Ga:** (+/0) 0.980, (0/−) 1.624, (−/2−) 2.093, (2−/3−) 2.535 eV.
- **V_N:** (3+/+) 0.657, (+/0) 3.156, (0/−) 3.276 eV. The last is not graded (within 0.3 eV of the CBM).
- **Ga_i:** (3+/+) 2.294 eV.
- **N_i:** (2+/+) 0.659, (+/0) 2.102, (0/−) 2.920 eV. (3+/2+) at −0.04 eV lies outside the gap.
- **C_N:** (+/0) 0.084 (not graded), (0/−) 0.905 eV.
- **Divacancy binding:** E_b = 3.165 eV (V_Ga 3−, V_N +, divacancy 2−).
- **Frenkel pairs:** N 9.304 eV (V_N +, N_i −), Ga 13.319 eV (V_Ga 3−, Ga_i +).

## Decision points

| DP | Kind | Wrong path (ablation) |
|---|---|---|
| dp1 | k-point metadata, folded eigenvalues | R-point band edges (gap 6.53, levels +0.66); Fermi energy as VBM (+0.33) |
| dp2 | host consistency, cascade | Γ-only host energy: E_b −4.71, Frenkel −9.42 eV |
| dp3 | screening model | ε∞ (2+ states appear, E_b +0.55, FP_Ga +0.92); no correction |
| dp4 | alignment | point charge only: V_Ga −0.10 to −0.13, FP_Ga −0.41, FP_N +0.15 |
| dp5 | sign convention in raw output | potentials read as electrostatic: E_b −0.59, FP_Ga −2.10, V_N neutral state lost |
| dp6 | purchases, site switching | i1 only: Ga_i (3+/+) +0.52 plus spurious (+/0); FP +1.04 / +1.11 eV |
| dp7 | stability (negative U) | adjacent (q/q−1) levels for V_N and Ga_i (inverted order) |
| dp8 | PAW-set metadata | Ga_d band edges: VBM +2.33 eV, gap −0.05 |

dp1, dp2, dp5, dp6 and dp8 rest on evidence found in the metadata or raw outputs, or bought under the budget. dp3, dp4 and dp7 are methodological.

## Deviations from the design (candidates.json / DETAILS.md §35)

- **k-point.** The design said "Γ-only". The defect runs use (½,½,½), so the band-edge decision point became R versus Γ plus the Ga_d trap (dp1, dp8), and the host-energy consistency was added (dp2).
- **Dielectric tensors.** HSE DFPT was not possible, so the builders computed PBE DFPT values at the archive's structure and say so. The prompt requires using the served tensors.
- **Targets.**
  - The N Frenkel "pinning level" was replaced by the N and Ga Frenkel-pair energies at E_CBM − 0.5 eV. A pair's net charge changes sign at a constituent's level, so the pinning level would duplicate a CTL.
  - The binding energy is kept.
  - The band gap was added.
  - Levels are graded as full ladders inside [0.3 eV, gap − 0.3 eV], with the charge pairs required.
  - Ga_i: the data give one (3+/+) level with U ≈ 0, not "(3+/2+)/(2+/+)", hence the 0.1 eV merge rule.
  - C_N: only (0/−) is graded, because (+/0) lies 0.08 eV above the VBM.
- **Dropped decision points.**
  - Perturbed-host states: no in-gap effect with this sampling; now a rubric item.
  - Spin states: NUPDOWN is fixed.
  - Thermodynamic versus optical levels: no vertical data.
  - "Charge states for pair energies": contrary to the stated convention, so not an inferential fork.
- **Budget.** 43 of 68 units (63 %), slightly above "half the catalog". About 22 units are irrelevant to the targets (antisites, extra perfect runs, positive divacancy states).

## Calibration (full table in `design.json` → `calibration`)

| Field | Reference | Robustness max dev | Nearest wrong path | Tolerance |
|---|---|---|---|---|
| gap | 3.5655 | 0 | Ga_d −0.052 | ±0.03 |
| levels | above | +0.051 (all-site VBM), +0.049 (5.5 Å region) | point charge only −0.125 (V_Ga +/0), N-site VBM +0.105 (not tested) | ±0.12, same pairs |
| E_b | 3.165 | +0.079 (5.5 Å) | point charge only +0.095 (passes E_b, fails elsewhere); ε∞ +0.55; wrong sign −0.59 | ±0.15 |
| FP_N | 9.304 | +0.020 | point charge only +0.149, ε∞ +0.18 | ±0.12 |
| FP_Ga | 13.319 | +0.101 (all-site VBM), +0.096 (5.5 Å) | N-site VBM +0.21 (passes), point charge only −0.405 | ±0.25 |

- **Robustness variants:** all-site VBM alignment, isotropic ε0, alignment sites beyond 5.5 Å or 4.0 Å. All pass.
- **Ablations:** all 10 fail.
- **Rounding:** 4 significant digits passes.

## Oracle and length

- **Purchases:** 38 of 43 units.
  - perfect, kpt_gam and the dielectric tensors;
  - V_Ga +1…−3 and V_N +3…−1;
  - Ga_i +3…0 at both sites;
  - N_i +3…−2 at both sites;
  - C_N 0 and −1;
  - divacancy −1…−3.
- **Accuracy:** the oracle reproduces the reference to 1e-7 eV.
- **Runtime:** 12 s in the gate container.
- **Expected agent length:** 50–80 tool calls.

## Results

- **Gate:** passed. 8 decision points, 10 ablations failing, 4 robustness variants and rounding passing, 17 hidden values with no leaks, 1.13 MB.
- **Harbor oracle:** reward 1.0 (`runs/build-gan-irradiation-defect-levels/oracle-r1-gan-irradiation-defect-levels`).

## Open risks

- **Dev-easy risk.** A careful agent that reads the k-point and PAW metadata, codes eFNV with the right sign and buys both interstitial sites passes.
- **Recall.** C_N (0/−) ≈ 0.9 eV is close to a well-known HSE literature value.
- **Narrow margins.**
  - The point-charge-only V_Ga (+/0) misses the level tolerance by only 0.005 eV, but that path fails robustly on FP_Ga.
  - N-only VBM alignment (+0.105 in levels) would pass the levels and fail nothing else.
- **Reference is protocol-specific.** eFNV with the builders' PBE tensor; the paper's levels may differ.

## Disk

- **Peak:** about 145 MB in the cache: the zip, the 115 MB unpacked runs and the potentials.
- **Deleted:** the raw downloads, after the gate and the Harbor run.
- **Kept:** `pool.json`, `dielectric.json`, `dfpt/out` and the 3.5 MB of potentials (2.4 MB cache).
