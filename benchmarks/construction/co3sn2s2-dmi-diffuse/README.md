# Construction of co3sn2s2-dmi-diffuse

This is a released-start task for the benchmark (see `benchmarks/tools/BUILDING.md`).

- **Source:** D. W. Tam et al., *Noncollinear kagome magnet driven by giant antisymmetric exchange*, Nat. Commun. (2026), DOI 10.1038/s41467-026-77000-6.
- **Data:** Materials Cloud record 0cwzy-zz096 (DOI 10.24435/materialscloud:xw-6f), CC BY 4.0.
- **Provenance, checksums and what is hidden:** the task's `SOURCE.md`.

## Rebuild

From this folder. All compute runs in `ai-theorist-benchmark-runtime:20261007` (numpy 1.21 / scipy 1.8, the versions the gate grades with).

1. `python3 -B fetch.py [CACHE]` downloads the six record files (18.9 MB) through the per-file API links and checks every md5. It then extracts the `.npy` members of both zips into `CACHE/x/`. The default CACHE is `~/.cache/ai-theorist/co3sn2s2-dmi-diffuse`.
2. Run `make_task.py` in the image. It takes about 1 minute:
   ```sh
   docker run --rm --cpuset-cpus <set> --user $(id -u):$(id -g) -e HOME=/tmp -e PYTHONDONTWRITEBYTECODE=1 \
       -e CO3SN2S2_CACHE=/cache -v <CACHE>:/cache:ro -v <repo>:/repo \
       -w /repo/benchmarks/construction/co3sn2s2-dmi-diffuse ai-theorist-benchmark-runtime:20261007 python3 -B make_task.py
   ```
   It writes:
   - `catalog.json` and `archive/` (the README, `readcombined.py` and the SpinW script, verbatim);
   - the lab `data.json` (5.8 MB);
   - `private/hidden.json` and `hidden_values.json`;
   - `solution/ablation_values.json`;
   - both `public.json` files;
   - `calibration.json` here.

   It also asserts that:
   - every served crop reproduces the full-map profile;
   - the DM matrix is compatible with the symmetry of its bond;
   - the LSWT gap equals the linearised Landau–Lifshitz k = 0 frequency.
3. From the repository root, `.venv/bin/python -m benchmarks.tools.check_task co3sn2s2-dmi-diffuse --cpuset <set>` writes `check.json`.
4. `python3 -B write_design.py` fills `design.json → calibration` from `calibration.json` and the gate's oracle state.

Modules:

- `dms.py`: array loading, crops, and the graded annular profile.
- `predict.py`: the reference prediction recipe, the robustness recipes and the wrong-path recipes.
- `hamiltonian.py`: R-3m operations, SpinW-style bond numbering, DM propagation, the classical ground state, LSWT and the Landau–Lifshitz check.

`solution/solve.py` is a stand-alone re-implementation; it imports nothing from here.

## What the data showed (design inputs)

- **Scattering geometry.**
  - The maps are in the (H 0 L) plane, with Qx = h·a* and Qy = l·c* in 1/Å. Bragg peaks were indexed at (003), (101), (−102), (10−2) and (104).
  - The four DMCpy channels are counts, monitor, vanadium normalisation and normcount. `readcombined.py` defines the normalized intensity.
- **What the diffuse scattering is.**
  - It is critical scattering around (003), peaked at T_C. This is the paper's own DMS analysis (SI: Lorentzian κ(T), critical exponent γ ≈ 0.6, (101) used as a background reference).
  - In the graded profiles, the excess at the first annulus is about 0 at 140 K, 0.25 at 170 K and 0.85 at 175–185 K, falling slowly to about 0.2 at 250 K. The 165 K wide scan shows no excess at (003).
  - The onset below T_C is therefore steep and asymmetric. The narrow scans alone leave the 140–178 K gap open, and the 165 K wide region closes it (dp1, dp2).
- **Polarisation and structure factors.**
  - (101) has the same |Q| as (003): 1.4365 versus 1.4359 1/Å, so the form factor is the same.
  - Per kagome triangle, ferromagnetic correlations have |F|² = 9 at (003) and 1 at (101). An in-plane q = 0 120° component has |F|² = 0 at (003) and 4 at (101).
  - At (003) only in-plane spin components scatter. The (003) critical signal therefore implies critical in-plane FM fluctuations.
  - The 180 K (101) map shows no critical peak, consistent with 1/9 of the (003) term and incompatible with AFM correlations (dp4, dp5).
- **Background.**
  - It rises smoothly with T at every reflection. For example, at (101) it is about 8.1 (1.5 K), 8.3 (165 K), 8.4 (180 K) and 8.7 (300 K), in 1e-6 units.
  - A (101) reference taken at 1.5 K or at 300 K fails (dp3).
  - For the (003) predictions, a T-independent background is absorbed by the joint fit, so it is not an ablation.
- **Bragg intensities are not usable for magnetism.** Integrated peaks change non-monotonically, even between 180 and 300 K, above T_C. For example, (−102) changes 5424 → 2448 and (104) changes 15177 → 7808. These peaks are served, as a budget trap and as a rubric item, but they are not graded.
- **Archive defects found:**
  - `3870_Qys.csv` is a byte-identical copy of `3870_Qxs.csv`.
  - `readcombined.py` writes lower edge + 0.015, which is the upper edge of the 0.015 bin, not its centre.
  - The component runs 3870, 3872 and 3873 are pixel-for-pixel contained in the combined 3871 and 3874 files (99.5% of shared pixels identical), so averaging them double-counts.
  - The narrow scans come on three different grids: offsets of 0, −0.287 and −8.56 pixels, and 175 K has a different shape.
  - The notebook's `makeNorm` maps 0/0 to 0, so uncovered pixels would enter means as zeros.

  None of these is graded. The CSVs and component runs are not served. The profile definition uses integer pixel offsets and pooled exposure, so the grid and centre conventions cannot change a graded number.
- **The spin Hamiltonian.**
  - The SpinW script has J1 = −56 and a J1dmi matrix whose DM vector is (0, −15, −6) meV, and Jc1 = −3 on bond 3. Bond 3 is the 4.6409 Å interlayer bond; bond 2 is the 4.6393 Å in-plane shell, 1.6 mÅ away and separated by SpinW's 1e-3 Å tolerance. S = 1.
  - (0, −15, −6) is Moriya-compatible only for bonds along a. On the R-3m stabiliser of a nearest-neighbour bond, D must be ⊥ to the bond, with in-plane and z components allowed.
  - With the matrix on Co(0, ½, ½) → Co(½, ½, ½), the classical ground state is a k = 0 umbrella: tilt 8.661° from c with a 120° in-plane pattern, energy −121.957 meV/Co. LSWT gives a zone-centre gap of 4.1505 meV, identical to the linearised Landau–Lifshitz frequency.
  - The model reproduces the IN8 magnon positions plotted by the script (for example H at 32 meV versus L, from 0.177 at L = −6 to 0.089 at L = −7.2) to about 0.02 r.l.u. This confirms the once-per-bond, S = 1 reading.
  - The measured canting is 17° (peer-review file) and the measured gap is 2–3 meV (Liu et al.); neither is reproduced exactly by the model.

## Deviations from the design (candidates.json / DETAILS.md section 15)

- **No diffuse-scattering simulator; the lab serves archived arrays.**
  - The archived DMS is critical scattering near T_C. It does not determine D/J1: no archived observable fixes D independently of the gap model.
  - A classical MC/SCGA lab tuned to these maps would therefore produce a graded D/J1 that is not uniquely recoverable.
  - The design's target "D/J1 within 30%" was replaced by two consequences of the authors' archived Hamiltonian: the classical canting angle and the zone-centre LSWT gap. D/J1 itself is read off the script and is a rubric item.
- **The "temperature window of the in-plane AFM component" does not exist in this source.** The authors find one canted state at all T < T_C. The AFM question is instead graded through the withheld (101) region at 180 K. The structure factors make this a discriminating prediction: AFM correlations give chi2r 2182, FM correlations about 2.
- **Withheld narrow scans: 170, 175 and 200 K.** The design suggested 175, 185 and 200 K. 170 K replaced 185 K because the steep sub-T_C onset is where data choice and model class matter; 185 K is interpolable.
- **Generator bond fixed in the instruction.**
  - SpinW's internal choice of generator bond could not be verified without MATLAB. The other orientation gives 7.75° and 3.31 meV.
  - The instruction therefore names the bond, Co(0, ½, ½) → Co(½, ½, ½), while leaving the symmetry propagation to the agent.
- **Decision points dropped:**
  - *Q-integration/resolution:* the profile definition makes it unique.
  - *Energy scale from T_C (mean field):* no T_C calculation is graded; dp9 (spin length) carries the energy-scale fork.
  - *Easy-axis anisotropy:* there is no single-ion term in the authors' model. It is covered by dp7 and rubric item 3.
- **Extra wide-scan positions** (half-integer L and a centring-forbidden one) let an agent test interlayer AFM stacking. They are real data and are not graded. They enlarge the catalog to 47 units against a budget of 16.

## Decision points (10)

| dp | Fork | Ablation result |
|---|---|---|
| dp1 | Narrow scans only for the sub-T_C predictions (no 165 K wide region) | 170 K chi2r 19.6 (global fit), 18.2 (independent fits), 13.3 (per-annulus interpolation); limit 12 |
| dp2 | Critical scattering mirrored about T_C | 170 K 20.4, 200 K 45.3 |
| dp3 | (101) background from 1.5 K or 300 K | 59.8 / 136; limit 15 |
| dp4 | (101) diffuse equal to (003), or c-axis polarisation only | 409 / 317 |
| dp5 | In-plane AFM (120°) correlations near T_C | 2182 |
| dp6 | DM matrix copied to every bond (no symmetry operation) | tilt 77.0°, gap 0 |
| dp7 | Out-of-plane DM only / in-plane DM only | tilt 0, gap 0 / tilt 8.18°, gap 3.70 meV (−10.9%) |
| dp8 | Matrix on the reversed (or other) +a bond | tilt 7.75°, gap 3.31 meV |
| dp9 | S = 0.15 from the ordered moment | gap 0.62 meV |
| dp10 | Bonds counted twice | gap 8.30 meV |

Where the difficulty comes from:

- dp1 is a budget and discriminating-data decision: the relevant 165 K data sit in the wide scans, not the narrow ones.
- dp2, dp3 and dp5 are regime and competing-explanation judgements that are visible only in purchased data.
- dp4 is a transfer to a withheld observable.
- dp6–dp10 are discovered from the script's SpinW semantics and cascade from the ground state into the gap.

## Calibration (full table in `design.json → calibration`)

| Field | Reference | Defensible spread (pass) | Nearest wrong path (fail) | Tolerance |
|---|---|---|---|---|
| 170 K chi2r | 4.92 | 3.58–10.68 (more data, per-annulus interpolation, independent fits) | 13.28 narrow-only interpolation; 18.2–19.6 narrow-only fits; 20.4 mirror | ≤ 12 |
| 175 K chi2r | 1.25 | 2.15–6.47 | (dp1/dp2 fail at 170 or 200 K) | ≤ 12 |
| 200 K chi2r | 1.75 | 2.19–3.18 | 45.3 mirror | ≤ 12 |
| (101) 180 K chi2r | 2.34 | 1.48–3.87 (165 K-only background, XY ratio, no critical term, unsmoothed) | 59.8 (1.5 K background); 136; 317; 409; 2182 | ≤ 15 |
| canting angle (°) | 8.661 | 8.66 (rounded); 9.12 if the interlayer bond is lost (passes, not an ablation) | 7.749 reversed bond; 0; 77.0; measured 17 | ±0.5 |
| gap (meV) | 4.1505 (Landau–Lifshitz 4.1505) | 4.151 in the 9-site cell; 4.384 without the interlayer bond (+5.6%, passes) | 3.697 (−10.9%); 3.312; 0.62; 8.30; 0 | ±6% |

**Tightest margins:**

- **170 K:** the independent-fit recipe (10.68) passes by 1.3, and narrow-only per-annulus interpolation (13.28) fails by 1.3. The noise in the visible data limits this separation.
- **Canting angle:** the reversed bond fails by 0.41°.
- **Gap:** the no-interlayer variant passes by 0.4 points. It is not an ablation.

Rounding to 4 significant digits passes.

## Oracle

The oracle makes 9 evaluations for 14 of 16 units (88%):

- narrow scans at 140, 178, 185, 192 and 220 K (10 units);
- wide (003) regions at 165 and 180 K;
- wide (101) regions at 165 and 300 K.

Its method:

- **(003) profiles:** a joint fit to every (003) profile, y_T(q) = b0 + b1(T − 200) + s·q + A_T/(1 + (q/κ_T)²). A_T and κ_T are interpolated between bracketing temperatures.
- **(101) prediction:** the background (linear in q, interpolated between 165 and 300 K) plus 1/9 of the (003) critical term.
- **Hamiltonian:** the SpinW bond numbering, R J Rᵀ propagation, 12 BFGS starts for the 9-site ground state, and LSWT at k = 0 in the primitive cell.

Scores: chi2r 4.92, 1.25, 1.75 and 2.34; tilt error −5e-6°; gap error −6e-8.

Results:

- **Gate (2026-10-20, cpuset 4,5,20,21):** passed. 10 decision points; all 15 ablations fail; 8 robustness variants and 4-significant-digit rounding pass; 8 rubric items; 34 hidden values; no leaks; 5.6 MB package; oracle 23 s of container time.
- **Harbor oracle:** reward 1.0 (`runs/build-co3sn2s2-dmi-diffuse/oracle-r1-co3sn2s2-dmi-diffuse`, 55 s, charge 14).

## Open risks

- **Evaluations below target.** The oracle makes 9 evaluations, below the 15–60 target. Each purchase is a whole map, and the length comes from local analysis. An agent exploring Bragg peaks or the extra positions spends more of the budget.
- **Noise-limited profile checks.** The 170 K margin is about 1.3 in chi2r on both sides. An unusual but defensible recipe could exceed 12. A strongly smoothed or parametric model is safer, and the threshold sits between the measured spreads.
- **Model check versus SpinW.** The canting and gap are graded for the model as the instruction defines it. If SpinW's own generator were the reversed bond, the authors' numbers would be 7.75° and 3.31 meV. The instruction removes this ambiguity explicitly.
- **Memorisation.** The paper's verdict (one canted state) is in the instruction's framing. Its measured 17° fails the canting check, and the measured 2–3 meV gap fails the gap check.
- **Dev-easy risk.** An agent that buys the 165 K (003) region, models the critical term jointly, applies kagome structure factors and implements SpinW's symmetry propagation and LSWT will pass. The LSWT step (canted local frames, bond orientation) is the likeliest failure point.

## Disk

- **Peak:** about 130 MB. That is 95 MB in `~/.cache/ai-theorist/co3sn2s2-dmi-diffuse` (19 MB raw, 77 MB unpacked), 5.5 MB of scratch (the SI and peer-review PDFs) and about 13 MB of gate and Harbor run directories. A refetch test briefly added another 96 MB.
- **After the build:** the raw downloads and the gate's lab image were deleted. `fetch.py` and the pinned md5s regenerate them.
