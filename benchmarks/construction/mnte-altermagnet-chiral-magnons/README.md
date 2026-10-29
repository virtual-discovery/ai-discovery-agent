# Construction of mnte-altermagnet-chiral-magnons

This is a released-start task for the benchmark (see `benchmarks/tools/BUILDING.md`).

- **Source:** Liu et al., *Observation of Switchable Chiral Magnons in an Altermagnet*, PRL 136, 236705 (2026), arXiv:2605.14124.
- **Data:** Zenodo 20049187 (`Figures.zip`, CC BY 4.0).
- **Provenance and what is hidden:** the task's `SOURCE.md`.

## Rebuild

All compute runs in `ai-theorist-benchmark-runtime:20261007`. Mount this folder with the repository's `benchmarks/` directory, so that `make_task.py` finds `../../tasks/<id>`, and the cache at `/cache`.

1. `python3 -B fetch.py CACHE` downloads `Figures.zip` (5.3 MB), checks the md5 and unpacks 60 members to `CACHE/x`.
2. `MNTE_CACHE=/cache python3 -B calib.py` writes `calibration.json`, in about 1 minute. It contains:
   - the cut-peak fits;
   - the IN20 truth;
   - the exchange-difference estimates;
   - the IN20 predictions of every variant and wrong path.
3. `MNTE_CACHE=/cache python3 -B make_task.py` writes:
   - `catalog.json` and `data.json`;
   - `private/hidden.json` and `hidden_values.json`;
   - `solution/ablation_values.json`;
   - both `public.json` files.
4. Optional: `MNTE_CACHE=/cache python3 -B verify_model.py` runs two checks:
   - the spin-wave code against the archived dispersion: within 0.014 meV;
   - the kinematics against the archived HYSPEC trajectory: within 1e-8 meV.
5. From the repository root, `.venv/bin/python -m benchmarks.tools.check_task mnte-altermagnet-chiral-magnons --cpuset 2,3,18,19` writes `check.json`.

Modules:

- `ill.py`: ILL TAS parser.
- `reduce.py`: IN20 binning and two-Gaussian fits.
- `sym.py`: P6₃/mmc operations and Mn–Mn bond classes.
- `lswt.py`: two-sublattice BdG/Colpa LSWT with single-ion anisotropy.
- `forward.py`: HYSPEC kinematics and neutron weights.

`solution/solve.py` is a standalone re-implementation.

## What the data showed (design inputs)

- **Archive contents.**
  - HYSPEC: the authors' reduced symmetric and chiral slices; the cuts g–j; the per-channel slices and cuts (FigS4: their mean equals Fig2c, and their difference equals Fig2e); and the PND at (1 0 −1).
  - IN20: raw Cryopad files. These are the PND files of the two field-cooled states and 47 constant-energy scans at 32 meV, in the I_x^{+-} and I_x^{-+} channels.
  - Model files: `Fig1c.mat` (dispersion and trajectory), `Fig2df.mat` (model slices) and `Fig2ghij.mat` (model cuts).
- **Field-cooled states.**
  - Red file 043731 is state A, −4.2 mT, scans `Fig3d`, measured 14–16 Sep.
  - Blue file 043806 is state B, +4.2 mT, scans `Fig3c`, measured 16–18 Sep.
  - Paper Fig. 3 confirms both assignments.
  - Flipper convention: I9 is the incident flipper and I10 the scattered one. This is verified because state B gives I^{-+} − I^{+-} < 0 at Ω₂, as in paper Fig. 3(e).
- **Lattice in the IN20 headers.**
  - The headers use AS = 4.17, BS = 4.15, CS = 6.87 Å. MnTe has a = 4.148 and c = 6.711 Å.
  - The nominal scans along (h, 0, 1.25) are therefore at l = 1.221. The paper quotes (h, 0, −1.22) in a frame with l reversed.
  - Predicting at the nominal l shifts both 32 meV crossings by about −0.025 in h.
- **Peaks and chirality in the withheld scans.**
  - The peaks are at QH = 1.3022 and 1.3960 (both states summed, binned at 0.005).
  - Chiral amplitudes: state A is (+, −) at 10 and 14 σ; state B is (−, +) at 9 and 12 σ.
  - The chiral fraction is 6–7%.
- **PND asymmetries (I^{++} − I^{−−})/(I^{++} + I^{−−}).**
  - y channel: A +0.089, B −0.098.
  - z channel: A −0.098, B +0.101. The z channel measures the b projection of n̄, and its sign is opposite to y.
  - HYSPEC: I^{+0} − I^{−0} < 0 at (1 0 −1).
  - The empirical rule "y-NMI negative ⇔ the upper branch has a negative chiral term" holds on both crystals. It is consistent with theory: both N(1 0 ±1) = ±√3 b_Te and the g-wave form factor are odd in l.
- **Bond classes.**
  - The 9.831 Å same-sublattice shell (√3a in plane, ±c) is the first that splits into two classes exchanged between the sublattices.
  - The 7.18 Å in-plane √3a shell cannot split: each pair contains ±r.
  - The 7.89 Å shell is a single class.
- **Splitting form factor.**
  - It is ∝ sin 2πl [sin 2π(2h+k) + sin 2π(k−h) − sin 2π(h+2k)], with its maximum at (⅓, 0, ¼).
  - The global maximum of |ω₂ − ω₁| is instead at Γ: the 4.05 meV easy-plane gap (dp2).
- **Prior inconsistency.**
  - The 2024 paper prints J10 + J11 = 0.0931 meV from the nodal-point fit, but its J10 = 0.0681 and J11 = −0.0221 sum to 0.046.
  - The latter reproduces the archived dispersion (`Fig1c.mat`) to 0.014 meV.
  - The catalog stages the consistent average, 0.0230 per bond.
- **HYSPEC energies.**
  - The archive has no per-pixel energy. The trajectory `Etraj(h)` along (h, 0, −1.22) is exactly quadratic.
  - With ω = 32° (paper), it fixes E_i = 37.000 meV and a reduction lattice of a = 4.100 Å, c = 6.843 Å.
  - The resulting grid matches the Fig. 2(g) energy axis.
- **Known model–data offset (not graded).** In the steep region of cut g, the simple crossing model sits about 0.01 in u above the data and the authors' model map. The authors' map is presumably resolution-convolved in Q; that cannot be reproduced without their resolution model. Cut h (near the band top) agrees to 0.004. The offset enters the exchange difference only through peak separations.

## Deviations from the design (candidates.json / DETAILS.md section 19)

- **The lab serves archived data, not an LSWT/polarization calculator.** A calculator would not bind the budget, because the agent can run LSWT locally.
- **Withheld data.** All IN20 inelastic scans are withheld, not "one domain state". The agent predicts them from the HYSPEC data plus the PND of both states. If one state's scans were sold, the peak positions and the reversal would be trivial.
- **Graded targets.**
  - (1) The energies of the chiral modes at withheld Q become the IN20 peak positions in instrument QH: the mean within 0.012 and the separation within 35%.
  - (2) The chirality sign of the reversed domain is graded for both states and both peaks.
  - (3) The direction of maximal splitting becomes a wave vector, scored by the symmetry form factor (≥ 0.97).
  - (4) The exchange difference uses the window [0.055, 0.12] meV instead of ±30% of a single value, because no archived value exists and the defensible estimates span 0.068–0.097.
- **Not graded.**
  - The domain-population magnitude: the archive's own PND and PINS estimates differ by up to 40% (9.7% vs 6.2% on IN20).
  - Resolution convolution: no resolution parameters are archived. It is a rubric item through the cut-g offset.
- **The prior Hamiltonian is staged.** Without J1–J3, D and the 9.83 Å average the data cannot fix the absolute band position. The "which shell splits" fork is therefore weakened: the shell is named. It survives as the class partition, the form factor (dp7) and rubric item 1.

## Decision points (8)

| dp | Fork | Ablation result |
|---|---|---|
| dp1 | Shell-averaged exchange (J_a = J_b) | ΔJ 0, IN20 separation 0.0045 |
| dp2 | Anisotropy gap taken as the maximal splitting | q = Γ, score 0 |
| dp3 | IN20 nominal l = 1.25 (lattice in the headers) | mean 1.3216 (1.3286 with h rescaled) vs 1.3491 ± 0.012 |
| dp4 | HYSPEC chirality pattern copied to both IN20 states (PND ignored) | state A wrong |
| dp5 | Pattern flipped for l > 0 without flipping the NMI sign | all four signs wrong |
| dp6 | z-channel PND used for the a* projection | all four signs wrong |
| dp7 | Maximal splitting taken on the measured line (0.375, 0, 0.22) | score 0.913 < 0.97 |
| dp8 | Ordered-pair exchange convention | ΔJ 0.0375, separation 0.0383 |

## Calibration (full table in `design.json → calibration`)

| Field | Truth / reference | Defensible spread (pass) | Nearest wrong path (fail) | Tolerance |
|---|---|---|---|---|
| ΔJ (meV) | oracle 0.0749 (no archived value) | 0.0677 (cut h) – 0.0821 (cut g); published 0.0902 | 0.0375 (ordered pairs); 0 | [0.055, 0.12] |
| IN20 peak mean (QH) | 1.3491 (withheld scans) | 1.3447 – 1.3523 | 1.3286 (nominal l, h rescaled); 1.3216 | ± 0.012 |
| IN20 separation (QH) | 0.0939 | 0.0673 (−28%) – 0.0874 | 0.0383 (−59%); 0.0045 | ± 35% |
| Chirality | A (+, −), B (−, +), 9–14 σ | oracle exact | copy pattern; flips | exact |
| q_max score | 1.0 at (⅓, 0, ¼) | (0, ⅓, ¼), (⅔, 0, ¾) | 0.913 measured line; 0 at Γ | ≥ 0.97 |

Tightest margins:

- **IN20 mean:** the nearest wrong path (nominal l with h rescaled) fails by 0.0085.
- **IN20 separation:** the cut-h-only estimate passes by 7 points. All HYSPEC-based estimates undershoot the withheld IN20 separation by 7–28%, a real tension between the two crystals and instruments.
- **ΔJ:** the cut-h estimate is 0.013 above the lower bound.
- **q:** the measured-line answer fails by 0.057.

Rounding to 4 significant digits passes.

## Oracle

The oracle buys cuts g, h, i and j and the three PND items: 7 evaluations, 7 of 8 units (88%), 10 s of container time.

Method:

1. Two-Gaussian fits to cuts g and h.
2. Crossing separations along the cut centre lines, with the catalog energy grid, interpolated over a ΔJ grid; the result is the mean of the two cuts.
3. The CS value read from the header gives l = 1.2211; the 32 meV crossings follow from that.
4. Chirality: the HYSPEC upper-branch sign from cut i, multiplied by the ratio of the HYSPEC and IN20 y-NMI signs.
5. A form-factor maximisation over the zone.

Results:

- **Gate (2026-10-26, cpuset 2,3,18,19):** passed. 8 decision points; all 9 ablations fail; 5 robustness variants pass; rounding passes; 9 rubric items; 5 hidden values, no leaks; 0.6 MB package.
- **Harbor oracle:** reward 1.0 (`runs/build-mnte-altermagnet-chiral-magnons/oracle-r1-*`, 58 s).

## Open risks

- **Few evaluations (7) and a small catalog (11 units).** The length comes from local analysis: LSWT, bond classes, kinematics, raw ILL parsing and the sign chain. The budget forces only a choice between the slices and the cuts.
- **Memorisation.** The published ΔJ (0.090) passes the ΔJ window, and an agent that recalls it predicts the IN20 separation well. The IN20 mean still needs the header lattice correction, and the chirality still needs the PND sign chain. The PRL (June 2026) states the field-sign/chirality relation in prose, but the states are labelled A/B without field signs.
- **Chirality is categorical.** Once an agent sees that A and B have opposite NMI, it has a 50% chance of guessing both states right. The joint check (all four signs) and the pattern within each state limit this to one bit.
- **Model–data offset in cut g** (above). An agent that refits the prior isotropic exchanges to absorb it could move the IN20 mean by a few 0.001. This is within the tolerance in our tests, which did not refit.

## Disk

- **Peak:** about 30 MB. That is 21 MB of cache (zip and unpacked members), 2 MB of scratch (paper HTML/PDF and figures) and under 1 MB of gate and Harbor run directories.
- **After the build:** the raw downloads were deleted. `fetch.py` and the md5 regenerate them.

## Fix round after the dev pilots (2026-10-26, cpuset 14,15,30,31)

**What the pilots did.** In `dev-mntemag-pilot1`, both models submitted almost the same failing answers, and both got all four chirality signs and the q_max check right:

| Model | \|J_a − J_b\| (meV) | IN20 peaks (QH) |
|---|---|---|
| Sonnet 5 | 0.014 | 1.3095 / 1.3241 |
| GPT-5.6 Terra | 0.015 | 1.309 / 1.325 |

**Shared errors.** Both models fell into two errors at once:

1. **Nominal l (dp3).** Both kept the nominal l = 1.25 from the scan headers. This moves the peak mean to about 1.317.
2. **Splitting underestimated by about a factor of 5.**
   - Sonnet read the chiral cut i as the difference of two fully chiral branches, M ≈ −(δu/2) dS/du. The measured chiral term is diluted by the domain imbalance (6–15% of S), so δu comes out at 0.004–0.006 instead of the resolved peak separation of 0.032.
   - Terra ranked ridge residuals over the whole symmetric slice instead of resolving the two peaks in the cuts.

**Reproduction.** `calib.py` reproduces the Sonnet path from the staged data (`calibration.json → dipole`): ΔJ = 0.011 meV, and IN20 peaks at 1.3108 / 1.3227 with the nominal l.

**Why this is not an ambiguity.** Paper Fig. 2(g) labels the two cut-g peaks, at u = 0.842 and 0.875, as Ω₂ and Ω₁. Neither the instruction nor the staged inputs invite the dipole reading: the prior and the convention are unambiguous, and both models used the pair-once convention.

**Changes:**

- **New ablation dp9** (`chiral_dipole_full_chirality`, with nominal and true l). Both variants fail the exchange and separation checks.
- **Slices sold as tiles.** Each slice is now sold as four tiles (u, v quadrants) at 1 unit each. The catalog grows from 11 to 15 units; the budget stays at 8.
- **Oracle unchanged:** 7 evaluations, 7 units.

**Why the budget cannot bind much harder.** The archive is small: two slices, four cuts and three PND scans. Adding more evaluations would require inventing purchasable items that the archive does not contain.

**Results:**

- **Gate:** passed. 9 decision points; all 11 ablations fail; robustness and rounding pass.
- **Harbor oracle:** reward 1.0 (`runs/build-mnte-altermagnet-chiral-magnons-fix`).
