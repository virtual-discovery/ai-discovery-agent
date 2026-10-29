# Construction of li3ps4-superionic-transport

Released-start task built from Gigli, Tisi, Grasselli and Ceriotti, Chem. Mater. 36, 1482 (2024), using Materials Cloud record pkn2d-h8190 (CC BY 4.0).

## Rebuild

1. Run `bash fetch.sh [DEST]`. It downloads `MaterialsCloudArchive-ChemMat2024.zip` (171 MB, md5 e838e8716fb42200c34f06119337e909), extracts the 47 members listed in `members.md5` and verifies them. The default DEST is `~/.cache/ai-theorist/li3ps4-superionic-transport`.
2. Run `python3 extract.py [DEST]`. It writes `pool.json` (118 KB): the archived arrays, keyed by nominal temperature, with time series decimated as described in SOURCE.md.
3. Run `python3 -B make_task.py`. It writes:
   - `environment/agent/catalog.json`;
   - `environment/lab/data.json`;
   - `private/hidden.json` and `hidden_values.json`;
   - `solution/ablation_values.json`;
   - the budget in both `public.json` files.

   It also prints the calibration numbers. It asserts that, for each functional, the lower heat-capacity peak, the largest jump of the γ-started Nernst–Einstein conductivity and the step of the PS4 angle fluctuation fall in the same temperature interval.
4. Run `python3 -B write_design.py`. It writes `design.json`, including the calibration table and the oracle errors from `oracle_errors.json`.
5. Run `.venv/bin/python -m benchmarks.tools.check_task li3ps4-superionic-transport --cpuset <cpus>`. It writes `check.json`.

## What the task is

The agent buys archived per-run analyses of the ML-MD simulations:

- conductivities: Nernst–Einstein for 1 unit, Green–Kubo for 2 units;
- PS4/P sublattice dynamics;
- heat-capacity intervals;
- the Li MSD of the PBE0 γ-started run at 750 K;
- the σ(τ) breakdown in either frame;
- Li pair correlations;
- order parameters.

The budget is 60 units; the catalogue is worth about 780 units.

The graded quantities are:

- the PBE0 and r2SCAN transition temperatures;
- σ of the three withheld runs, each the first simulated temperature above its reference's transition;
- the Haven ratio of the withheld PBE0 run;
- the PBE0 Green–Kubo activation energies below the transition and from the transition to 1000 K;
- the dynamical Li–S fraction of σ;
- the Li-hop to PS4-flip rate ratio at 800 K;
- the paddle-wheel verdict.

## Deviations from the design (candidates.json / DETAILS.md section 44)

- **No trajectories.** The archive holds figure-level analyses and MD inputs, not trajectories, so nothing could be preprocessed from trajectories. The lab re-serves the archived analyses per run. The only time series are the 750 K Li MSD (decimated to 1,015 points), the σ(τ) breakdown and the order parameters.
- **Heat capacity served for all functionals instead of withheld.** The difficulty is a confounder instead: for every reference the largest Cp peak is the melting of the PS4 sublattice, which coincides with the onset of P diffusion. The superionic transition is the lower peak.
- **Transition tolerance ±30 K (design ±25 K).** The grid spacing is 25 K. The tolerance accepts either bracketing temperature or the midpoint, and still rejects the experimental α transition of 746 K for r2SCAN (41.5 K off), the melting peaks, and transferring the PBEsol transition.
- **Activation energies graded for σ (Green–Kubo) rather than tracer diffusion.** The paper tabulates Nernst–Einstein values (PBE0: 1.43 and 0.269 eV), and the memorised 0.269 eV fails the superionic band. Below the transition the band is wide (±0.35 eV), because the resolved window spans only 600–725 K and the fit spread is 1.42–1.76 eV.
- **Withheld conductivities graded at ±0.15 dex (design ±0.25 dex).** Reporting σ_NE as σ is only 0.26–0.29 dex off, so ±0.25 dex would let that wrong path pass. Defensible predictions fall within ±0.11 dex.
- **Haven ratio graded at the withheld PBE0 γ 750 K run, ±0.12 (design ±0.15).** The α-started value at the same temperature (0.367) must fail.
- **Rate ratio graded at 800 K (design 600 K).** At 600 K the PBE0 γ-started runs are below their transition: there are no flips, D_ω is librational, and the ratio is ill-posed. At 800 K the system is superionic and below melting. Tolerance ±1 dex.
- **Added Li–S cross fraction (frame decision, dp7).** In the barycentric frame this term equals the kinematic recoil of the matrix: 0.275 from the data versus 0.2727 from the masses alone.
- **β-started runs not served above 625 K.** At 650–1150 K the archived r2SCAN-β and PBE0-β series coincide to 0.1 % (z ≈ 0.01), so at least one of them is mislabelled. The catalogue says so; nothing graded depends on β.

## Decision points (9) and how difficulty arises

| dp | Fork | Kind |
|---|---|---|
| dp1 | Nernst–Einstein vs Green–Kubo, with a cascade into the superionic E_a and the Haven ratio | convention |
| dp2 | Green–Kubo noise floor (about 0.03 S/m, apparent H_R > 1) below about 600 K | regime judgement, found in the data |
| dp3 | Largest Cp peak is melting, not the transition; cascades into the withheld predictions | confounder; discriminating purchase (D_P, angle step) |
| dp4 | Transition depends on the reference; the budget must be split over three references | found in the data under budget |
| dp5 | Each withheld run sits just above its jump; interpolating across it fails | discriminating purchase at that temperature; cascade |
| dp6 | α-started runs remain a different, more conductive state (lower H_R) | hysteresis found in the data |
| dp7 | Barycentric Li–S term is recoil, not coupling | confounder |
| dp8 | Librations (fast) vs flips (slow) for the paddle-wheel test | discriminating observable |
| dp9 | One Arrhenius law through the transition | model class |

## Calibration

Values are reference / robustness / nearest wrong path / tolerance; the full table is in `design.json` under `calibration`.

| Field | Reference | Defensible spread | Nearest wrong path | Tolerance |
|---|---|---|---|---|
| Tc PBE0, r2SCAN (K) | 737.5, 787.5 | 725–750, 775–800 | 746 for r2SCAN (−41.5); melting 875/925; PBEsol 612.5 | ±30 K |
| σ PBEsol/γ/625 (S/m) | 145.4 | +0.036 to +0.110 dex | α branch +0.278; tracer −0.289; interpolation −0.503 | ±0.15 dex |
| σ PBE0/γ/750 | 232.3 | −0.017 to +0.048 dex | α branch +0.247; tracer −0.275; interpolation −0.273 | ±0.15 dex |
| σ r2SCAN/γ/800 | 242.2 | −0.016 to +0.059 dex | α branch +0.307; tracer −0.262; interpolation −0.515 | ±0.15 dex |
| H_R PBE0/γ/750 | 0.530 | 0.529–0.637 | α branch 0.367; 1.0 | ±0.12 |
| E_a below (eV) | 1.610 | 1.42–1.76 (oracle 1.47) | single Arrhenius 0.646; noise floor 0.160 | ±0.35 |
| E_a superionic (eV) | 0.184 | 0.166–0.230 | tracer 0.277 | ±0.06 |
| Li–S fraction | 0.050 | 0–0.057 | barycentric 0.275 | ±0.10 |
| log10 hop/flip | 2.68 | 2.57–3.10 | libration −0.29 | ±1 dex |

The reference definitions:

- **Tc:** midpoint of the lower Cp-peak interval.
- **σ and H_R:** the archived values of the withheld runs.
- **E_a:** block-error-weighted σT fits of the archived γ-started Green–Kubo series, over the resolved points at 600–725 K and over 750–1000 K.
- **Li–S fraction:** matrix-frame Li–S term over the total, averaged over τ = 2–8 ps.
- **Rate ratio:**
  - D_Li comes from the archived PBE0 γ 800 K Nernst–Einstein conductivity and the ideal γ cell volume (14,759 Å³);
  - the hop rate is 6D/(2.5 Å)²;
  - the flip rate is 4D_ω/φ², with φ = 70.5° and D_ω = 2.48e-4 rad²/ps.

## Oracle

The oracle makes 40 evaluations for 50 of 60 units (83 %). Its steps:

1. Bisects the PS4 angle step for each reference: 15 purchases.
2. Buys the heat capacity of each transition interval and of the PBE0 and r2SCAN melting intervals: 5.
3. Buys D_ω at 800 K: 1.
4. Buys Green–Kubo and Nernst–Einstein points on each superionic branch: 12 requests.
5. Buys three resolved low-temperature Green–Kubo points for PBE0: 3.
6. Buys the 800 K Nernst–Einstein point: 1.
7. Buys the σ(τ) breakdown in both frames: 2.

Its graded errors are:

- Tc: 0 K;
- σ: +0.036, +0.048 and +0.027 dex;
- H_R: −0.001;
- E_a: −0.138 eV below the transition and −0.018 eV above it;
- Li–S fraction and rate ratio: 0.

Results:

- **Gate (2026-10-19):** passed; all 10 ablations fail, 4 robustness variants pass, and rounding to 4 significant digits passes.
- **Harbor oracle:** reward 1.0, runtime 47 s.

## Open risks

- **Memorisation.** The paper's text gives Tc = 750 K for PBE0, which falls inside the tolerance. No r2SCAN Tc appears as a number in the text. The tabulated Nernst–Einstein activation energies differ from the graded Green–Kubo ones, and the superionic value fails.
- **Labels inferred rather than archived.**
  - The Fig. 7 series (Cp, angle fluctuation, D_ω, D_P) has no initial-phase label in the archive. Its steps coincide with the γ-started conductivity jump for every reference (asserted), but the label is inferred.
  - The condition of the σ(τ) breakdown run is not stated.
  - The labels of the Li MSD, order parameters and pair correlations come from the paper's text.
- **Rate-ratio conventions.** The ratio depends on the hop length, flip angle and the 2D-vs-3D factor (spread about 0.5 dex). ±1 dex absorbs this spread but not the libration path.
- **Tight margin.** The superionic E_a margin to the Nernst–Einstein path is 0.024 eV, and the defensible window-to-875 K value of 0.230 eV sits 0.014 eV inside the band.
- **Dev-easy risk.** Current models probably know dp1 (Nernst–Einstein vs Green–Kubo), dp7 (frames) and dp8 (librations vs flips). The difficulty should come from dp2, dp3, dp5 and dp6, which require discovering the noise floor, the melting peak, the jump location and the branch dependence under a binding budget.
