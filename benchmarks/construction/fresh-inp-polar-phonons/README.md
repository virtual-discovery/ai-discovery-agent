# Construction of fresh-inp-polar-phonons

Fresh-start the benchmark task (see `benchmarks/tools/BUILDING.md`): the lab runs new Quantum ESPRESSO 6.7
pw.x/ph.x calculations for zincblende InP that the agent requests; the hidden truth comes from reference
calculations made through the same `task_lab.run` code path.

## Rebuild

1. `python3 benchmarks/construction/fresh-inp-polar-phonons/fetch_pseudos.py --install` downloads the SSSP 1.3.0
   PBE-precision archive (Materials Cloud rcyfm-68h65) into `~/.cache/ai-theorist/fresh-inp-polar-phonons/`,
   checks the archive md5s, extracts In and P, checks them against the SSSP metadata md5 and the pinned sha256
   (`pseudo_checksums.json`) and copies them into `environment/lab/`. (`environment/lab/SSSP_LICENSE.txt` is the
   archive's `LICENSE.txt`.)
2. `bash benchmarks/construction/fresh-inp-polar-phonons/reference_runs.sh` runs every reference calculation with
   `run_request.sh` (runtime image, `--cpuset-cpus 12,13,28,29`, `mpirun -np 2`); results land in
   `~/.cache/.../runs/<name>/result.json`. About 2.5 h on 2 cores. The parsed results are kept in
   `reference_results.json` so step 3 can be rerun without QE.
3. `python3 -B benchmarks/construction/fresh-inp-polar-phonons/make_task.py` writes `private/hidden.json`,
   `private/hidden_values.json`, `solution/ablation_values.json` and `calibration.json`, and prints the calibration.
   It uses the oracle's own analysis functions (`solution/solve.py`).
4. `.venv/bin/python -m benchmarks.tools.check_task fresh-inp-polar-phonons --cpuset 12,13,28,29` writes `check.json`.

## Lab

- `environment/lab/qe_driver.py` writes the inputs (ibrav 0 fcc primitive cell, optional z scaling for
  `strain_zz`, nbnd 13 so the band edges on the mesh are reported, fixed occupations, conv_thr 1e-12,
  ecutrho = 8 ecutwfc; ph.x at q = 0 with epsil, tr2_ph 1e-15, trans true for `dfpt` and false for
  `dielectric`), runs `mpirun -np 2 pw.x/ph.x -npool 2`, parses pw.x output and the ph.x dyn file (raw force
  constants, eps, Z*, raw frequencies) or ph.x output (eps, Z* for `dielectric`), and deletes the scratch.
- `task_lab.validate` checks ranges and charges from the cost table in `public.json`; charges track measured wall
  time (about 1 unit per minute on 2 cores, minimum 1).

Measured wall times (2 cores of this host, shared load 4-16): scf 13-60 s; dielectric k4 1.0 min, k6 2.0, k8 3.8-6.1,
k10 6.3-8.5, k12 11.9 (not offered); dfpt k4 2.2, k6 4.1-6.0, k8 Gamma-centred 4.6 (dfpt k8 is not offered:
the shifted 8x8x8 phonon run is estimated at 8-12 min under load); strained dfpt k4 3.1,
k6 8.3-8.6 min (the slowest offered request); dfpt at 65 Ry k6 7.9 min (not offered; 65 Ry is offered for scf,
dielectric k<=6 and dfpt k4).

## Key measurements (a = 5.96 A unless stated)

| mesh | eps_inf | Z*(In) raw | Z*(P) raw | gap on mesh (eV) |
|---|---|---|---|---|
| 4 Gamma-centred (dfpt, design survey) | 139.8 | 3.31 | | 0.43 |
| 6 Gamma-centred | 46.73 | 2.811 | -3.317 | 0.44 |
| 8 Gamma-centred | 25.82 | 2.688 | -2.851 | 0.44 |
| 12 Gamma-centred | 15.51 | 2.628 | -2.657 | 0.44 |
| 4 shifted | 12.23 | 2.628 | -2.712 | 1.63 |
| 6 shifted | 11.700 | 2.612 | -2.621 | 1.25 |
| 8 shifted | 11.691 | 2.610 | -2.611 | 1.02 |
| 10 shifted | 11.733 | 2.609 | -2.609 | 0.88 |
| 12 shifted | 11.770 | 2.609 | -2.609 | 0.79 |

- TO is k-insensitive (284.57 k6 shifted, 284.64 k4 shifted, 284.63 k8 Gamma-centred at a0); 65 Ry changes TO by
  0.04 cm^-1, eps by 0.003 and Z* by 0.002 (cutoff converged).
- EOS (6x6x6 shifted, 50 Ry, 5.88-6.04 A): Birch-Murnaghan a0 = 5.9624 A, B0 = 59.3 GPa; zero-pressure 5.9616 A;
  65 Ry 5.9619 A.
- PBE gap at Gamma (Gamma-centred 4x4x4 SCF): 1.16 eV at 5.70 A, 0.43 at 5.96, 0.26 at 6.02, 0.15 at 6.06,
  0.04 at 6.10 A (the lab's upper limit; the gap nearly closes on expansion).
- Strained cell (a0, strain_zz = -0.02, 6x6x6 shifted): TO 290.56 (z singlet), 292.88 (xy doublet);
  Z*_zz/Z*_xx(In) 2.529/2.611; eps_zz/eps_xx 11.309/11.419.

## Reference values (private/hidden.json)

a0 = 5.9624 A; TO 284.53, LO 321.06 cm^-1; Z*(In) 2.6116; eps_inf 11.737; eps_0 14.944; dlnZ/dlnV 0.596;
dln eps/dlnV 1.778; gamma_TO 1.325; splitting grows under compression (d ln(LO-TO)/d ln V = -0.32);
strained LO q||z 326.29, q||x 330.19 cm^-1 (difference 3.90).
TO and force constants from 6x6x6-shifted DFPT; eps_inf and Z* from the mean of the 8x8x8- and 10x10x10-shifted
responses; derivatives from centred +-1% steps in a.

## Build notes and deviations from the design

- **Agent-chosen numerics, no post-processing endpoint.** labkit charges every evaluation at least 1 unit, so the
  design's uncharged dynmat endpoint was dropped; the lab returns raw force constants, raw Z* and eps, and the agent
  applies sum rules and the nonanalytic term itself.
- **Added a `dielectric` kind** (ph.x electric field only, trans false) and offered shifted and Gamma-centred meshes
  at every k. This makes k convergence of eps/Z* a budget decision separate from the phonons (the TO converges on
  coarse meshes), and makes the Gamma-centred trap a silent one at 8x8x8 and above (eps 26 at 8, 15.5 at 12).
- **ASR decision point dropped.** At converged shifted meshes raw Z*(In) and the sum-rule value differ by < 0.005
  (violation 0.001-0.009) and the acoustic frequencies are 5 cm^-1, inside the method spread, so it cannot be graded.
  The sum-rule violation remains as evidence for dp1 (0.17 at 8x8x8 Gamma-centred) and dp2 (0.08 at 4x4x4 shifted).
- **No withheld configurations.** The design's hidden strain (-1%) and hidden lattice runs were replaced by grading
  the agent's own calculations at specified configurations (strain_zz = -0.02, larger than the design's -1% because
  the LO anisotropy at -1% is only ~2 cm^-1) against reference runs of the same code path.
- **Splitting trend graded by sign only.** d ln(LO-TO)/d ln V is a near-cancellation (2 dlnZ - dln eps - 1 = -1.59
  against the TO stiffening); defensible estimates range from -0.14 to -0.37 because the eps derivative still drifts
  with k (1.62, 1.73, 1.83 for 6/8/10 shifted), so only the sign ("grows") is graded. The informative wrong path is
  the Z*-only argument (Z* falls, so "shrinks").
- **k drift of eps_inf.** Shifted meshes undersample the small-gap region around Gamma, so eps_inf drifts up slowly
  beyond 8x8x8 (0.04 per step); the reference uses the 8/10 mean and the tolerance (0.20) covers 6-12 shifted.
- **A coarse mesh used consistently at all volumes is not a wrong path** for the derivatives (4x4x4 at both side
  volumes gives dlnZ 0.63, dln eps 1.56, both inside tolerance), so only the mixed-mesh path is an ablation (dp3).

## Calibration

`calibration.json` (full) and `design.json` -> `calibration`. Tightest margins (reference error of the nearest
failing ablation vs tolerance; robustness spread):

| field | tolerance | max robust error | nearest failing ablation |
|---|---|---|---|
| eps_inf | 0.20 | 0.022 | 0.51 (k4 shifted) |
| born_charge_In | 0.020 | 0.0067 | 0.060 (k4 shifted) |
| eps_0 | 0.35 | 0.022 | 0.66 (k4 shifted) |
| dln_eps_inf_dlnV | 20 % | 6.3 % (one-sided) | 74 % (cheaper side mesh) |
| strained LO anisotropy | 0.8 cm^-1 | 0.015 | 1.87 (isotropic tensors) |
| omega_LO | 2.0 cm^-1 | 0.22 | 7.5 (eps_0 screening) |
| strained LO q_z, q_x | 2.0 cm^-1 | 0.15 | 2.04 (isotropic tensors; also fails the anisotropy) |

## Oracle

`solution/solve.py`: 5 SCF (EOS) + 1 Gamma-centred SCF (gap) + dielectric 4/6/8 shifted at a0 + DFPT 6x6x6 shifted
at a0 and at a0(1 +- 1%) + strained DFPT 6x6x6 shifted: 13 evaluations, 33 of 40 units, about 30-35 min of compute.
Gate and Harbor results: see the end of this file.

The dp1 ablation uses the 8x8x8 Gamma-centred response (eps 26.47, Z* raw 2.692/-2.857); the agent can buy the same
tensors as `dielectric` k8 kshift 0 (identical numbers), and the TO is k-insensitive (284.63 vs 284.57 cm^-1).

## Open risks

- Strained DFPT k6 is the slowest request (8.3-8.6 min measured on 2 dedicated cores under host load); on a slower
  or oversubscribed host it could exceed 10 minutes (the driver timeout is 50 minutes, so it would still finish).
- The truly k-converged eps_inf of PBE InP may lie a few tenths above the shifted-mesh values offered here; the
  graded quantity is defined as converged with respect to the lab's parameters, and the tolerance covers 6-12 shifted.
- Models that already default to shifted meshes for semiconductors bypass dp1; dp2/dp3 still require a paid k ladder.

## Results

- **Gate** (`check_task --cpuset 12,13,28,29`, 2026-10-19T21:22Z, full rerun on the final package): passed, no issues;
  7 decision points, 11 ablations all failing, 9 robustness variations all passing, 4-significant-digit rounding passes,
  grading deterministic; oracle 13 evaluations, 33 of 40 units, 0 failed; 1975 s container time.
- **Harbor oracle** (`runs/build-fresh-inp-polar-phonons/oracle-r1-fresh-inp-polar-phonons`): reward 1.0; all 14 numerical
  checks true; lab pinned by the coordinator's daemon to cores 2,18,3,19; 35 min wall (20:48-21:23Z).
- **Disk:** peak about 90 MB in `~/.cache/ai-theorist/fresh-inp-polar-phonons/` (61 MB archive, deleted after the
  gate); 25 MB of QE outputs from the reference runs kept there; task package 2.7 MB.
