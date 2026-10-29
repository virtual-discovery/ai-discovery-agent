# Construction of metal-resistivity-te-vs-anharmonicity

This is a released-start task for the benchmark (see `benchmarks/tools/BUILDING.md`).

- **Source:** Wang, Yin, Goudreault, Côté, Hellman and Poncé (2026), DOI 10.1103/yk57-bt6t.
- **Archive:** Materials Cloud record zxf0b-rs406 (DOI 10.24435/materialscloud:4h-r1, CC BY 4.0).
- **Provenance, checksums and what is hidden:** the task's `SOURCE.md`.

## Rebuild

From this folder:

1. **Fetch.** `python3 -B fetch.py [DEST]` downloads the 237.6 MB tarball, `README.md` and `Acknowledgements.md` through the archive's per-file links and checks their md5. It then extracts and verifies the 150 members in `members.md5`; the three 0.2–0.3 GB `inv_tau_freq_0.fmt` files are skipped. DEST defaults to `~/.cache/ai-theorist/metal-resistivity-te-vs-anharmonicity`.
2. **Extract.** `bash run_in_image.sh extract.py /data` writes `pool.json` (2.1 MB): 36 EPW runs and 32 force-constant sets, parsed without changing any number.
3. **Generate.** `bash run_in_image.sh make_task.py` (about 15 s) writes:
   - the catalog (it quotes `/data/README.md`, so step 1 must have run);
   - `lab/data.json` and `lab/phonons.py`, copied from here;
   - the private targets and `hidden_values.json`;
   - `solution/ablation_values.json`;
   - the budget in both `public.json` files;
   - `calibration.json`.

   `model.py` holds the estimators, and `phonons.py` the numpy phonon code shared with the lab.
4. **Gate.** From the repository root, `.venv/bin/python -m benchmarks.tools.check_task metal-resistivity-te-vs-anharmonicity --cpuset 2,3,18,19`.
5. **Design.** `python3 -B record_oracle.py`, then `python3 -B write_design.py` (writes `design.json`).

`run_in_image.sh` runs a script in `ai-theorist-benchmark-runtime:20261007` on cpuset 2,3,18,19 (override with CPUSET).

## What the archive shows

- **36 EPW iterative-BTE runs.**
  - Pb: harmonic (one run, 10–610 K), TE, TE+anh (10–610 K) and anharmonic at the 0 K lattice (210/410/610 K).
  - Al: harmonic at 50³/60³/70³ fine meshes (100–900 K), TE and TE+anh at 481/682/900 K.
  - Nb: harmonic at 40³/50³/60³ (100–2700 K), TE at 1000/1500/2000 K, TE+anh at 1000/2000 K.
  - All iterative runs converge to below 1e-6.
- **Nb TE+anh 1500 K has no results.** Its `ibte.out` is a copy of `ibte.in`. The design's Nb 1500 K target therefore has no truth.
- **Mesh mismatch in the Al and Nb harmonic references.** The TE and TE+anh runs use 50³ meshes; the harmonic runs exist at three meshes. For Al at 481 K the mesh changes Δ_TE by 2 points (8.94% against 10.91% with 70³). The archived `1-DFPT/ibte.in` corresponds to the 50³ output only.
- **The effects are large and not separable (Pb).**
  - ρ_TE/ρ_harm is 1.25, 1.40, 1.60 and 1.92 at 310–610 K.
  - The anharmonic correction is −18/−25/−29% at the 0 K lattice (210/410/610 K) but −18/−30/−48% at the expanded lattice.
  - TE+anh is within 2% of harmonic above 400 K for Pb, and within 1.5% for Al.
  - Nb does not cancel: TE+anh/harmonic is 0.865 at 1000 K and 0.890 at 2000 K.
- **10 K is grid-limited.** Pb gives 1.57 μΩ cm (0 K lattice) against 2.95 (10 K lattice, 0.1% larger) and 2.86 (TDEP). These are 100× Bloch–Grüneisen values.
- **Phonon code.** `phonons.py` (WS-weighted Fourier interpolation, simple ASR) reproduces the archived matdyn.x paths of Pb, Al and Nb to 5×10⁻⁵ cm⁻¹. The Pb 410/610 K DFPT force constants violate the ASR by 2% of the on-site term (−0.76 meV acoustic modes at Γ without the ASR); EPW and matdyn both impose it.
- **Archive quirk (harmless).** The Pb TDEP 210 K `updated_qe.fc.xml` carries the 410 K geometry in its header. TDEP itself ran in the 210 K cell, and frequencies do not depend on the header alat. Both values are served.
- **Surrogate tests (central to the design).** For two force-constant sets at the same lattice, the EPW ratio was compared with phonon kernels. Simple moments (⟨ω⁻²⟩, ⟨ω⁻¹⟩, ω_log, ⟨ω²⟩) miss by 5–17% (Pb up to +16%, Nb +24% at 1000 K). The λ-like ⟨ω⁻²⟩ is also grid- and cutoff-sensitive (near-Γ modes). The transport-weighted kernel Σ|q|²(x/sinh x)²/ω², with q folded into the first BZ, reproduces:

  | Pairs | Surrogate / EPW |
  |---|---|
  | Pb, T-lattice, 110–610 K | 1.000–1.03 |
  | Pb, 0 K lattice | 1.001–1.011 |
  | Al | 1.005–1.011 |
  | Nb | 1.054 and 1.081 |
  | Pb, 10 K | 0.86 |

  It is grid-converged at n = 12–24. A longitudinally projected q² kernel fits Nb better (1.03, 1.01) but Pb worse (1.05–1.15).
- **Thermal expansion goes beyond phonons.** Across lattices, thermal expansion includes electronic and deformation-potential changes. A phonon-only q² estimate of ρ_TE/ρ_harm is 7–10% low for Pb (the residual grows with T), 3–6% high for Al, and within 2% for Nb.

## Deviations from the design (candidates.json / DETAILS.md section 8)

- **Targets.**
  - The design's (1), Pb TE+anh at 310 and 510 K (±5%), became Pb TE+anh at 510 K (±5%), with TE+anh at 410–610 K withheld. With 410 and 610 K served, both points are trivial interpolations, within 1–2.5%.
  - The design's (2), Nb TE+anh at 1500 K, is not archived. It became Nb TE+anh at 2000 K (±5.5%), with the 1000 K TE+anh run served for calibration.
  - The design's (3), Al TE+anh at 682 K, is not discriminating: every route lands within ±4.5%, and only SERTA and TE-only fail. It was replaced by Al Δ_TE at 481 K (±1 pp), which tests the mesh-matched reference.
  - The design's (4), the Pb 0 K-lattice anharmonic change at 410 K (±3 pp), is kept, and 610 K is added. All 0 K-lattice runs are withheld, so the change must come from a surrogate validated on released expanded-lattice pairs.
  - Added: Pb ρ_TE at 610 K (±6%), the superlinear regime near melting, with the 610 K TE run withheld.
- **Lab.** The lab serves archived transport runs, and computes phonons (grid, polarisation vectors, DOS, moments, dispersions) from any archived force-constant set.
  - No electronic-structure request: the ibte.out DOS column is adaptive-smearing noise (Al 1.06–5.2 states/Ry across T and meshes), and the archive has no Wannier Hamiltonians for N(E_F) or v_F.
  - The spectral decompositions (Pb 610 K) are not served, because the TE one would leak the withheld TE 610 K scattering rate.
- **Decision points.**
  - Kept from the design, as dp1–dp8 here: SERTA vs IBTE (dp7), the 10 K regime (dp8), the surrogate class (dp4), which force constants represent TE+anh and non-separability (dp2, dp3), electronic changes with expansion (dp6, which matters for Pb rather than Nb), and Pb T-dependence (part of dp3 and dp6).
  - Added: the mesh-matched reference (dp1), from the archive's mesh variants, and the Nb transfer (dp5).
  - Unit conversion is a stated convention and is not a separate fork.

## Decision points (8)

| dp | Fork | Kind |
|---|---|---|
| dp1 | Harmonic reference at the TE runs' 50³ meshes, not the "converged" 70³ run | metadata inconsistency, visible after purchase |
| dp2 | TE+anh from combined physics, not by combining separate effects | model class / cascade |
| dp3 | The anharmonic correction depends on the lattice; it cannot be transferred from the expanded lattice | discriminating evidence (0 K-lattice runs withheld) |
| dp4 | Transport-weighted q² kernel validated on released pairs vs the printed ⟨ω⁻²⟩ | model class, found by validation |
| dp5 | Nb: no cancellation, and its own calibration (d-band coupling) | transfer where the Pb/Al model breaks down |
| dp6 | Pb TE near melting: superlinear, with a growing electronic residual | regime judgement |
| dp7 | IBTE, not SERTA (printed first) | data reading |
| dp8 | 10 K grid-limited / Bloch–Grüneisen: excluded from fits and calibration | data quality / regime |

dp4 cascades into the anharmonic and TE+anh fields, and dp8 into the surrogate calibration.

## Calibration (full table in `design.json` → `calibration`)

| Field | Reference | Oracle | Defensible routes (pass) | Nearest wrong paths (fail) | Tolerance |
|---|---|---|---|---|---|
| Al Δ_TE 481 K (%) | 8.945 | 8.945 | 8.915 (linear interpolation) | 10.91 (70³ ref., +1.97 pp), 11.09 (60³) | ±1.0 pp |
| Pb ρ_TE 610 K (μΩ cm) | 85.44 | −1.8% | −1.3% to −4.6%: q² + residual (−1.8), uniform kernel + residual (−1.3), longitudinal q² (−2.7), cubic fit (−2.2), quadratic ratio (−3.2), constant residual (−3.8), quadratic fit (−4.6) | linear extrapolation −9.4%, phonon-only −12.2% | ±6% |
| Pb Δ_anh,0K 410 K (%) | −25.19 | −25.48 | −24.4 to −25.5 (q² calibrated or uncalibrated, longitudinal q², calibration on 110–310 K) | expanded-lattice transfer −29.25 (4.1 pp), ⟨ω⁻²⟩ −29.68 (4.5 pp) | ±3 pp |
| Pb Δ_anh,0K 610 K (%) | −29.46 | −29.59 | −28.7 to −29.6 | expanded-lattice transfer −45.1 (15.7 pp) | ±3 pp |
| Pb ρ_TE+anh 510 K (μΩ cm) | 36.72 | +1.7% | −3.6% to +3.9%: q² (+1.7), longitudinal q² (+3.9), calibration on 110–310 K (+2.1), ratio-to-harmonic extrapolation (+0.8), linear fit 110–310 K (−3.6) | calibration incl. 10 K +5.9%, SERTA +7.0%, ⟨ω⁻²⟩ +9.3%, linear fit incl. 10 K −9.8%, multiplicative +17.9%, additive +34.7% | ±5% |
| Nb ρ_TE+anh 2000 K (μΩ cm) | 91.76 | +2.5% | −2.9% to +4.3%: q² Nb-calibrated (+2.5), longitudinal q² (−2.2), ratio to harmonic (−2.9), ρ ∝ T (+4.3) | Pb-calibrated surrogate +6.75%, SERTA +7.9%, uncalibrated +8.1%, constant TE+anh/TE +11.6%, cancellation +12.3% | ±5.5% |

Not ablations, because they pass:

- SERTA on the Al field (9.52) and on the anharmonic fields;
- calibration including 10 K on the anharmonic fields (−22.36/−26.64; that path fails at 510 K);
- the uncalibrated ω_log⁻² moment (2.5/2.8 pp, grey zone);
- the calibrated ω_log⁻² moment for Pb at 510 K (+4.8%, grey zone).

Rounding every number to 4 significant digits passes.

## Oracle

The oracle makes 28 evaluations for 41 of 46 units (89%):

- 13 transport purchases:
  - Al 50³ harmonic and TE 481 K (it checks that the meshes match);
  - Pb harmonic, TE 210–510 K and TE+anh 210/310 K;
  - Nb TE 1000/2000 K, TE+anh 1000 K and the 50³ harmonic run;
- 15 phonon grids (n = 16):
  - Pb DFPT at a(0) and a(210–610 K), TDEP 210/310/510 K, and TDEP at a(0) for 410/610 K;
  - Nb DFPT and TDEP at 1000/2000 K.

It calibrates the q² surrogate on the 210/310 K Pb pairs and leaves 10 K out. It extrapolates the Pb TE residual linearly from 410/510 K, and calibrates Nb on its own 1000 K pair. Its graded errors equal the construction model's to all digits.

The catalog is worth 86 units: 27 served transport runs × 2 and 32 force-constant sets × 1. The budget is 46 (53% of the catalog). The expected agent length is 45–75 tool calls.

## Results

- **Gate (2026-10-20, cpuset 2,3,18,19):** passed.
  - 8 decision points; all 15 ablations fail.
  - 13 robustness variants and 4-significant-digit rounding pass.
  - 8 rubric items, 14 hidden values (6 graded plus the withheld archived resistivities), no leaks.
  - Package 2.1 MB; the oracle took about 90 s in the gate container.
- **Harbor oracle (`runs/build-metal-resistivity-te-vs-anharmonicity/oracle-r1-metal-resistivity-te-vs-anharmonicity`):** reward 1.0, solver charge 41, runtime 66 s.

## Open risks

- **Dev-easy risk.** A careful agent that buys the TE+anh series can pass Pb 510 K through near-cancellation, since TE+anh ≈ harmonic. It can also pass Nb through ratio-to-harmonic, and Al by matching meshes as the convention states. The strongest discriminators are the 0 K-lattice anharmonic changes, which need a validated transport-weighted surrogate, and the Nb transfer.
- **Thin margins.**
  - Pb 510 K with calibration including 10 K fails by 0.9 points.
  - Al with the 70³ reference fails by 0.97 pp.
  - The Pb anharmonic change at 410 K from expanded-lattice transfer fails by 1.1 pp, and from ⟨ω⁻²⟩ by 1.5 pp.
  - Nb with the Pb calibration fails by 1.25 points.
  - Grey-zone surrogates (ω_log⁻²) pass by 0.2–0.5 points: luck for agents who use them unvalidated.
- **Truth defined by archived numerics.** EPW grids, adaptive smearing and the Fermi window. The targets are the authors' values, not converged physical predictions. The Al reference requires interpolating the archived 100 K temperature grid (method spread 0.03 pp).
- **The surrogate is an approximation with known limits.** The oracle's Nb estimate (+2.5%) depends on one calibration point; a second Nb point does not exist in the archive.
- **Memorisation.** The paper is from 2026. Its headline statements (opposite signs; near-cancellation for Pb and Al) are rubric-level, and no graded number is a headline value.

## Disk

- **Peak:** about 760 MB. The cache held 749 MB: the 238 MB tarball and the 522 MB extraction without the spectral files, done during exploration. The scratch subfolder held 11 MB (exploration scripts and phonon grids).
- **After the Harbor run:** both were deleted, together with the gate's lab image. `pool.json` (2.1 MB) is kept here, and `fetch.py` with `members.md5` refetches the rest.
