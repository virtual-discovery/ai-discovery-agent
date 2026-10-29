# Construction of halfheusler-strain-hinge

This is a released-start task, built from S. K. Das, I. C. Fulga, R. Dhawan, H. C. Kandpal, J. van den Brink and J. I. Facio, *Electronic structure and hinge states of strained half-Heusler compounds LiSbZn and LiBiZn*, Phys. Rev. B 112, 125132 (2025) (arXiv:2504.20190). It uses Zenodo record 16927281 (CC BY 4.0), file `strained_hH_zenodo.zip`. Provenance and checksums are in the task's `SOURCE.md`.

## Rebuild

1. Run `python3 -B fetch.py [DEST]`. It downloads the 2.7 MB zip, then checks its md5 and sha256 and the md5 of all 9 members (26 MB extracted). The default DEST is `~/.cache/ai-theorist/halfheusler-strain-hinge`.
2. Run `bash run_in_image.sh make_task.py`, which runs in the runtime image on cpuset 4,5,20,21. It does the following:
   - converts the four `+hamdata` files with `hamdata.py` into `lab/data/*.npz`; the values are unchanged and the files total 0.9 MB;
   - copies `tb.py` into the lab;
   - stages the code file, the four header blocks and the record description;
   - writes `catalog.json` and the budget in both `public.json` files;
   - computes every reference value with the lab's own code (`task_lab.run`, `tb.Model`);
   - asserts the archive facts below;
   - writes `private/hidden.json`, `hidden_values.json`, `solution/ablation_values.json`, `refs.json` and `calibration.json`.

   Notes on running it:
   - Expensive stages are cached in `refs_partial.json`. Delete that file to recompute from scratch, which takes about 40 min on 4 CPUs.
   - `--reuse` rebuilds every output from `refs.json` in about 1 min.
3. Run `.venv/bin/python -m benchmarks.tools.check_task halfheusler-strain-hinge --cpuset 4,5,20,21`, then `python3 -B record_oracle.py` (writes `oracle_run.json`) and `python3 -B write_design.py` (writes `design.json`).

`run_in_image.sh` mounts the repository, the archive cache (read-only, at `/data`) and, if `HH_SCRATCH` is set, a scratch folder.

## What the task is

The lab evaluates the four archived 8-band Wannier Hamiltonians with numpy and scipy (`tb.py`):

| Request | Cost | Returns |
|---|---|---|
| `bulk_bands` | 1 | bands and Zn-4s weights at ≤ 300 k |
| `band_edges` | 1 (mesh ≤ 24) or 2 (≤ 48) | the top of band n and the bottom of band n+1 on a mesh, and the smallest direct gap |
| `fermi_level` | 1 (mesh ≤ 32) or 2 (≤ 64) | the energy holding x electrons per cell, and the band fillings |
| `slab` | 1 | slab states finite along one basis row, with face weights |
| `prism` | 1, 2 or 3 (N ≤ 12, ≤ 20, ≤ 28) | prism states with corner, facet and bulk weights; optional onsite potentials |

How slabs and prisms are built:

- They take any unimodular integer basis in the model's own lattice vectors.
- Sites are selected by position: the Zn site is regrouped so that its reduced coordinates relative to X lie in [0, 1). This makes a prism a physical object that does not depend on how the file groups atoms into cells.
- Each result reports the basis vectors and the Zn reduced coordinates.

The budget is 30 units. The hopping blocks are not staged, so every calculation goes through the service.

## What the archive shows (and the task is built on)

- **The "strained" files are not strained.**
  - Both `+hamdata_*_eq_strained_setup.txt` files are body-centred tetragonal cells with c = √2·a_t = c₀ and the volume of the cubic cell: δ = 0 to 1e-5 %.
  - The Γ8 quartet is degenerate to 0.1 meV.
  - The LiBiZn strained-setup bands agree with `orig` to ≤ 7 meV along Γ–X, Γ–L and Γ–K (after rotating the frame).
  - The kwant pickles `*_strain.ham` are the same hoppings: verified exactly equal.
  - "eq_strained_setup" therefore means the equilibrium structure in the setup used for strain, not a strained structure. The candidate design misread it, and the paper's strain-dependent results cannot be reproduced from the archive.
- **Frames and bases differ.**
  - The strained-setup Cartesian frame is the orig frame rotated by −45° about z.
  - Its lattice vectors are (a/2)(0,1,1), (a/2)(0,−1,1) and (a/2)(1,0,−1) in cubic axes.
  - Its file groups each X with a Zn displaced by (a/4)(1,1,1) + (a/2)(−1,0,1).
  - The orig prism (axis a1, sides a2, a3) is [[1,0,0],[1,1,1],[1,0,1]] in the strained-setup basis. The C2z-equivalent basis [[0,1,0],[0,0,−1],[−1,0,−1]] gives the identical level.
  - The orig integers give a (100)/(111) prism. Its gap holds 12 facet-dominated states, and its most corner-weighted state lies at 205.5 meV instead of 175.2 meV.
- **Energy zero.**
  - The FPLO zero lies 46.7 meV above the VBM in LiSbZn orig and 288.6 meV above it in the LiSbZn strained setup.
  - The two LiSbZn calculations of the same geometry also differ in the fundamental gap (343.8 against 251.7 meV) and in the depth of the L valley below Γ6 (8 against 61 meV).
  - This is a calculation difference that mimics a strain effect. LiBiZn shows no such difference.
- **Band edges.**
  - LiSbZn has its VBM at Γ (Γ8) and its CBM about 0.003 (fractional) from L, 1.2 meV below L itself. The gap is indirect: 343.8 meV, against 352.0 meV at Γ.
  - LiBiZn is inverted: the Zn-4s Γ6 lies 1.72 eV below Γ8.
  - LiBiZn is a compensated semimetal: band 7 near L lies at −470 meV, against +400.5 meV at Γ8, so the gap is −870.6 meV.
  - With 6 electrons per cell (Li 1, Zn 2, X-p 3), E_F lies at −401.9 meV relative to Γ8 (96³ mesh, offset 0.5), i.e. at the FPLO zero to within 1.5 meV.
  - E_F against mesh:

    | Mesh | E_F − Γ8 (meV) |
    |---|---|
    | 8 | −545.6 |
    | 12 | −386.5 |
    | 16 | −379.2 |
    | 24 | −411.2 |
    | 48 | −405.5 |
    | 64 (offset 0.5) | −401.9 |
- **Hinge state (LiSbZn, prism P_N, identity orig basis).**
  - At k = 0 the gap holds one Kramers doublet. It is localized at corner cNN, where both facets are Zn-terminated (Zn at +1/4 along both sides), and carries 0.68 of its weight on Zn-4s.
  - Its level converges slowly with N:

    | N | 6 | 8 | 10 | 12 | 16 | 20 | 24 | 28 |
    |---|---|---|---|---|---|---|---|---|
    | E (meV from model zero) | 221.1 | 203.2 | 194.1 | 189.0 | 184.2 | 182.4 | 181.6 | 181.3 |

    The geometric extrapolation is 181.14 meV, i.e. 227.8 meV above the VBM. The strained-setup model gives 175.2 meV above its VBM.
  - The hinge band has its minimum at k = π, 218.3 meV above the VBM, and is converged at N ≥ 12.
  - Its two spin-split branches peak at k = 1.00, the upper one at 773.4 meV. That is inside the projected gap but far above the global CBM (343.8 meV).
  - The band never touches the valence band.
- **Removal by a local potential.**
  - V on all orbitals of the 2 × 2 cNN cells moves the k = 0 level to the VBM at V = −0.380 eV (N = 20). N = 12 and N = 16 give −0.378 and −0.379.
  - Strong attractive V pulls new bound states out of the conduction band. Potentials of either sign on the X/X corner also create in-gap states there.
  - A potential on the X/X corner (c00) leaves the cNN level unchanged, within 4 meV for |V| ≤ 1.5 eV.
- **Conflict with the paper.**
  - The paper places the LiSbZn hinge mode at the Bi/Sb hinge with about 95 % Bi/Sb weight.
  - In the authors' kwant code, `process_hamdata` stores ⟨row, 0|H|col, R⟩ under the key R, and `build_template` assigns it to `HoppingKind(R, row, col)`. That is the hopping from col at x to row at x + R, i.e. ⟨row|H|col⟩ at displacement −R. Every hopping vector is therefore reversed relative to the site positions.
  - The bulk spectrum is unchanged (time reversal). The connectivity, however, is that of the inverted twin: each Zn behaves as if it sat at X − (a/4)(1,1,1).
  - With R → −R the state moves to corner c00, which the authors draw as the Sb/Sb hinge. This was verified with the same evaluator: the level is identical, at the opposite corner.
  - The Zn-4s weight (0.62–0.68 for all N and k) is not reproduced as 5 % by either convention. The task grades the self-consistent geometry of the archived hopping vectors and states the conflict in the rubric.

## Deviations from the design (candidates.json / DETAILS.md section 29)

- **No strain targets.** The archive has no strained Hamiltonian, so graded strain effects (Γ8 splitting, Weyl or Dirac phases, strain shifts of the hinge level) were replaced by the strain of each archived model, which is 0. Recognizing this is dp1, and the LiSbZn calculation difference is the confounder.
- **Strong Z2 not graded.** It is trivially 0 for normal-ordered LiSbZn and undefined for cubic LiBiZn, whose Γ8 touching remains. "Insulating versus semimetallic" became the fundamental gap of three models, which is negative for LiBiZn, plus the LiBiZn charge-neutral E_F.
- **Hinge targets made quantitative and unique:**
  - the termination of the hinge;
  - the converged k = 0 energy relative to the VBM, for two LiSbZn calculations;
  - the hinge-band window;
  - the Zn-4s weight;
  - the potential that removes the state.

  "Number of hinge channels" was dropped: one Kramers pair at one hinge, which cannot be checked independently of the window. The LiBiZn hybridization and the visibility argument are rubric items.
- **No Wilson-loop request.** Without Z2 targets it would only serve rubric statements. Slabs give the facet spectra.
- **Prism geometry.** The paper's prism (axis a1 and 21 × 21 cells of the orig file) is defined in the instruction in cubic terms with position-based site selection, so that it can be rebuilt in the rotated file.

## Decision points (9)

| dp | Fork | Kind |
|---|---|---|
| dp1 | "strained_setup" files are cubic (lattice vectors, Γ8, bands) | archive mislabel; metadata |
| dp2 | LiSbZn CBM next to L, not Γ6 | full-zone search |
| dp3 | LiBiZn is a compensated semimetal; E_F from 6 electrons on a converged mesh | regime and convergence |
| dp4 | the FPLO zero sits at an arbitrary point in the gap and differs between calculations | energy reference; cascades into the hinge energies and the removal potential |
| dp5 | the hinge is Zn/Zn (Zn-4s character), against the paper's Bi/Sb | evidence against recall |
| dp6 | the k = 0 hinge level converges slowly with N | finite size |
| dp7 | the hinge band: minimum at π, two branches, top far above the CBM | k-resolution |
| dp8 | the same prism in the rotated file needs a different integer basis | convention of the data |
| dp9 | trivial origin from a local potential at the right hinge | discriminating calculation |

## Calibration

The full table is in `design.json` under `calibration`.

| Field | Reference | Defensible variants (pass) | Nearest wrong paths (fail) | Tolerance |
|---|---|---|---|---|
| strain (%) | 0, 0, 0, 0 | — | ±2 (file names) | ±0.25 |
| gap LiSbZn orig (meV) | 343.76 | 344.99 (L grid point) | 352.02 (Γ) | ±3.5 |
| gap LiSbZn strained setup (meV) | 251.71 | 253.38 (L) | 313.00 (Γ) | ±3.5 |
| gap LiBiZn (meV) | −870.6 | −866.9 (L) | 0 (zero-gap picture) | ±15 |
| E_F LiBiZn − Γ8 (meV) | −401.9 | −405.5 (48³), −401.9 (64³) | −386.5 (12³), −545.6 (8³), 0 (Γ8) | ±12 |
| hinge E(k = 0), orig (meV) | 227.8 | 228.0 (N = 28), 229.1 (N = 20) | 235.7 (N = 12), 240.8 (N = 10), 181.1 (FPLO zero) | ±6 |
| hinge E(k = 0), strained setup (meV) | 175.2 | 175.4 (N = 28), 176.6 (N = 20) | 183.6 (N = 12), −113.4 (FPLO zero), 205.5 (orig integers) | ±6 |
| band window (meV) | [218.3, 773.4] | [218.9, 773.4] (N = 12), [218.9, 771.6] (25-point grid) | top 734.3 (5-point grid), 684.4 (one branch), 343.8 (clipped); bottom 227.8 (at k = 0) | bottom ±6, top ±15 |
| termination | Zn/Zn | — | X/X (paper) | exact |
| Zn-4s weight | 0.681 | 0.665 (N = 16) | 0.05 (paper) | ±0.08 |
| removal V (eV) | −0.380 | −0.379 (N = 16, linear or root) | −0.313 (to the FPLO zero), null (wrong corner) | ±0.035 |

The tightest margins are these:

- the hinge energy at N = 12 fails by 1.9 meV (orig) and 2.4 meV (strained setup);
- E_F at 12³ fails by 3.4 meV;
- the window bottom taken at k = 0 fails by 3.5 meV;
- the LiSbZn Γ gap fails by 4.7 meV.

Memorized paper values fail: δ = ±2 %, the Bi/Sb hinge, "95 % on Bi/Sb", and "ε = 0.2 eV" read as the hinge energy (it is relative to the FPLO zero).

## Oracle

The oracle makes 19 evaluations for 25 of 30 units (83 %):

- 2 Γ spectra;
- band edges for 3 models, each a mesh-24 scan plus a local 6³ refinement;
- the LiBiZn E_F on a 64³ mesh;
- prisms at k = 0 with N = 12, 20 and 28, plus k = π at N = 20, and the strained-setup prism at N = 28 in its own basis, found by matching rotated lattice vectors and the Zn reduced coordinates;
- 16 + 8 k-points of the hinge band at N = 12;
- three potentials at N = 12, with a quadratic root.

The oracle uses no literature value. The slowest request is an N = 28 prism, at 21 s on 2 CPUs.

The lab threads are pinned (`OMP_NUM_THREADS=2`). The worst case allowed per request, an N = 28 prism with 2 k-points and more than 150 states in the window, takes about 7 min.

Results:

- **Gate (2026-10-25, cpuset 4,5,20,21):** passed. All 18 ablations fail, the 10 robustness variants pass, rounding to 4 significant digits passes, and grading is deterministic. The task has 9 decision points and 9 rubric items, the package is 1.0 MB with 10 hidden values, and the gate container took 143 s.
- **Harbor oracle (`runs/build-halfheusler-strain-hinge/oracle-r1-halfheusler-strain-hinge`):** reward 1.0, total runtime 2 min 36 s, 25 of 30 units.

## Open risks

- **Dev-easy risk.**
  - dp1 and dp4 can be defused by reading the staged headers carefully.
  - The strongest discriminators should be these:
    - dp5: the geometry has to be worked out, against the paper;
    - dp8: the rotated basis;
    - dp6 and dp7: convergence and k-resolution under a budget where an N = 28 prism costs 3 units;
    - dp9: the potential must go at the right corner, and it must be attractive.
- **Paper conflict (dp5).** The graded termination and Zn weight contradict the paper. I traced the cause to the authors' kwant `HoppingKind` convention, which I verified numerically. The 95 % Bi/Sb figure is not reproduced under either convention, so a reviewer may question it. The geometry used here follows directly from the archived hopping vectors and Wannier centres.
- **Thin margins:** the N = 12 hinge energy (about 2 meV) and E_F at 12³ (3.4 meV).
- **Wannier fidelity.** No DFT bands are archived, so the fidelity of the 8-band models cannot be checked; this is a rubric item.

## Disk

- **Peak:** under 0.1 GB of files (2.7 MB zip, 26 MB extracted, about 12 MB of scratch: paper PDF and npz copies), plus the gate's and Harbor's lab images, which are thin layers on the runtime image.
- **After the build:** the raw downloads, the scratch files, the gate image and the Harbor images were deleted. `fetch.py` (with its checksums), `refs.json` and `refs_partial.json` are kept.
