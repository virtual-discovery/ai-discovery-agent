# Shared runbook — the stripped round-8 GW-BSE chain for nonmagnetic band insulators

This is the pipeline every card in this folder reuses. It is the α-RuCl3 round-8 chain (QE 7.4.1 → BerkeleyGW 4.0) with the DFT+U, Hund and
spin–orbit stages removed. Everything below was exercised on RuCl3 (16 atoms, spinor, 440 bands, 16–64 k) and is cheaper here: 8-atom scalar
cells with 30–60 occupied bands.

## 1. Rules that are not negotiable

1. **Freeze before you compute.** Rename `acceptance_draft.json` → `acceptance.json`, fill every `TO_FILL` field, record `sha256sum`, and
   commit it before any production stage (S2 in every card). Never edit it afterwards; a change is a new card.
2. **Reproduce the Li first (G0).** Every card names the closest published calculation. The LITCHECK gate runs before production and can stop
   the campaign. We were the outlier once (RuCl3: a published GW-BSE got 1.23 eV where our chain got 1.48; the fault was ours).
3. **No fitted parameters.** No U, no scissors tuned to experiment, no strain or geometry adjusted to hit a number. The one allowed freedom is
   the single global shift defined in the card, and it is printed on every figure.
4. **Report components separately.** E_QP and E_b, or d_S and F_in, are always reported next to each other; a gap error and under-binding
   compensate in the absolute energy (round-8 "compensation stated" lesson).
5. **Gates on observables, not residuals.** Every gate is scored by a frozen scorer that reads a results JSON; it must be possible for each
   gate to pass and to fail, and the record states which outcomes were possible before the numbers were read.
6. **Stop at the verdict.** Extra chains need a new card.

## 2. Software

| component | version | notes |
|---|---|---|
| Quantum ESPRESSO | 7.4.1 (pw.x, pw2bgw.x, projwfc.x) | CPU build, gfortran + OpenMPI + OpenBLAS; no Hubbard patch needed for these cards |
| BerkeleyGW | 4.0, complex flavour, HDF5 | GPU: NVIDIA HPC SDK, OpenACC, `env/arch_nvhpc_cc80.mk` (cc80 = A100; use cc90 for H100), HDF5 1.14.5 parallel, ScaLAPACK 2.2.2: `env/setup_gpu_vast.sh` |
| Pseudopotentials | PseudoDojo v0.4 NC PBE, "standard" | scalar-relativistic for the primary chains; fully-relativistic copies only for the SOC covariates. Download from pseudo-dojo.org; record file sha256 in the card |
| Python | numpy, h5py, spglib, matplotlib | `env/fit_scissors.py`, `env/pick_states.py` are ours |

Build traps seen on rented GPU hosts (all in `env/setup_gpu_vast.sh` comments): `make -j` returns before sub-makes finish → wait for `pgrep -x make` to
empty before collecting `bin/`; `rm -rf bin` before `make` lets the first `ln -sf` create a FILE named bin → `mkdir -p bin` right after; a
toolchain tarball built with `-fast` on Intel dies with "Illegal instruction" on AMD EPYC → rebuild from source or build with a portable `-tp`.

## 3. Stage list (one chain = one strain point / one stacking / one system)

| stage | tool | template | output | notes |
|---|---|---|---|---|
| S-A SCF | pw.x | `templates/scf.in` | charge density | fixed occupations; converge ecutwfc to 5 meV on the KS gap (S0 ladder) |
| S-B NSCF coarse | pw.x | `templates/nscf.in` | WFN bands | full uniform grid, `nosym`, `noinv`; bands ≈ 10× occupied for epsilon/sigma |
| S-C NSCF shifted | pw.x | `templates/nscfq.in` | WFNq | q0 = (0.001, 0.001, 0) crystal; occupied + 8 bands |
| S-D export | pw2bgw.x, wfn2hdf.x | `templates/pw2bgw_WFN.in`, `pw2bgw_WFNq.in` | WFN.h5, WFNq.h5, RHO, vxc.dat | vxc diag range = the sigma band range |
| S-E epsilon | epsilon.cplx.x | `templates/epsilon.inp` | epsmat.h5, eps0mat.h5 | 10 Ry, ladder to 14 Ry once; slabs: `cell_slab_truncation` |
| S-F sigma | sigma.cplx.x | `templates/sigma.inp` | eqp1.dat | G0W0, `exact_static_ch 1`; bands VBM−3…CBM+3 on every coarse k |
| S-G scissors | `env/fit_scissors.py` | — | evs/ev0/evdel/ecs/ec0/ecdel | two-branch linear fit of eqp1 vs mean field; NEVER interpolate eqp_co from a coarser grid |
| S-H kernel | kernel.cplx.x | `templates/kernel.inp` | bsemat.h5 | nv × nc per card; window ladder ×1.5 is a gate |
| S-I absorption | absorption.cplx.x | `templates/absorption_b1.inp` (b2, b3 copies) | absorption_eh.dat, eigenvalues.dat, eigenvectors.h5 | momentum operator, one polarization per run; exact readout co = fi; `write_eigenvectors -1` |
| S-J readout | `env/pick_states.py` + card scorer | — | results/*.json | leading sticks (not envelope maxima), overlaps, polarization ratios |
| S-K real space | plotxct.cplx.x, `env/run_realspace.sh` | `templates/plotxct.inp` | conditional densities | §5 rules |

Reference cost on RuCl3 (for scaling): 16-k spinor chain with 440 bands, epsilon 16 q + sigma 8 bands + kernel 16v16c + absorption ≈ 2.75 h on
2× H100; 64-q epsilon ≈ 4.5 h per run on 2× A100; 64-k 16v16c kernel 1.5 h. Scalar 8-atom cells are ≈ 4–8× cheaper per k-point.

## 4. Trap list (each cost us hours or dollars on RuCl3)

- **Gauge.** `eigenvectors.h5` is tied to the exact WFN it was computed with. Regenerating the NSCF and reusing stored eigenvectors gives real-space
  maps uncorrelated with the originals while summary numbers look fine. Keep WFN + eigenvectors together or recompute absorption.
- **plotxct WFNq_fi** must be the UNSHIFTED fine WFN; the shifted file gives incoherent volumes.
- **Velocity operator.** The shifted NSCF for `use_velocity` needs nbnd well above the BSE window and `diago_full_acc = .true.`; otherwise the
  finite-difference velocity blows up by 100–800×. Check the velocity/momentum peak ratio ≈ 1 before using velocity spectra.
- **Envelope vs member.** Gate on the leading stick(s), not on the maximum of a broadened envelope; report both.
- **eqp_co interpolation** from a coarser grid produced chaotic spectra; the scissors fit is stable.
- **Energy-shift exclusions** (if ever used): ±1e5 eV, not ±1e3, and gate at the exciton level, not on the dielectric head.
- **Stale directories.** "copy unless exists" reuses a failed attempt's NSCF dir; make stages resumable by a state file, not by existence.
- **Interrupted drivers.** A record assembled from partial outputs reads as complete. Audit the pre-registered stage list against files on disk
  and state "N of M planned stages ran" (EXECUTION_PROTOCOL.md §4).

## 5. Real-space rules (plotxct)

1. `supercell_size` = the fine k-grid exactly (Born–von Kármán period). A smaller box is a window: density outside is never computed and `np.roll`
   wraps the floor, producing fake rectangular support.
2. Translate the hole by whole unit cells to the box centre (same site by translation symmetry).
3. `OMP_NUM_THREADS=1` (unpinned threads made a 4-minute job run > 2 h).
4. Render on a log scale over 4 decades; bound excitons span > 5 orders of magnitude.
5. Containment: face density / max ≤ 1e−5 before calling a state bound; otherwise report "box-limited" and do not score the character gate.
6. Chain/layer assignment of grid points is a deliverable (script + figure), not an assumption.

## 6. Cloud GPU rules (Vast.ai or equivalent)

- Read the real price from the raw offer listing (`dph_total`) and write it into `vast_run.json` before creating; abort if price × estimated hours
  > 60 % of the remaining credit. We lost $92 once to a mis-read $14.79/h offer.
- Launch with a dead-man deadline (MAX_HOURS = estimate × 1.5 + 1); run every job detached with the log on disk; never `stop`, only `destroy`.
- Destroy in the SAME pipeline as the pull + checksum verification; match the exact completion string, never a generic "FAILED"/"error" word on
  an accumulating log.
- Never `pkill -f <pattern>` from a shell whose own command line contains the pattern; kill by PID.
- Heavy files (WFN, epsmat, bsemat, eigenvectors) stay on the instance unless the card needs them locally; pull readouts, eqp tables, logs,
  and record sha256 of what was left behind.

## 7. Record format (every card)

`results/*.json` (scorer inputs), `RESULT.md` with a verdict box, the gate table with every number and its source file, figure captions that
stand alone, provenance (inputs, versions, instance ids), `MANIFEST.sha256`, `make_figures.py` that regenerates every figure from `results/`.
Every data file gets a figure with a full caption. The first paragraph of RESULT.md names the limitations the card declares out of scope.
