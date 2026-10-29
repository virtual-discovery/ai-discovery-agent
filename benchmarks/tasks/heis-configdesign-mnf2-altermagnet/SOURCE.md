# Source and transformations

- **Source:** B. Li, S. von Malottki, G.-M. Rignanese, *Magnetic-configuration design for reliable Heisenberg exchange parameters*, Materials Cloud Archive (2026), record 3p6fv-8hd59, DOI 10.24435/materialscloud:zg-1w.
- **Archive:** https://archive.materialscloud.org/records/3p6fv-8hd59
- **License:** CC BY 4.0 (all files).
- **File used:** `MnF2.tar.gz` (51,865,693 bytes, md5 4311cefc120956738f6d26440354f3fd), fetched with the archive API's per-file link.
- **Not used:** `Fe.tar.gz` and the archive README (the README is quoted only for its documented units).

## Content used

From `MnF2/spin-spiral/` we used:

- `spin-spiral-energies.csv`: 101 flat spirals, namely q = 0 and 20 q-points each along [001], [111], [113], [121] and [214], with t = 0.05 … 1.00;
- for each spiral, the VASP `INCAR`, `KPOINTS`, `POSCAR`, `OUTCAR` and `OSZICAR`.

The calculations are VASP 6.4.2 with PAW_PBE (Mn_pv, F), no +U, `I_CONSTRAINED_M = 4`, `LSPIRAL`, `ENCUT = 1100 eV`, `ENINI = 700 eV`, and a 9×9×13 k-mesh.

We did not use `fitted_J.csv`, `Figure7_rmse.csv` or the authors' compiled `__pycache__` fitting code.

No VASP POTCAR is staged or redistributed: POTCARs are skipped on extraction. We read only their TITEL lines locally, to confirm the functional.

## Construction

The scripts are in `benchmarks/construction/heis-configdesign-mnf2-altermagnet/`:

- **`fetch.py`** downloads the archive, checks its md5 and unpacks it without POTCARs.
- **`extract.py`** builds `pool.json`. For each spiral it records:
  - the OUTCAR energy(sigma→0) of the 6-atom cell;
  - VASP's "choose ENMAX > X" value;
  - the verbatim OUTCAR header warnings;
  - the INCAR tags and the k-mesh;
  - the final constrained moments `M_int` from the OSZICAR;
  - the number of electronic steps.

  It asserts two things:
  - every CSV energy equals 1000 × (E_cell(q) − E_cell(0)) to within 1.5 µeV. The CSV energies are therefore meV per 6-atom primitive cell, although the archive README documents "meV / atom".
  - the CSV moments equal the OSZICAR `M_int`.
- **`reference.py`** fits the reference model: 10 shells, sublattice-symmetrized, with an intercept, fitted to the 98 spirals whose basis is complete (the three [214] spirals with t ≥ 0.90 have X > ENCUT).
- **`tn_mc.py`** computes the reference Néel temperature of that model with the C Monte Carlo in `mc_tn.c`: classical heat bath plus over-relaxation, and Binder-ratio crossings for L = 12, 16 and 20.
- **`make_task.py`** writes the catalog, the lab data, the hidden data and the ablation values.

Lab responses are the released values unchanged. The catalog quotes the archive's documented units verbatim.

## What is hidden

The 40 spirals along [111] and [113] are never served. Their released energies (the CSV values, meV per cell, relative to q = 0) grade the predictions. This includes the ground-state energy: two of the hidden spirals, (1,1,1) and (1,1,3), realize the Néel state.

Reference values derived from the released data:

- J1, J2, J3 and the sublattice-mean exchange of the 7.0 Å shell, from the reference fit;
- the classical Monte Carlo T_N of the reference model.

The tolerances sit between the defensible analyses in `ablations.py` (ROBUSTNESS) and the plausible wrong paths (ABLATIONS). See `design.json` → `calibration`.

## Scope

This is a retrospective reconstruction from released DFT data with the authors' protocol, not new calculations.

The altermagnetic framing is our task design, not a claim of the source. It follows from the rutile symmetry: uniform spirals contain only the sublattice sum of each split shell.

Grading uses the released PBE values, never neutron or literature exchange constants, nor the experimental T_N.
