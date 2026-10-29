# Source and transformations

- **Source:** J.-C. Chen, J. Leveillee, C. G. Van de Walle, F. Giustino, *Design of high-mobility p-type GaN via the piezomobility tensor*, Materials Cloud Archive (2025-08-04), record qfdfa-cs517, DOI 10.24435/materialscloud:zy-qw. The archive accompanies the paper of the same title (submitted 2025).
- **License:** CC BY 4.0 (all files).
- **Files used:**
  - `strain_GaN.tar.gz` (12,934,938 bytes, md5 9a8747242f4759d235c5552aaea0ba52);
  - `README.txt` (973 bytes, md5 e72a76b7c622cd39d0d3a91ea75811e6).

  Per-member md5 checksums of the 148 members read (`scf.in`, `epw2.in`, `gan.eig` and `epw2.out` of the 37 states) are in `benchmarks/construction/gan-hole-piezomobility/members.md5`.

## Construction

The scripts are in `benchmarks/construction/gan-hole-piezomobility/`.

- **`fetch.sh`** downloads both files from the archive's per-file links, checks them and unpacks the tarball.
- **`extract.py`** reads, for each state:
  - from `scf.in`: the cell in bohr and the crystal coordinates;
  - from `epw2.in`: `fermi_energy`, which is the GW valence-band maximum the authors used;
  - from `gan.eig`: the GW quasiparticle energies on the 10×10×10 grid. Bands 47–56 are staged: the top six valence and the lowest four conduction bands, unchanged at 6 decimals.
  - from `epw2.out`: the zero-field SERTA table, the final table of the zero-field iterative BTE solve, the iteration count, the last "Max error" and the "did not converge" flag.

  The magnetic-field solve that follows in each `epw2.out` (Hall factor, B = 1e-10 T) is not used. Its tensors differ from the zero-field ones by about 0.1%. The script asserts that each `fermi_energy` equals the archived valence-band maximum at Γ.
- **`make_task.py`** writes:
  - `catalog.json`: labels, cells, positions, the k-grid, methods and costs;
  - `lab/data.json`: per-state SERTA and IBTE records of the purchasable states;
  - `lab/bands.json`: band energies of all 37 states;
  - the private targets;
  - the ablation and robustness values.
- **`write_design.py`** writes `design.json`, including the calibration table.

Mobility tensors are served as printed by EPW (3×3, cm²/Vs, at 100–500 K). Grading uses the symmetrized 300 K IBTE tensor. No pseudopotential, Wannier or phonon data are staged.

## What is hidden

The mobilities (SERTA and IBTE) of five states are never served: `compress_3%_e3`, `expand_3%_e3`, `expand_3%_e2`, `expand_3%_e4` and `compress_3%_e6`. Their bands are served. Their archived 300 K IBTE tensors grade the predictions:

- `compress_3%_e3`: a regime window on μ_zz;
- the other four: ±10% on the diagonal elements and ±15% on the strained off-diagonal element.

Two design answers are derived from the archived 2% states and graded:

- the 2% state with the largest μ_zz: `compress_2%_e3`;
- the 2% state with the largest in-plane principal mobility: e6 shear, either sign.

Both answers must be reported as Voigt strains computed from the cells. The archive's shear labels apply ε_ij = X%, an engineering strain of 2X%.

## Scope

This is a retrospective reconstruction from released GW+EPW data, not new calculations. Each state was relaxed and computed independently. Mirror-equivalent shear pairs therefore differ numerically by up to 3% on the diagonal elements (6% for μ_zz at 1% e4).
