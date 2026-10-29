# Source and transformations

- **Source:** A. Tellez-Mora, X. He, E. Bousquet, L. Wirtz, A. H. Romero, *Systematic determination of a material's magnetic ground state from first principles*, Materials Cloud Archive record rbrcp-0tc34 (2024), DOI 10.24435/materialscloud:5m-2t. Paper: npj Comput. Mater. 10, 20 (2024), DOI 10.1038/s41524-024-01202-z; its Supplementary Information (section 1.3) describes the CuO calculations.
- **License:** CC BY 4.0 (all files of the record).
- **File used:** `CuO.tar.gz` (6,482,521 bytes; md5 b57a770ecd719fe8e9713b727fb1290b; sha256 6d2507643f20c18873e8a4de656c9bddebdaf1bcbfa83dcdd63b4dec80650a71).
- **Members used:** `README`, `{ligand,noligand}/{FM,AFM}/exchange.out`, the four `dispersion_*` files, and `ligand/AFM/CuO_AFM.xsf` (full-precision AFM cell for the path check). Member sha256 values are in `benchmarks/construction/tb2j-groundstate-cuo-multiferroic/fetch.py`.
- **Not used:** the `TB2J.yaml` files, which duplicate `exchange.out` (they agree to 5e-5 meV). Also not used: every other material of the record. FeP, FeCl2, FePS3 and NiO feed the sibling task `tb2j-reference-selfconsistency-set`.

## Transformations

The scripts are in `benchmarks/construction/tb2j-groundstate-cuo-multiferroic/`:

- `fetch.py` downloads and verifies the archive.
- `cuomodel.py` parses `exchange.out`: the header lines, cell, atoms (charges and moments) and every isotropic pair (i, j, R, J_iso, vector, distance). The `[Testing!]` DMI and J_ani lines are dropped, because the archive README says only the isotropic constants are valid.
- `make_task.py` performs these steps:
  - splits each calculation's pair list into ten distance windows (lab data);
  - stages the conventional cell and atoms of the FM-reference calculations;
  - computes the graded targets and the ablation values from the pair lists.
- No values are altered. Pairs are served exactly as archived (J to 1e-4 meV, vectors to 1e-3 Å). The `dispersion_*` files are served unchanged, with their header line.

## Validation against the authors' own output

- **FM/ligand spectrum.** The archived h(k) eigenvalues of the FM/ligand calculation are reproduced from its pair list to 2.4e-5 eV, using the paper's bosonic matrix (Toth–Lake, S = 1, values in eV). The path nodes are not archived: the paper gives only Y = (1/2, 0, 1/2) in its own axes, and A and X were fitted.
- **AFM/ligand spectrum.** The archived spectrum is reproduced from the archived pair list, as listed in its 2×1×2 cell and its own moments, to 1.0e-5 eV on Y–Γ–X and to 1.2e-4 eV on A–Y (the fitted node). The authors evaluated it at the FM path's Cartesian k-vectors in the AFM calculation's frame, which is rotated by about 2°.
- **Paper k_min values.** The FM/ligand Luttinger–Tisza minimum, (0.4911, 0, −0.4736) in the staged cell, equals the SI value k_min = (0.491, 0.0, 0.474), whose c axis has the opposite sign. The uncorrected value, (0, 0.5042, 0), matches the SI value (0.0, 0.5, 0.0).
- **Unreproducible spectra.** The uncorrected-calculation spectra cannot be reproduced from their archived pair lists. The FM/no_ligand spectrum is within 0.4 %; the AFM/no_ligand spectrum cannot be matched by any collinear configuration. They are served as archived, and nothing graded depends on them.

## What is hidden

No graded number is served. The graded values are derived from the archived pair lists by `make_task.py`:

- the k_min vectors of the two FM-reference calculations (full pair lists);
- the stability verdicts (lowest bosonic-matrix eigenvalue on an 8×8×8 grid of each calculation's own zone);
- the propagation vector of the converged AFM/ligand state (Fourier analysis of its moments; the state is also the Luttinger–Tisza optimum of its own pair list);
- its chain coupling, −8 × the mean TB2J value of the 32 chain bonds;
- the mean and maximum LSWT magnon energies for S = 1/2: positive eigenvalues of g·h on an 8×8×8 grid of the 2×1×2 cell's zone, with the maximum refined.

The tolerances admit the defensible alternatives listed in `solution/ablations.py` (exchange cutoffs of 8–17 Å, coarse grids, single-bond chain values) and reject the realistic wrong paths.

## Scope

This is a retrospective analysis of released calculations, not new DFT. The AFM-reference calculations used structures slightly different from the FM one: the ligand-corrected one has its origin shifted by c/2 and cell lengths differing by up to 0.4 %, and the uncorrected one has two inequivalent Cu sites. Each calculation's own cell, atoms and moments are served with its pairs.
