# Source and transformations

- **Sources** (author lists as in the archive records):
  - N. Biniskos, F. J. dos Santos, K. Schmalzl, S. Raymond, M. dos Santos Dias, J. Persson, N. Marzari, S. Blügel, S. Lounis, T. Brückel, *Complex magnetic structure and spin waves of the noncollinear antiferromagnet Mn5Si3*, Phys. Rev. B 105, 104404 (2022); Materials Cloud Archive record paxx3-a1k09, DOI 10.24435/materialscloud:85-ys.
  - N. Biniskos, F. J. dos Santos, M. dos Santos Dias, S. Raymond, K. Schmalzl, P. Steffens, J. Persson, N. Marzari, S. Blügel, S. Lounis, T. Brückel, *An overview of the spin dynamics of antiferromagnetic Mn5Si3*, APL Mater. 11, 081103 (2023); Materials Cloud Archive record 4tgt6-1p378, DOI 10.24435/materialscloud:p1-pr.
- **License:** CC BY 4.0 (both records).
- **Files used:**
  - `paxx3-a1k09/Mn5Si3.zip` (74,339,438 bytes, md5 bd4cfd1a5a489c22560ff0bf50f5432b): `DFT/inputcard` (atomic positions), `Spin-wave calculations/exp data/Fig3/data/*` (10 K energy scans near (1,2,0), gap table, (2,1,0) energy scan), `Spin-wave calculations/exp data/Fig4/*/E=*meV.txt` (10 K constant-energy scans).
  - `4tgt6-1p378/Mn5Si3.zip` (51,784,614 bytes, md5 1c5a20826c98f71ab10322b50b9f2315): `Fig 2/inputfiles/` (lattice and exchange tables `pair_MnSi_rad4.5noncol.txt`, `pairMnSi_AFMz.txt`, `pair_MnSi_rad17.txt`, `pair_MnSi_rad4.5noncol_temp.txt`), `Fig 2/INS_data_Mn5Si3/*_10K_dispersion.txt`, `Fig 2/spirit_output/spinconfig_MnSinoncol_final_h=0.00.ovf`, `Fig 1c/data/H-T phase diagram_new.txt`, `Fig 3a`, `Fig 3c`, `Fig 4a`, `Fig 4b`, `Fig 4c`, `Fig 5b`, `Fig 5c` raw triple-axis scans.

## Construction

The scripts are in `benchmarks/construction/mn5si3-noncollinear-spinwaves/` (`fetch.py`, `extract.py`, `make_task.py`, `calibrate.py`):

- Atomic positions come from the KKR input card, mapped onto the orthorhombic cell; the 16 sites of the authors' spin model are labelled Mn2a_1–8 and Mn1_1–8, the other four Mn2 sites Mn2b.
- Exchange tables: site indices are replaced by labels. Bond vectors are converted to Å: for the 4.5 Å and 17 Å noncollinear tables from the Cartesian bond vectors archived in the files, for the collinear AFMz table and the Mn1 couplings from the listed cell offsets. Every vector was checked to connect the two listed sites. The authors' own spin-wave run used a copy of the 4.5 Å table whose cell offsets do not match these bond vectors; the sold tables use the archived bond vectors.
- Neutron data are served as archived; the raw triple-axis files are reduced to Q, E, counts, errors, monitor, temperature and field.

## What is hidden

- The (0,3,Ql) dispersion table (13 energies, both branches) and all (0,3,Ql) constant-energy scans.
- The 5–20 meV energy scan at (2,1,0) and the 9.84 meV point derived from it.
- The authors' relaxed zero-field ground state (used for the canting reference).

The zone-centre excitation and canting windows and the branch-position tolerance were set from the reference route, defensible variants and wrong paths in `calibrate.py` (table in `design.json`).

## Scope

Retrospective analysis of released data with no new calculation or measurement. The Mn1 couplings are the authors' model choices, not KKR results, and they are labelled as such.
