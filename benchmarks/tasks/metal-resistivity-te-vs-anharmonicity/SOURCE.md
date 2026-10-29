# Source and transformations

- **Paper:** A. Wang, J. Yin, F. A. Goudreault, M. Côté, O. Hellman, S. Poncé, *Opposite impact of thermal expansion and phonon anharmonicity on the phonon-limited resistivity of elemental metals from first principles* (2026), DOI 10.1103/yk57-bt6t.
- **Archive:** Materials Cloud Archive record zxf0b-rs406, DOI 10.24435/materialscloud:4h-r1, https://archive.materialscloud.org/records/zxf0b-rs406.
- **License:** CC BY 4.0, which applies to all files used.
- **Files used:**

  | File | Bytes | md5 |
  |---|---|---|
  | `data_for_paper.tar.gz` | 237,561,636 | 83486659df1249b85d0025c2231d06b6 |
  | `README.md` | 15,213 | 089442c4cc8bfd798f062c178048897c |
  | `Acknowledgements.md` | 610 | 968797ce577c1c08e4ab2b2195796e48 |

  150 members of the tarball are used, each listed with its md5 in `benchmarks/construction/metal-resistivity-te-vs-anharmonicity/members.md5`:
  - every EPW `ibte.in`/`ibte.out*` and `scf.in` under `*/4-epw/` (the `spec/` and `comp_EPW_interpolation/` subfolders are excluded);
  - every `2-ph/ph-*K/espresso.xml`;
  - every TDEP `updated_qe.fc.xml` and the TDEP `infile.ucposcar`;
  - the READMEs.

## Construction

The scripts are in `benchmarks/construction/metal-resistivity-te-vs-anharmonicity/`.

- **`fetch.py`** downloads and verifies the files above.
- **`extract.py`** parses them into `pool.json`, without changing any number:
  - 36 EPW runs: input parameters, lattice parameter, fine meshes, and per temperature the Fermi level, printed DOS, the SERTA tensor, the IBTE iteration history and the final IBTE tensor;
  - 32 force-constant sets.
- **`phonons.py`** Fourier-interpolates the force constants with Wigner–Seitz weights and the simple acoustic sum rule. It reproduces the archived matdyn.x frequencies to 5×10⁻⁵ cm⁻¹, the print precision. The lab uses the same file.
- **`model.py`** holds the estimators.
- **`make_task.py`** writes the catalog, lab data, private targets and ablation values.

Transformations:

- Resistivity is 10⁶/(trace σ/3) of the final IBTE tensor, in μΩ cm.
- The Al harmonic reference at 481 K is a quadratic interpolation of the archived 400, 500 and 600 K values of the 50³ run.
- The catalog quotes the archive README verbatim.

## What is hidden

These archived runs are withheld (never served):

- `Pb/4-epw/2-with_TE/610K`;
- `Pb/4-epw/3-anharm/{210,410,610}K`;
- `Pb/4-epw/4-TE_and_anharm/{410,510,610}K`;
- `Nb/4-epw/3-TE_and_anharm/2000K`.

The graded values are derived from them and from released runs:

- Δ_TE of Al at 481 K, using the released runs;
- ρ_TE of Pb at 610 K;
- the 0 K-lattice anharmonic change of Pb at 410 and 610 K;
- ρ_TE+anh of Pb at 510 K;
- ρ_TE+anh of Nb at 2000 K.

`Nb/4-epw/3-TE_and_anharm/1500K` is listed but not purchasable: its archived `ibte.out` repeats the input and contains no results.

## Scope

This is a retrospective analysis of released calculations; no new EPW calculation is run. The lab computes new harmonic phonons from the archived force constants.
