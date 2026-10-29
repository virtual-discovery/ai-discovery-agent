# Source and transformations

- **Source:** K. Aicher, T. Berger, U. Aschauer, O. Diwald, *Ubiquitous shallow trap states and lattice hydrogenation of ZnO particles*, Materials Cloud Archive (2025), record c91vp-p9h03, DOI 10.24435/materialscloud:f8-mc. Companion paper: J. Mater. Chem. C (2025), DOI 10.1039/D5TC03036F.
- **License:** CC BY 4.0 (archive). The paper's text is not used.
- **File used:** `ZnO_ubiquitous.tar.gz` (727,535,698 bytes, md5 d21e5cfb405c798dcd23ad41c9678614). It was streamed once, and only members under 12 MiB were extracted; vasprun.xml files were skipped.

## Construction

The scripts are in `benchmarks/construction/zno-hydrogen-shallow-donor/`.

**Fetch and extract.**

- `fetch.py` streams the tarball, verifies the md5 and extracts the small members.
- `extract.py` parses every relaxation and hybrid single point under `chem_pot/` and `supercell_3_3_2/` into `pool.json`. From the OUTCAR it reads the final energies (TOTEN and σ→0), NELECT, ISPIN, ENCUT, LHFCALC/AEXX, the Fermi energy, the k-points and weights, the last eigenvalue/occupation block, the last "average (electrostatic) potential at core" block and the site moments. It also reads the INCAR, KPOINTS, POTCAR titles (POTCAR.spec only) and the POSCAR/CONTCAR structures.
- Three PBE relaxation OUTCARs exceed 12 MiB and were not extracted: H_O q_0, H_O q_+2 and Zn_i-H_O q_+1. For those runs only the stdout `log` is parsed (energies, magnetization, ionic steps), and the lab says so. Their hybrid single points are complete.

**Make the task.** `make_task.py` writes:

- the catalog (run ids are the archive folder paths);
- the lab data (one record per run, plus eigenvalues within −5…+4 eV of each run's Fermi energy);
- the private targets and the ablation values.

**Check against the authors.** The authors' notebook protocol (`analyze.ipynb` and `defect_analysis/`) was re-implemented in numpy. It uses the hybrid energies, chemical potentials at the Zn-rich, maximal-μ_H vertex, the Lany–Zunger 2/3 Madelung term with ε = 10.4, and core-potential alignment beyond 5 Å. It reproduces the paper's formation energies at E_F = 3.85 eV to about 0.01 eV: V_O 1.32, OH_O 1.34/1.42, H_O 1.71, H_Zn 4.38, Zn_i 3.95, Zn_i–H_O 4.94 and H_Zn–OH_O 3.41 eV.

**What is not staged.** No VASP POTCAR, vasprun.xml, notebook, analysis code or experimental file is staged.

## How the truth differs from the authors' headline numbers

- **Band filling.** The neutral OH_O and H_O cells (and Zn_i⁰/Zn_i⁺) hold their extra electron in conduction-band states at several k-points. The reference removes this band filling, using the Lany–Zunger/Persson form with the OUTCAR k-weights and the alignment-shifted host CBM. The authors did not; their OH_O(+/0) level lies at 3.78 eV, above their own 3.25 eV gap.
- **Fermi level and concentrations.** These come from charge neutrality with parabolic bands (m_e* = 0.24, m_h* = 0.59, Fermi–Dirac), as stated in the instruction. The authors instead used py-sc-fermi on the 4-k-point supercell DOS. Their self-consistent E_F of 3.85 eV lies above the CBM, and is neither staged nor graded.

## What is hidden

Every number derived from the released energies:

- the gap;
- Δμ_O and Δμ_H;
- the OH_O, H_O and V_O transition levels;
- the 1000 K Fermi level;
- log₁₀ of the V_O, OH_O⁺ and OH_O⁰ concentrations;
- the categorical answers (the EPR centre OH_O⁰, the BC∥ hydrogen site, and the instability of V_O⁺).

The lab serves only raw run data.
