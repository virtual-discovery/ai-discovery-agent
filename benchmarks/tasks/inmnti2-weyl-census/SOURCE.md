# Source and transformations

- **Study:** D. Grassano, L. Binci, N. Marzari, *Type-I antiferromagnetic Weyl semimetal InMnTi2*, Phys. Rev. Research 6, 013140 (2024), DOI 10.1103/PhysRevResearch.6.013140; arXiv:2208.11412.
- **Archive:** Materials Cloud Archive record hkprp-xw234, *Prediction of a novel type-I antiferromagnetic Weyl semimetal*, DOI 10.24435/materialscloud:ph-3c, https://archive.materialscloud.org/records/hkprp-xw234.
- **License:** CC BY 4.0 (all files).
- **Files used:**
  - `upload.zip`: 160,482,239 bytes, md5 27d1aa76a5cc384765108beb2ae05289, sha256 a360dae1e68dc8fc0d060d89088d00b6ba68aa820b510fb0a0cc5715057dfa91;
  - `files_description.md`: md5 37d5be45010d2ae698b72c26c1dcd574.
- **Members:** the 96 members used are listed with their md5 in `benchmarks/construction/inmnti2-weyl-census/members.md5`. They include, for example, `hubbard/wannier/InMnTi2_hr.dat` (md5 1a14e924…), `standard/wannier/InMnTi2_hr.dat` (md5 5b52edd7…), `hubbard/nscfXYZ/nscf.out` (md5 4cd3cf8a…) and `hubbard/wannier/kubo/300/InMnTi2-kubo_A_xy.dat` (md5 28ecad3b…). The two Wannier90 `.wout` logs (579 MB and 49 MB) enter only as excerpts (header and final state). The surface spectral functions, Fermi-arc maps, bulk band and node-dispersion outputs of WannierTools are not used.

## Construction

The scripts are in `benchmarks/construction/inmnti2-weyl-census/`:

- `fetch.py` downloads and verifies the archive and extracts the members.
- `make_task.py` runs in the runtime image and does the following:
  - converts the two `InMnTi2_hr.dat` files to `environment/lab/data/{hubbard,standard}.npz` (R, degeneracies, H as complex64; the conversion changes band energies by less than 0.001 meV);
  - parses the band energies of the four non-self-consistent runs (as printed by pw.x) into `environment/lab/data/archive.json`;
  - keeps the rows with ω ≤ 1.5 eV of every Kubo output file;
  - copies the WannierTools node and chirality outputs, the latter with the per-node progress lines removed;
  - copies 21 small archive files unchanged into `environment/agent/archive/`, and writes excerpts of `hubbard/relax/relax.out` and of the two `.wout` logs;
  - computes the references and the ablation values with the lab's own evaluator (`tbmodel.py`) on the lab's npz files.
- No file is altered in content. No pseudopotentials are in the archive or the task.

## What is hidden

Graded values are derived quantities stored only in `environment/lab/private/hidden.json`:

- **Node counts:** 24 nodes between DFT bands 52/53 in each calculation, from a 40³ (DFT+U) or 36³ (PBE) gap scan, gap minimisation and Berry flux on the archived Wannier models.
- **Orbit sizes:** the partition of the DFT+U nodes under the 8 operations of the archived magnetic state (`hubbard/scf/scf.out`).
- **Node energy, pair separation, k_z offset:** the crossings of DFT bands 52/53 on the three archived DFT+U lines (`hubbard/nscfXYZ`), relative to the scf Fermi energy. Chirality comes from the Wannier model.
- **PBE node energy:** the mean energy of the 24 PBE Wannier nodes relative to the scf Fermi energy. The PBE model agrees with its DFT lines to 1.7 meV.
- **Hall conductivities:** σ_xy(ω = 0) of the densest archived Kubo mesh; σ_yz = σ_zx = 0 by symmetry.

The lab never returns these values directly. The paper's node energy (96 meV below E_F) and pair separation (0.37 Å⁻¹) agree with the DFT-line references within their tolerances. Its k_z shift ("0.0047 Å⁻¹") and the statement that all nodes are equivalent do not.
