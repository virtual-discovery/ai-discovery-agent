# Source and transformations

- **Study:** K. Nakazawa, T. Yamaguchi, A. Yamakage, *Orbital-related gyrotropic responses in Cu2WSe4 and chirality indicator*, Phys. Rev. Materials 10, 064603 (2026); arXiv:2512.16387.
- **Archive:** Zenodo record 20131063, *Dataset for "Orbital-related gyrotropic responses in Cu2WSe4 and chirality indicator"*, DOI 10.5281/zenodo.20131063, https://zenodo.org/records/20131063.
- **License:** CC BY 4.0 (all files).
- **File used:** `CWS_PRM.zip` (24,610,495 bytes, md5 9a50ca688ae85f597776e7d6e5931765). All 34 members are listed with their md5 in `benchmarks/construction/cu2wse4-orbital-edelstein/members.md5`, for example:
  - `Cu2WSe4_soc.HWR`: md5 3dbab966f32b7d1b848e9fdd07d973fb
  - `Cu2WSe4_wosoc.HWR`: md5 a7afc29a389172bf3ca94f785e887388
  - `Fig2_orb.d`: md5 d2b668ff60c3d830221d8854485ecb34
  - `Fig2_spin.d`: md5 550371de189f2038c7424c8aa20782c1
  - `Fig4c_xyz.d`: md5 1139902fb07a5ff541f808a738da3f72

## Construction

The scripts are in `benchmarks/construction/cu2wse4-orbital-edelstein/`:

- `fetch.py` downloads and verifies the zip and its members.
- `make_model.py` parses the two OpenMX `.HWR` Wannier Hamiltonians (Hartree, Bohr, Wigner–Seitz degeneracies). The lab copy stores H(R) as integers in units of 1e-9 Hartree; the file has 12 decimals, so the rounding error is at most 5e-10 Hartree (14 neV).
- `cws_tb.py` is the numpy evaluator. It is copied unchanged into the lab and computes:
  - bands and band velocities;
  - Berry curvature;
  - the Berry-phase orbital moment, the intra-atomic −μ_B⟨L⟩ and the spin moment.

  It uses Wannier centres at the atomic positions, and treats band pairs within 5 meV as degenerate.
- `make_cache.py` stores 1 meV energy histograms of the Fermi-surface summands on Γ-centred N × N × N/2 meshes (N = 16 … 128).
- `make_task.py` writes:
  - the staged files (archive files unchanged, the two band files gzip-compressed, and the first 9 lines of each `.HWR` file);
  - the compact lab caches;
  - the catalog, the hidden targets and the ablation values.

**Validation against the archive** (construction README):

- The Wannier bands match the archived DFT bands (Fig1b) to 2.5 meV near the gap.
- The evaluator reproduces the archived Edelstein curves (Fig2, in 10⁻⁵ μ_B per formula unit) to 2–3% rms over −2…2 eV, with no energy shift:
  - the spin curve to a scale factor of 1.002;
  - the orbital curves (with and without SOC) to a factor of −2 (×1.01). The archived orbital curves carry a factor −2 relative to the paper's Eq. (3); the evaluator's sign and normalisation equal the direct wave-packet moment (q/2)⟨r × v⟩ (checked on a model molecule).
- The archived Berry-curvature-dipole curves (Fig4c,d) are reproduced in shape with a factor of about 1.9, a convention factor of 2 in the definition of σ_ijl.

## What is hidden

Graded values are stored only in `environment/lab/private/hidden.json` and come from the validated evaluator on its densest mesh (128 × 128 × 64) at the stated conventions:

- the orbital-to-spin ratios at VBM − 0.10 eV and CBM + 0.10 eV;
- M_x^orb at VBM − 0.10 eV;
- the SOC-to-no-SOC orbital ratio at p = 5 × 10²⁰ cm⁻³.

The symmetry fields come from the D2d point group derived from Fig1a.cif.

The archived response data Fig2, Fig3 and Fig4 are neither staged nor served. Archive values are not graded directly because of the orbital normalisation issue above.
