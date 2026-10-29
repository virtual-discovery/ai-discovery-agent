# Source and transformations

- **Paper:** A. S. Hauck, M. Jin, B. R. Tuttle, *Nonradiative carrier capture by point defects in irradiated GaN*, J. Appl. Phys. 139, 245702 (2026), DOI 10.1063/5.0325688.
- **Data:** A. Hauck, M. Jin, B. Tuttle, *Structures and energies of GaN native point defects, divacancy, and carbon substitutional from hybrid functional calculations*, Zenodo (2026), record 21893631, DOI 10.5281/zenodo.21893631, https://zenodo.org/records/21893631.
- **License:** CC BY 4.0 (Zenodo record metadata, `cc-by-4.0`). Every staged file derived from the archive carries this license.
- **File fetched:** `gan_defects.zip`, 25,125,074 bytes, md5 85d4a96b1bf0bcf032c89bf92dbb021f, sha256 5112c84e0d6cd8f5041c67da3fb1ece63a087dcd1d191641e5de1664e929fe25.
- **Members used:** all 67 calculation folders (POSCAR, CONTCAR, OSZICAR, OUTCAR; 269 files). Their sha256 values are listed in `benchmarks/construction/gan-irradiation-defect-levels/members.sha256`.

## What the archive contains

VASP 6.3.1 runs of 96-atom orthorhombic wurtzite supercells (9.584 × 11.067 × 10.377 Å) with a screened hybrid functional (AEXX = 0.305, HFSCREEN = 0.2), ENCUT = 400 eV, ISPIN = 2, Gaussian smearing of 0.05 eV:

- `perfect`, `perfect/kpt_spec` and `perfect/kpt_gam`: defect-free cells with the Ga PAW set (3d in the core, 384 electrons), at k = (½,½,½) (`perfect`, `kpt_spec`) or at Γ (`kpt_gam`);
- `perfect/WAV_ref` and `perfect/kpt_gam/d10`: Γ-point single points with the Ga_d PAW set (864 electrons);
- the defect folders V_Ga (q = −3…+1), V_N (−1…+3), Ga_i at two sites (−1…+3), N_i at two sites (−3…+3), Ga_N (−2…+6), N_Ga (−3…+3), V_Ga–V_N (−4…+4) and C_N (−1…+1). They are all at k = (½,½,½), with the Ga PAW set and fixed NUPDOWN (0 or 1).

## Construction

The scripts are in `benchmarks/construction/gan-irradiation-defect-levels/`:

- **`fetch.py`** downloads the zip through Zenodo's per-file link, verifies md5 and sha256 and unpacks it.
- **`extract.py`** parses every OUTCAR, OSZICAR, POSCAR and CONTCAR into `pool.json`:
  - the INCAR tags printed in the OUTCAR, the PAW titles (no POTCAR content), k-point, NELECT, NUPDOWN;
  - ionic-step count and convergence flag;
  - energy(sigma→0), TOTEN and E-fermi of the last step;
  - the initial (POSCAR) and final (CONTCAR) structures;
  - the final eigenvalues and occupations, core potentials and LORBIT site moments.
- **`dfpt/run_dfpt.sh`** computes the dielectric tensors served by the lab (see below), and **`make_dielectric.py`** writes `dielectric.json` from its output.
- **`make_task.py`** writes the catalog, the lab data, the private references and the ablation values; **`reference.py`** and **`model.py`** hold the analysis.

The lab serves each calculation's parsed record unchanged, one purchase at a time. The catalog lists the folder names and quotes the archive's own description.

## Dielectric tensors (not from the archive)

The archive contains no dielectric data. The task builders computed the electronic (ε∞) and static (ε0) tensors of wurtzite GaN by density-functional perturbation theory, as follows:

- Quantum ESPRESSO 6.7 (pw.x, ph.x, dynmat.x), PBE;
- pslibrary 1.0.0 ultrasoft potentials `Ga.pbe-dn-rrkjus_psl.1.0.0.UPF` (Ga 3d in valence; sha256 4451a9da…) and `N.pbe-n-rrkjus_psl.1.0.0.UPF` (sha256 311bcb45…);
- 45/360 Ry cutoffs and an 8×8×5 k-mesh;
- the lattice (a = 3.1948 Å, c = 5.1886 Å) and internal parameter (u = 0.3775) of the archive's relaxed `perfect` supercell;
- ε0 adds the Γ-point polar-phonon contribution (dynmat.x, crystal acoustic sum rule).

The lab returns both tensors, and the catalog states that they come from the builders. They are PBE values, not hybrid-functional values.

## What is hidden

Nothing in the archive is withheld: every calculation can be purchased (budget 43 of 68 units). The graded quantities are derived values that the lab never serves:

- the band gap of the defect calculations' Hamiltonian;
- the transition levels of V_Ga, V_N, Ga_i, N_i and C_N;
- the V_Ga–V_N binding energy and the Ga and N Frenkel-pair energies at E_F = E_CBM − 0.5 eV.

The reference protocol, its alternatives and the wrong paths are documented in the construction README and in `design.json`.

## Scope

This is a retrospective analysis of released calculations, not new defect calculations. The reference values follow the builders' analysis protocol (eFNV corrections with the builders' PBE dielectric tensor), not numbers printed in the paper; the paper's own levels may differ.
