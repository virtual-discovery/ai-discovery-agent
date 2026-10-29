# Source and transformations

- **Source:** Y. L. Müller, A. U. Kurbay, X. Zhang, E. Kioupakis, A. R. Natarajan, *Thermodynamic and electronic properties of rutile Sn1−xGexO2 alloys from first principles*, arXiv:2601.12184 (2026). Data: Materials Cloud Archive record 1qt4h-tr585, DOI 10.24435/materialscloud:e2-xp.
- **License:** CC BY 4.0 (archive record).
- **Files used:**
  - `data.tar.xz` (8,574,684 bytes, md5 ba7dcda1497afecf1863f2975a73eca0, sha256 98aef23f3b7adb9d7c707bd5cefe727d431f813946a4706be562c4b1365eaa26). Members used: `data_CE_{PBE,r2SCAN}.json`, `data_SQS_{PBE,r2SCAN}.json`, `data_{001,100,101}_{PBE,r2SCAN}.json`, `Alloy_PBE_Values/*/{OUTCAR,CONTCAR,INCAR,KPOINTS}`, `{GeO2,SnO2}_PBE_values/{OUTCAR,CONTCAR,INCAR,KPOINTS}`.
  - `README.txt` (4,813 bytes, md5 c44362bc1f34d2d6b172f9f2a982a24e).

## Construction

The scripts are in `benchmarks/construction/sngeo2-coherent-miscibility/`:

- `fetch.py` downloads and verifies the two files.
- `extract.py` parses every POSCAR string and OUTCAR into `pool.json.gz`, with structures as lattice, species and fractional coordinates. It asserts the archive relations used:
  - the strain-scan minima equal the SQS end-member energies;
  - the subregular fit of the released SQS values reproduces the paper's Table 1 (α = 0.915 / 1.071 eV, β = −0.117 / −0.219 eV per primitive cell);
  - the HSE-setting PBE totals give half of the released x = 0.5 values per cation, within 2 meV.
- `make_task.py` writes the catalog, the lab data, the private targets and the ablation values. `model.py` holds the thermodynamics.
- Energies, structures and settings are served unchanged. No POTCAR content is staged; only the POTCAR labels printed in the OUTCAR are.

## What is hidden

Nothing graded is served. The targets are derived quantities, recomputed by the task builders from the full released data:

- the PBE mixing enthalpy at x = 0.5;
- the coherency strain energies at x = 0.5 on the three released planes, and the habit plane;
- the incoherent critical temperature (graded as a window that admits mean-field and cluster-expansion Monte Carlo estimates);
- the coherent spinodal critical temperatures (PBE, and r2SCAN after excluding three r2SCAN SQS that did not reach their equilibrium volume);
- the 700 K coherent binodal and spinodal (PBE).

The paper's headline numbers (incoherent Tc above 2300 K, coherent Tc of about 900 K) are not graded as such. The r2SCAN coherent Tc differs from the paper's 890 K because the three unrelaxed SQS are excluded; see `benchmarks/construction/sngeo2-coherent-miscibility/README.md`.

## Scope

This is a retrospective reconstruction from released DFT data with the paper's free-energy model, not new calculations. The HSE band-gap data in the archive are not used.
