# Source and transformations

Fresh-start task: every number the agent sees comes from new Quantum ESPRESSO 6.7 calculations (pw.x, ph.x) run by
the lab on request. No released results are served.

- **Paired study (context only, nothing used):** T. Esswein and N. A. Spaldin, *Ferroelectric, quantum paraelectric,
  or paraelectric? Calculating the evolution from BaTiO3 to SrTiO3 to KTaO3 using a single-particle quantum
  mechanical description of the ions*, Phys. Rev. Research 4, 033020 (2022); Materials Cloud Archive record
  4bdqc-a9g25 (DOI 10.24435/materialscloud:r3-df, CC BY 4.0). Its inputs were not downloaded or used; it motivates
  the question (shallow PBE polar wells, quantum paraelectricity).
- **Pseudopotentials:** SSSP 1.3.0 PBE precision (Prandini, Marrazzo, Castelli, Mounet, Marzari, npj Comput. Mater. 4,
  72 (2018)); Materials Cloud Archive record rcyfm-68h65, DOI 10.24435/materialscloud:f3-ym, CC BY 4.0. Archive
  `SSSP_1.3.0_PBE_precision.tar.gz` (md5 fde94756886f32ada7bf597547557eb5), metadata
  `SSSP_1.3.0_PBE_precision.json` (md5 1692c5c9ce89e1c7c783f8f0eee0cbfa), `LICENSE.txt` (md5
  c0f9df05fb4cd4921f4caa25606815ef, vendored as `environment/lab/SSSP_LICENSE.txt`). Vendored unchanged into
  `environment/lab/` (lab image only):
  - `Sr_pbe_v1.uspp.F.UPF` (GBRV 1.2, 4s4p5s valence; GBRV license GPL v3), sha256
    292f7f1da5cae4b18da817ab2e82943ad90e74f267a06a6172db8dbcb6edf0ce
  - `ti_pbe_v1.4.uspp.F.UPF` (GBRV 1.4, 3s3p3d4s valence; GPL v3), sha256
    747afa52fa17dc061e9eff72bcea534e81ffa5a4a5d33af3eb8fc9f1b3aee580
  - `O.pbe-n-kjpaw_psl.0.1.UPF` (PSlibrary 0.3.1 PAW; GPL v2 or later), sha256
    7c4b6ed541f83d0afdf5c1d3a8c611340f073f3e2b11e95b7963bb2ee26929aa
  Refetch: `benchmarks/construction/fresh-srtio3-soft-modes/fetch_pseudos.py --install`.

## What is hidden

The graded reference values (`environment/lab/private/hidden.json`) come from construction calculations made through
the lab's own `task_lab.run` (`benchmarks/construction/fresh-srtio3-soft-modes/reference_runs.jsonl`), analysed with
the estimators of `solution/solve.py` on converged settings. The lab never serves them; the agent recovers them from
its own calculations.
