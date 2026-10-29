# Source and transformations

- **Source:** Z. Liasi, R. Zerdoumi, F. Thelen, G. Arruda de Oliveira, R. Zehl, N. Pukhareva, L. H. Morais, H. H. Kristoffersen, A. Ludwig, W. Schuhmann, J. Rossmeisl, *Work Function and High-Coverage Adsorption Energy as Hydrogen-Evolution Descriptors on Ag-Au-Pd-Pt Alloys*, arXiv:2608.28347 (2026).
- **Dataset:** Zenodo record 22031607, DOI 10.5281/zenodo.22031607, https://zenodo.org/records/22031607.
- **License:** CC BY 4.0.

## Files used

| file | bytes | md5 |
|---|---|---|
| `01-dataset.tar.gz` | 220,296,618 | eff9518f36a306a6db04347c5b427032 |
| `03-inference.tar.gz` | 105,625,609 | 77eacf62b4e5d107a018649fc55b9ea9 |
| `04-activity.tar.gz` | 331,894,430 | 1541a46cecb4c58b3aea0852a233079a |
| `02-models.tar.gz` (read only for the held-out surrogate metrics) | 107,299,438 | 4f184f3bd3798bf60da36aa27ee0a624 |

## Construction

The scripts are in `benchmarks/construction/her-alloy-workfunction-descriptor/`.

- `fetch.sh` downloads the archives and verifies the checksums.
- `extract.py` reads, for each library:
  - `04-activity/experimental-data/<lib>/`: the EDX compositions and the SECCM currents at -100 to -400 mV;
  - `03-inference/output/<lib>/`: the DimeNet++ and linear-combination surrogate outputs per measurement area.

  The 9,000 dilute site energies per area are stored as a histogram (0.005 eV bins). This reproduces the authors' production fits: Lib-1 model A, fitted centre, R2_log 0.9759 vs 0.976 in `runs/summary-all.csv`. The surrogate work functions are reduced to their mean, as the authors do.

  From `01-dataset/data/` it also reads:
  - the DFT adsorption energies and work functions (`wf_top_ev`);
  - the coverage paths (`coverage-isotherms.json`).

  Lib-3 MDI 80 has no EDX composition (`#DIV/0!`) and is dropped, as in the authors' fits.
- `make_task.py` writes the catalog, the descriptor file, the lab data, the private targets and the ablation values.
  - Coverage paths are served as the single-site (1/9 ML) energies plus the differential energies along the filling path. The single-site energy of the hollow added at step n is reconstructed as dG_diff(n) - dG_rep(n); this is exact by the authors' definition.
  - The authors' alpha/omega slopes and `repulsion.json` are not served.

## What is hidden

- The SECCM currents of 301 Lib-3 measurement areas (`spots_to_predict`) at all four potentials. Only 40 randomly chosen Lib-3 areas can be bought, and at most 5 of them.
- Reference values derived from the released data:
  - the volcano centre, from a joint fit to all Lib-1/Lib-2 rows;
  - the equimolar coverage shift, omega*(8/9) from all 10 decorations;
  - alpha_wf, from the combined model fitted to all Lib-1/Lib-2 rows.

The authors' fitted parameters, leave-one-library-out results and controls (`04-activity/runs/`) are not staged.

## Scope

This is a retrospective re-analysis of released measurements and calculations; no new data are produced. Tolerances: see `design.json` -> `calibration`.
