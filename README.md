# DSIIE: dispersion-based sliced inverse interval estimation

Research code for the article

> Y. Chen and H.-M. Wu, *Dispersion-based sliced inverse interval estimation for sufficient dimension reduction* (manuscript).

The repository holds the following:

- R implementations of the eight estimators compared in the article: SIR, SAVE, pHd, DR, SIME, csMAVE, SIIE and DSIIE.
- The scripts, random seeds and processed data that reproduce every table and figure.
- The replication-level simulation records and the out-of-fold predictions of the applications.

It is research code, not an R package.

## Estimators and conventions

| Method | Implementation | Response resolution |
|---|---|---|
| SIR, SAVE | `fit_SIR_eq`, `fit_SAVE_eq` (`R/methods_v8.R`); equal-count slices, divisor-*n* moments | *HK* slices |
| pHd | `dr::dr(method = "phdy")` | no slicing |
| DR | `fit_DR_LW` (`R/methods_corrected.R`); sliced kernel of Li and Wang (2007, Sec. 4) | *HK* slices |
| SIME | `fit_SIME` (`R/methods_v8.R`); weighted covariance of slice spatial medians (Dong, Yu and Zhu, 2015). Medians use the modified Weiszfeld algorithm of Vardi and Zhang (2000), and convergence is certified by the subgradient optimality condition | *HK* slices |
| csMAVE | `MAVE::mave(method = "csMAVE")` | – |
| SIIE | `siie_q` = DSIIE with *K* = 1 (called SIIE-Q in the predecessor paper) | *H* slices |
| DSIIE | `dsiie` (`R/methods_base.R`) | *H* primary slices × *K* fine cells |

The slices are equal-count blocks of the response order statistics. Ties are broken by observation index, and any remainder goes to the earliest blocks. Predictors are whitened with the divisor-*n* covariance. Directional widths use type-7 sample quantiles with α = 0.20 (endpoints 0.10 and 0.90) unless stated otherwise.

`R/methods_base.R` holds the base definitions: utilities, DSIIE, SIIE, the comparators, the metrics and the data generators. They are extracted verbatim from the authors' simulation script; some of its comments are in Chinese. `R/methods_corrected.R` and `R/methods_v8.R` add the corrected DR kernel, SIR and SAVE with equal-count slicing, SIME, and the common eight-method dispatch `fit_one8`.

## Requirements

- R ≥ 4.6 (results were produced with R 4.6.1 on Windows 11).
- CRAN packages:

  | Package | Version used |
  |---|---|
  | dr | 3.0.11 |
  | MAVE | 1.3.12 |
  | MASS | 7.3-65 |
  | mgcv | 1.9-4 (figures only) |
  | faraway | 1.0.9 |
  | gamlss.data | 6.0-7 |
  | ISLR2 | 1.3.2 |

- `parallel` is used for the simulations and the applications.

```r
install.packages(c("dr", "MAVE", "MASS", "mgcv", "faraway", "gamlss.data", "ISLR2"))
```

## Reproducing the results

Run every command from the repository root. The number of worker processes is the first argument of the simulation script; for the applications it is set by the environment variable `DSIIE_CORES`. Both default to the number of cores minus one.

| Step | Command | Time (30 cores) |
|---|---|---|
| Method checks | `Rscript checks/check_methods.R` | < 1 min |
| Simulations, 29 cells, 6050 replications | `Rscript simulation/run_simulation.R 30` | 18 min |
| Simulation summaries | `Rscript simulation/summarize_simulation.R` | < 1 min |
| Applications, 18 data sets, 8 methods, nested CV | `DSIIE_CORES=30 Rscript applications/run_applications.R` | 3.8 h |
| Application summaries | `Rscript applications/summarize_applications.R` | < 1 min |
| Table bodies (LaTeX) | `Rscript tables/make_sim_tables.R` and `Rscript tables/make_app_tables.R` | < 1 min |
| Figures | `Rscript figures/sim_figure.R` and `Rscript figures/coord_figures.R` | < 1 min |
| Quick reproduction check | `Rscript checks/reproduce_subset.R 30` | < 1 min |

- **Simulations.** A subset of cells can be run by naming them, for example `Rscript simulation/run_simulation.R 8 SIM-C_M7_n400_p10_H4_K2_a020_c005`. The cell names and data seeds are in `simulation/design_seeds.csv`.
- **Applications.** `applications/run_applications.R Mussels RD06` runs only the named data sets.
- **Quick reproduction check.** It copies the repository to a temporary directory and recomputes one simulation cell and one application data set there. It also regenerates every summary file and table body from the archived results, then compares everything with the archived versions. `checks/reproduce_subset.log` records a run in which every difference was exactly zero.

## Applications protocol

The analysis protocol is in `applications/candidate_registry.md`. It covers:

- the candidate data sets and the exclusions made before analysis;
- the fold designs (random, grouped or blocked);
- the tuning grids;
- the NRMSE metric.

It was fixed before any eight-method result was computed. `applications/apps_data.R` defines each data set, and `applications/feasibility_screen.csv` records the pre-analysis feasibility screen.

NRMSE for repetition *r* is {Σ(y − ŷ_OOF)² / Σ(y − ȳ)²}^{1/2}, with ȳ the full-sample mean, used only for evaluation. It is reported as the mean and standard deviation over the repetitions.

Two departures from a clean run are recorded:

1. **SIME on the mussel data.** In the first application run, 64 of these fits stopped because the spatial-median convergence check failed. In 5–8-point slices the minimizer can be a data point, which the Weiszfeld iterates approach only sublinearly. The check now tests the nearest data point against the subgradient optimality condition (`spatial_median` in `R/methods_v8.R`), and the mussel data were rerun with all eight methods. The seven other methods, and every previously converged SIME fit, reproduced exactly. The original predictions are kept by the authors.
2. **csMAVE on the yacht data (RD14).** `MAVE::mave` returns NaN directions on some training samples, so 79 of the 100 outer folds have at least one failed fit. These are reported as failures and are not imputed.

## Data

- **R packages.** The mussel (`dr`), plasma (`gamlss.data`), body-fat (`faraway`), baseball and college (`ISLR2`) data are loaded from R packages.
- **Processed CSVs.** The other thirteen data sets are in `data/processed/`. They were prepared under the predecessor's documented preprocessing. `data/sources.csv` gives, for each one:
  - the original source and citation;
  - the DOI or package;
  - the licence (CC BY 4.0, GPL or MIT);
  - the access date;
  - the SHA-256 checksum recorded when the file was prepared (`sha256_matches` confirms the copy here is identical).

  Please cite the original sources when you use these data.

## Results

| Path | Contents |
|---|---|
| `results/simulation/cell_*_v8.rds` | replication-level records (recovery `r_ave`, `r_min`, status, timing) |
| `results/simulation/summary_v8.csv`, `paired_v8.csv` | means, Monte Carlo SEs and paired DSIIE differences |
| `results/simulation/reproduction_vs_corrected.csv` | agreement with the authors' earlier corrected run (≤ 1.3e-13) |
| `results/applications/oof/oof_*.csv.gz` | out-of-fold predictions for both tuning grids, with selected tuning values, *H*, fit status and failure messages |
| `results/applications/*.csv` | NRMSE summaries, paired differences, tuning-boundary rates and the screening ledger |
| `tables/*.tex`, `figures/*.pdf` | the table bodies and figures in the article |

## Licence

The code is distributed under the GNU General Public License, version 2 or (at your option) any later version; see `LICENSE`. This matches the licence of the authors' `dataSDA` package. The data in `data/processed/` keep their original licences, listed in `data/sources.csv`.

## Contact

Han-Ming Wu, Department of Statistics, National Chengchi University, Taipei, Taiwan (wuhm@g.nccu.edu.tw).
