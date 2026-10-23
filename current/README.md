# DSIIE current computational release

Version **DSIIE-analysis-R2-20261024-public1**. Research code and data for Yvonne Chen and Han-Ming Wu, *Dispersion-based sliced inverse interval estimation for sufficient dimension reduction*, prepared for the Scandinavian Journal of Statistics.

This is the maintained computational release. The repository's older top-level directories reproduce the baseline at commit `432135f`; their original files are preserved. Use this `current/` directory for the revised article, including the R1 response-resolution studies, candidate models, R2 matched-resolution M5 comparison, current summaries, and equal-unit-axis figures. No manuscript or confidential review documents are distributed in this release.

The eight methods are SIR, SAVE, pHd, DR, SIME, csMAVE, SIIE and DSIIE. SIIE is the predecessor's SIIE-Q. SIIE with more slices is a response-resolution diagnostic of the same estimator, not a ninth method. The release retains all reported outcomes, including cases where DSIIE loses to the fine-resolution SIIE comparator. Option B was not run.

## Start here

Run commands from `DSIIE-Github/current` (or `DSIIE/current` after cloning).

```r
install.packages(c("dr", "MAVE", "MASS", "mgcv"))
```

R 4.6.1 was used, with dr 3.0.11, MAVE 1.3.12, MASS 7.3-65 and mgcv 1.9-4. R >= 4.6 is required by the SHA-256 manifest checks. Exported application inputs are bundled, so application execution does not require loading the source-data packages. To re-export or independently check source-data identities, also install faraway 1.0.9, gamlss.data 6.0-7 and ISLR2 1.3.2; dr supplies the mussel data.

| Command | Work performed |
|---|---|
| `Rscript --vanilla inspect_archives.R` | Verify every manifest hash; count simulations, diagnostic fits, candidate results and application failures. No fits. |
| `Rscript --vanilla data/check_analysis_ready.R` | Verify all 18 bundled analysis inputs and all 13 original processed-file hashes. No fits. |
| `Rscript --vanilla data/check_analysis_ready.R --source` | Additionally compare every bundled input exactly with its original processed/package definition. Requires the source-data packages and versions above. |
| `Rscript --vanilla regenerate_outputs.R` | Recompute all 24 summary/diagnostic CSVs and table bodies in a temporary copy; compare with their archived versions. Recreate figures, including 16 descriptive full-data fits, and check equal axis units. |
| `Rscript --vanilla rerun_selected.R` | In a temporary copy, refit all 600 Option A fits, 1,600 fits in one M7 simulation cell and one eight-method Mussels outer split with nested tuning; compare with the archives. |
| `Rscript --vanilla check_validation_failures.R` | Five disposable-copy controls: corrupt/missing manifested file, failed refit, child exit 0 without new outputs, and failed validation child must return failure. No fits. |
| `Rscript --vanilla validate_package.R` | Copy the release to a fresh directory and run inspection, input checking, regeneration and selected refits in fresh R sessions. Any child failure gives a nonzero exit. |

These commands retain the release archives. An optional log-file argument to `validate_package.R` should point outside the release directory; generated files inside a manifested release intentionally invalidate its manifest. `regenerate_outputs.R --in-place` is an explicit developer command that overwrites generated outputs. After an intentional release edit, rebuild the manifest with `Rscript --vanilla provenance/write_manifest.R`.

Validation uses the executing machine and installed packages. It is not a claim of independent-computer reproduction. New simulation/model studies and a full application rerun are not part of the validation commands.

## Code, seeds and records

| Path | Contents |
|---|---|
| `external/DSIIE_simulation_all.R` | Frozen primary definitions; `integration_work/simulation/lib_sim.R` evaluates only Parts 1-6, excluding the original drivers. |
| `integration_work/simulation/methods_corrected.R` | Correct Li-Wang directional-regression kernel and safe-fit handling. |
| `integration_work/author_revision/methods_v8.R` | Equal-count SIR/SAVE, spatial-median SIME and the eight-method dispatch. |
| `integration_work/author_revision/simulation/` | 29 baseline cells, 6,050 replications, 48,400 method fits; exact seeds and settings are stored in the replicate records. |
| `integration_work/author_revision/realdata/` | Definitions, frozen eligibility registry, nested-CV runner, all 18 application archives and summaries. |
| `integration_work/author_revision/R1/diagnostic/` | Frozen diagnostic protocol, seed ledger, 1,650 raw fits and identity checks. |
| `integration_work/author_revision/R1/candidates/` | All exploration, confirmation and sensitivity records for the three candidates, including unsuccessful inclusion candidates. |
| `integration_work/author_revision/R2_SJoS_20261024/optionA/` | Protocol, 200 existing M5 seeds, 600 matched-resolution fits, both recovery measures and fine-cell/archive identity checks. |
| `integration_work/author_revision/tables/` | Current table generators and table bodies, with leading zeros. The R1 resolution body is retained as a historical output. |
| `integration_work/author_revision/figures/`, `DSIIE_main_assets/` | Figure generators, coordinate diagnostics and output PDFs; all 16 application panels have equal physical x/y units. |
| `provenance/` | Historical frozen runners, source hash records and public-release transformations. Historical runners retain their original workstation paths and are records, not entrypoints. |

Simulated predictor/response samples are generated from the frozen generator and exact seed/design records. Replicate result archives store fitted recovery measures, timing and failures, not redundant copies of every simulated predictor matrix. No simulated samples or recorded outcomes were changed in preparing this release.

The estimator files retain their archived bytes and MD5s. Public-release edits concern relative paths, input bundling, validation, documentation and number formatting. The same frozen scientific estimators and tuning rules are used. Dataset loading from bundled RDS inputs preserves the binary doubles and attributes of the original definitions.

## Data and attribution

All **18 analysis-ready datasets** are in `data/analysis_ready/`, as exact RDS objects and CSVs written with 17 significant digits. `sources.csv` records response, predictor names, preprocessing, sample size, source package/version, licence and both hashes. `identity_checks.csv` records the export's exact round-trip checks. The original 13 processed CSVs remain byte-identical in `data/processed/`; their earlier source and checksum registry is `data/sources.csv`.

The five additionally bundled datasets are Plasma (`gamlss.data::plasma`), Bodyfat (`faraway::fat`), Mussels (`dr::mussels`), Hitters and College (`ISLR2`). `data/export_analysis_ready.R` reconstructs all inputs from the original definitions and verifies exact equality. `DSIIE_USE_SOURCE_DATA=1` restores those source definitions instead of reading the bundled RDS inputs. The source-data packages are needed only for that mode.

NRMSE is the square root of total out-of-fold squared error divided by the full-sample response sum of squares, summarized over repetitions. Nested training-fold tuning, original grouped/blocked splits, both prediction tuning grids and the recorded csMAVE failures are preserved. The full 18-dataset results are retained, not only the six featured applications.

The code licence is GPL version 2 or later (`LICENSE`). Dataset licences and source citations remain as recorded in the two data registries; bundling does not change them. Cite the original data sources when reusing datasets. In particular, the aquatic-toxicity source is Cassotti et al. (2014), *Alternatives to Laboratory Animals*, 42(1), 31-41, doi:10.1177/026119291404200106.

## Explicit full refits

Full refits are optional and can take hours. Run them in a copy or set an output directory for the baseline runners. To set options portably, start R in `current/`:

```r
Sys.setenv(DSIIE_WORK = normalizePath(".", winslash = "/"), DSIIE_CORES = "4",
           DSIIE_OUT = "rerun/simulation")
system2(file.path(R.home("bin"), "Rscript"),
        c("--vanilla", "integration_work/author_revision/run_sim_v8.R"))
```

On Windows the executable is `Rscript.exe`. The baseline runner accepts optional numeric cell positions; position 12 is the M7 contamination 0.05 cell used by selected verification. The application runner is `integration_work/author_revision/realdata/run_apps_v8.R`, with optional dataset IDs and a separate `DSIIE_OUT`. R1 diagnostic/candidate and Option A runners write into their own analysis directories; run them only in a disposable copy. Defaults are relative to `current/`; optional environment overrides are `DSIIE_WORK`, `DSIIE_SIM_SRC`, `DSIIE_DATA`, `DSIIE_ANALYSIS_READY`, `DSIIE_CORES` and `DSIIE_OUT`.

Contact: Han-Ming Wu, National Chengchi University, wuhm@g.nccu.edu.tw.
