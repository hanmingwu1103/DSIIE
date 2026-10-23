# Application inputs

`processed/` preserves the 13 predecessor CSV files and hashes. The complete current analysis uses 18 inputs; `analysis_ready/` includes all 18 as RDS objects and 17-digit CSV files, with source/licence/preprocessing metadata and export checks.

Use RDS for estimator execution: it preserves the exact original matrix attributes and binary-double values. CSV values round-trip exactly but intentionally omit R-specific attributes. Row order is the original analysis order and must be preserved for folds, grouped observations and time blocks.

Run `Rscript --vanilla data/check_analysis_ready.R` from `current/` to verify the bundle. Add `--source` to compare it against the original package/processed definitions; this requires the data-source package versions listed in the release README. `export_analysis_ready.R` is a deterministic reconstruction utility, not a new preprocessing specification.
