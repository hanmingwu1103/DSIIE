# Option A execution record

Times are UTC (Taipei = UTC + 8). The author selected Option A and approved the plan on 24 October 2026 (Taipei).

1. **16:31.** Read the R1 candidate protocol, runner and archive. Confirmed that the C1 rows of `R1/candidates/design_seeds_confirm.csv` are 200 replications with seeds 9210001–9210200, M5, n = 400, p = 10, H = 4, K = 2, alpha = .20, c = 0, and that the archived C1 records have d = 2. Confirmed that `siie_q(X, Y, d, S, a)` is `dsiie(X, Y, d, H = S, K = 1, a)`.
2. **16:32:45.** Wrote `PROTOCOL.md`, `run_optionA.R` and `summarize_optionA.R` and recorded their SHA-256 in `FROZEN_SHA256.txt`, before any SIIE fit with 8 slices on an M5 sample. The runner also checks the seed ledger, the archive and the source MD5 values before fitting.
3. **16:32:50.** Ran `run_optionA.R`: 600 fits (200 samples × 3 configurations) in 0.06 min, with 0 failures. Source MD5 values equal those of the R1 confirmatory run.
   - Archive reuse: regenerated SIIE (4 slices) and DSIIE r_ave and r_min equal the archived C1 values exactly, in all 200 replications.
   - Matched fine objects: in all 200 replications the eight fine cells have equal counts (50 each), and the maximum absolute differences in weights, centres, scatters and whitening matrices are exactly 0. The centres of both fits also equal the column means of the whitened rows of the recomputed cell memberships (difference 0). Cell sizes, Y ranges and identity diagnostics are in `check_matched_fine_objects.csv`; observation memberships are reconstructed from the archived seeds and runner, not stored as indices in that CSV.
4. **16:33.** Ran `summarize_optionA.R` (`optionA_summary.csv`, `optionA_paired.csv`).

## Results (R = 200, no failures)

| Configuration | mean r_ave (MCSE) | mean r_min (MCSE) |
|---|---|---|
| SIIE (4 slices) | .8725 (.0041) | .7420 (.0108) |
| SIIE (8 slices) | .9284 (.0020) | .8788 (.0041) |
| DSIIE (4 primary slices, 2 cells each) | .8931 (.0032) | .7979 (.0075) |

| Paired contrast | r_ave | r_min | DSIIE larger (r_ave; r_min) |
|---|---|---|---|
| DSIIE − SIIE (4 slices) | +.0206 ± .0048 (gain) | +.0558 ± .0140 (gain) | 154; 160 of 200 |
| DSIIE − SIIE (8 slices) | −.0353 ± .0052 (loss) | −.0810 ± .0120 (loss) | 31; 29 of 200 |
| SIIE (8 slices) − SIIE (4 slices) | +.0559 ± .0074 | +.1368 ± .0199 | 176; 176 of 200 (first larger) |

The half-width is 1.96 sd / sqrt(200), with sd the (R − 1)-denominator standard deviation of the 200 paired differences.

Under the frozen reporting rule, DSIIE − SIIE (8 slices) is a **loss** in both metrics: the M5 gain over SIIE with 4 slices does not persist at the same fine response resolution. All results are reported in the manuscript (Table 3, Sections 4.4 and 4.6, Discussion, abstract).

Not done: no other model, sample or configuration was run, and Option B was not started.
