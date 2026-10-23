# Ledger of additional candidates (simulation models and data sets)

All candidates considered for the R1 revision, evaluated or not. The rule for featuring a simulation
candidate was fixed in `candidates/PROTOCOL.md` before any fit. Under that rule, DSIIE's paired 95%
interval for r_ave must lie above zero against at least five of the seven comparators, including SIIE.
At most one comparator may beat DSIIE significantly. For d = 2, DSIIE's mean r_min must be among the
two largest. Only confirmatory runs (fresh seeds, R = 200) count as evidence.

## Simulation models

| ID | Model | Rationale | Exploration (R=50, seeds 91x0001..) | Confirmatory (R=200, seeds 92x0001..) | Decision |
|---|---|---|---|---|---|
| C1 | M5, adaptation of Xia's (2007) quadratic location–scale example as used in Wu (2026), d=2, (n,p)=(400,10), H=4, K=2, alpha=.20 | Quadratic plus heteroscedastic effects | DSIIE better than 6/7; worse than csMAVE | DSIIE 0.893 (r_min 0.798); better than SIR, SAVE, pHd, DR, SIME, SIIE (SIIE +0.021±0.005; r_min +0.056±0.014; DR +0.052±0.010); worse than csMAVE (−0.073±0.006) | **Featured** (Table 1, Table 2, Sect. 4.6). Rule met |
| C1 sens | M5 at (200,10), (800,10), (400,20); seeds 93x0001.. | Sensitivity of C1 | – | Same order in all three: DSIIE better than the six other methods, worse than csMAVE (e.g. vs SIIE +0.024, +0.011, +0.026) | Reported in text as a qualification |
| C2 | TAQ3 two-index quantile-width design (historical student model), d=2, H=5, K=2 | Signal in interior quantile widths | DSIIE better than SIR, SIME only | DSIIE 0.726; better than SIR, SIME; n.s. vs SAVE, DR; worse than pHd (−0.141), csMAVE (−0.261), SIIE (−0.025±0.013) | Not featured; adverse outcome reported in Sect. 4.6 |
| C3 | M2t: M2 with t_3/sqrt(3) nuisance coordinates (new), d=1 | Heavy-tailed nuisance; quantile widths vs moments | DSIIE better than 4/7 | DSIIE 0.985; better than SIR, SAVE, pHd, DR, SIME; n.s. vs SIIE (−0.001±0.004); worse than csMAVE (−0.015±0.002) | Not featured (SIIE not beaten); reported in Sect. 4.6 |

Not evaluated: M8 (heteroscedastic nuisance tied to the central coordinate; violates the complement
condition and has no specific DSIIE rationale). Historical TAM2/TAM3 endpoint-reach designs were not
evaluated either: their mechanism is the extreme-endpoint reach targeted by SIIE-MM, not the
interquantile widths of SIIE-Q/DSIIE. No M7 oracle-whitening ablation was run; it is not required
by the plan.

Attribution clarification from the follow-up source check: the current M5 design agrees with Wu's
predecessor specification, whereas Xia's accessible primary preprint gives a related example with
different predictor and noise specifications. The manuscript and this ledger therefore use
"adaptation" rather than the unverified exact designation "Model I of Xia (2007)". This changes
neither the model that was fitted nor the outcomes; the frozen protocols retain their historical
wording. The baseline plus the three sensitivity settings give four M5 configurations in total.

## Data sets

No new data set was evaluated. The existing candidate pool, covering RD1R, the student-only sets
and the documented exclusions, had yielded the 18 analysed regressions. The approved plan allowed
relevant additional data sets with fair evaluation, performance-qualified inclusion and transparent
candidate records. No suitable additional candidate was identified, and the revision did not extend
the search beyond this pool. Retaining the existing applications was a discretionary scope decision,
not a prohibition on further exploration. No performance claim is made about unevaluated data sets.

Editorial alternatives considered with existing analyses:

- **Replace College with Concrete strength (DSIIE rank 3, csMAVE better by 0.12).** Not done. The
  contrary evidence for Concrete strength is already quoted in Section 5.2 and Table 7, so the
  default six were kept.
- **Replace the Plasma figure with Real estate.** Not done. Plasma is kept as the weak-signal
  contrast (plan default), and no corrected Real estate figure exists.

The original 18-dataset summaries are unchanged and are not combined with any new evaluation.
