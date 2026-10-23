# Application candidate registry and analysis specification (frozen 23 Oct 2026)

Fixed **before** any eight-method application result was computed. The only earlier outcome information is the archived student report (19 datasets, uncorrected methods) and the corrected two-dataset analysis (Bodyfat, Plasma; no SIME). Definitions live in `apps_data.R`; the feasibility screen is in `feasibility_screen.csv`.

## Candidate pool
The pool is the union of:
- the corrected predecessor registry RD1R (RD01–RD20);
- the student-only datasets: Mussels, Hitters, Boston, Swiss, College, and the student AIS version.

## Excluded before analysis
| Candidate | Reason |
|---|---|
| RD03 soil evaporation (n=46, p=10) | Infeasible: inner training samples give H < 3. |
| RD04 wheat protein (n=50, p=6) | Infeasible: inner H < 3. |
| RD05 Tecator (p=100 > inner n) | Infeasible: singular covariance in training samples. |
| Swiss fertility (n=47, p=5) | Infeasible: inner H < 3. |
| RD12 penguins (p=3) | With the common working dimension d=2 the reduction from p=3 is not meaningful. |
| RD18 forest fires | 48% of responses tied at zero (log1p area). Equal-count slices within the tie are an arbitrary index partition, so slice-based comparisons are not interpretable. |
| Boston housing | Contains the race-derived variable `B`, whose construction is ethically and statistically problematic. The data set is not used. |
| Student AIS version (dr::ais, 202×10, LBM) | Definitional leakage: body-fat percentage and weight jointly determine lean body mass. It is replaced by RD02. |

## Eligible set analysed (18 datasets, all 8 methods)
- **Random design:** RD02, RD06 (raw response, as in the corrected analysis), RD07, RD08, RD09 (13 named measurements), RD10, RD17, RD20, Mussels, Hitters, College.
- **Grouped design**, with group keys:
  - RD11: the 7 mixture-composition predictors;
  - RD13: angle, chord, velocity and thickness;
  - RD14: the 5 hull-geometry predictors;
  - RD15: identical descriptor vectors;
  - RD16: identical latitude/longitude;
  - RD19: identical descriptor vectors.
- **Blocked design:** RD01 ozone, using 10 contiguous time blocks.

Generalization targets:
- grouped: new mixtures, airfoil configurations, hulls, compounds or locations;
- blocked: new periods.

## Estimation and evaluation (all datasets, all methods)
**Methods and slicing**
- Methods, in order: SIR, SAVE, pHd, DR, SIME, csMAVE, SIIE, DSIIE.
- Settings: d=2, K=2, α=0.20.
- H = min{5, ⌊n_fit/(K(p+1))⌋}, recomputed for every fit; H ≥ 3 is required.
- SIR, SAVE, DR and SIME use HK equal-count slices; SIIE uses H slices; DSIIE uses H×K cells; pHd is response-based.

**Folds**
- Outer:
  - random: 10 repetitions of 10-fold, response-stratified, seed 1000r+7 (the corrected-analysis rule);
  - grouped: 10 repetitions of grouped 10-fold, groups in random order with greedy size balancing, seed 2000r+11;
  - blocked: one repetition, 10 contiguous blocks.
- Inner: the same design restricted to the training rows, 5 folds, seed 100000+100r+k; blocked uses 5 contiguous blocks.

**Prediction**
- Gaussian kernel ridge regression on the whitened training coordinates, with γ = multiplier × (median squared training distance)⁻¹.
- The ridge objective is Σ(y−f)² + λ‖f‖² on centred responses (λ not scaled by n).
- Fully nested: SDR, whitening and bandwidth are refitted on every inner training part.

**Tuning grids**
- Primary grid: multipliers {0.25, 0.5, 1, 2, 4} × λ ∈ 10^(−4:0).
- Wider grid: 2^(−6:2) × 10^(−4:2). History: added in the earlier corrected analysis after boundary selection was observed. It is applied here to all datasets as a fixed sensitivity analysis.
- Both grids are evaluated on the same inner fits. Ties go to the first value in grid order.

**Metric**
- Per repetition: NRMSE_r = {Σ(y_i − ŷ_i)² / Σ(y_i − ȳ)²}^{1/2} over all out-of-fold predictions, with ȳ the full-sample mean (used only for evaluation).
- Reported as the mean and SD over repetitions.
- Failed fits are recorded, never imputed. A method–repetition with any missing prediction is reported as incomplete.

## Selection of manuscript applications
The article reports at least five of these datasets, chosen after the complete eight-method analysis. Datasets where DSIIE has the lowest mean NRMSE or is competitive with the leaders are emphasised, with scientific variety and full disclosure.
- The full ledger for all 18 datasets is kept in the reproducibility materials.
- Coordinate figures: Plasma (primary) and possibly Bodyfat, per the revision brief.
