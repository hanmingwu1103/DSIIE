# Additional simulation candidates: protocol and ledger rules

Written 2026-10-23T12:15Z, before any candidate fit. This is a focused, finite
candidate set chosen for a statistical reason connected to DSIIE. It is not a
search over many models.

## Candidates (all three, no others)

| ID | Model | d | Reason for candidacy | n | p | H | K | alpha |
|---|---|---:|---|---:|---:|---:|---:|---:|
| C1 | M5: historical model (Xia 2007, Model I), X ~ N(0, 0.5^|i-j|), Y = 2(b1'X)^2 + exp(b2'X) e | 2 | Quadratic location plus heteroscedastic scale; both directions act partly through within-slice dispersion. Earlier favourable student numbers are not evidence. | 400 | 10 | 4 | 2 | .20 |
| C2 | TAQ3: historical quantile-width design (G uniform on 1..5, Y = G + U(-.2,.2); A1, A2 have a body of half-width a_G and 10% tail mass at +-[2.9,3]) | 2 | The signal is in interior directional quantile widths with a fixed envelope, the quantity DSIIE summarizes. | 400 | 10 | 5 | 2 | .20 |
| C3 | M2t (new): M2 with nuisance coordinates Z2..Zp independent t_3/sqrt(3); Y = Z1^2; X = R Z | 1 | Symmetric dispersion signal with heavy-tailed nuisance coordinates; interquantile widths are less sensitive to heavy tails than second moments. Outside Condition 1 (complement not spherical); central subspace span(R e1). | 400 | 10 | 4 | 2 | .20 |

C2 uses H = 5 because its response has five latent levels; this is fixed here and
not varied.

## Methods and metrics

All eight v8 methods through `one_rep_v8()` (`methods_v8.R`), same sample for all
methods within a replication; d is supplied. Metrics r_ave and, for d = 2, r_min.
Paired differences DSIIE minus each comparator with 1.96 sd/sqrt(R) half-widths,
for r_ave and (d = 2) r_min. Failures recorded by `safe_fit`, never retried.

## Stages

1. Exploration: R = 50, seeds 9110001.. (C1), 9120001.. (C2), 9130001.. (C3).
2. Confirmatory: R = 200 on fresh seeds 9210001.. (C1), 9220001.. (C2),
   9230001.. (C3), run for all three candidates whatever the exploration showed,
   with no change of design, tuning or methods. Only confirmatory results may be
   reported as evidence.

## Rule for featuring a candidate in the article (fixed now)

Feature a candidate only if, in the confirmatory run, DSIIE's paired 95% interval
for r_ave lies entirely above zero against at least five of the seven comparators
including SIIE, no comparator exceeds DSIIE by an interval entirely above zero
in more than one case, and, for d = 2, DSIIE's mean r_min is within the two
largest. State exactly which comparators are improved upon and which are not.
Candidates that do not meet the rule stay in the ledger and are not featured.

## Sensitivity stage for C1 (added 2026-10-23T12:25Z, after the confirmatory run, before these fits)

Only C1 met the rule. Sensitivity settings, all eight methods, R = 200 fresh seeds
each, all other choices as in the confirmatory run: (n, p) = (200, 10) seeds
9310001..; (800, 10) seeds 9320001..; (400, 20) seeds 9330001... All three are
reported whatever their outcome; they qualify, but cannot by themselves create,
the featured claim.
