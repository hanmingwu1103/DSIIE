# Frozen application definitions for the author revision (fixed 23 Oct 2026,
# before any eight-method application result was computed).
#
# Sources: the corrected predecessor registry (SIIE-Rcode RD1R processed files,
# with recorded SHA-256 checksums) and the student-only datasets of the earlier
# exploratory report. Plasma and Bodyfat keep the definitions of the corrected
# DSIIE analysis (raw betaplasma; 13 named body measurements). The fold design
# is "random" (exchangeable units), "grouped" (repeated units kept intact in
# outer and inner folds) or "blocked" (contiguous time blocks).
RD1R <- Sys.getenv("DSIIE_DATA", "data/processed")

rd1r <- function(id) {
  df <- read.csv(file.path(RD1R, paste0(id, ".csv")), check.names = FALSE)
  list(X = as.matrix(df[, setdiff(names(df), "y")]), y = df$y)
}
canon_key <- function(M) apply(M, 1, function(r) paste(formatC(r, digits = 16, format = "g"), collapse = "|"))

APPS <- list(
  RD01 = list(title = "Los Angeles ozone", resp = "daily maximum ozone (raw)", design = "blocked",
              load = function() rd1r("RD01")),
  RD02 = list(title = "AIS male lean body mass", resp = "lean body mass (kg)", design = "random",
              load = function() rd1r("RD02")),
  RD06 = list(title = "Plasma beta-carotene", resp = "plasma beta-carotene (ng/ml, raw)", design = "random",
              load = function() {
                e <- new.env(); data("plasma", package = "gamlss.data", envir = e); d <- e$plasma
                xn <- c("age", "bmi", "calories", "fat", "fiber", "alcohol", "cholesterol", "betadiet", "retdiet")
                list(X = as.matrix(d[, xn]), y = as.numeric(d$betaplasma)) }),
  RD07 = list(title = "Automobile price", resp = "log price", design = "random",
              load = function() rd1r("RD07")),
  RD08 = list(title = "Prostate cancer", resp = "log PSA", design = "random",
              load = function() rd1r("RD08")),
  RD09 = list(title = "Body fat", resp = "percent body fat (Brozek)", design = "random",
              load = function() {
                e <- new.env(); data("fat", package = "faraway", envir = e); d <- e$fat
                xn <- c("age", "weight", "height", "neck", "chest", "abdom", "hip", "thigh",
                        "knee", "ankle", "biceps", "forearm", "wrist")
                list(X = as.matrix(d[, xn]), y = d$brozek) }),
  RD10 = list(title = "Biomass heating value", resp = "higher heating value", design = "random",
              load = function() rd1r("RD10")),
  RD11 = list(title = "Concrete compressive strength", resp = "compressive strength (MPa)", design = "grouped",
              load = function() rd1r("RD11"), group = function(X) canon_key(X[, 1:7])),
  RD13 = list(title = "Airfoil self-noise", resp = "scaled sound pressure (dB)", design = "grouped",
              load = function() rd1r("RD13"), group = function(X) canon_key(X[, c("angle", "chord", "velocity", "thickness")])),
  RD14 = list(title = "Yacht hydrodynamics", resp = "residuary resistance", design = "grouped",
              load = function() rd1r("RD14"), group = function(X) canon_key(X[, 1:5])),
  RD15 = list(title = "QSAR aquatic toxicity", resp = "LC50 (-log mol/L)", design = "grouped",
              load = function() rd1r("RD15"), group = function(X) canon_key(X)),
  RD16 = list(title = "Real-estate valuation", resp = "price per unit area", design = "grouped",
              load = function() rd1r("RD16"), group = function(X) canon_key(X[, c("lat", "long")])),
  RD17 = list(title = "Concrete slump test", resp = "28-day compressive strength (MPa)", design = "random",
              load = function() rd1r("RD17")),
  RD19 = list(title = "QSAR fish toxicity", resp = "LC50 (-log mol/L)", design = "grouped",
              load = function() rd1r("RD19"), group = function(X) canon_key(X)),
  RD20 = list(title = "Computer hardware", resp = "log relative performance", design = "random",
              load = function() rd1r("RD20")),
  Mussels = list(title = "Horse mussels", resp = "log muscle mass", design = "random",
              load = function() {
                e <- new.env(); data("mussels", package = "dr", envir = e); d <- e$mussels
                d <- d[complete.cases(d) & apply(d[, c("H", "L", "W", "S", "M")] > 0, 1, all), ]
                list(X = as.matrix(log(d[, c("H", "L", "W", "S")])), y = log(d$M)) }),
  Hitters = list(title = "Baseball salaries", resp = "log salary", design = "random",
              load = function() {
                e <- new.env(); data("Hitters", package = "ISLR2", envir = e); d <- e$Hitters
                d <- d[!is.na(d$Salary), ]
                xn <- c("AtBat", "Hits", "HmRun", "Runs", "RBI", "Walks", "Years", "CAtBat", "CHits",
                        "CHmRun", "CRuns", "CRBI", "CWalks", "PutOuts", "Assists", "Errors")
                list(X = as.matrix(d[, xn]), y = log(d$Salary)) }),
  College = list(title = "US colleges", resp = "graduation rate (%)", design = "random",
              load = function() {
                e <- new.env(); data("College", package = "ISLR2", envir = e); d <- e$College
                xn <- setdiff(names(d), c("Private", "Grad.Rate"))
                list(X = as.matrix(d[, xn]), y = d$Grad.Rate) })
)

## ---- fold designs ------------------------------------------------------------
make_folds_strat <- function(y, K_fold = 10, seed = 1) {        # corrected-analysis rule
  set.seed(seed); n <- length(y); ord <- order(y); fold <- integer(n)
  grp <- rep_len(seq_len(K_fold), n)
  fold[ord] <- ave(grp, ceiling(seq_along(grp) / K_fold), FUN = sample)
  fold
}
grouped_folds <- function(g, K, seed) {   # random group order, greedy size balancing
  set.seed(seed); ug <- unique(g); ug <- ug[sample.int(length(ug))]
  sz <- table(g)[ug]; f <- integer(0); load <- numeric(K)
  for (i in seq_along(ug)) { k <- which.min(load); f[i] <- k; load[k] <- load[k] + sz[i] }
  unname(setNames(f, ug)[g])
}
outer_folds <- function(app, dat, r) {
  n <- length(dat$y)
  switch(app$design,
    random  = make_folds_strat(dat$y, 10, 1000 * r + 7),
    grouped = grouped_folds(app$group(dat$X), 10, 2000 * r + 11),
    blocked = rank_split(n, 10))                          # rows are in time order
}
inner_folds <- function(app, dat, tr, r, k) {
  m <- length(tr)
  switch(app$design,
    random  = { set.seed(100000L + 100L * r + k); sample(rep(seq_len(5), length.out = m)) },
    grouped = grouped_folds(app$group(dat$X[tr, , drop = FALSE]), 5, 100000L + 100L * r + k),
    blocked = rank_split(m, 5))                           # contiguous blocks of the training rows
}
n_repeats <- function(app) if (app$design == "blocked") 1L else 10L

# Public release: prefer the exact, versioned analysis inputs bundled with the
# repository. DSIIE_USE_SOURCE_DATA=1 restores the original definitions above
# for export/identity verification. RDS preserves the original binary doubles.
ANALYSIS_READY <- Sys.getenv("DSIIE_ANALYSIS_READY", "data/analysis_ready")
if (Sys.getenv("DSIIE_USE_SOURCE_DATA", "0") != "1") {
  for (id in names(APPS)) {
    cached <- file.path(ANALYSIS_READY, paste0(id, ".rds"))
    if (file.exists(cached)) APPS[[id]]$load <- local({ path <- cached; function() readRDS(path) })
  }
}
