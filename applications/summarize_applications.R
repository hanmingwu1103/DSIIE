# NRMSE summaries from the held-out (out-of-fold) records of run_apps_v8.R.
# NRMSE_r = sqrt( sum_i (y_i - yhat_i)^2 / sum_i (y_i - ybar)^2 ), computed per
# repetition over all n out-of-fold predictions (ybar = full-sample mean, used
# for evaluation only); reported as mean and SD over repetitions.
setwd(Sys.getenv("DSIIE_ROOT", getwd()))
IN <- "results/applications/oof"
OUTD <- "results/applications"
M8 <- c("SIR", "SAVE", "pHd", "DR", "SIME", "csMAVE", "SIIE", "DSIIE")
fs <- list.files(IN, pattern = "^oof_.*[.]csv([.]gz)?$")
summ <- list(); pair <- list(); ledger <- list(); tune <- list()
for (f in fs) {
  x <- read.csv(file.path(IN, f), stringsAsFactors = FALSE); id <- x$dataset[1]
  n <- length(unique(x$row_id)); R <- length(unique(x$repeat_))
  for (grid in c("primary", "wide")) {
    yh <- x[[paste0("yhat_", grid)]]
    per <- do.call(rbind, lapply(split(seq_len(nrow(x)), list(x$method, x$repeat_)), function(i) {
      d <- x[i, ]; complete <- all(d$fit_status == "ok") && length(unique(d$row_id)) == n &&
        nrow(d) == n && all(is.finite(yh[i]))
      data.frame(method = d$method[1], repeat_ = d$repeat_[1], complete = complete,
                 nrmse = if (complete) sqrt(sum((d$yobs - yh[i])^2) / sum((d$yobs - mean(d$yobs))^2)) else NA_real_)
    }))
    s <- do.call(rbind, lapply(M8, function(m) { q <- per[per$method == m, ]
      data.frame(dataset = id, grid = grid, method = m, repeats = nrow(q), complete = sum(q$complete),
                 nrmse = if (all(q$complete)) mean(q$nrmse) else NA_real_,
                 sd = if (all(q$complete) && nrow(q) > 1) sd(q$nrmse) else NA_real_) }))
    s$rank <- rank(s$nrmse, na.last = "keep", ties.method = "min"); summ[[length(summ) + 1]] <- s
    w <- reshape(per[, c("method", "repeat_", "nrmse")], idvar = "repeat_", timevar = "method", direction = "wide")
    pair[[length(pair) + 1]] <- do.call(rbind, lapply(setdiff(M8, "DSIIE"), function(m) {
      dd <- w[["nrmse.DSIIE"]] - w[[paste0("nrmse.", m)]]
      data.frame(dataset = id, grid = grid, comparator = m, mean_diff = mean(dd),
                 sd_diff = if (length(dd) > 1) sd(dd) else NA_real_, repeats_DSIIE_lower = sum(dd < 0), R = length(dd)) }))
  }
  first <- x[!duplicated(x[, c("method", "repeat_", "fold")]), ]
  tune[[id]] <- do.call(rbind, lapply(M8, function(m) { q <- first[first$method == m, ]
    data.frame(dataset = id, method = m, failed_splits = sum(q$fit_status != "ok"),
               H_values = paste(sort(unique(c(q$H_outer, unlist(strsplit(as.character(q$H_inner), ";"))))), collapse = ","),
               prim_gamma_low = mean(q$gamma_mult_primary == 0.25, na.rm = TRUE),
               prim_lambda_high = mean(q$lambda_primary == 1, na.rm = TRUE),
               wide_gamma_low = mean(q$gamma_mult_wide == 2^-6, na.rm = TRUE),
               wide_lambda_high = mean(q$lambda_wide == 100, na.rm = TRUE),
               median_sec = median(q$seconds)) }))
}
S <- do.call(rbind, summ); P <- do.call(rbind, pair); Tn <- do.call(rbind, tune)
write.csv(S, file.path(OUTD, "nrmse_summary_v8.csv"), row.names = FALSE)
write.csv(P, file.path(OUTD, "nrmse_paired_v8.csv"), row.names = FALSE)
write.csv(Tn, file.path(OUTD, "tuning_failures_v8.csv"), row.names = FALSE)
# screening ledger: one line per dataset
L <- do.call(rbind, lapply(split(S[S$grid == "primary", ], S$dataset[S$grid == "primary"]), function(s) {
  w <- S[S$grid == "wide" & S$dataset == s$dataset[1], ]
  data.frame(dataset = s$dataset[1], best_primary = s$method[which.min(s$nrmse)],
             DSIIE_rank_primary = s$rank[s$method == "DSIIE"], DSIIE_rank_wide = w$rank[w$method == "DSIIE"],
             DSIIE_nrmse = s$nrmse[s$method == "DSIIE"], best_nrmse = min(s$nrmse),
             gap_to_best = s$nrmse[s$method == "DSIIE"] - min(s$nrmse),
             SIIE_minus_DSIIE = s$nrmse[s$method == "SIIE"] - s$nrmse[s$method == "DSIIE"],
             DSIIE_sd = s$sd[s$method == "DSIIE"], incomplete = sum(s$complete < s$repeats)) }))
write.csv(L, file.path(OUTD, "screening_ledger_v8.csv"), row.names = FALSE)
options(width = 200); print(L, digits = 4, row.names = FALSE)
cat("total failed method-splits:", sum(Tn$failed_splits), "\n")
