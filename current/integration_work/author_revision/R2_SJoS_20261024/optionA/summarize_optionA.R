# Summaries for Option A (PROTOCOL.md). Run from LaTeX/Working (or set DSIIE_WORK).
setwd(Sys.getenv("DSIIE_WORK", "."))
OUT <- "integration_work/author_revision/R2_SJoS_20261024/optionA"
x <- readRDS(file.path(OUT, "optionA_records.rds"))
SET <- c("SIIE_4", "SIIE_8", "DSIIE")
S <- do.call(rbind, lapply(SET, function(s) { d <- x[x$setting == s, ]; ok <- d$fit_status == "success"
  data.frame(setting = s, R = nrow(d), R_ok = sum(ok), fails = sum(!ok),
             r_ave = mean(d$r_ave[ok]), r_ave_mcse = sd(d$r_ave[ok]) / sqrt(sum(ok)),
             r_min = mean(d$r_min[ok]), r_min_mcse = sd(d$r_min[ok]) / sqrt(sum(ok))) }))
pair <- function(a, b, label) {
  A <- x[x$setting == a & x$fit_status == "success", ]; B <- x[x$setting == b & x$fit_status == "success", ]
  j <- merge(A[, c("data_seed", "r_ave", "r_min")], B[, c("data_seed", "r_ave", "r_min")],
             by = "data_seed", suffixes = c(".a", ".b"))
  do.call(rbind, lapply(c("r_ave", "r_min"), function(m) {
    del <- j[[paste0(m, ".a")]] - j[[paste0(m, ".b")]]; hw <- 1.96 * sd(del) / sqrt(length(del))
    data.frame(contrast = label, metric = m, R_pairs = length(del), mean_diff = mean(del), half_width = hw,
               lower = mean(del) - hw, upper = mean(del) + hw, n_first_larger = sum(del > 0),
               n_ties = sum(del == 0),
               verdict = if (mean(del) - hw > 0) "gain" else if (mean(del) + hw < 0) "loss" else "no detectable difference") }))
}
P <- rbind(pair("DSIIE", "SIIE_4", "DSIIE - SIIE (4 slices)"),
           pair("DSIIE", "SIIE_8", "DSIIE - SIIE (8 slices)"),
           pair("SIIE_8", "SIIE_4", "SIIE (8 slices) - SIIE (4 slices)"))
write.csv(S, file.path(OUT, "optionA_summary.csv"), row.names = FALSE)
write.csv(P, file.path(OUT, "optionA_paired.csv"), row.names = FALSE)
options(width = 200); print(S, digits = 4, row.names = FALSE); print(P, digits = 4, row.names = FALSE)
cat("failures:", sum(x$fit_status != "success"), "\n")
