# Summaries for the candidate models; Usage: Rscript summarize_candidates.R <stage>
setwd(Sys.getenv("DSIIE_WORK", "."))
OUT <- "integration_work/author_revision/R1/candidates"
stage <- commandArgs(trailingOnly = TRUE)[1]
x <- readRDS(file.path(OUT, sprintf("records_%s.rds", stage)))
M8 <- c("SIR", "SAVE", "pHd", "DR", "SIME", "csMAVE", "SIIE", "DSIIE")
S <- do.call(rbind, lapply(split(x, list(x$candidate, x$method_id), drop = TRUE), function(d) {
  ok <- d$fit_status == "success"
  data.frame(candidate = d$candidate[1], model = d$model_id[1], method = d$method_id[1], R = nrow(d),
             fails = sum(!ok), r_ave = mean(d$r_ave[ok]), mcse = sd(d$r_ave[ok]) / sqrt(sum(ok)),
             r_min = if (d$d[1] > 1) mean(d$r_min[ok]) else NA_real_,
             r_min_mcse = if (d$d[1] > 1) sd(d$r_min[ok]) / sqrt(sum(ok)) else NA_real_) }))
S <- S[order(S$candidate, match(S$method, M8)), ]
P <- do.call(rbind, lapply(split(x, x$candidate), function(d) {
  ds <- d[d$method_id == "DSIIE", ]
  do.call(rbind, lapply(setdiff(M8, "DSIIE"), function(m) {
    j <- merge(ds[, c("data_seed", "r_ave", "r_min")], d[d$method_id == m, c("data_seed", "r_ave", "r_min")],
               by = "data_seed", suffixes = c(".D", ".o"))
    f <- function(a, b) { del <- a - b; del <- del[is.finite(del)]
      if (!length(del)) return(c(NA, NA, NA)); c(mean(del), 1.96 * sd(del) / sqrt(length(del)), length(del)) }
    ra <- f(j$r_ave.D, j$r_ave.o); rm <- f(j$r_min.D, j$r_min.o)
    data.frame(candidate = d$candidate[1], comparator = m, R_pairs = ra[3], diff_r_ave = ra[1], hw_r_ave = ra[2],
               diff_r_min = rm[1], hw_r_min = rm[2]) })) }))
P$sig <- ifelse(P$diff_r_ave - P$hw_r_ave > 0, "DSIIE better", ifelse(P$diff_r_ave + P$hw_r_ave < 0, "DSIIE worse", "n.s."))
rule <- do.call(rbind, lapply(split(P, P$candidate), function(p) {
  s <- S[S$candidate == p$candidate[1], ]
  better <- p$comparator[p$sig == "DSIIE better"]; worse <- p$comparator[p$sig == "DSIIE worse"]
  rmin_ok <- if (all(is.na(s$r_min))) TRUE else rank(-s$r_min)[s$method == "DSIIE"] <= 2
  data.frame(candidate = p$candidate[1], n_better = length(better), better = paste(better, collapse = ","),
             n_worse = length(worse), worse = paste(worse, collapse = ","), SIIE_better = "SIIE" %in% better,
             rmin_top2 = rmin_ok,
             feature = length(better) >= 5 && "SIIE" %in% better && length(worse) <= 1 && rmin_ok) }))
write.csv(S, file.path(OUT, sprintf("summary_%s.csv", stage)), row.names = FALSE)
write.csv(P, file.path(OUT, sprintf("paired_%s.csv", stage)), row.names = FALSE)
write.csv(rule, file.path(OUT, sprintf("rule_%s.csv", stage)), row.names = FALSE)
options(width = 200); print(S, digits = 3, row.names = FALSE); print(P, digits = 3, row.names = FALSE); print(rule, row.names = FALSE)
cat("failures:", sum(x$fit_status != "success"), "\n")
