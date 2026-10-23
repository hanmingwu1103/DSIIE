# Summaries of the corrected simulation records: replication counts, failures,
# mean recovery with Monte Carlo standard errors, weakest-direction recovery,
# common-seed paired differences, and reproduction of archived values.
setwd(Sys.getenv("DSIIE_WORK", "."))
IN <- "integration_work/author_revision/simulation/results_v8"; PREV <- "integration_work/simulation/results_corrected"
fs <- list.files(IN, pattern = "_v8[.]rds$")
recs <- do.call(rbind, lapply(fs, function(f) { x <- readRDS(file.path(IN, f)); x$file <- f; x }))
recs$cell <- sub("_v8[.]rds$", "", sub("^cell_", "", recs$file))

summ <- do.call(rbind, lapply(split(recs, list(recs$cell, recs$method_id), drop = TRUE), function(d) {
  ok <- d$fit_status == "success"
  data.frame(cell = d$cell[1], model = d$model_id[1], n = d$n[1], p = d$p[1], H = d$H[1],
             K = d$K[1], alpha = d$alpha[1], c = d$c[1], method = d$method_id[1],
             R = nrow(d), fails = sum(!ok),
             r_ave = mean(d$r_ave[ok]), sd = sd(d$r_ave[ok]), mcse = sd(d$r_ave[ok]) / sqrt(sum(ok)),
             r_min = if (d$d[1] > 1) mean(d$r_min[ok]) else NA_real_,
             r_min_mcse = if (d$d[1] > 1) sd(d$r_min[ok]) / sqrt(sum(ok)) else NA_real_,
             sec_median = median(d$seconds))
}))
summ <- summ[order(summ$cell, summ$method), ]
write.csv(summ, "integration_work/author_revision/simulation/summary_v8.csv", row.names = FALSE)

paired <- do.call(rbind, lapply(split(recs, recs$cell), function(d) {
  ds <- d[d$method_id == "DSIIE", c("data_seed", "r_ave", "r_min")]
  do.call(rbind, lapply(setdiff(unique(d$method_id), "DSIIE"), function(m) {
    o <- d[d$method_id == m, c("data_seed", "r_ave", "r_min")]
    j <- merge(ds, o, by = "data_seed", suffixes = c(".D", ".o"))
    j <- j[is.finite(j$r_ave.D) & is.finite(j$r_ave.o), ]
    del <- j$r_ave.D - j$r_ave.o
    data.frame(cell = d$cell[1], comparator = m, R_pairs = nrow(j), mean_diff = mean(del),
               half_width = 1.96 * sd(del) / sqrt(nrow(j)), prop_DSIIE_higher = mean(del > 0))
  }))
}))
write.csv(paired, "integration_work/author_revision/simulation/paired_v8.csv", row.names = FALSE)

# reproduction against the previous corrected run (same seeds): methods whose
# definition is unchanged must agree; SIIE was labelled SIIE-Q(H) there.
map <- c(SIR = "SIR", SAVE = "SAVE", pHd = "pHd", DR = "DR", csMAVE = "csMAVE", SIIE = "SIIE-Q(H)", DSIIE = "DSIIE")
repro <- do.call(rbind, lapply(fs, function(f) {
  x <- readRDS(file.path(IN, f)); pf <- file.path(PREV, sub("_v8", "_corrected", f))
  if (!file.exists(pf)) return(NULL); a <- readRDS(pf)
  do.call(rbind, lapply(names(map), function(m) {
    u <- merge(x[x$method_id == m, c("data_seed", "r_ave")], a[a$method_id == map[[m]], c("data_seed", "r_ave")], by = "data_seed")
    if (!nrow(u)) return(NULL)
    data.frame(file = f, method = m, pairs = nrow(u), max_abs_diff = max(abs(u$r_ave.x - u$r_ave.y)))
  }))
}))
write.csv(repro, "integration_work/author_revision/simulation/reproduction_vs_corrected.csv", row.names = FALSE)
cat("Total method fits:", nrow(recs), " failures:", sum(recs$fit_status != "success"), "\n")
print(aggregate(max_abs_diff ~ method, repro, max))
print(table(recs$method_id, recs$fit_status))
