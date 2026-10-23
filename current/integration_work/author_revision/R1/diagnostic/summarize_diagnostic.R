# Summaries and checks for the response-resolution diagnostic (PROTOCOL.md).
setwd(Sys.getenv("DSIIE_WORK", "."))
A <- "integration_work/author_revision"; OUT <- file.path(A, "R1/diagnostic")
x <- readRDS(file.path(OUT, "diagnostic_records.rds"))
idc <- read.csv(file.path(OUT, "identity_refits.csv"), stringsAsFactors = FALSE)
CELLS <- unique(x$cell)
arch <- do.call(rbind, lapply(CELLS, function(cl) {
  a <- readRDS(file.path(A, "simulation/results_v8", sprintf("cell_%s_v8.rds", cl))); a$cell <- cl
  a[, c("cell", "data_seed", "method_id", "r_ave")] }))

# 1. matched fine-cell objects
chk <- aggregate(cbind(max_diff_centre, max_diff_scatter, max_diff_whitening) ~ cell, x, max)
# 2. reuse of archived SIIE (H slices) and DSIIE, and identity of five further methods
key <- c(SIIE_H = "SIIE", DSIIE = "DSIIE")
re <- merge(transform(x[x$setting %in% names(key), c("cell", "data_seed", "setting", "r_ave")],
                      method_id = key[setting]), arch, by = c("cell", "data_seed", "method_id"))
ri <- merge(idc, arch, by = c("cell", "data_seed", "method_id"))
reuse <- rbind(
  aggregate(cbind(pairs = 1, max_abs_diff = abs(r_ave.x - r_ave.y)) ~ cell + method_id, re,
            function(v) if (all(v == 1)) length(v) else max(v)),
  aggregate(cbind(pairs = 1, max_abs_diff = abs(r_ave_refit - r_ave)) ~ cell + method_id, ri,
            function(v) if (all(v == 1)) length(v) else max(v)))
write.csv(chk, file.path(OUT, "check_matched_objects.csv"), row.names = FALSE)
write.csv(reuse, file.path(OUT, "check_reuse_identity.csv"), row.names = FALSE)

# 3. summaries
S <- do.call(rbind, lapply(split(x, list(x$cell, x$setting), drop = TRUE), function(d) {
  ok <- d$fit_status == "success"
  data.frame(cell = d$cell[1], setting = d$setting[1], R = nrow(d), fails = sum(!ok),
             r_ave = mean(d$r_ave[ok]), mcse = sd(d$r_ave[ok]) / sqrt(sum(ok))) }))
w <- reshape(x[, c("cell", "data_seed", "setting", "r_ave")], idvar = c("cell", "data_seed"),
             timevar = "setting", direction = "wide")
pd <- function(a, b, lab) do.call(rbind, lapply(split(w, w$cell), function(d) {
  del <- d[[paste0("r_ave.", a)]] - d[[paste0("r_ave.", b)]]; del <- del[is.finite(del)]
  data.frame(cell = d$cell[1], contrast = lab, R_pairs = length(del), mean_diff = mean(del),
             half_width = 1.96 * sd(del) / sqrt(length(del)), n_first_larger = sum(del > 0),
             n_ties = sum(del == 0)) }))
P <- rbind(pd("DSIIE", "SIIE_HK", "DSIIE - SIIE (HK slices)"),
           pd("DSIIE", "SIIE_H", "DSIIE - SIIE (H slices)"),
           pd("SIIE_HK", "SIIE_H", "SIIE (HK slices) - SIIE (H slices)"))
P$verdict <- ifelse(P$mean_diff - P$half_width > 0, "positive", ifelse(P$mean_diff + P$half_width < 0, "negative", "not detectable"))
write.csv(S, file.path(OUT, "diagnostic_summary.csv"), row.names = FALSE)
write.csv(P, file.path(OUT, "diagnostic_paired.csv"), row.names = FALSE)
options(width = 160); print(chk); print(reuse); print(S, digits = 4); print(P, digits = 4)
