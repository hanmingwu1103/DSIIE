# Refit Option A, one simulation cell, and one application outer split.
# Work in a fresh copy; clear expected outputs and require child exit 0 before
# reading any new records. The caller's archived records are never overwritten.
src <- normalizePath(Sys.getenv("DSIIE_WORK", "."), winslash = "/"); setwd(src)
source("R/validation_helpers.R")
run <- copy_package(src, "dsiie_rerun_")
set_package_environment(run)
Sys.setenv(DSIIE_CORES = Sys.getenv("DSIIE_CORES", "4"))
A <- "integration_work/author_revision"
ok <- logical()

cat("== Option A rerun\n")
rel <- file.path(A, "R2_SJoS_20261024/optionA")
arch_A <- readRDS(file.path(src, rel, "optionA_records.rds"))
fresh <- file.path(run, rel, c("optionA_records.rds", "optionA_records.csv", "check_archive_reuse.csv", "check_matched_fine_objects.csv"))
clear_expected_outputs(fresh)
run_r_checked(file.path(run, rel, "run_optionA.R"), tail_lines = 8)
require_files(fresh)
new_A <- readRDS(fresh[1])
keys <- c("replication", "data_seed", "setting")
stopifnot(nrow(new_A) == 600, same_rows(arch_A, new_A, keys), all(new_A$fit_status == "success"))
j <- merge(arch_A, new_A, by = keys, suffixes = c(".a", ".n"))
dA <- max(abs(j$r_ave.a - j$r_ave.n), abs(j$r_min.a - j$r_min.n))
stopifnot(is.finite(dA), dA == 0)
ok["optionA"] <- TRUE
cat(sprintf("Option A: %d fits, maximum recovery difference = %.3g; fresh outputs verified\n", nrow(new_A), dA))

cat("\n== Baseline simulation cell rerun\n")
cell <- "cell_SIM-C_M7_n400_p10_H4_K2_a020_c005_v8.rds"
outdir <- file.path(run, "fresh_simulation_output")
stopifnot(!dir.exists(outdir)); dir.create(outdir)
Sys.setenv(DSIIE_OUT = outdir, DSIIE_SEED_ARCH = file.path(run, A, "simulation/results_v8"))
run_r_checked(file.path(run, A, "run_sim_v8.R"), "12", tail_lines = 5)
require_files(file.path(outdir, cell))
new_s <- readRDS(file.path(outdir, cell)); old_s <- readRDS(file.path(src, A, "simulation/results_v8", cell))
keys <- c("data_seed", "method_id")
stopifnot(nrow(new_s) == 1600, same_rows(old_s, new_s, keys), all(new_s$fit_status == "success"))
j <- merge(old_s, new_s, by = keys, suffixes = c(".a", ".n"))
# The archived one-direction M7 cell reports r_ave and deliberately leaves
# r_min missing. Require exactly the same missingness, then compare every
# defined measure; do not treat omitted redundant metrics as failed fits.
stopifnot(all(is.finite(j$r_ave.a)), all(is.finite(j$r_ave.n)),
  identical(is.na(j$r_min.a), is.na(j$r_min.n)), identical(is.nan(j$r_min.a), is.nan(j$r_min.n)))
dS <- max(c(0, abs(j$r_ave.a - j$r_ave.n), abs(j$r_min.a - j$r_min.n)), na.rm = TRUE)
stopifnot(is.finite(dS), dS == 0)
ok["simulation_cell"] <- TRUE
cat(sprintf("Simulation cell: %d fits, maximum recovery difference = %.3g\n", nrow(new_s), dS))
Sys.unsetenv(c("DSIIE_OUT", "DSIIE_SEED_ARCH"))

cat("\n== Mussels: repetition 1, fold 1; eight methods, nested tuning\n")
setwd(run)
L <- readLines(file.path(A, "realdata/run_apps_v8.R"), encoding = "UTF-8")
stop_at <- grep('^source\\("integration_work/author_revision/realdata/apps_data.R"\\)', L)[1]
stopifnot(is.finite(stop_at), stop_at > 1)
eval(parse(text = L[seq_len(stop_at - 1)], encoding = "UTF-8"), envir = globalenv())
worker_init(); x <- run_split(list(id = "Mussels", r = 1, k = 1))
old <- read.csv(file.path(src, A, "realdata/results_v8/oof_Mussels.csv"), stringsAsFactors = FALSE, colClasses = c(H_inner = "character"))
old <- old[old$repeat_ == 1 & old$fold == 1, ]
keys <- c("method", "row_id")
stopifnot(same_rows(old, x, keys), all(x$fit_status == "ok"), length(unique(x$method)) == 8)
m <- merge(old, x, by = keys, suffixes = c(".a", ".n"))
dP <- max(abs(m$yhat_primary.a - m$yhat_primary.n), abs(m$yhat_wide.a - m$yhat_wide.n))
stopifnot(is.finite(dP), dP <= 1e-10)
ok["application_split"] <- TRUE
cat(sprintf("Application: %d predictions, maximum prediction difference = %.3g\n", nrow(m), dP))
print(ok)
stopifnot(all(ok))
cat("RERUN CHECK PASSED\n")
