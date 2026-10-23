# Clean-copy reproduction check, run from the repository root:
#   Rscript checks/reproduce_subset.R [ncores]
# Recomputes one simulation cell and one application data set into a temporary
# directory and compares them with the archived results; then regenerates the
# summaries and table bodies from the archived results and compares them with
# the archived copies. Prints the largest absolute discrepancies.
args <- commandArgs(trailingOnly = TRUE)
nc <- if (length(args)) args[1] else "4"
ROOT <- getwd(); TMP <- file.path(tempdir(), "dsiie_repro"); unlink(TMP, recursive = TRUE)
dir.create(TMP); file.copy(list.files(ROOT, full.names = TRUE), TMP, recursive = TRUE)
setwd(TMP); Rs <- file.path(R.home("bin"), "Rscript")
run <- function(...) { st <- system2(Rs, c(...), stdout = TRUE, stderr = TRUE); cat(tail(st, 3), sep = "\n") }
cmp_num <- function(a, b) { k <- intersect(names(a)[sapply(a, is.numeric)], names(b))
  k <- setdiff(k, c("seconds", "sec", "sec_median", "median_sec"))
  max(sapply(k, function(v) max(abs(a[[v]] - b[[v]]), na.rm = TRUE))) }

# 1. one simulation cell
cell <- "SIM-C_M7_n400_p10_H4_K2_a020_c005"
old <- readRDS(file.path(ROOT, "results/simulation", sprintf("cell_%s_v8.rds", cell)))
run("simulation/run_simulation.R", nc, cell)
new <- readRDS(file.path("results/simulation", sprintf("cell_%s_v8.rds", cell)))
cat(sprintf("simulation cell %s: %d records, max |difference| = %.3g\n", cell, nrow(new), cmp_num(old, new)))

# 2. one application data set (all eight methods, 10 x 10 folds)
Sys.setenv(DSIIE_CORES = nc)
run("applications/run_applications.R", "Mussels")
a <- read.csv(file.path(ROOT, "results/applications/oof/oof_Mussels.csv.gz"))
b <- read.csv("results/applications/oof/oof_Mussels.csv.gz")
key <- function(x) paste(x$method, x$repeat_, x$fold, x$row_id); b <- b[match(key(a), key(b)), ]
cat(sprintf("application Mussels: %d predictions, max |difference| = %.3g\n", nrow(b),
            max(abs(c(a$yhat_primary - b$yhat_primary, a$yhat_wide - b$yhat_wide)))))

# 3. summaries and tables from the archived results
file.copy(file.path(ROOT, "results/applications/oof/oof_Mussels.csv.gz"), "results/applications/oof", overwrite = TRUE)
run("simulation/summarize_simulation.R"); run("applications/summarize_applications.R")
for (f in c("results/simulation/summary_v8.csv", "results/simulation/paired_v8.csv",
            "results/applications/nrmse_summary_v8.csv", "results/applications/nrmse_paired_v8.csv"))
  cat(sprintf("%s: max |difference| = %.3g\n", f, cmp_num(read.csv(file.path(ROOT, f)), read.csv(f))))
unlink(list.files("tables", pattern = "[.]tex$", full.names = TRUE))   # must be regenerated
run("tables/make_sim_tables.R"); run("tables/make_app_tables.R")
for (f in list.files(file.path(ROOT, "tables"), pattern = "[.]tex$"))
  cat(sprintf("tables/%s identical: %s\n", f, identical(readLines(file.path(ROOT, "tables", f)), readLines(file.path("tables", f)))))
setwd(ROOT)
