# Inspect the archived results without fitting anything.
# Usage (from the package root): Rscript inspect_archives.R
# Checks every file against MANIFEST.sha256, then counts records, settings and
# failures in the simulation, diagnostic, candidate, Option A and application archives.
root <- normalizePath(Sys.getenv("DSIIE_WORK", "."), winslash = "/"); setwd(root)
A <- "integration_work/author_revision"
cat("Package root:", root, "\n\n")

# ---- manifest -------------------------------------------------------------------
man <- read.table("MANIFEST.sha256", header = FALSE, sep = "", quote = "", comment.char = "#",
                  col.names = c("sha256", "file"), stringsAsFactors = FALSE)
man$file <- sub("^[*]", "", man$file)
miss <- man$file[!file.exists(man$file)]
now <- tools::sha256sum(man$file[file.exists(man$file)])
bad <- names(now)[now != man$sha256[file.exists(man$file)]]
cat(sprintf("MANIFEST: %d files listed, %d missing, %d with a different SHA-256\n", nrow(man), length(miss), length(bad)))
if (length(miss)) print(miss); if (length(bad)) print(bad)
if (length(miss) || length(bad)) stop("Manifest verification failed", call. = FALSE)

# ---- baseline simulation (29 cells) -----------------------------------------------
fs <- list.files(file.path(A, "simulation/results_v8"), pattern = "_v8[.]rds$", full.names = TRUE)
sim <- do.call(rbind, lapply(fs, readRDS))
cat(sprintf("\nBaseline simulation: %d cells, %d replications, %d method fits, %d failures\n",
            length(fs), nrow(unique(sim[, c("design_id", "model_id", "n", "p", "H", "K", "alpha", "c", "data_seed")])),
            nrow(sim), sum(sim$fit_status != "success")))
print(table(sim$method_id))

# ---- R1 response-resolution diagnostic --------------------------------------------
dg <- readRDS(file.path(A, "R1/diagnostic/diagnostic_records.rds"))
cat(sprintf("\nR1 diagnostic: %d fits, %d failures; settings %s\n", nrow(dg), sum(dg$fit_status != "success"),
            paste(names(table(dg$setting)), table(dg$setting), collapse = ", ")))
mo <- read.csv(file.path(A, "R1/diagnostic/check_matched_objects.csv"))
cat("  maximum differences between matched fine-cell objects (SIIE with HK slices vs DSIIE):\n")
print(sapply(mo[grepl("^max_diff", names(mo))], max, na.rm = TRUE))

# ---- R1 candidate models --------------------------------------------------------------
for (st in c("explore", "confirm", "sens")) {
  x <- readRDS(file.path(A, sprintf("R1/candidates/records_%s.rds", st)))
  cat(sprintf("R1 candidates %-7s: %5d fits, %d failures, candidates %s\n", st, nrow(x),
              sum(x$fit_status != "success"), paste(unique(x$candidate), collapse = ",")))
}

# ---- R2 Option A ----------------------------------------------------------------------
oa <- readRDS(file.path(A, "R2_SJoS_20261024/optionA/optionA_records.rds"))
cat(sprintf("\nR2 Option A: %d fits, %d failures, seeds %d-%d\n", nrow(oa), sum(oa$fit_status != "success"),
            min(oa$data_seed), max(oa$data_seed)))
print(read.csv(file.path(A, "R2_SJoS_20261024/optionA/check_archive_reuse.csv")))

# ---- applications ---------------------------------------------------------------------
ff <- list.files(file.path(A, "realdata/results_v8"), pattern = "^oof_.*[.]csv$", full.names = TRUE)
st <- do.call(rbind, lapply(ff, function(f) { x <- read.csv(f, stringsAsFactors = FALSE)
  unique(x[, c("dataset", "method", "repeat_", "fold", "fit_status", "failure")]) }))
cat(sprintf("\nApplications: %d data sets, %d outer method-folds\n", length(ff), nrow(st)))
print(table(st$method, st$fit_status))
cat("Failure messages:\n"); print(table(st$failure[st$fit_status != "ok"]))
