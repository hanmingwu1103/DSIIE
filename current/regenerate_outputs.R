# Regenerate every summary, table body and figure from the archived records,
# without refitting any estimator except the 16 descriptive full-data fits that
# the two application figures display. Usage (from the package root):
#   Rscript regenerate_outputs.R            # works in a temporary copy (default)
#   Rscript regenerate_outputs.R --in-place # overwrites the package's outputs
# Regenerated files are compared with the package copies: CSV files numerically
# (tolerance 1e-12) and byte-wise, table bodies byte-wise against their versioned
# release copies, and figures through their diagnostics and the
# equal-unit axis check (PDF bytes contain creation metadata and are not compared).
src <- normalizePath(Sys.getenv("DSIIE_WORK", "."), winslash = "/")
in_place <- "--in-place" %in% commandArgs(trailingOnly = TRUE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
if (in_place) run <- src else {
  run <- file.path(tempdir(), "dsiie_regenerate"); unlink(run, recursive = TRUE); dir.create(run)
  fs <- list.files(src, recursive = TRUE)
  for (f in fs) { dir.create(dirname(file.path(run, f)), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(src, f), file.path(run, f), copy.date = TRUE) }
}
cat("Source package:", src, "\nWorking copy:  ", run, "\n")
A <- "integration_work/author_revision"
steps <- list(
  c("simulation summaries", file.path(A, "simulation/summarize_v8.R")),
  c("application summaries", file.path(A, "realdata/summarize_apps_v8.R")),
  c("R1 diagnostic summaries", file.path(A, "R1/diagnostic/summarize_diagnostic.R")),
  c("R1 candidate summaries (explore)", file.path(A, "R1/candidates/summarize_candidates.R"), "explore"),
  c("R1 candidate summaries (confirm)", file.path(A, "R1/candidates/summarize_candidates.R"), "confirm"),
  c("R1 candidate summaries (sens)", file.path(A, "R1/candidates/summarize_candidates.R"), "sens"),
  c("Option A summaries", file.path(A, "R2_SJoS_20261024/optionA/summarize_optionA.R")),
  c("simulation table bodies", file.path(A, "tables/make_sim_tables_v8.R")),
  c("application table bodies", file.path(A, "tables/make_app_tables_v8.R")),
  c("simulation figure", file.path(A, "figures/sim_figure_v8.R")),
  c("application figures (16 descriptive fits)", file.path(A, "figures/coord_figures_v8.R")))
Sys.setenv(DSIIE_WORK = run, DSIIE_SIM_SRC = file.path(run, "external/DSIIE_simulation_all.R"),
           DSIIE_DATA = file.path(run, "data/processed"),
           DSIIE_ANALYSIS_READY = file.path(run, "data/analysis_ready"))
Sys.unsetenv("DSIIE_USE_SOURCE_DATA")
status <- character(0)
for (s in steps) {
  t0 <- Sys.time()
  out <- system2(rscript, c(shQuote(file.path(run, s[2])), if (length(s) > 2) s[3]), stdout = TRUE, stderr = TRUE)
  ok <- is.null(attr(out, "status")) || attr(out, "status") == 0
  status[s[1]] <- if (ok) "ok" else "FAILED"
  cat(sprintf("%-45s %s (%.1f s)\n", s[1], status[s[1]], as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (!ok) cat(tail(out, 15), sep = "\n")
}

# ---- comparisons ---------------------------------------------------------------------
cmp_csv <- function(f) {
  a <- read.csv(file.path(src, f), stringsAsFactors = FALSE); b <- read.csv(file.path(run, f), stringsAsFactors = FALSE)
  bytes <- unname(tools::md5sum(file.path(src, f)) == tools::md5sum(file.path(run, f)))
  same <- identical(dim(a), dim(b)) && identical(names(a), names(b)) && all(vapply(seq_along(a), function(j) {
    x <- a[[j]]; y <- b[[j]]
    if (is.numeric(x) && is.numeric(y)) isTRUE(all(is.na(x) == is.na(y))) && isTRUE(max(c(0, abs(x - y)), na.rm = TRUE) <= 1e-12)
    else identical(as.character(x), as.character(y)) }, TRUE))
  data.frame(file = f, byte_identical = bytes, numerically_equal = same)
}
csvs <- c(file.path(A, "simulation", c("summary_v8.csv", "paired_v8.csv", "reproduction_vs_corrected.csv")),
          file.path(A, "realdata", c("nrmse_summary_v8.csv", "nrmse_paired_v8.csv", "tuning_failures_v8.csv",
                                     "screening_ledger_v8.csv")),
          file.path(A, "R1/diagnostic", c("diagnostic_summary.csv", "diagnostic_paired.csv",
                                          "check_matched_objects.csv", "check_reuse_identity.csv")),
          file.path(A, "R1/candidates", paste0(c(outer(c("summary", "paired", "rule"),
                                                       c("explore", "confirm", "sens"), paste, sep = "_")), ".csv")),
          file.path(A, "R2_SJoS_20261024/optionA", c("optionA_summary.csv", "optionA_paired.csv")),
          file.path(A, "figures", c("coord_diagnostics.csv", "coord_axis_scale.csv")))
stopifnot(all(file.exists(file.path(src, csvs))))
res <- do.call(rbind, lapply(csvs, cmp_csv))
tex <- list.files(file.path(src, A, "tables"), pattern = "[.]tex$")
tab <- do.call(rbind, lapply(tex, function(f) {
  a <- readLines(file.path(src, A, "tables", f), warn = FALSE); b <- readLines(file.path(run, A, "tables", f), warn = FALSE)
  data.frame(file = file.path(A, "tables", f), byte_identical = identical(a, b)) }))
ax <- read.csv(file.path(run, A, "figures/coord_axis_scale.csv"))
cat("\nCSV comparisons:\n"); print(res, row.names = FALSE)
cat("\nTable bodies compared with versioned release copies; R1 resolution body is retained for provenance:\n"); print(tab, row.names = FALSE)
cat(sprintf("\nFigure axis check: %d panels, max |x/y units - 1| = %.3g\n", nrow(ax),
            max(abs(ax$x_inches_per_unit / ax$y_inches_per_unit - 1))))
cat("Regenerated figure files:", paste(list.files(file.path(run, "DSIIE_main_assets"), pattern = "[.]pdf$"), collapse = ", "), "\n")
ok <- all(status == "ok") && all(res$numerically_equal) && all(tab$byte_identical) &&
  all(res$byte_identical) && nrow(ax) == 16L &&
  max(abs(ax$x_inches_per_unit / ax$y_inches_per_unit - 1)) < 1e-10
cat("\nREGENERATION CHECK:", if (ok) "PASSED" else "FAILED", "\n")
if (!ok) quit(status = 1)
