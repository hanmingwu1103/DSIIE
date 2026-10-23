# Shared validation utilities. Failures always stop the calling R session.
package_root <- function() normalizePath(Sys.getenv("DSIIE_WORK", "."), winslash = "/", mustWork = TRUE)
rscript_path <- function() file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
copy_package <- function(src, prefix) {
  out <- tempfile(prefix); dir.create(out, recursive = TRUE)
  fs <- list.files(src, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  fs <- fs[!file.info(file.path(src, fs))$isdir]
  for (f in fs) {
    dest <- file.path(out, f)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file.path(src, f), dest, copy.date = TRUE)) stop("Copy failed: ", f)
  }
  normalizePath(out, winslash = "/", mustWork = TRUE)
}
set_package_environment <- function(root) {
  Sys.setenv(DSIIE_WORK = root, DSIIE_SIM_SRC = file.path(root, "external/DSIIE_simulation_all.R"),
    DSIIE_DATA = file.path(root, "data/processed"), DSIIE_ANALYSIS_READY = file.path(root, "data/analysis_ready"))
  Sys.unsetenv(c("DSIIE_OUT", "DSIIE_SEED_ARCH", "DSIIE_USE_SOURCE_DATA", "RD_TEST"))
}
run_r_checked <- function(script, args = character(), tail_lines = Inf) {
  output <- suppressWarnings(system2(rscript_path(), c("--vanilla", shQuote(script), args), stdout = TRUE, stderr = TRUE))
  code <- attr(output, "status"); if (is.null(code)) code <- 0L
  cat(if (is.finite(tail_lines)) tail(output, tail_lines) else output, sep = "\n")
  if (code != 0L) stop("Child R process failed (exit ", code, "): ", script, call. = FALSE)
  invisible(output)
}
require_files <- function(files) {
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Expected fresh output missing: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}
clear_expected_outputs <- function(files) {
  existing <- files[file.exists(files)]
  if (length(existing) && any(!file.remove(existing))) stop("Could not remove copied output before refit")
  if (any(file.exists(files))) stop("Stale outputs remain before refit")
  invisible(TRUE)
}
same_rows <- function(a, b, keys) {
  if (nrow(a) != nrow(b) || anyDuplicated(a[, keys, drop = FALSE]) || anyDuplicated(b[, keys, drop = FALSE])) return(FALSE)
  # Drop column names before do.call: an application key named "method" must
  # be an ordering vector, not bind to order()'s formal method argument.
  aa <- a[do.call(order, unname(as.list(a[keys]))), keys, drop = FALSE]
  bb <- b[do.call(order, unname(as.list(b[keys]))), keys, drop = FALSE]
  rownames(aa) <- rownames(bb) <- NULL
  isTRUE(all.equal(aa, bb, check.attributes = FALSE, tolerance = 0))
}
