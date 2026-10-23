# Validate a clean copy in fresh R processes. Any child failure exits nonzero.
src <- normalizePath(Sys.getenv("DSIIE_WORK", "."), winslash = "/"); setwd(src)
source("R/validation_helpers.R")
args <- commandArgs(trailingOnly = TRUE)
logfile <- if (length(args)) normalizePath(args[1], winslash = "/", mustWork = FALSE) else tempfile("validation_", fileext = ".log")
clean <- copy_package(src, "dsiie_clean_")
set_package_environment(clean)
txt <- c(paste("Version: DSIIE-analysis-R2-20261024-public1"), paste("Clean copy:", clean), R.version.string,
  paste("Validation started:", format(Sys.time(), tz = "UTC", usetz = TRUE)))
status <- logical()
for (ep in c("inspect_archives.R", "data/check_analysis_ready.R", "regenerate_outputs.R", "rerun_selected.R")) {
  t0 <- Sys.time()
  out <- suppressWarnings(system2(rscript_path(), c("--vanilla", shQuote(file.path(clean, ep))), stdout = TRUE, stderr = TRUE))
  code <- attr(out, "status"); if (is.null(code)) code <- 0L
  status[ep] <- code == 0L
  txt <- c(txt, sprintf("==== %s: exit %d; %.1f seconds", ep, code, as.numeric(difftime(Sys.time(), t0, units = "secs"))), out)
  cat(tail(txt, min(5, length(out) + 1)), sep = "\n")
  writeLines(txt, logfile)
  if (code != 0L) break
}
passed <- length(status) == 4L && all(status)
txt <- c(txt, if (passed) "CLEAN VALIDATION PASSED" else "CLEAN VALIDATION FAILED")
writeLines(txt, logfile)
cat(tail(txt, 1), "\nLog:", logfile, "\n")
if (!passed) quit(save = "no", status = 1L)
