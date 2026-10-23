# Negative controls: deliberately damaged disposable copies must fail.
# No fits are run. The original files and archived results remain untouched.
root <- normalizePath(Sys.getenv("DSIIE_WORK", "."), winslash = "/"); setwd(root)
source("R/validation_helpers.R")
invoke <- function(dir, script) {
  set_package_environment(dir)
  out <- suppressWarnings(system2(rscript_path(), c("--vanilla", shQuote(file.path(dir, script))), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status"); if (is.null(status)) status <- 0L
  list(status = status, output = out)
}
checks <- logical()
bad_manifest <- copy_package(root, "dsiie_bad_manifest_")
cat("\n# Deliberate validation negative control\n", file = file.path(bad_manifest, "external/DSIIE_simulation_all.R"), append = TRUE)
r <- invoke(bad_manifest, "inspect_archives.R")
checks["modified_manifest_file_fails"] <- r$status != 0L && any(grepl("Manifest verification failed", r$output, fixed = TRUE))
cat("Altered file: child exit", r$status, "\n")

missing_manifest <- copy_package(root, "dsiie_missing_manifest_")
stopifnot(file.remove(file.path(missing_manifest, "external/DSIIE_simulation_all.R")))
r <- invoke(missing_manifest, "inspect_archives.R")
checks["missing_manifest_file_fails"] <- r$status != 0L && any(grepl("Manifest verification failed", r$output, fixed = TRUE))
cat("Missing listed file: child exit", r$status, "\n")

bad_refit <- copy_package(root, "dsiie_bad_refit_")
option <- "integration_work/author_revision/R2_SJoS_20261024/optionA"
stopifnot(file.exists(file.path(bad_refit, option, "optionA_records.rds")))
writeLines('stop("Deliberate Option A child failure")', file.path(bad_refit, option, "run_optionA.R"))
r <- invoke(bad_refit, "rerun_selected.R")
checks["failed_refit_rejects_stale_archive"] <- r$status != 0L && any(grepl("Child R process failed", r$output, fixed = TRUE)) &&
  !any(grepl("RERUN CHECK PASSED", r$output, fixed = TRUE))
cat("Failing Option A child with copied archive present: parent exit", r$status, "\n")

empty_refit <- copy_package(root, "dsiie_empty_refit_")
writeLines('quit(save = "no", status = 0L)', file.path(empty_refit, option, "run_optionA.R"))
r <- invoke(empty_refit, "rerun_selected.R")
checks["successful_child_without_output_fails"] <- r$status != 0L && any(grepl("Expected fresh output missing", r$output, fixed = TRUE)) &&
  !any(grepl("RERUN CHECK PASSED", r$output, fixed = TRUE))
cat("Option A child exit 0 without outputs: parent exit", r$status, "\n")

bad_wrapper <- copy_package(root, "dsiie_bad_wrapper_")
writeLines('stop("Deliberate archive-check child failure")', file.path(bad_wrapper, "inspect_archives.R"))
r <- invoke(bad_wrapper, "validate_package.R")
checks["wrapper_propagates_child_failure"] <- r$status != 0L && any(grepl("CLEAN VALIDATION FAILED", r$output, fixed = TRUE))
cat("Failing first validation child: wrapper exit", r$status, "\n")
print(checks)
if (!all(checks)) stop("A negative control did not fail as required")
cat("NEGATIVE CONTROLS PASSED\n")
