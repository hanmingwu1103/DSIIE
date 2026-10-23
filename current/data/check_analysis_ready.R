# Validate exported inputs, original processed-file hashes, and source identities.
root <- normalizePath(Sys.getenv("DSIIE_WORK", "."), winslash = "/"); setwd(root)
meta <- read.csv("data/analysis_ready/sources.csv", stringsAsFactors = FALSE)
stopifnot(nrow(meta) == 18, !anyDuplicated(meta$id))
source("integration_work/author_revision/realdata/apps_data.R")
for (i in seq_len(nrow(meta))) {
  id <- meta$id[i]; rds <- file.path("data/analysis_ready", paste0(id, ".rds")); csv <- sub("[.]rds$", ".csv", rds)
  stopifnot(unname(tools::sha256sum(rds)) == meta$rds_sha256[i], unname(tools::sha256sum(csv)) == meta$csv_sha256[i])
  d <- readRDS(rds); z <- read.csv(csv, check.names = FALSE)
  stopifnot(nrow(d$X) == meta$n[i], ncol(d$X) == meta$p[i],
    identical(unname(as.matrix(z)), unname(cbind(y = d$y, d$X))), identical(APPS[[id]]$load(), d))
}
orig <- read.csv("data/sources.csv", stringsAsFactors = FALSE)
for (i in seq_len(nrow(orig))) stopifnot(unname(tools::sha256sum(file.path("data/processed", paste0(orig$id[i], ".csv")))) == orig$processed_sha256[i])
if ("--source" %in% commandArgs(trailingOnly = TRUE)) {
  Sys.setenv(DSIIE_USE_SOURCE_DATA = "1")
  source("integration_work/author_revision/realdata/apps_data.R")
  for (id in meta$id) stopifnot(identical(APPS[[id]]$load(), readRDS(file.path("data/analysis_ready", paste0(id, ".rds")))))
  cat("Original package/processed definitions: all 18 inputs exactly equal.\n")
}
cat("DATA CHECK PASSED: 18 datasets; 13 original processed-file hashes unchanged.\n")
