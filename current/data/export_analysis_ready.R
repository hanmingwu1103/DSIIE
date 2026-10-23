# Recreate the exact 18 analysis inputs from the frozen source definitions.
# RDS files retain binary-double precision and attributes. CSV files use 17
# significant digits for an exact numeric round-trip, independent of write.table.
root <- normalizePath(Sys.getenv("DSIIE_WORK", "."), winslash = "/"); setwd(root)
Sys.setenv(DSIIE_USE_SOURCE_DATA = "1")
source("integration_work/author_revision/realdata/apps_data.R")
dir.create("data/analysis_ready", recursive = TRUE, showWarnings = FALSE)
old <- read.csv("data/sources.csv", stringsAsFactors = FALSE)
extra <- list(
  RD06 = c("gamlss.data", "plasma", "Raw betaplasma; age,bmi,calories,fat,fiber,alcohol,cholesterol,betadiet,retdiet."),
  RD09 = c("faraway", "fat", "Brozek body-fat percentage; 13 measurements age,weight,height,neck,chest,abdom,hip,thigh,knee,ankle,biceps,forearm,wrist."),
  Mussels = c("dr", "mussels", "Keep complete positive H,L,W,S,M; log H,L,W,S predictors and log M response."),
  Hitters = c("ISLR2", "Hitters", "Exclude missing Salary; 16 numeric predictors; log Salary response."),
  College = c("ISLR2", "College", "Grad.Rate response; numeric predictors excluding Private and Grad.Rate."))
meta <- list(); checks <- list()
for (id in names(APPS)) {
  d <- APPS[[id]]$load(); stopifnot(is.matrix(d$X), nrow(d$X) == length(d$y), all(is.finite(d$X)), all(is.finite(d$y)))
  rds <- file.path("data/analysis_ready", paste0(id, ".rds")); csv <- sub("[.]rds$", ".csv", rds)
  saveRDS(d, rds, version = 3)
  z <- cbind(y = d$y, d$X)
  header <- paste(sprintf('"%s"', gsub('"', '""', colnames(z), fixed = TRUE)), collapse = ",")
  rows <- apply(z, 1, function(r) paste(sprintf("%.17g", r), collapse = ","))
  writeLines(c(header, rows), csv, useBytes = TRUE)
  back <- read.csv(csv, check.names = FALSE)
  stopifnot(identical(readRDS(rds), d), identical(unname(as.matrix(back)), unname(z)))
  if (id %in% names(extra)) {
    x <- extra[[id]]; pkg <- packageDescription(x[1])
    source <- paste0(x[1], "::", x[2]); version <- pkg$Version; licence <- pkg$License; transform <- x[3]
    cite <- paste(capture.output(print(citation(x[1]))), collapse = " ")
  } else {
    s <- old[old$id == id, ]; stopifnot(nrow(s) == 1)
    source <- s$obtain_source; version <- s$version; licence <- s$license; cite <- s$original_citation
    transform <- "Exact copy of the predecessor processed input, read without new preprocessing; response is y."
  }
  meta[[id]] <- data.frame(id, title = APPS[[id]]$title, response = APPS[[id]]$resp, fold_design = APPS[[id]]$design,
    n = nrow(d$X), p = ncol(d$X), source, source_version = version, source_licence = licence, citation = cite,
    preprocessing = transform, predictors = paste(colnames(d$X), collapse = ";"),
    csv_sha256 = unname(tools::sha256sum(csv)), rds_sha256 = unname(tools::sha256sum(rds)), stringsAsFactors = FALSE)
  checks[[id]] <- data.frame(id, n = nrow(d$X), p = ncol(d$X), rds_exact = TRUE, csv_exact = TRUE)
}
write.csv(do.call(rbind, meta), "data/analysis_ready/sources.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, checks), "data/analysis_ready/identity_checks.csv", row.names = FALSE)
cat("Exported and exactly verified", length(meta), "analysis-ready datasets in RDS and CSV form.\n")
print(do.call(rbind, checks), row.names = FALSE)
