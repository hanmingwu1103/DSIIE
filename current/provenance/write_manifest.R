# Rebuild the manifest after an intentional release edit, from current/.
root <- normalizePath(Sys.getenv("DSIIE_WORK", "."), winslash = "/"); setwd(root)
fs <- sort(list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE))
fs <- fs[!file.info(fs)$isdir & fs != "MANIFEST.sha256"]
stopifnot(!any(grepl("\\s", fs)))
hashes <- tools::sha256sum(fs)
writeLines(c(sprintf("# DSIIE-analysis-R2-20261024-public1: SHA-256 of %d files", length(fs)),
  sprintf("%s  %s", unname(hashes), fs)), "MANIFEST.sha256")
cat(length(fs), "files recorded in MANIFEST.sha256\n")
