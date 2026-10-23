# Controlled loader for the method and generator definitions in the
# author-designated primary simulation script. Only PART 1--6 (utilities,
# DSIIE, comparators, metrics, generators, single-replication engine) are
# evaluated; directory creation and all run/plot drivers are skipped.
SIM_SRC <- Sys.getenv("DSIIE_SIM_SRC", "external/DSIIE_simulation_all.R")
local({
  L <- readLines(SIM_SRC, encoding = "UTF-8", warn = FALSE)
  start <- grep("^suppressPackageStartupMessages", L)[1]
  stop_ <- grep("PART 7", L)[1] - 2L
  L <- L[start:stop_]
  L <- L[!grepl("^dir\\.create", L)]
  eval(parse(text = L, encoding = "UTF-8"), envir = globalenv())
})
SIM_SRC_MD5 <- unname(tools::md5sum(SIM_SRC))
