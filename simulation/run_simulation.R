# Eight-method simulation study (SIR, SAVE, pHd, DR, SIME, csMAVE, SIIE, DSIIE).
# Every replication regenerates its sample from the data seed in
# design_seeds.csv, so all methods see identical samples.
#
# Usage, from the repository root:
#   Rscript simulation/run_simulation.R [ncores] [cell ...]
# Without cell names all 29 cells are run (about 20 minutes on 30 cores).
# Output: results/simulation/cell_<cell>_v8.rds (replication-level records);
# existing files are overwritten.
args <- commandArgs(trailingOnly = TRUE)
ncores <- if (length(args)) as.integer(args[1]) else max(1L, parallel::detectCores() - 1L)
library(parallel)
design <- read.csv("simulation/design_seeds.csv", stringsAsFactors = FALSE)
cells <- if (length(args) > 1) args[-1] else unique(design$cell)
stopifnot(all(cells %in% design$cell))
dir.create("results/simulation", showWarnings = FALSE, recursive = TRUE)

Sys.setenv(DSIIE_ROOT = getwd())
init <- function() {
  setwd(Sys.getenv("DSIIE_ROOT"))
  source("R/methods_base.R"); source("R/methods_corrected.R"); source("R/methods_v8.R"); NULL
}
cl <- makeCluster(ncores); clusterExport(cl, "init"); clusterEvalQ(cl, init())
for (cell in cells) {
  d <- design[design$cell == cell, ]
  t0 <- Sys.time()
  recs <- parLapplyLB(cl, split(d, seq_len(nrow(d))), function(r)
    one_rep_v8(r$model_id, r$design_id, r$replication, r$data_seed, r$n, r$p, r$H, r$K, r$alpha, r$c))
  out <- do.call(rbind, recs)
  saveRDS(out, file.path("results/simulation", sprintf("cell_%s_v8.rds", cell)))
  cat(sprintf("%s: %d fits, %d failures, %.1f min\n", cell, nrow(out), sum(out$fit_status != "success"),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
