# Targeted corrected rerun of the approved simulation cells.
# Data seeds are taken from the archived records (data_seed column), so every
# replication regenerates exactly the archived sample (verified in
# validate_reconstruction.R). Output goes to a fresh directory; archived
# records are never read as results, only for their seed lists and for a
# reproduction comparison of the unchanged methods.
setwd(Sys.getenv("DSIIE_WORK", "."))
library(parallel)
ARCH <- Sys.getenv("DSIIE_SEED_ARCH", "integration_work/author_revision/simulation/results_v8")
OUT  <- Sys.getenv("DSIIE_OUT", "integration_work/author_revision/simulation/results_v8")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

FULL <- c("SIR", "SAVE", "pHd", "DR", "SIME", "csMAVE", "SIIE", "DSIIE")
INTV <- FULL  # every reported cell now carries all eight methods

cells <- list(
  list("cell_SIM-PRIMARY_M1_n400_p10_H4_K2_a050_c000_fast.rds", FULL),
  list("cell_SIM-PRIMARY_M2_n400_p10_H4_K2_a020_c000_fast.rds", FULL),
  list("cell_SIM-PRIMARY_M4_n400_p10_H4_K2_a020_c000_fast.rds", FULL),
  list("cell_SIM-PRIMARY_M6_n400_p10_H4_K2_a020_c000_fast.rds", FULL),
  list("cell_SIM-RANK_M2_n400_p50_H4_K2_a020_c000.rds", FULL),
  list("cell_SIM-P50_M2_n800_p50_H4_K2_a020_c000.rds", FULL))
for (nn in c(200, 400, 800, 1600))
  cells[[length(cells) + 1]] <- list(sprintf("cell_SIM-N_M1_n%d_p10_H4_K2_a050_c000.rds", nn), FULL)
for (cc in c("000", "005", "010"))
  cells[[length(cells) + 1]] <- list(sprintf("cell_SIM-C_M7_n400_p10_H4_K2_a020_c%s.rds", cc), FULL)
for (hh in 2:5)
  cells[[length(cells) + 1]] <- list(sprintf("cell_SIM-H_M3_n400_p10_H%d_K2_a020_c000.rds", hh), INTV)
for (kk in 1:4)
  cells[[length(cells) + 1]] <- list(sprintf("cell_SIM-K_M1_n400_p10_H4_K%d_a050_c000.rds", kk), INTV)
for (mm in c("M1", "M2")) for (aa in c("010", "020", "030", "050"))
  cells[[length(cells) + 1]] <- list(sprintf("cell_SIM-A_%s_n400_p10_H4_K2_a%s_c000.rds", mm, aa), INTV)

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) cells <- cells[as.integer(args)]

cl <- makeCluster(as.integer(Sys.getenv("DSIIE_CORES", 30)))
clusterEvalQ(cl, {
  setwd(Sys.getenv("DSIIE_WORK", "."))
  source("integration_work/simulation/lib_sim.R")
  source("integration_work/simulation/methods_corrected.R")
  source("integration_work/author_revision/methods_v8.R")
  NULL
})

for (cell in cells) {
  f <- cell[[1]]; meths <- cell[[2]]
  out_f <- file.path(OUT, sub("[.]rds$", "_v8.rds", f))
  if (file.exists(out_f)) { cat("exists, skipping:", out_f, "\n"); next }
  arch <- readRDS(file.path(ARCH, if (file.exists(file.path(ARCH, f))) f else sub("[.]rds$", "_v8.rds", f)))
  design <- unique(arch[, c("replication", "data_seed", "model_id", "design_id",
                            "n", "p", "H", "K", "alpha", "c")])
  stopifnot(!anyDuplicated(design$replication))
  t0 <- Sys.time()
  recs <- parLapplyLB(cl, split(design, seq_len(nrow(design))), function(r, meths)
    one_rep_v8(r$model_id, r$design_id, r$replication, r$data_seed,
                      r$n, r$p, r$H, r$K, r$alpha, r$c, meths), meths = meths)
  out <- do.call(rbind, recs)
  attr(out, "source_md5") <- c(
    primary = unname(tools::md5sum(Sys.getenv("DSIIE_SIM_SRC", "external/DSIIE_simulation_all.R"))),
    corrected = unname(tools::md5sum("integration_work/simulation/methods_corrected.R")),
    v8 = unname(tools::md5sum("integration_work/author_revision/methods_v8.R")))
  attr(out, "R_version") <- R.version.string
  saveRDS(out, out_f)
  cat(sprintf("%s: %d rows, %d failures, %.1f min\n", f, nrow(out),
              sum(out$fit_status != "success"),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
