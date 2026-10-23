# Candidate simulation models (PROTOCOL.md). Usage: Rscript run_sensitivity.R <stage>
# with stage = explore (R = 50) or confirm (R = 200). Run from LaTeX/Working.
setwd(Sys.getenv("DSIIE_WORK", "."))
library(parallel)
A <- "integration_work/author_revision"; OUT <- file.path(A, "R1/candidates")
stage <- commandArgs(trailingOnly = TRUE)[1]; stopifnot(stage == "sens")
R <- 200L; base <- 9300000L
CAND <- data.frame(id = c("C1n200", "C1n800", "C1p20"), model_id = "M5",
                   n = c(200L, 800L, 400L), p = c(10L, 10L, 20L), H = 4L, K = 2L, alpha = 0.20, c = 0,
                   stream = c(10000L, 20000L, 30000L), stringsAsFactors = FALSE)
design <- do.call(rbind, lapply(seq_len(nrow(CAND)), function(i) {
  cc <- CAND[i, ]; data.frame(cc[rep(1, R), c("id", "model_id", "n", "p", "H", "K", "alpha", "c")],
                              replication = seq_len(R), data_seed = base + cc$stream + seq_len(R)) }))
write.csv(design, file.path(OUT, sprintf("design_seeds_%s.csv", stage)), row.names = FALSE)

cl <- makeCluster(as.integer(Sys.getenv("DSIIE_CORES", 24)))
clusterEvalQ(cl, {
  setwd(Sys.getenv("DSIIE_WORK", "."))
  source("integration_work/simulation/lib_sim.R")
  source("integration_work/simulation/methods_corrected.R")
  source("integration_work/author_revision/methods_v8.R")
  # new candidate C3; all other models use the archived generators unchanged
  gen_M2t <- function(n, p) {
    Rp <- Rp_reflector(p)
    Z <- cbind(rnorm(n), matrix(rt(n * (p - 1), df = 3) / sqrt(3), n, p - 1))
    list(X = Z %*% t(Rp), Y = Z[, 1]^2, d = 1L, B_true = matrix(Rp[, 1], p, 1))
  }
  gen_model_base <- gen_model
  gen_model <- function(name, n, p, c = 0) if (name == "M2t") gen_M2t(n, p) else gen_model_base(name, n, p, c)
  NULL })
t0 <- Sys.time()
recs <- parLapplyLB(cl, split(design, seq_len(nrow(design))), function(r) {
  x <- one_rep_v8(r$model_id, r$id, r$replication, r$data_seed, r$n, r$p, r$H, r$K, r$alpha, r$c)
  x$candidate <- r$id; x })
stopCluster(cl)
out <- do.call(rbind, recs)
attr(out, "source_md5") <- c(
  primary = unname(tools::md5sum(Sys.getenv("DSIIE_SIM_SRC", "external/DSIIE_simulation_all.R"))),
  corrected = unname(tools::md5sum("integration_work/simulation/methods_corrected.R")),
  v8 = unname(tools::md5sum(file.path(A, "methods_v8.R"))),
  runner = unname(tools::md5sum(file.path(OUT, "run_sensitivity.R"))))
saveRDS(out, file.path(OUT, sprintf("records_%s.rds", stage)))
cat(sprintf("%s: %d fits, %d failures, %.1f min\n", stage, nrow(out), sum(out$fit_status != "success"),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
