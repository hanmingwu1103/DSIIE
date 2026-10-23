# Response-resolution diagnostic (see PROTOCOL.md, frozen before this run).
# SIIE with H slices, SIIE with HK slices and DSIIE on the archived samples of
# three v8 cells; identity checks of the matched fine-cell objects and of the
# archived records. Run from LaTeX/Working.
setwd(Sys.getenv("DSIIE_WORK", "."))
library(parallel)
A <- "integration_work/author_revision"
ARCH <- file.path(A, "simulation/results_v8")
OUT <- file.path(A, "R1/diagnostic")
CELLS <- c("SIM-PRIMARY_M1_n400_p10_H4_K2_a050_c000_fast",
           "SIM-RANK_M2_n400_p50_H4_K2_a020_c000",
           "SIM-P50_M2_n800_p50_H4_K2_a020_c000")

one_rep <- function(r) {
  set.seed(r$data_seed)
  dat <- gen_model(r$model_id, r$n, r$p, c = r$c)
  P_true <- tcrossprod(qr.Q(qr(dat$B_true)))
  H <- r$H; K <- r$K; al <- r$alpha
  fits <- list(
    "SIIE_H"  = safe_fit({ f <- siie_q(dat$X, dat$Y, 1, H, al); f$fail <- FALSE; f }),
    "SIIE_HK" = safe_fit({ f <- siie_q(dat$X, dat$Y, 1, H * K, al); f$fail <- FALSE; f }),
    "DSIIE"   = safe_fit({ f <- dsiie(dat$X, dat$Y, 1, H, K, al); f$fail <- FALSE; f }))
  # matched fine-cell objects: SIIE slice j vs DSIIE cell (h,k), j = (h-1)K + k
  dc <- dO <- dW <- NA_real_
  if (!fits$SIIE_HK$fail && !fits$DSIIE$fail) {
    s <- fits$SIIE_HK; g <- fits$DSIIE; dc <- dO <- 0
    for (h in seq_len(H)) for (k in seq_len(K)) {
      j <- (h - 1) * K + k
      dc <- max(dc, abs(s$c_hk[[j]][[1]] - g$c_hk[[h]][[k]]))
      dO <- max(dO, abs(s$Om_hk[[j]][[1]] - g$Om_hk[[h]][[k]]))
      stopifnot(s$n_hk[[j]][1] == g$n_hk[[h]][k])
    }
    dW <- max(abs(s$Sx_isqrt - g$Sx_isqrt))
  }
  main <- do.call(rbind, lapply(names(fits), function(nm) { f <- fits[[nm]]
    data.frame(cell = r$cell, replication = r$replication, data_seed = r$data_seed, setting = nm,
               r_ave = if (f$fail) NA_real_ else r_ave(P_true, f$projector),
               fit_status = if (f$fail) "fail" else "success",
               fail_msg = if (f$fail) as.character(f$msg) else "", seconds = f$sec,
               max_diff_centre = dc, max_diff_scatter = dO, max_diff_whitening = dW,
               stringsAsFactors = FALSE) }))
  idc <- do.call(rbind, lapply(c("SIR", "SAVE", "pHd", "DR", "SIME"), function(nm) {
    f <- fit_one8(nm, dat$X, dat$Y, 1, H, K, al)
    data.frame(cell = r$cell, replication = r$replication, data_seed = r$data_seed, method_id = nm,
               r_ave_refit = if (isTRUE(f$fail)) NA_real_ else r_ave(P_true, f$projector)) }))
  list(main = main, idc = idc)
}

design <- do.call(rbind, lapply(CELLS, function(cl) {
  a <- readRDS(file.path(ARCH, sprintf("cell_%s_v8.rds", cl)))
  d <- unique(a[, c("replication", "data_seed", "model_id", "n", "p", "H", "K", "alpha", "c")])
  stopifnot(!anyDuplicated(d$replication), !anyDuplicated(d$data_seed)); d$cell <- cl; d }))
write.csv(design, file.path(OUT, "design_seeds.csv"), row.names = FALSE)

cl <- makeCluster(as.integer(Sys.getenv("DSIIE_CORES", 24)))
clusterExport(cl, "A")
clusterEvalQ(cl, {
  setwd(Sys.getenv("DSIIE_WORK", "."))
  source("integration_work/simulation/lib_sim.R")
  source("integration_work/simulation/methods_corrected.R")
  source(file.path(A, "methods_v8.R")); NULL })
t0 <- Sys.time()
res <- parLapplyLB(cl, split(design, seq_len(nrow(design))), one_rep)
stopCluster(cl)
main <- do.call(rbind, lapply(res, `[[`, "main")); idc <- do.call(rbind, lapply(res, `[[`, "idc"))
attr(main, "source_md5") <- c(
  primary = unname(tools::md5sum(Sys.getenv("DSIIE_SIM_SRC", "external/DSIIE_simulation_all.R"))),
  corrected = unname(tools::md5sum("integration_work/simulation/methods_corrected.R")),
  v8 = unname(tools::md5sum(file.path(A, "methods_v8.R"))),
  runner = unname(tools::md5sum(file.path(OUT, "run_diagnostic.R"))))
attr(main, "R_version") <- R.version.string
saveRDS(main, file.path(OUT, "diagnostic_records.rds"))
write.csv(main, file.path(OUT, "diagnostic_records.csv"), row.names = FALSE)
write.csv(idc, file.path(OUT, "identity_refits.csv"), row.names = FALSE)
writeLines(c(capture.output(sessionInfo()), "", capture.output(print(attr(main, "source_md5")))),
           file.path(OUT, "sessionInfo.txt"))
cat(sprintf("%d fits, %d failures, %.1f min\n", nrow(main), sum(main$fit_status != "success"),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
