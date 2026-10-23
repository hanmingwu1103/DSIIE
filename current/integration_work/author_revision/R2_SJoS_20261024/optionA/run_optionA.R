# Option A: matched-resolution M5 comparison (PROTOCOL.md, frozen before this run).
# SIIE (4 slices), SIIE (8 slices) and DSIIE (H = 4, K = 2) on the 200 R1
# confirmatory M5 samples, with archive-reuse and matched-fine-object checks.
# Run from LaTeX/Working (or set DSIIE_WORK); DSIIE_SIM_SRC locates the primary
# estimator script loaded by lib_sim.R.
setwd(Sys.getenv("DSIIE_WORK", "."))
A <- "integration_work/author_revision"
OUT <- file.path(A, "R2_SJoS_20261024/optionA")
CAND <- file.path(A, "R1/candidates")
source("integration_work/simulation/lib_sim.R")
source("integration_work/simulation/methods_corrected.R")
source(file.path(A, "methods_v8.R"))

# ---- pre-run checks of the seed ledger, archive and code version -------------
ledger <- read.csv(file.path(CAND, "design_seeds_confirm.csv"), stringsAsFactors = FALSE)
ledger <- ledger[ledger$id == "C1", ]
stopifnot(nrow(ledger) == 200, identical(ledger$replication, 1:200),
          identical(as.integer(ledger$data_seed), 9210000L + 1:200),
          all(ledger$model_id == "M5"), all(ledger$n == 400), all(ledger$p == 10),
          all(ledger$H == 4), all(ledger$K == 2), all(ledger$alpha == 0.2), all(ledger$c == 0))
arch <- readRDS(file.path(CAND, "records_confirm.rds"))
md5_arch <- attr(arch, "source_md5")
md5_now <- c(primary = unname(tools::md5sum(SIM_SRC)),
             corrected = unname(tools::md5sum("integration_work/simulation/methods_corrected.R")),
             v8 = unname(tools::md5sum(file.path(A, "methods_v8.R"))))
stopifnot(identical(md5_now, md5_arch[names(md5_now)]))
arch <- arch[arch$candidate == "C1" & arch$method_id %in% c("SIIE", "DSIIE"), ]
stopifnot(nrow(arch) == 400, setequal(arch$data_seed, ledger$data_seed), all(arch$d == 2))

H <- 4L; K <- 2L; al <- 0.20; d <- 2L

one_rep <- function(r) {
  set.seed(r$data_seed)
  dat <- gen_model(r$model_id, r$n, r$p, c = r$c)
  stopifnot(dat$d == d)
  P_true <- tcrossprod(qr.Q(qr(dat$B_true)))
  fits <- list(
    "SIIE_4" = safe_fit({ f <- siie_q(dat$X, dat$Y, d, H, al); f$fail <- FALSE; f }),
    "SIIE_8" = safe_fit({ f <- siie_q(dat$X, dat$Y, d, H * K, al); f$fail <- FALSE; f }),
    "DSIIE"  = safe_fit({ f <- dsiie(dat$X, dat$Y, d, H, K, al); f$fail <- FALSE; f }))
  # matched fine-cell objects, SIIE slice j versus DSIIE cell (h, k), j = (h - 1)K + k
  chk <- data.frame(data_seed = r$data_seed, replication = r$replication, n_cells = NA_integer_,
                    counts_equal = NA, max_diff_weight = NA_real_, max_diff_centre = NA_real_,
                    max_diff_scatter = NA_real_, max_diff_whitening = NA_real_,
                    max_diff_centre_vs_membership = NA_real_, cell_sizes = "", y_breaks = "")
  if (!fits$SIIE_8$fail && !fits$DSIIE$fail) {
    s <- fits$SIIE_8; g <- fits$DSIIE
    n <- nrow(dat$X); ord <- order(dat$Y); lab <- rank_split(n, H * K)
    Z <- scale(dat$X, center = TRUE, scale = FALSE) %*% g$Sx_isqrt
    cnt <- TRUE; dw <- dc <- dO <- dm <- 0; sizes <- integer(0); brk <- character(0)
    for (h in seq_len(H)) for (k in seq_len(K)) {
      j <- (h - 1L) * K + k; idx <- ord[lab == j]
      cnt <- cnt && s$n_hk[[j]][1] == g$n_hk[[h]][k] && g$n_hk[[h]][k] == length(idx)
      dw <- max(dw, abs(s$weights$om[j] - g$weights$om[h] * g$weights$q[[h]][k]))
      dc <- max(dc, abs(s$c_hk[[j]][[1]] - g$c_hk[[h]][[k]]))
      dO <- max(dO, abs(s$Om_hk[[j]][[1]] - g$Om_hk[[h]][[k]]))
      m <- colMeans(Z[idx, , drop = FALSE])
      dm <- max(dm, abs(m - s$c_hk[[j]][[1]]), abs(m - g$c_hk[[h]][[k]]))
      sizes <- c(sizes, length(idx)); brk <- c(brk, sprintf("%.6g:%.6g", min(dat$Y[idx]), max(dat$Y[idx])))
    }
    chk$n_cells <- H * K; chk$counts_equal <- cnt; chk$max_diff_weight <- dw
    chk$max_diff_centre <- dc; chk$max_diff_scatter <- dO
    chk$max_diff_whitening <- max(abs(s$Sx_isqrt - g$Sx_isqrt))
    chk$max_diff_centre_vs_membership <- dm
    chk$cell_sizes <- paste(sizes, collapse = ";"); chk$y_breaks <- paste(brk, collapse = ";")
  }
  main <- do.call(rbind, lapply(names(fits), function(nm) { f <- fits[[nm]]
    data.frame(replication = r$replication, data_seed = r$data_seed, setting = nm,
               r_ave = if (f$fail) NA_real_ else r_ave(P_true, f$projector),
               r_min = if (f$fail) NA_real_ else r_min(dat$B_true, f$projector),
               fit_status = if (f$fail) "fail" else "success",
               fail_msg = if (f$fail) as.character(f$msg %||% "NULL projector") else "",
               seconds = f$sec, stringsAsFactors = FALSE) }))
  list(main = main, chk = chk)
}

t0 <- Sys.time()
cat("start", format(t0, tz = "UTC", usetz = TRUE), "\n")
res <- lapply(split(ledger, seq_len(nrow(ledger))), one_rep)
main <- do.call(rbind, lapply(res, `[[`, "main")); chk <- do.call(rbind, lapply(res, `[[`, "chk"))
t1 <- Sys.time()

# ---- archive reuse check -------------------------------------------------------
key <- c(SIIE_4 = "SIIE", DSIIE = "DSIIE")
reuse <- do.call(rbind, lapply(names(key), function(nm) {
  a <- arch[arch$method_id == key[[nm]], ]; m <- main[main$setting == nm, ]
  j <- merge(m, a[, c("data_seed", "r_ave", "r_min", "fit_status")], by = "data_seed", suffixes = c("", ".arch"))
  stopifnot(nrow(j) == 200)
  data.frame(setting = nm, archived_method = key[[nm]], R = nrow(j),
             status_agree = sum(j$fit_status == j$fit_status.arch),
             max_abs_diff_r_ave = max(abs(j$r_ave - j$r_ave.arch)),
             max_abs_diff_r_min = max(abs(j$r_min - j$r_min.arch)),
             exact_r_ave = sum(j$r_ave == j$r_ave.arch), exact_r_min = sum(j$r_min == j$r_min.arch)) }))

attr(main, "source_md5") <- c(md5_now, runner = unname(tools::md5sum(file.path(OUT, "run_optionA.R"))))
attr(main, "R_version") <- R.version.string
saveRDS(main, file.path(OUT, "optionA_records.rds"))
write.csv(main, file.path(OUT, "optionA_records.csv"), row.names = FALSE)
write.csv(chk, file.path(OUT, "check_matched_fine_objects.csv"), row.names = FALSE)
write.csv(reuse, file.path(OUT, "check_archive_reuse.csv"), row.names = FALSE)
writeLines(c(capture.output(sessionInfo()), "", capture.output(print(attr(main, "source_md5"))),
             sprintf("SIM_SRC = %s", SIM_SRC),
             sprintf("start %s; end %s", format(t0, tz = "UTC", usetz = TRUE), format(t1, tz = "UTC", usetz = TRUE))),
           file.path(OUT, "sessionInfo.txt"))
options(width = 160); print(reuse)
cat(sprintf("matched objects: counts equal in %d/%d; max diffs weight %.3g centre %.3g scatter %.3g whitening %.3g membership-centre %.3g\n",
            sum(chk$counts_equal, na.rm = TRUE), nrow(chk), max(chk$max_diff_weight), max(chk$max_diff_centre),
            max(chk$max_diff_scatter), max(chk$max_diff_whitening), max(chk$max_diff_centre_vs_membership)))
cat(sprintf("%d fits, %d failures, %.2f min\n", nrow(main), sum(main$fit_status != "success"),
            as.numeric(difftime(t1, t0, units = "mins"))))
