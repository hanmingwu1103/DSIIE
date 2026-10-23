# Eight-method repeated cross-validation for the frozen eligible application set.
# Specification: candidate_registry.md. Both tuning grids are evaluated on the
# same nested fits (the primary grid is a subset of the wider grid); the
# primary choice is the first minimizer within the primary subset in the
# original grid order, the wider choice the first minimizer over the wider grid.
setwd(Sys.getenv("DSIIE_WORK", "."))
library(parallel)
OUT <- Sys.getenv("DSIIE_OUT", "integration_work/author_revision/realdata/results_v8")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
args <- commandArgs(trailingOnly = TRUE)

worker_init <- function() {
  setwd(Sys.getenv("DSIIE_WORK", "."))
  source("integration_work/simulation/lib_sim.R")
  source("integration_work/simulation/methods_corrected.R")
  source("integration_work/author_revision/methods_v8.R")
  source("integration_work/author_revision/realdata/apps_data.R")
  D_DIM <<- 2L; KK <<- 2L; ALPHA <<- 0.20
  GM_W <<- 2^(-6:2); LAM_W <<- 10^(-4:2)
  GRID <<- expand.grid(lam = LAM_W, gm = GM_W)[, c("gm", "lam")]   # gm outer, lam inner
  PRIM <<- which(GRID$gm %in% c(0.25, 0.5, 1, 2, 4) & GRID$lam %in% 10^(-4:0))
  H_rule <<- function(n, p, K) min(5L, floor(n / (K * (p + 1))))
  sdr_basis <<- function(nm, X, y) {
    H <- H_rule(nrow(X), ncol(X), KK); if (H < 3) stop("inadmissible H = ", H)
    f <- fit_one8(nm, X, y, D_DIM, H, KK, ALPHA)
    if (isTRUE(f$fail)) stop(f$msg %||% "SDR fit failed")
    B <- as.matrix(f$beta); if (ncol(B) < D_DIM || any(!is.finite(B))) stop("bad basis")
    list(beta = B, H = H)
  }
  score_map <<- function(Xtr, beta) {
    mu <- colMeans(Xtr); S <- sweep(Xtr, 2, mu) %*% beta
    sv <- svd(S); if (min(sv$d) <= 1e-10 * max(sv$d)) stop("degenerate scores")
    Tm <- sv$v %*% diag(sqrt(nrow(S)) / sv$d, length(sv$d))
    function(Xnew) sweep(Xnew, 2, mu) %*% beta %*% Tm
  }
  rbf <<- function(A, B, g) exp(-g * pmax(outer(rowSums(A^2), rowSums(B^2), "+") - 2 * A %*% t(B), 0))
  g0f <<- function(S) { D2 <- as.numeric(dist(S))^2; D2 <- D2[D2 > 0]; 1 / median(D2) }
  krr_sse <<- function(Str, ytr, Ste, yte) {
    g0 <- g0f(Str); yb <- mean(ytr); out <- numeric(nrow(GRID))
    for (gm in GM_W) {
      Kt <- rbf(Str, Str, gm * g0); Ks <- rbf(Ste, Str, gm * g0)
      for (lam in LAM_W) {
        a <- solve(Kt + lam * diag(nrow(Kt)), ytr - yb)
        out[GRID$gm == gm & GRID$lam == lam] <- sum((yte - (as.numeric(Ks %*% a) + yb))^2)
      }
    }
    out
  }
  NULL
}

run_split <- function(task) {
  id <- task$id; a <- APPS[[id]]; dat <- a$load(); r <- task$r; k <- task$k
  fo <- outer_folds(a, dat, r); te <- which(fo == k); tr <- which(fo != k)
  Xtr <- dat$X[tr, , drop = FALSE]; ytr <- dat$y[tr]
  fi <- inner_folds(a, dat, tr, r, k)
  rows <- list()
  for (nm in METHODS8) {
    t0 <- proc.time()[[3]]
    res <- tryCatch({
      sse <- 0; Hin <- integer(0)
      for (j in 1:5) {
        itr <- which(fi != j); iva <- which(fi == j)
        b <- sdr_basis(nm, Xtr[itr, , drop = FALSE], ytr[itr]); Hin <- c(Hin, b$H)
        sm <- score_map(Xtr[itr, , drop = FALSE], b$beta)
        sse <- sse + krr_sse(sm(Xtr[itr, , drop = FALSE]), ytr[itr], sm(Xtr[iva, , drop = FALSE]), ytr[iva])
      }
      ip <- PRIM[which.min(sse[PRIM])]; iw <- which.min(sse)
      b <- sdr_basis(nm, Xtr, ytr); sm <- score_map(Xtr, b$beta)
      Str <- sm(Xtr); Ste <- sm(dat$X[te, , drop = FALSE]); g0 <- g0f(Str)
      pred <- function(i) { g <- GRID$gm[i] * g0
        a_ <- solve(rbf(Str, Str, g) + GRID$lam[i] * diag(nrow(Str)), ytr - mean(ytr))
        as.numeric(rbf(Ste, Str, g) %*% a_) + mean(ytr) }
      list(yp = pred(ip), yw = pred(iw), gp = GRID$gm[ip], lp = GRID$lam[ip], gw = GRID$gm[iw],
           lw = GRID$lam[iw], H = b$H, Hin = paste(Hin, collapse = ";"), status = "ok", msg = "")
    }, error = function(e) list(yp = rep(NA_real_, length(te)), yw = rep(NA_real_, length(te)),
                                gp = NA, lp = NA, gw = NA, lw = NA, H = NA, Hin = "",
                                status = "failed", msg = conditionMessage(e)))
    rows[[nm]] <- data.frame(dataset = id, method = nm, repeat_ = r, fold = k, row_id = te,
      yobs = dat$y[te], yhat_primary = res$yp, yhat_wide = res$yw,
      gamma_mult_primary = res$gp, lambda_primary = res$lp, gamma_mult_wide = res$gw,
      lambda_wide = res$lw, H_outer = res$H, H_inner = res$Hin, fit_status = res$status,
      failure = res$msg, seconds = proc.time()[[3]] - t0, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

source("integration_work/author_revision/realdata/apps_data.R")
ids <- if (length(args)) args else names(APPS)
if (Sys.getenv("RD_TEST") == "1") {
  worker_init(); t0 <- Sys.time()
  x <- run_split(list(id = ids[1], r = 1, k = 1))
  print(unique(x[, c("method", "H_outer", "H_inner", "gamma_mult_primary", "lambda_primary", "fit_status", "failure")]))
  print(Sys.time() - t0); quit(save = "no")
}
size <- sapply(ids, function(id) nrow(APPS[[id]]$load()$X))
ids <- ids[order(-size)]                                      # large data sets first
tasks <- list()
for (id in ids) for (r in seq_len(n_repeats(APPS[[id]]))) for (k in 1:10)
  tasks[[length(tasks) + 1]] <- list(id = id, r = r, k = k)
cl <- makeCluster(as.integer(Sys.getenv("DSIIE_CORES", 30))); clusterExport(cl, "worker_init"); clusterEvalQ(cl, worker_init())
clusterExport(cl, "run_split")
t0 <- Sys.time()
res <- parLapplyLB(cl, tasks, run_split)
stopCluster(cl)
all <- do.call(rbind, res)
for (id in unique(all$dataset))
  write.csv(all[all$dataset == id, ], file.path(OUT, sprintf("oof_%s.csv", id)), row.names = FALSE)
cat(sprintf("done: %d rows, %d failed method-splits, %.1f min\n", nrow(all),
            sum(!duplicated(all[all$fit_status != "ok", c("dataset", "method", "repeat_", "fold")])),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(c(capture.output(sessionInfo()),
  sprintf("md5 %s %s", c("methods_v8.R", "apps_data.R", "run_apps_v8.R"),
          tools::md5sum(c("integration_work/author_revision/methods_v8.R",
                          "integration_work/author_revision/realdata/apps_data.R",
                          "integration_work/author_revision/realdata/run_apps_v8.R")))),
  file.path(OUT, "provenance.txt"))
