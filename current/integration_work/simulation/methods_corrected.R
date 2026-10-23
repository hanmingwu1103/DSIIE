# Corrections layered on top of the primary simulation definitions (lib_sim.R).
#
# 1. Directional regression. The primary script calls dr::dr(method = "dr"),
#    which does not exist in the dr package (dispatch error), so every archived
#    "DR" value came from its fallback kernel sum_h w_h {2(V_h - I) + m_h m_h'}^2.
#    That is not the Li and Wang (2007) estimator. fit_DR_LW implements the
#    sliced sample version of Li and Wang (2007, Sect. 4):
#      F = 2 sum_h w_h (W_h - I)^2 + 2 (sum_h w_h m_h m_h')^2
#          + 2 (sum_h w_h m_h'm_h) (sum_h w_h m_h m_h'),
#    with W_h = E_n(ZZ' | slice h) and m_h = E_n(Z | slice h), Z standardized
#    with the divisor-n covariance and equal-count slices (rank_split).
# 2. Every method is fitted inside its own tryCatch so that one failure does
#    not discard a whole replication; failures are recorded with messages.

fit_DR_LW <- function(X, Y, d, H) {
  n <- nrow(X); p <- ncol(X); Si <- sym_inv_sqrt(cov_n(X))
  Z <- scale(X, center = TRUE, scale = FALSE) %*% Si
  ord <- order(Y); Zo <- Z[ord, , drop = FALSE]; idx <- rank_split(n, H)
  A1 <- matrix(0, p, p); A2 <- matrix(0, p, p); a3 <- 0
  for (h in seq_len(H)) {
    Zh <- Zo[idx == h, , drop = FALSE]; wh <- nrow(Zh) / n
    m <- colMeans(Zh); W <- crossprod(Zh) / nrow(Zh)
    D <- W - diag(p)
    A1 <- A1 + wh * (D %*% D)
    A2 <- A2 + wh * tcrossprod(m)
    a3 <- a3 + wh * sum(m^2)
  }
  Fm <- 2 * A1 + 2 * (A2 %*% A2) + 2 * a3 * A2
  eg <- eigen((Fm + t(Fm)) / 2, symmetric = TRUE)
  beta <- Si %*% eg$vectors[, seq_len(d), drop = FALSE]
  list(projector = proj_of(beta), beta = beta, fail = FALSE)
}

safe_fit <- function(expr) {
  t0 <- proc.time()[[3]]
  r <- tryCatch(expr, error = function(e) list(projector = NULL, fail = TRUE,
                                               msg = conditionMessage(e)))
  if (is.null(r$projector)) r$fail <- TRUE
  r$sec <- proc.time()[[3]] - t0
  r
}

fit_one <- function(nm, X, Y, d, H, K, alpha) {
  switch(nm,
    "SIR"        = safe_fit(fit_dr(X, Y, d, H * K, "sir")),
    "SAVE"       = safe_fit(fit_dr(X, Y, d, H * K, "save")),
    "SIR(H)"     = safe_fit(fit_dr(X, Y, d, H, "sir")),
    "SAVE(H)"    = safe_fit(fit_dr(X, Y, d, H, "save")),
    "pHd"        = safe_fit(fit_dr(X, Y, d, H, "phdy")),
    "DR"         = safe_fit(fit_DR_LW(X, Y, d, H * K)),
    "DR(H)"      = safe_fit(fit_DR_LW(X, Y, d, H)),
    "CUME"       = safe_fit(fit_CUME(X, Y, d)),
    "csMAVE"     = safe_fit(fit_csMAVE(X, Y, d)),
    "SIIE-Q(H)"  = safe_fit({ r <- siie_q(X, Y, d, H, alpha); list(projector = r$projector, fail = FALSE) }),
    "SIIE-Q(HK)" = safe_fit({ r <- siie_q(X, Y, d, H * K, alpha); list(projector = r$projector, fail = FALSE) }),
    "DSIIE"      = safe_fit({ r <- dsiie(X, Y, d, H, K, alpha); list(projector = r$projector, fail = FALSE) }),
    stop("unknown method ", nm))
}

one_rep_corrected <- function(model_id, design_id, rep, seed, n, p, H, K, alpha,
                              c_val, methods) {
  set.seed(seed)
  dat <- gen_model(model_id, n, p, c = c_val)
  d <- dat$d
  P_true <- tcrossprod(qr.Q(qr(dat$B_true)))
  do.call(rbind, lapply(methods, function(nm) {
    fit <- fit_one(nm, dat$X, dat$Y, d, H, K, alpha)
    fail <- isTRUE(fit$fail)
    data.frame(model_id = model_id, design_id = design_id, replication = rep,
               data_seed = seed, method_id = nm, n = n, p = p, H = H, K = K,
               alpha = alpha, c = c_val, d = d,
               r_ave = if (fail) NA_real_ else r_ave(P_true, fit$projector),
               r_min = if (fail || d < 2) NA_real_ else r_min(dat$B_true, fit$projector),
               fit_status = if (fail) "fail" else "success",
               fail_msg = if (fail) as.character(fit$msg %||% "NULL projector") else "",
               seconds = fit$sec, stringsAsFactors = FALSE)
  }))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
