# Eight-method dispatch for the author revision:
#   SIR, SAVE, pHd, DR, SIME, csMAVE, SIIE, DSIIE.
# Builds on lib_sim.R (DSIIE/SIIE kernels, csMAVE, CUME, metrics, generators)
# and methods_corrected.R (Li-Wang DR, safe_fit).
#
# Slicing convention. SIR, SAVE, DR and SIME use the same deterministic
# equal-count partition into HK slices: observations ordered by Y (ties broken
# by original index, R's stable order()), block sizes differing by at most one,
# remainder to the earliest blocks (rank_split). SIIE uses H such slices and
# DSIIE H primary slices of K fine cells built by the same rule. All moment
# kernels use divisor-n global whitening and divisor-n_h slice moments.
# pHd (response-based principal Hessian directions) uses no slicing and is
# computed with dr::dr(method = "phdy").

# ---- standardization ---------------------------------------------------------
std_n <- function(X) {
  Si <- sym_inv_sqrt(cov_n(X))
  list(Z = scale(X, center = TRUE, scale = FALSE) %*% Si, Si = Si)
}
basis_out <- function(Si, M, d) {
  eg <- eigen((M + t(M)) / 2, symmetric = TRUE)
  beta <- Si %*% eg$vectors[, seq_len(d), drop = FALSE]
  list(projector = proj_of(beta), beta = beta, eigenvalues = eg$values, fail = FALSE)
}
slice_list <- function(Y, H) {
  ord <- order(Y); lab <- rank_split(length(Y), H)
  lapply(seq_len(H), function(h) ord[lab == h])
}

# ---- SIR and SAVE with the common slicing rule ---------------------------------
fit_SIR_eq <- function(X, Y, d, H) {
  s <- std_n(X); n <- nrow(X); M <- 0
  for (idx in slice_list(Y, H)) M <- M + (length(idx) / n) * tcrossprod(colMeans(s$Z[idx, , drop = FALSE]))
  basis_out(s$Si, M, d)
}
fit_SAVE_eq <- function(X, Y, d, H) {
  s <- std_n(X); n <- nrow(X); p <- ncol(X); M <- matrix(0, p, p)
  for (idx in slice_list(Y, H)) {
    Zh <- s$Z[idx, , drop = FALSE]; V <- cov_n(Zh); D <- diag(p) - V
    M <- M + (length(idx) / n) * (D %*% D)
  }
  basis_out(s$Si, M, d)
}

# ---- spatial (geometric, L1) median: Vardi-Zhang modified Weiszfeld ----------
# Minimizes f(m) = sum_i ||x_i - m||. At a data point y of multiplicity eta,
# y is optimal iff ||R(y)|| <= eta, R(y) = sum_{x_i != y} (x_i - y)/||x_i - y||.
# Returns the median, the iteration count and a convergence flag based on the
# subgradient optimality condition, never silently accepting nonconvergence.
spatial_median <- function(X, tol = 1e-10, maxit = 10000L) {
  X <- as.matrix(X); n <- nrow(X)
  if (n == 1) return(list(m = X[1, ], converged = TRUE, iter = 0L, grad = 0))
  sc <- max(1, max(abs(X)))
  m <- apply(X, 2, median)                      # coordinatewise start
  for (it in seq_len(maxit)) {
    D <- sweep(X, 2, m); r <- sqrt(rowSums(D^2))
    at <- r <= 1e-12 * sc; eta <- sum(at)
    w <- ifelse(at, 0, 1 / r)
    Tt <- colSums(X[!at, , drop = FALSE] * w[!at]) / sum(w[!at])
    if (eta == 0) m_new <- Tt else {
      Rv <- colSums(D[!at, , drop = FALSE] * w[!at]); nr <- sqrt(sum(Rv^2))
      if (nr <= eta) return(list(m = m, converged = TRUE, iter = it, grad = 0))
      m_new <- (1 - eta / nr) * Tt + (eta / nr) * m
    }
    step <- sqrt(sum((m_new - m)^2)); m <- m_new
    if (step <= tol * sc) break
  }
  # optimality check: norm of the gradient (or subgradient slack) relative to n
  D <- sweep(X, 2, m); r <- sqrt(rowSums(D^2)); at <- r <= 1e-12 * sc
  g <- if (any(!at)) colSums(D[!at, , drop = FALSE] / r[!at]) else 0
  gn <- max(0, sqrt(sum(g^2)) - sum(at))
  if (gn > 1e-6 * n) {
    # The iterates converge sublinearly when the minimizer is a data point. A
    # data point x_j is the minimizer iff || sum_{x_i != x_j} (x_j - x_i)/|x_j - x_i| ||
    # <= #{i : x_i = x_j}; test the data point nearest to the current iterate.
    j <- which.min(r); Dj <- sweep(X, 2, X[j, ]); rj <- sqrt(rowSums(Dj^2)); atj <- rj <= 1e-12 * sc
    gj <- max(0, sqrt(sum(colSums(Dj[!atj, , drop = FALSE] / rj[!atj])^2)) - sum(atj))
    if (gj <= 1e-6 * n) return(list(m = X[j, ], converged = TRUE, iter = it, grad = gj))
  }
  list(m = m, converged = gn <= 1e-6 * n, iter = it, grad = gn)
}

# ---- SIME (Dong, Yu and Zhu, 2015) -------------------------------------------
# Weighted, centred covariance of the slice spatial medians of the standardized
# predictor; leading eigenvectors back-transformed by Sigma^{-1/2}.
fit_SIME <- function(X, Y, d, H) {
  s <- std_n(X); n <- nrow(X); p <- ncol(X)
  sl <- slice_list(Y, H); w <- sapply(sl, length) / n
  med <- matrix(0, H, p); ok <- TRUE
  for (h in seq_len(H)) {
    sm <- spatial_median(s$Z[sl[[h]], , drop = FALSE])
    if (!sm$converged) ok <- FALSE
    med[h, ] <- sm$m
  }
  if (!ok) stop("SIME: spatial median did not converge")
  mbar <- colSums(med * w); M <- matrix(0, p, p)
  for (h in seq_len(H)) M <- M + w[h] * tcrossprod(med[h, ] - mbar)
  basis_out(s$Si, M, d)
}

METHODS8 <- c("SIR", "SAVE", "pHd", "DR", "SIME", "csMAVE", "SIIE", "DSIIE")

# Common dispatch; returns list(projector, beta, eigenvalues, fail, msg, sec).
fit_one8 <- function(nm, X, Y, d, H, K, alpha) {
  switch(nm,
    "SIR"    = safe_fit(fit_SIR_eq(X, Y, d, H * K)),
    "SAVE"   = safe_fit(fit_SAVE_eq(X, Y, d, H * K)),
    "pHd"    = safe_fit(fit_dr(X, Y, d, H, "phdy")),
    "DR"     = safe_fit(fit_DR_LW(X, Y, d, H * K)),
    "SIME"   = safe_fit(fit_SIME(X, Y, d, H * K)),
    "csMAVE" = safe_fit(fit_csMAVE(X, Y, d)),
    "SIIE"   = safe_fit({ r <- siie_q(X, Y, d, H, alpha)
                          list(projector = r$projector, beta = r$beta, eigenvalues = r$eigenvalues, fail = FALSE) }),
    "DSIIE"  = safe_fit({ r <- dsiie(X, Y, d, H, K, alpha)
                          list(projector = r$projector, beta = r$beta, eigenvalues = r$eigenvalues, fail = FALSE) }),
    stop("unknown method ", nm))
}

one_rep_v8 <- function(model_id, design_id, rep, seed, n, p, H, K, alpha, c_val,
                       methods = METHODS8) {
  set.seed(seed)
  dat <- gen_model(model_id, n, p, c = c_val)
  d <- dat$d; P_true <- tcrossprod(qr.Q(qr(dat$B_true)))
  do.call(rbind, lapply(methods, function(nm) {
    fit <- fit_one8(nm, dat$X, dat$Y, d, H, K, alpha); fail <- isTRUE(fit$fail)
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
