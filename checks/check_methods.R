# Substantive checks of the eight-method dispatch (author revision).
setwd(Sys.getenv("DSIIE_ROOT", getwd()))
source("R/methods_base.R")
source("R/methods_corrected.R")
source("R/methods_v8.R")
rep_ <- function(label, val, tol) cat(sprintf("%-66s %s (%.2e)\n", label, if (val <= tol) "agree" else "DIFFER", val))
set.seed(20261024)
f_obj <- function(m, X) sum(sqrt(rowSums(sweep(X, 2, m)^2)))

## (1) Spatial median against an independent numerical minimizer of sum ||x_i - m||.
worst <- 0
for (k in 1:40) {
  p <- sample(2:8, 1); n <- sample(5:60, 1)
  X <- matrix(rt(n * p, df = 3), n, p)
  a <- spatial_median(X)
  b <- optim(colMeans(X), f_obj, X = X, method = "BFGS", control = list(reltol = 1e-14, maxit = 5000))
  b2 <- optim(b$par, f_obj, X = X, method = "Nelder-Mead", control = list(reltol = 1e-15, maxit = 20000))
  worst <- max(worst, f_obj(a$m, X) - min(b$value, b2$value))
  stopifnot(a$converged)
}
rep_("spatial median: objective not above independent minimizer (40 cases)", max(0, worst), 1e-7)

## (2) Coincident observations: median at a data point (majority) and duplicates.
X <- rbind(matrix(0, 6, 3), diag(3), -diag(3)[1:2, ])
a <- spatial_median(X)
rep_("coincident majority point: median equals the repeated point", sqrt(sum(a$m^2)), 1e-10)
X2 <- rbind(matrix(rnorm(20 * 4), 20), matrix(rep(c(1, 2, 3, 4), each = 5), 5))
a2 <- spatial_median(X2); b2 <- optim(colMeans(X2), f_obj, X = X2, method = "Nelder-Mead",
                                     control = list(reltol = 1e-15, maxit = 50000))
rep_("duplicated block of 5 rows: finite, converged, objective optimal",
     if (all(is.finite(a2$m)) && a2$converged) max(0, f_obj(a2$m, X2) - b2$value) else Inf, 1e-7)
# median located exactly at a data point that is optimal but not a majority
X3 <- rbind(c(0, 0), c(1, 0), c(-1, 0), c(0, 1), c(0, -1))
rep_("symmetric configuration: median at the central data point", sqrt(sum(spatial_median(X3)$m^2)), 1e-10)

## (3) Orthogonal equivariance and translation equivariance.
X <- matrix(rnorm(50 * 5), 50); Q <- qr.Q(qr(matrix(rnorm(25), 5))); b <- rnorm(5)
rep_("spatial median: m(XQ' + b) = Q m(X) + b",
     max(abs(spatial_median(sweep(X %*% t(Q), 2, b, "+"))$m - (Q %*% spatial_median(X)$m + b))), 1e-7)

## (4) SIME: affine equivariance of the projector (whitening/back-transformation).
X <- matrix(rnorm(400 * 6), 400) %*% matrix(runif(36), 6); Y <- X[, 1] + X[, 2] * rnorm(400)
A <- matrix(rnorm(36), 6); s1 <- fit_SIME(X, Y, 2, 8); s2 <- fit_SIME(X %*% t(A), Y, 2, 8)
rep_("SIME: projector of X A' equals P of A^{-T} beta", max(abs(proj_of(solve(t(A), s1$beta)) - s2$projector)), 1e-6)

## (5) SIME is not coordinatewise-median substitution: kernels differ on skewed data.
Xs <- matrix(rexp(300 * 4), 300); Ys <- Xs[, 1] + rnorm(300)
sz <- std_n(Xs); sl <- slice_list(Ys, 6)
cw <- t(sapply(sl, function(i) apply(sz$Z[i, ], 2, median)))
sp <- t(sapply(sl, function(i) spatial_median(sz$Z[i, ])$m))
cat(sprintf("%-66s %.3e (should be clearly > 0)\n", "spatial vs coordinatewise slice medians, max difference", max(abs(cw - sp))))

## (6) SIME with duplicated observations: finite output.
Xd <- rbind(X, X[1:50, ]); Yd <- c(Y, Y[1:50])
sd_ <- fit_SIME(Xd, Yd, 2, 8)
rep_("SIME with 50 duplicated rows: finite basis", if (all(is.finite(sd_$beta))) 0 else Inf, 0)

## (7) SIR with the common rule agrees with dr::dr when slices coincide (n divisible).
X <- matrix(rnorm(400 * 6), 400); Y <- X[, 1] + 0.3 * rnorm(400)
rep_("SIR (equal-count) vs dr::dr sir, 8 slices of 50", max(abs(fit_SIR_eq(X, Y, 2, 8)$projector -
     fit_dr(X, Y, 2, 8, "sir")$projector)), 1e-8)
cat(sprintf("%-66s %.3e\n", "SAVE (equal-count, divisor n_h) vs dr::dr save, max |P diff|",
            max(abs(fit_SAVE_eq(X, Y, 2, 8)$projector - fit_dr(X, Y, 2, 8, "save")$projector))))

## (8) Previously validated: SIIE = DSIIE at K = 1; DR expanded = pairwise (see check_methods.R).
rep_("SIIE kernel = DSIIE kernel at K = 1", max(abs(siie_q(X, Y, 1, 5, .2)$M - dsiie(X, Y, 1, 5, 1, .2)$M)), 1e-14)

## (9) Dispatch returns eight finite fits on a model sample.
set.seed(1); dat <- gen_model("M2", 400, 10)
ok <- sapply(METHODS8, function(m) !isTRUE(fit_one8(m, dat$X, dat$Y, 1, 4, 2, .2)$fail))
cat("dispatch, M2 sample, all eight methods succeed:", all(ok), "\n")
writeLines(capture.output(sessionInfo()), "checks/sessionInfo_checks.txt")
