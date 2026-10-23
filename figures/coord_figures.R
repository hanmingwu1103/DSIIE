# Descriptive full-data two-coordinate displays (2 x 4, eight methods) for the
# representative applications. Adapted from the predecessor's construction
# (SIIE-Rcode/audit_v9/05c-twocoordinate-figures-v9.R) with a response-
# independent orientation:
#   1. fit each method to the complete data (d = 2, K = 2, alpha = .20, H rule);
#   2. Q = qr.Q(qr(beta)); S = centred X Q; rotate S onto its principal axes (SVD);
#   3. scale each coordinate to unit variance; sign chosen so that each coordinate
#      has nonnegative third central moment (no use of the response);
#   4. contours: thin-plate smooth y ~ s(V1, V2, k = 30) per panel, drawn at
#      common levels within a data set, masked where no observation lies within
#      a radius of 0.1 of the scaled range (mgcv::exclude.too.far);
#   5. point colour = response rank on one sequential scale per data set.
args <- commandArgs(trailingOnly = TRUE)
ids <- if (length(args)) args else c("RD06", "RD09")
setwd(Sys.getenv("DSIIE_ROOT", getwd()))
source("R/methods_base.R")
source("R/methods_corrected.R")
source("R/methods_v8.R")
source("applications/apps_data.R")
suppressPackageStartupMessages(library(mgcv))
ramp <- colorRampPalette(c("#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#104281", "#0d366b"))(100)
LAB <- c(RD06 = "Plasma beta-carotene (ng/ml)", RD09 = "Percent body fat")
diag_rows <- list()

for (id in ids) {
  dat <- APPS[[id]]$load(); X <- dat$X; y <- dat$y; n <- nrow(X); p <- ncol(X)
  H <- min(5L, floor(n / (2 * (p + 1))))
  coords <- list(); lam <- list()
  for (m in METHODS8) {
    f <- fit_one8(m, X, y, 2, H, 2, .2); stopifnot(!isTRUE(f$fail))
    Q <- qr.Q(qr(f$beta)); S <- scale(X, TRUE, FALSE) %*% Q; S <- scale(S, TRUE, FALSE)
    D <- S %*% svd(S)$v
    D <- sweep(D, 2, apply(D, 2, sd), "/")
    for (j in 1:2) if (sum(D[, j]^3) < 0) D[, j] <- -D[, j]
    coords[[m]] <- D
    ev <- f$eigenvalues
    lam[[m]] <- if (is.null(ev)) NA_real_ else ev[2] / ev[1]
    diag_rows[[length(diag_rows) + 1]] <- data.frame(dataset = id, method = m, H = H,
      lambda2_over_lambda1 = lam[[m]], stringsAsFactors = FALSE)
  }
  # smooths and common contour levels
  surf <- list(); allz <- c()
  for (m in METHODS8) {
    D <- coords[[m]]; df <- data.frame(y = y, V1 = D[, 1], V2 = D[, 2])
    g <- gam(y ~ s(V1, V2, bs = "tp", k = 30), data = df)
    g1 <- seq(min(D[, 1]), max(D[, 1]), length.out = 70); g2 <- seq(min(D[, 2]), max(D[, 2]), length.out = 70)
    gr <- expand.grid(V1 = g1, V2 = g2)
    z <- as.numeric(predict(g, gr)); far <- exclude.too.far(gr$V1, gr$V2, D[, 1], D[, 2], dist = 0.1)
    z[far] <- NA; surf[[m]] <- list(g1 = g1, g2 = g2, z = matrix(z, 70))
    allz <- c(allz, fitted(g))
  }
  lev <- pretty(range(allz), 8)
  col <- ramp[pmax(1, ceiling(100 * rank(y, ties.method = "average") / n))]
  cairo_pdf(sprintf("figures/fig_app_%s.pdf", id), width = 7.2, height = 4.4)
  layout(rbind(matrix(1:8, 2, 4, byrow = TRUE), rep(9, 4)), heights = c(1, 1, 0.24))
  par(mar = c(2.9, 2.9, 1.6, 0.5), mgp = c(1.6, 0.45, 0), tcl = -0.25, cex.axis = 0.75, cex.lab = 0.8)
  for (m in METHODS8) {
    D <- coords[[m]]; s <- surf[[m]]
    plot(D, type = "n", xlab = "First coordinate", ylab = "Second coordinate",
         main = m, font.main = 1, cex.main = 0.95, asp = NA)
    contour(s$g1, s$g2, s$z, levels = lev, add = TRUE, col = "grey45", lwd = 0.6,
            labcex = 0.5, drawlabels = TRUE)
    points(D, pch = 21, bg = col, col = "grey30", lwd = 0.25, cex = 0.62)
  }
  par(mar = c(2.3, 8, 0.4, 8))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rasterImage(as.raster(matrix(ramp, nrow = 1)), 0, 0.35, 1, 0.9, interpolate = TRUE)
  rect(0, 0.35, 1, 0.9, border = "grey30", lwd = 0.5)
  qs <- c(0, .25, .5, .75, 1)
  axis(1, at = qs, labels = formatC(quantile(y, qs, type = 7), format = "fg", digits = 3),
       cex.axis = 0.7, mgp = c(1, 0.3, 0), tcl = -0.2)
  mtext(sprintf("%s: colour by rank (labels: response quantiles)", LAB[[id]]), side = 1,
        line = 1.25, cex = 0.62)
  dev.off()
  cat("written", id, "H =", H, "\n")
}
write.csv(do.call(rbind, diag_rows), "figures/coord_diagnostics.csv", row.names = FALSE)
