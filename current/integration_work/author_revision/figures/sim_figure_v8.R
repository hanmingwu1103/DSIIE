# M1 sample-size trajectory, eight methods. Categorical colours follow the
# pre-validated reference order (slots 1-8); every series also has its own
# marker and line type, so identity does not depend on colour (print/CVD).
setwd(Sys.getenv("DSIIE_WORK", "."))
S <- read.csv("integration_work/author_revision/simulation/summary_v8.csv", stringsAsFactors = FALSE)
M8 <- c("SIR", "SAVE", "pHd", "DR", "SIME", "csMAVE", "SIIE", "DSIIE")
COL <- c("#eda100", "#e87ba4", "#008300", "#6250d6", "#1baf7a", "#e34948", "#eb6834", "#2a78d6")
PCH <- c(1, 2, 5, 6, 0, 4, 17, 19); LTY <- c(4, 5, 3, 6, 2, 1, 2, 1); LWD <- c(rep(1.3, 7), 2)
ns <- c(200, 400, 800, 1600)
cairo_pdf("DSIIE_main_assets/fig_sim_M1_samplesize.pdf", width = 6.2, height = 3.9)
par(mar = c(4, 4.2, 1, 7.6), mgp = c(2.4, .6, 0), tcl = -.3)
plot(NA, xlim = range(ns), ylim = c(0, 1), log = "x", axes = FALSE,
     xlab = "Sample size n (log scale)", ylab = expression(bar(r)[ave]))
abline(h = seq(0, 1, .2), col = "#e4e3df", lwd = .6)
axis(1, at = ns, col = "#52514e", col.axis = "#52514e", lwd = .6)
axis(2, at = seq(0, 1, .2), las = 1, col = "#52514e", col.axis = "#52514e", lwd = .6)
for (i in seq_along(M8)) {
  v <- t(sapply(ns, function(n) { r <- S[S$model == "M1" & grepl("^SIM-N", S$cell) & S$n == n & S$method == M8[i], ]
    stopifnot(nrow(r) == 1); c(r$r_ave, r$mcse) }))
  lines(ns, v[, 1], col = COL[i], lty = LTY[i], lwd = LWD[i])
  arrows(ns, v[, 1] - 1.96 * v[, 2], ns, v[, 1] + 1.96 * v[, 2], angle = 90, code = 3, length = .02, col = COL[i], lwd = .7)
  points(ns, v[, 1], col = COL[i], pch = PCH[i], cex = .95)
}
par(xpd = NA)
legend(2150, 1.0, legend = M8, col = COL, pch = PCH, lty = LTY, lwd = LWD, bty = "n", cex = .8, text.col = "#0b0b0b")
dev.off()
