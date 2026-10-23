# LaTeX table bodies for the simulation section (eight methods), generated
# from summary_v8.csv / paired_v8.csv. Bold = largest mean in a row at the displayed precision.
setwd(Sys.getenv("DSIIE_ROOT", getwd()))
S <- read.csv("results/simulation/summary_v8.csv", stringsAsFactors = FALSE)
P <- read.csv("results/simulation/paired_v8.csv", stringsAsFactors = FALSE)
M8 <- c("SIR", "SAVE", "pHd", "DR", "SIME", "csMAVE", "SIIE", "DSIIE")
dir.create("tables", showWarnings = FALSE)
f3 <- function(x) sub("^0[.]", ".", sprintf("%.3f", x))
row8 <- function(cell, label, what = "r_ave") {
  s <- S[S$cell == cell, ]; v <- s[[what]][match(M8, s$method)]; stopifnot(length(v) == 8, !anyNA(v))
  cells <- f3(v); b <- which(round(v, 3) == max(round(v, 3))); cells[b] <- paste0("\\textbf{", cells[b], "}")
  sprintf("%s & %s \\\\", label, paste(cells, collapse = " & "))
}
c1 <- function(x) sprintf("cell_%s", x)
core <- c(
  row8("SIM-PRIMARY_M1_n400_p10_H4_K2_a050_c000_fast", "M1 & $(400,10)$"),
  row8("SIM-PRIMARY_M2_n400_p10_H4_K2_a020_c000_fast", "M2 & $(400,10)$"),
  row8("SIM-RANK_M2_n400_p50_H4_K2_a020_c000", "M2 & $(400,50)$"),
  row8("SIM-P50_M2_n800_p50_H4_K2_a020_c000", "M2 & $(800,50)$"),
  row8("SIM-C_M7_n400_p10_H4_K2_a020_c000", "M7, $c=0$ & $(400,10)$"),
  row8("SIM-C_M7_n400_p10_H4_K2_a020_c005", "M7, $c=.05$ & $(400,10)$"),
  row8("SIM-C_M7_n400_p10_H4_K2_a020_c010", "M7, $c=.10$ & $(400,10)$"),
  row8("SIM-PRIMARY_M6_n400_p10_H4_K2_a020_c000_fast", "M6 & $(400,10)$"),
  row8("SIM-PRIMARY_M6_n400_p10_H4_K2_a020_c000_fast", "M6, $r_{\\min}$ & $(400,10)$", "r_min"),
  row8("SIM-PRIMARY_M4_n400_p10_H4_K2_a020_c000_fast", "M4 & $(400,10)$"),
  row8("SIM-PRIMARY_M4_n400_p10_H4_K2_a020_c000_fast", "M4, $r_{\\min}$ & $(400,10)$", "r_min"))
core_cells <- unique(sub("^(SIM-[^ ]*).*$", "\\1", c("SIM-PRIMARY_M1_n400_p10_H4_K2_a050_c000_fast",
  "SIM-PRIMARY_M2_n400_p10_H4_K2_a020_c000_fast", "SIM-RANK_M2_n400_p50_H4_K2_a020_c000",
  "SIM-P50_M2_n800_p50_H4_K2_a020_c000", "SIM-C_M7_n400_p10_H4_K2_a020_c000",
  "SIM-C_M7_n400_p10_H4_K2_a020_c005", "SIM-C_M7_n400_p10_H4_K2_a020_c010",
  "SIM-PRIMARY_M6_n400_p10_H4_K2_a020_c000_fast", "SIM-PRIMARY_M4_n400_p10_H4_K2_a020_c000_fast")))
sc <- S[S$cell %in% core_cells, ]
writeLines(c(core, sprintf("%% max MC se r_ave %.4f ; r_min %.4f", max(sc$mcse), max(sc$r_min_mcse, na.rm = TRUE))),
           "tables/tab_sim_core_v8.tex")

sens <- c(
  sapply(2:5, function(h) row8(sprintf("SIM-H_M3_n400_p10_H%d_K2_a020_c000", h), sprintf("M3 & $H=%d$", h))),
  sapply(1:4, function(k) row8(sprintf("SIM-K_M1_n400_p10_H4_K%d_a050_c000", k), sprintf("M1 & $K=%d$", k))),
  sapply(c("010", "020", "030", "050"), function(a) row8(sprintf("SIM-A_M1_n400_p10_H4_K2_a%s_c000", a),
         sprintf("M1 & $\\alpha=.%s$", sub("0$", "", substr(a, 2, 3))))),
  sapply(c("010", "020", "030", "050"), function(a) row8(sprintf("SIM-A_M2_n400_p10_H4_K2_a%s_c000", a),
         sprintf("M2 & $\\alpha=.%s$", sub("0$", "", substr(a, 2, 3))))))
ss <- S[grepl("^SIM-[HKA]_", S$cell), ]
writeLines(c(sens, sprintf("%% max MC se %.4f", max(ss$mcse))), "tables/tab_sim_sens_v8.tex")

f3s <- function(x) { s <- sprintf("%+.3f", x); s <- sub("([+-])0[.]", "\\1.", s); sub("^[+]", "", s) }
prow <- function(cell, label, comps = c("SIR", "DR", "SIME", "csMAVE", "SIIE")) {
  p <- P[P$cell == cell, ]; p <- p[match(comps, p$comparator), ]; stopifnot(nrow(p) == length(comps))
  sprintf("%s & %s \\\\", label, paste(sprintf("$%s\\pm%s$", f3s(p$mean_diff),
          sub("^0[.]", ".", sprintf("%.3f", p$half_width))), collapse = " & "))
}
pr <- c(sapply(c(200, 400, 800, 1600), function(n) prow(sprintf("SIM-N_M1_n%d_p10_H4_K2_a050_c000", n), sprintf("M1 & $(%d,10)$", n))),
  prow("SIM-PRIMARY_M2_n400_p10_H4_K2_a020_c000_fast", "M2 & $(400,10)$"),
  prow("SIM-RANK_M2_n400_p50_H4_K2_a020_c000", "M2 & $(400,50)$"),
  prow("SIM-P50_M2_n800_p50_H4_K2_a020_c000", "M2 & $(800,50)$"),
  prow("SIM-C_M7_n400_p10_H4_K2_a020_c000", "M7, $c=0$ & $(400,10)$"),
  prow("SIM-C_M7_n400_p10_H4_K2_a020_c005", "M7, $c=.05$ & $(400,10)$"),
  prow("SIM-C_M7_n400_p10_H4_K2_a020_c010", "M7, $c=.10$ & $(400,10)$"))
writeLines(pr, "tables/tab_sim_paired_v8.tex")
cat(core, sep = "\n"); cat(sens, sep = "\n"); cat(pr, sep = "\n")
