# LaTeX table bodies for the simulation section (eight methods), generated
# from summary_v8.csv / paired_v8.csv. Bold = largest mean in a row at the displayed precision.
setwd(Sys.getenv("DSIIE_WORK", "."))
S <- read.csv("integration_work/author_revision/simulation/summary_v8.csv", stringsAsFactors = FALSE)
P <- read.csv("integration_work/author_revision/simulation/paired_v8.csv", stringsAsFactors = FALSE)
M8 <- c("SIR", "SAVE", "pHd", "DR", "SIME", "csMAVE", "SIIE", "DSIIE")
# R1 revision: the confirmatory run of the additional model M5 (R1/candidates,
# candidate C1) enters the core and paired tables under its own cell identifier.
R1C <- "integration_work/author_revision/R1/candidates"; M5CELL <- "R1-C1_M5_n400_p10_H4_K2_a020_c000"
s5 <- read.csv(file.path(R1C, "summary_confirm.csv"), stringsAsFactors = FALSE); s5 <- s5[s5$candidate == "C1", ]
stopifnot(nrow(s5) == 8, all(s5$fails == 0))
S <- rbind(S, data.frame(cell = M5CELL, model = "M5", n = 400, p = 10, H = 4, K = 2, alpha = 0.2, c = 0,
  method = s5$method, R = s5$R, fails = s5$fails, r_ave = s5$r_ave, sd = s5$mcse * sqrt(s5$R), mcse = s5$mcse,
  r_min = s5$r_min, r_min_mcse = s5$r_min_mcse, sec_median = NA))
p5 <- read.csv(file.path(R1C, "paired_confirm.csv"), stringsAsFactors = FALSE); p5 <- p5[p5$candidate == "C1", ]
P <- rbind(P, data.frame(cell = M5CELL, comparator = p5$comparator, R_pairs = p5$R_pairs,
  mean_diff = p5$diff_r_ave, half_width = p5$hw_r_ave, prop_DSIIE_higher = NA))
dir.create("integration_work/author_revision/tables", showWarnings = FALSE)
f3 <- function(x) sprintf("%.3f", x)
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
  row8("SIM-C_M7_n400_p10_H4_K2_a020_c000", "M7, $c=0$ (M3) & $(400,10)$"),
  row8("SIM-C_M7_n400_p10_H4_K2_a020_c005", "M7, $c=0.05$ & $(400,10)$"),
  row8("SIM-C_M7_n400_p10_H4_K2_a020_c010", "M7, $c=0.10$ & $(400,10)$"),
  row8("SIM-PRIMARY_M6_n400_p10_H4_K2_a020_c000_fast", "M6 & $(400,10)$"),
  row8("SIM-PRIMARY_M6_n400_p10_H4_K2_a020_c000_fast", "M6, $r_{\\min}$ & $(400,10)$", "r_min"),
  row8("SIM-PRIMARY_M4_n400_p10_H4_K2_a020_c000_fast", "M4 & $(400,10)$"),
  row8("SIM-PRIMARY_M4_n400_p10_H4_K2_a020_c000_fast", "M4, $r_{\\min}$ & $(400,10)$", "r_min"),
  row8(M5CELL, "M5 & $(400,10)$"),
  row8(M5CELL, "M5, $r_{\\min}$ & $(400,10)$", "r_min"))
core_cells <- c(unique(sub("^(SIM-[^ ]*).*$", "\\1", c("SIM-PRIMARY_M1_n400_p10_H4_K2_a050_c000_fast",
  "SIM-PRIMARY_M2_n400_p10_H4_K2_a020_c000_fast", "SIM-RANK_M2_n400_p50_H4_K2_a020_c000",
  "SIM-P50_M2_n800_p50_H4_K2_a020_c000", "SIM-C_M7_n400_p10_H4_K2_a020_c000",
  "SIM-C_M7_n400_p10_H4_K2_a020_c005", "SIM-C_M7_n400_p10_H4_K2_a020_c010",
  "SIM-PRIMARY_M6_n400_p10_H4_K2_a020_c000_fast", "SIM-PRIMARY_M4_n400_p10_H4_K2_a020_c000_fast"))), M5CELL)
sc <- S[S$cell %in% core_cells, ]
writeLines(c(core, sprintf("%% max MC se r_ave %.4f ; r_min %.4f", max(sc$mcse), max(sc$r_min_mcse, na.rm = TRUE))),
           "integration_work/author_revision/tables/tab_sim_core_v8.tex")

sens <- c(
  sapply(2:5, function(h) row8(sprintf("SIM-H_M3_n400_p10_H%d_K2_a020_c000", h), sprintf("M3 & $H=%d$", h))),
  sapply(1:4, function(k) row8(sprintf("SIM-K_M1_n400_p10_H4_K%d_a050_c000", k), sprintf("M1 & $K=%d$", k))),
  sapply(c("010", "020", "030", "050"), function(a) row8(sprintf("SIM-A_M1_n400_p10_H4_K2_a%s_c000", a),
         sprintf("M1 & $\\alpha=0.%s$", sub("0$", "", substr(a, 2, 3))))),
  sapply(c("010", "020", "030", "050"), function(a) row8(sprintf("SIM-A_M2_n400_p10_H4_K2_a%s_c000", a),
         sprintf("M2 & $\\alpha=0.%s$", sub("0$", "", substr(a, 2, 3))))))
ss <- S[grepl("^SIM-[HKA]_", S$cell), ]
writeLines(c(sens, sprintf("%% max MC se %.4f", max(ss$mcse))), "integration_work/author_revision/tables/tab_sim_sens_v8.tex")

f3s <- function(x) sprintf("%.3f", x)
prow <- function(cell, label, comps = c("SIR", "DR", "SIME", "csMAVE", "SIIE")) {
  p <- P[P$cell == cell, ]; p <- p[match(comps, p$comparator), ]; stopifnot(nrow(p) == length(comps))
  sprintf("%s & %s \\\\", label, paste(sprintf("$%s\\pm%s$", f3s(p$mean_diff),
          f3(p$half_width)), collapse = " & "))
}
pr <- c(sapply(c(200, 400, 800, 1600), function(n) prow(sprintf("SIM-N_M1_n%d_p10_H4_K2_a050_c000", n), sprintf("M1 & $(%d,10)$", n))),
  prow("SIM-PRIMARY_M2_n400_p10_H4_K2_a020_c000_fast", "M2 & $(400,10)$"),
  prow("SIM-RANK_M2_n400_p50_H4_K2_a020_c000", "M2 & $(400,50)$"),
  prow("SIM-P50_M2_n800_p50_H4_K2_a020_c000", "M2 & $(800,50)$"),
  prow("SIM-C_M7_n400_p10_H4_K2_a020_c000", "M7, $c=0$ & $(400,10)$"),
  prow("SIM-C_M7_n400_p10_H4_K2_a020_c005", "M7, $c=0.05$ & $(400,10)$"),
  prow("SIM-C_M7_n400_p10_H4_K2_a020_c010", "M7, $c=0.10$ & $(400,10)$"),
  prow(M5CELL, "M5 & $(400,10)$"))
writeLines(pr, "integration_work/author_revision/tables/tab_sim_paired_v8.tex")
cat(core, sep = "\n"); cat(sens, sep = "\n"); cat(pr, sep = "\n")

# R1 revision: response-resolution diagnostic (R1/diagnostic). Columns: SIIE (H
# slices), SIIE (HK slices), DSIIE (H primary slices, K cells each), then the
# paired differences DSIIE minus each SIIE setting; failures are reported.
R1D <- "integration_work/author_revision/R1/diagnostic"
DS <- read.csv(file.path(R1D, "diagnostic_summary.csv"), stringsAsFactors = FALSE)
DP <- read.csv(file.path(R1D, "diagnostic_paired.csv"), stringsAsFactors = FALSE)
pm <- function(m, h) sprintf("$%s\\pm%s$", f3s(m), f3(h))
drow <- function(cell, label) {
  s <- DS[DS$cell == cell, ]; s <- s[match(c("SIIE_H", "SIIE_HK", "DSIIE"), s$setting), ]; stopifnot(nrow(s) == 3)
  v <- sprintf("%s {\\scriptsize(%s)}", f3(s$r_ave), f3(s$mcse))
  b <- which(round(s$r_ave, 3) == max(round(s$r_ave, 3))); v[b] <- sub("^([.0-9]+)", "\\\\textbf{\\1}", v[b])
  p <- DP[DP$cell == cell, ]
  q1 <- p[p$contrast == "DSIIE - SIIE (HK slices)", ]; q2 <- p[p$contrast == "DSIIE - SIIE (H slices)", ]
  stopifnot(sum(s$fails) == 0)  # the caption states that no fit failed
  sprintf("%s & %s & %s & %s \\\\", label, paste(v, collapse = " & "),
          pm(q1$mean_diff, q1$half_width), pm(q2$mean_diff, q2$half_width))
}
res <- c(drow("SIM-PRIMARY_M1_n400_p10_H4_K2_a050_c000_fast", "M1 & $(400,10)$"),
         drow("SIM-RANK_M2_n400_p50_H4_K2_a020_c000", "M2 & $(400,50)$"),
         drow("SIM-P50_M2_n800_p50_H4_K2_a020_c000", "M2 & $(800,50)$"))
writeLines(res, "integration_work/author_revision/tables/tab_sim_resolution_R1.tex"); cat(res, sep = "\n")

# R2 revision, Option A (R2_SJoS_20261024/optionA): matched-resolution M5 on the
# 200 confirmatory samples, both metrics; SIIE_4 and DSIIE equal the archived C1 values.
R2A <- "integration_work/author_revision/R2_SJoS_20261024/optionA"
AS <- read.csv(file.path(R2A, "optionA_summary.csv"), stringsAsFactors = FALSE)
AP <- read.csv(file.path(R2A, "optionA_paired.csv"), stringsAsFactors = FALSE)
stopifnot(all(AS$fails == 0), all(AS$R_ok == 200))
arow <- function(metric, label) {
  s <- AS[match(c("SIIE_4", "SIIE_8", "DSIIE"), AS$setting), ]
  m <- s[[metric]]; se <- s[[paste0(metric, "_mcse")]]
  v <- sprintf("%s {\\scriptsize(%s)}", f3(m), f3(se))
  b <- which(round(m, 3) == max(round(m, 3))); v[b] <- sub("^([.0-9]+)", "\\\\textbf{\\1}", v[b])
  q1 <- AP[AP$contrast == "DSIIE - SIIE (8 slices)" & AP$metric == metric, ]
  q2 <- AP[AP$contrast == "DSIIE - SIIE (4 slices)" & AP$metric == metric, ]
  stopifnot(nrow(q1) == 1, nrow(q2) == 1)
  sprintf("%s & %s & %s & %s \\\\", label, paste(v, collapse = " & "),
          pm(q1$mean_diff, q1$half_width), pm(q2$mean_diff, q2$half_width))
}
res2 <- c(res, arow("r_ave", "M5 & $(400,10)$"), arow("r_min", "M5, $r_{\\min}$ & $(400,10)$"))
writeLines(res2, "integration_work/author_revision/tables/tab_sim_resolution_R2.tex"); cat(res2, sep = "\n")
