# LaTeX table bodies for the application section, generated from
# nrmse_summary_v8.csv and feasibility_screen.csv. Bold = smallest mean NRMSE in
# a row at the displayed precision; a method-dataset with an incomplete
# repetition is shown as a dash.
setwd(Sys.getenv("DSIIE_ROOT", getwd()))
RD <- "results/applications"
S <- read.csv(file.path(RD, "nrmse_summary_v8.csv"), stringsAsFactors = FALSE)
F <- read.csv("applications/feasibility_screen.csv", stringsAsFactors = FALSE)
Tn <- read.csv(file.path(RD, "tuning_failures_v8.csv"), stringsAsFactors = FALSE)
M8 <- c("SIR", "SAVE", "pHd", "DR", "SIME", "csMAVE", "SIIE", "DSIIE")
OUT <- "tables"
f3 <- function(x) sprintf("%.3f", x)

SEL <- c("Mussels", "RD06", "RD09", "RD16", "RD15", "College")
LAB <- c(Mussels = "Mussels", RD06 = "Plasma", RD09 = "Bodyfat", RD16 = "Real estate",
         RD15 = "Aquatic toxicity", College = "College")
RESP <- c(Mussels = "log muscle mass", RD06 = "plasma $\\beta$-carotene", RD09 = "body fat (\\%)",
          RD16 = "price per unit area", RD15 = "LC50 ($-\\log$ mol/L)", College = "graduation rate (\\%)")
SRC <- c(Mussels = "\\citet{CookWeisberg1999}", RD06 = "\\citet{Nierenberg1989}", RD09 = "\\citet{Johnson1996}",
         RD16 = "\\citet{YehHsu2018}", RD15 = "\\citet{Ballabio2019}", College = "\\citet{James2021}")
DES <- c(random = "random", grouped = "grouped", blocked = "blocked")

means <- function(id, grid) { s <- S[S$dataset == id & S$grid == grid, ]; s$nrmse[match(M8, s$method)] }
sds <- function(id) { s <- S[S$dataset == id & S$grid == "primary", ]; s$sd[match(M8, s$method)] }
boldrow <- function(v) {
  out <- ifelse(is.na(v), "--", f3(v)); r <- round(v, 3)
  b <- which(!is.na(r) & r == min(r, na.rm = TRUE)); out[b] <- paste0("\\textbf{", out[b], "}"); out
}

# description table
desc <- sapply(SEL, function(id) { f <- F[F$id == id, ]
  sprintf("%s & %s & %d & %d & %s & %s \\\\", LAB[id], RESP[id], f$n, f$p,
          if (f$design == "grouped") sprintf("grouped (%d)", f$groups) else DES[f$design], SRC[id]) })
writeLines(desc, file.path(OUT, "tab_app_data_v8.tex"))

# main NRMSE table: primary-grid mean with SD beneath, then the wider-grid mean
main <- unlist(lapply(SEL, function(id) {
  mp <- means(id, "primary"); sp <- sds(id); mw <- means(id, "wide")
  c(sprintf("%s & primary & %s \\\\", LAB[id], paste(boldrow(mp), collapse = " & ")),
    sprintf(" & (SD) & %s \\\\", paste(sprintf("{\\scriptsize(%s)}", f3(sp)), collapse = " & ")),
    sprintf(" & wider & %s \\\\[2pt]", paste(boldrow(mw), collapse = " & ")))
}))
main[length(main)] <- sub("\\\\\\\\\\[2pt\\]$", "\\\\\\\\", main[length(main)])
writeLines(main, file.path(OUT, "tab_app_nrmse_v8.tex"))

# appendix: primary-grid means for all 18 analysed data sets
ALL <- F$id[F$feasible]
TITLE <- c(RD01 = "Ozone", RD02 = "AIS lean body mass", RD06 = "Plasma", RD07 = "Automobile price",
           RD08 = "Prostate cancer", RD09 = "Bodyfat", RD10 = "Biomass heating", RD11 = "Concrete strength",
           RD13 = "Airfoil noise", RD14 = "Yacht resistance", RD15 = "Aquatic toxicity",
           RD16 = "Real estate", RD17 = "Concrete slump", RD19 = "Fish toxicity", RD20 = "Computer hardware",
           Mussels = "Mussels", Hitters = "Baseball salaries", College = "College")
stopifnot(all(ALL %in% names(TITLE)))
allrows <- sapply(ALL, function(id) { f <- F[F$id == id, ]
  sprintf("%s & %d & %d & %s & %s \\\\", TITLE[id], f$n, f$p, substr(DES[f$design], 1, 1),
          paste(boldrow(means(id, "primary")), collapse = " & ")) })
writeLines(allrows, file.path(OUT, "tab_app_all_v8.tex"))

# numbers quoted in the text
P <- read.csv(file.path(RD, "nrmse_paired_v8.csv"), stringsAsFactors = FALSE)
sp <- S[S$grid == "primary", ]
rk <- table(factor(sp$rank[sp$method == "DSIIE"], levels = 1:8))
best <- table(sp$method[!is.na(sp$rank) & sp$rank == 1])
cat("DSIIE primary rank distribution over", length(ALL), "data sets:\n"); print(rk)
cat("best method counts (primary):\n"); print(best)
cat("DSIIE rank (wide):\n"); print(table(S$rank[S$grid == "wide" & S$method == "DSIIE"]))
cat("DSIIE minus SIIE (primary) over all data sets, number negative:",
    sum(P$mean_diff[P$grid == "primary" & P$comparator == "SIIE"] < 0), "of", length(ALL), "\n")
cat("median DSIIE - best gap:", median(sapply(ALL, function(id) { m <- means(id, "primary")
  m[8] - min(m, na.rm = TRUE) })), "\n")
x <- P[P$grid == "primary" & P$dataset %in% SEL, c("dataset", "comparator", "mean_diff", "sd_diff", "repeats_DSIIE_lower")]
x$se <- x$sd_diff / sqrt(10); print(x[order(x$dataset), ], digits = 3, row.names = FALSE)
cat("boundary rates by grid (all methods, all data sets): primary gamma-low", round(mean(Tn$prim_gamma_low), 3),
    "lambda-high", round(mean(Tn$prim_lambda_high), 3), "; wider gamma-low", round(mean(Tn$wide_gamma_low), 3),
    "lambda-high", round(mean(Tn$wide_lambda_high), 3), "\n")
for (id in SEL) cat(id, "boundary primary gamma-low", round(mean(Tn$prim_gamma_low[Tn$dataset == id]), 2),
                    "lambda-high", round(mean(Tn$prim_lambda_high[Tn$dataset == id]), 2), "\n")
