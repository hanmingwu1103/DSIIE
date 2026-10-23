# Feasibility screen of the frozen candidate set (no method is fitted here).
setwd(Sys.getenv("DSIIE_WORK", "."))
source("integration_work/simulation/lib_sim.R")           # rank_split, cov_n
source("integration_work/author_revision/realdata/apps_data.R")
H_rule <- function(n, p, K = 2) min(5L, floor(n / (K * (p + 1))))
out <- do.call(rbind, lapply(names(APPS), function(id) {
  a <- APPS[[id]]; dat <- a$load(); X <- dat$X; n <- nrow(X); p <- ncol(X)
  ev <- eigen(cov_n(X), symmetric = TRUE, only.values = TRUE)$values
  ng <- if (a$design == "grouped") length(unique(a$group(X))) else NA
  min_outer <- Inf; min_inner <- Inf; min_cell <- Inf
  for (r in seq_len(n_repeats(a))) {
    fo <- outer_folds(a, dat, r)
    for (k in 1:10) {
      tr <- which(fo != k); min_outer <- min(min_outer, length(tr))
      fi <- inner_folds(a, dat, tr, r, k)
      for (j in 1:5) {
        m <- sum(fi != j); min_inner <- min(min_inner, m)
        H <- H_rule(m, p); if (H >= 3) min_cell <- min(min_cell, floor(m / (2 * H)))
      }
    }
  }
  Hin <- H_rule(min_inner, p)
  data.frame(id = id, title = a$title, design = a$design, n = n, p = p,
             rank = sum(ev > 1e-10 * ev[1]), cond = signif(ev[1] / ev[p], 3),
             dup_rows = sum(duplicated(X)), groups = ng, ties_y = sum(duplicated(dat$y)),
             min_outer_train = min_outer, min_inner_train = min_inner,
             H_inner_min = Hin, min_fine_cell = min_cell,
             feasible = Hin >= 3 && min_cell > p)
}))
options(width = 220); print(out, row.names = FALSE)
write.csv(out, "integration_work/author_revision/realdata/feasibility_screen.csv", row.names = FALSE)
