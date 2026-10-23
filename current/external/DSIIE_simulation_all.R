################################################################################
#
#  DSIIE_simulation_all.R
#
#  整合版：一個檔完成所有模擬、填表、畫圖。
#
#  對應 DSIIE_Simulation_Study_Standalone.tex 的 Section 7 "Results template"：
#    Table 3  (tab:primary)   — §7.1  Primary recovery, n=400, p=10
#    Table 4  (tab:p50)       — §7.2  Moderate-dim recovery, n=800, p=50
#    Table 5  (tab:paired)    — §7.3  Paired DSIIE vs SIIE-Q comparisons
#    Table 6  (tab:rmin)      — §7.4  Minimum canonical correlation (M4–M6)
#    Table 7  (tab:null)      — §7.5  Null calibration & diagnostics
#    Table 8  (tab:runtime)   — §7.5  Runtime
#    Figure 1 (fig:trajectories) — §7.6  Sample-size trajectories
#    Figure 2 (fig:sensitivity)  — §7.6  Sensitivity panels (H, K, α, c)
#
#  使用方式：
#    1. 在 R 裡 setwd() 到此檔所在目錄
#    2. source("DSIIE_simulation_all.R")
#       → 會在同目錄下產生 results/, figures/, tex/ 三個資料夾
#    3. 要放大 replications → 修改下面 PART 0 的 R_* 常數
#
#  需要套件：MASS, dr, MAVE
#
################################################################################


########################################################################
#  PART 0 ── 設定：replications 與方法清單
#  （想跑完整規格就改這裡的數字）
########################################################################

R_PRIMARY <- 300L     # spec: 300  (Table 3, 5, 6 用)
R_P50     <- 150L     # spec: 150  (Table 4 用)
R_TRAJ    <- 200L     # spec: 200  (Figure 1 用)
R_SENS_H  <- 200L     # spec: 200  (Figure 2 Panel A 用)
R_SENS_K  <- 200L     # spec: 200  (Figure 2 Panel B 用)
R_SENS_A  <- 200L     # spec: 200  (Figure 2 Panel C 用)
R_SENS_C  <- 200L     # spec: 200  (Figure 2 Panel D 用)
R_RANK    <- 100L     # spec: 100  (Table 4 附加)
R_NULL    <- 500L     # spec: 500  (Table 7 用)
B_NULL    <- 199L     # spec: 199  (Table 7：每次 permutation 數)
R_RUNTIME <- 50L     # spec: 50   (Table 8 用)
R_CSMAVE  <- 300L     # csMAVE 太慢，PRIMARY 只跑這麼多

suppressPackageStartupMessages({ library(MASS); library(dr); library(MAVE) })

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("tex",     showWarnings = FALSE)


########################################################################
#  PART 1 ── 基礎工具函數
########################################################################

# 對稱正定矩陣的對稱平方根逆（不是 Cholesky）
sym_inv_sqrt <- function(S, tol = 1e-12) {
  eg <- eigen(S, symmetric = TRUE)
  if (any(eg$values <= tol)) stop("sym_inv_sqrt: matrix not positive definite")
  eg$vectors %*% diag(1 / sqrt(eg$values)) %*% t(eg$vectors)
}

# divisor-n 的樣本共變異數（Algorithm 1 指定用 n 而非 n-1）
cov_n <- function(X) { Xc <- scale(X, center=TRUE, scale=FALSE); crossprod(Xc)/nrow(Xc) }

# type-7 quantile（R 預設）
q7 <- function(x, p) as.numeric(quantile(x, probs=p, type=7, names=FALSE))

# 投影矩陣 P = Q Q^T，Q 為 span(B) 的正交基
proj_of <- function(B) {
  if (is.null(B) || length(B)==0) return(NULL)
  B <- as.matrix(B); qrB <- qr(B)
  Q <- qr.Q(qrB)[, seq_len(qrB$rank), drop=FALSE]; Q %*% t(Q)
}

# 把 1..n 等分成 H 段（每段大小差至多 1）
rank_split <- function(n, H) {
  base <- n %/% H; rem <- n - base*H
  rep.int(seq_len(H), rep(base,H) + c(rep(1L,rem), rep(0L,H-rem)))
}


########################################################################
#  PART 2 ── DSIIE 演算法（對應主論文 Algorithm 1）
########################################################################

dsiie <- function(X, Y, d, H, K, alpha = 0.20) {
  n <- nrow(X); p <- ncol(X)
  stopifnot(H>=2, K>=1, n>=H*K, d>=1, d<p, alpha>0, alpha<1)
  
  # Step 1: 白化
  Sx <- cov_n(X)
  Sx_isqrt <- sym_inv_sqrt(Sx)
  Z <- scale(X, center=TRUE, scale=FALSE) %*% Sx_isqrt
  
  # Step 3: 依 Y 排序後切 HK 個 fine cell
  ord <- order(Y); Zo <- Z[ord,,drop=FALSE]
  fine_idx <- rank_split(n, H*K)
  
  c_hk <- Om_hk <- n_hk <- vector("list", H)
  for (h in seq_len(H)) {
    c_hk[[h]] <- vector("list", K); Om_hk[[h]] <- vector("list", K); n_hk[[h]] <- integer(K)
    for (k in seq_len(K)) {
      rows <- which(fine_idx == (h-1L)*K + k)
      Zhk <- Zo[rows,,drop=FALSE]; nhk <- nrow(Zhk); n_hk[[h]][k] <- nhk
      
      # Step 6-7: fine cell 中心與 cell 內共變異數
      chk <- colMeans(Zhk)
      Whk <- sweep(Zhk, 2L, chk, "-")
      Gam <- (crossprod(Whk)/nhk + t(crossprod(Whk)/nhk))/2
      
      # Step 8: eigendecomp
      eg <- eigen(Gam, symmetric=TRUE); Vv <- eg$vectors
      
      # sign convention
      for (j2 in seq_len(p)) {
        idx <- which(abs(Vv[,j2])>1e-12)[1]
        if (!is.na(idx) && Vv[idx,j2]<0) Vv[,j2] <- -Vv[,j2]
      }
      
      # Step 9-10: 沿每個 eigenvector 的 type-7 interquantile 半寬
      S <- Whk %*% Vv
      r <- 0.5*(apply(S,2L,q7,p=1-alpha/2) - apply(S,2L,q7,p=alpha/2))
      
      # Step 12: Omega_hk = (1/3) sum r_j^2 v_j v_j^T
      Om <- (Vv %*% diag(r^2, p, p) %*% t(Vv))/3; Om <- (Om+t(Om))/2
      c_hk[[h]][[k]] <- chk; Om_hk[[h]][[k]] <- Om
    }
  }
  
  # Step 4: 實際計數權重
  om_h <- sapply(n_hk, sum)/n
  q_kh <- lapply(n_hk, function(v) v/sum(v))
  
  # Step 15-16: Psi_h（slice-level dispersion）
  c_h <- Psi_h <- vector("list", H)
  for (h in seq_len(H)) {
    ch <- Reduce("+", lapply(seq_len(K), function(k) q_kh[[h]][k]*c_hk[[h]][[k]]))
    ctr <- Reduce("+", lapply(seq_len(K), function(k)
      q_kh[[h]][k] * tcrossprod(c_hk[[h]][[k]]-ch)))
    omavg <- Reduce("+", lapply(seq_len(K), function(k) q_kh[[h]][k]*Om_hk[[h]][[k]]))
    c_h[[h]] <- ch; Psi_h[[h]] <- (ctr + omavg + t(ctr + omavg))/2
  }
  
  # Step 17-18: kernel
  c_bar <- Reduce("+", lapply(seq_len(H), function(h) om_h[h]*c_h[[h]]))
  Psi_bar <- Reduce("+", lapply(seq_len(H), function(h) om_h[h]*Psi_h[[h]]))
  M <- matrix(0,p,p)
  for (h in seq_len(H)) {
    M <- M + om_h[h]*tcrossprod(c_h[[h]]-c_bar) +
      om_h[h]*((Psi_h[[h]]-Psi_bar) %*% (Psi_h[[h]]-Psi_bar))
  }
  M <- (M+t(M))/2
  
  # Step 19-20: 取 top-d eigenvectors，回原座標
  eg <- eigen(M, symmetric=TRUE)
  beta <- Sx_isqrt %*% eg$vectors[,seq_len(d),drop=FALSE]
  
  list(beta=beta, projector=proj_of(beta), M=M, eigenvalues=eg$values,
       Sx_isqrt=Sx_isqrt, c_hk=c_hk, Om_hk=Om_hk, n_hk=n_hk,
       c_h=c_h, Psi_h=Psi_h, weights=list(om=om_h, q=q_kh))
}

# SIIE-Q(H) = DSIIE 的 K=1 特例
siie_q <- function(X, Y, d, H, alpha=0.20) dsiie(X,Y,d,H,K=1L,alpha=alpha)


########################################################################
#  PART 3 ── 比較方法（SIR, SAVE, pHd, DR, CUME, csMAVE）
########################################################################

fit_dr <- function(X, Y, d, H, method) {
  df <- data.frame(Y=Y, X)
  m <- try(dr(Y~., data=df, method=method, nslices=H), silent=TRUE)
  if (inherits(m,"try-error")) return(list(projector=NULL, beta=NULL, fail=TRUE))
  B <- m$evectors[, seq_len(d), drop=FALSE]
  list(projector=proj_of(as.matrix(B)), beta=as.matrix(B), fail=FALSE)
}

fit_DR <- function(X, Y, d, H) {
  df <- data.frame(Y=Y, X)
  m <- try(dr(Y~., data=df, method="dr", nslices=H), silent=TRUE)
  if (inherits(m,"try-error")) {
    # fallback: 手寫 Li & Wang 2007 PSD kernel
    n <- nrow(X); p <- ncol(X); Sx <- cov_n(X); Si <- sym_inv_sqrt(Sx)
    Z <- scale(X,center=TRUE,scale=FALSE) %*% Si
    ord <- order(Y); Zo <- Z[ord,,drop=FALSE]; idx <- rank_split(n,H)
    M <- matrix(0,p,p)
    for (h in seq_len(H)) {
      Zh <- Zo[idx==h,,drop=FALSE]; nh <- nrow(Zh); wh <- nh/n
      mu <- colMeans(Zh); V <- crossprod(sweep(Zh,2,mu,"-"))/nh
      A <- 2*(V-diag(p))+tcrossprod(mu); M <- M+wh*(A%*%A)
    }
    M <- (M+t(M))/2; eg <- eigen(M,symmetric=TRUE)
    beta <- Si %*% eg$vectors[,seq_len(d),drop=FALSE]
    return(list(projector=proj_of(beta), beta=beta, fail=FALSE))
  }
  B <- m$evectors[,seq_len(d),drop=FALSE]
  list(projector=proj_of(as.matrix(B)), beta=as.matrix(B), fail=FALSE)
}

fit_CUME <- function(X, Y, d) {
  n <- nrow(X); p <- ncol(X); Sx <- cov_n(X); Si <- sym_inv_sqrt(Sx)
  Z <- scale(X,center=TRUE,scale=FALSE) %*% Si
  ord <- order(Y); Zo <- Z[ord,,drop=FALSE]; Yo <- Y[ord]
  csum <- apply(Zo,2,cumsum)/n
  tie_end <- which(c(diff(Yo)>0, TRUE))
  grp_end <- rep(tie_end, times=diff(c(0,tie_end)))
  m_at_j <- csum[grp_end,,drop=FALSE]
  M <- crossprod(m_at_j)/n; M <- (M+t(M))/2
  eg <- eigen(M,symmetric=TRUE)
  beta <- Si %*% eg$vectors[,seq_len(d),drop=FALSE]
  list(projector=proj_of(beta), beta=beta, fail=FALSE)
}

fit_csMAVE <- function(X, Y, d) {
  m <- try(suppressWarnings(mave(Y ~ X, method="csMAVE")), silent=TRUE)
  if (inherits(m,"try-error")) return(list(projector=NULL, beta=NULL, fail=TRUE))
  B <- try(coef(m, dim=d), silent=TRUE)
  if (inherits(B,"try-error") || is.null(B)) return(list(projector=NULL, beta=NULL, fail=TRUE))
  B <- as.matrix(B); list(projector=proj_of(B), beta=B, fail=FALSE)
}

# 統一介面：一次跑指定的 methods
fit_all <- function(X, Y, d, H, K, alpha, methods) {
  out <- vector("list", length(methods)); names(out) <- methods
  for (nm in methods) {
    out[[nm]] <- switch(nm,
                        "SIR"        = fit_dr(X,Y,d,H*K,"sir"),
                        "SAVE"       = fit_dr(X,Y,d,H*K,"save"),
                        "pHd"        = fit_dr(X,Y,d,H,  "phdy"),
                        "DR"         = fit_DR(X,Y,d,H*K),
                        "CUME"       = fit_CUME(X,Y,d),
                        "csMAVE"     = fit_csMAVE(X,Y,d),
                        "SIIE-Q(H)"  = { r <- try(siie_q(X,Y,d,H,alpha),silent=TRUE)
                        if(inherits(r,"try-error")) list(projector=NULL,fail=TRUE)
                        else list(projector=r$projector, beta=r$beta, fail=FALSE) },
                        "SIIE-Q(HK)" = { r <- try(siie_q(X,Y,d,H*K,alpha),silent=TRUE)
                        if(inherits(r,"try-error")) list(projector=NULL,fail=TRUE)
                        else list(projector=r$projector, beta=r$beta, fail=FALSE) },
                        "DSIIE"      = { r <- try(dsiie(X,Y,d,H,K,alpha),silent=TRUE)
                        if(inherits(r,"try-error")) list(projector=NULL,fail=TRUE)
                        else list(projector=r$projector, beta=r$beta, fail=FALSE) })
  }
  out
}


########################################################################
#  PART 4 ── 評估準則（r_ave, r_min, projector error）
########################################################################

r_ave <- function(P_true, P_hat) {
  if (is.null(P_hat)) return(NA_real_)
  d <- as.integer(round(sum(diag(P_true))))
  sqrt(max(0, sum(diag(P_true %*% P_hat))/d))
}

r_min <- function(B_true, P_hat) {
  if (is.null(P_hat)) return(NA_real_)
  B_true <- qr.Q(qr(B_true))
  M <- t(B_true) %*% P_hat %*% B_true
  ev <- eigen((M+t(M))/2, symmetric=TRUE)$values
  sqrt(max(0, min(ev)))
}


########################################################################
#  PART 5 ── 資料生成模型 M0–M8
########################################################################

Rp_reflector <- function(p) {
  u <- seq_len(p)/sqrt(sum(seq_len(p)^2)); diag(p) - 2*tcrossprod(u)
}

gen_M0 <- function(n,p) list(X=matrix(rnorm(n*p),n,p), Y=rnorm(n), d=0L, B_true=NULL)

gen_M1 <- function(n,p) {
  Rp <- Rp_reflector(p)
  G <- sample.int(4,n,replace=TRUE); L <- sample.int(2,n,replace=TRUE)
  a <- c(0.2,0.2,5,5)[G]; delt <- a/sqrt(1+a^2); s <- c(-1,1)[L]
  A <- s*delt*abs(rnorm(n)) + sqrt(1-delt^2)*rnorm(n)
  Zmat <- cbind(A, matrix(rnorm(n*(p-1)),n,p-1))
  X <- Zmat %*% t(Rp); Y <- (2*(G-1)+L-1) + runif(n)
  list(X=X, Y=Y, d=1L, B_true=matrix(Rp[,1],p,1))
}

gen_M2 <- function(n,p) {
  Rp <- Rp_reflector(p); Z <- matrix(rnorm(n*p),n,p)
  X <- Z %*% t(Rp); Y <- Z[,1]^2
  list(X=X, Y=Y, d=1L, B_true=matrix(Rp[,1],p,1))
}

gen_M3 <- function(n,p) {
  Rp <- Rp_reflector(p); Z <- matrix(rnorm(n*p),n,p)
  X <- Z %*% t(Rp); Y <- Z[,1]*rnorm(n)
  list(X=X, Y=Y, d=1L, B_true=matrix(Rp[,1],p,1))
}

gen_M4 <- function(n,p) {
  Rp <- Rp_reflector(p); Z <- matrix(rnorm(n*p),n,p)
  X <- Z %*% t(Rp); Y <- Z[,1] + 0.75*Z[,2]*rnorm(n)
  list(X=X, Y=Y, d=2L, B_true=Rp[,1:2,drop=FALSE])
}

gen_M5 <- function(n,p) {
  Sigma <- outer(seq_len(p),seq_len(p), function(i,j) 0.5^abs(i-j))
  X <- mvrnorm(n,rep(0,p),Sigma)
  b1 <- rep(0,p); b1[c(1,2,p)] <- c(1,2,2); b1 <- b1/3
  b2 <- rep(0,p); b2[c(3,4)] <- c(3,4); b2 <- b2/5
  Y <- 2*(X%*%b1)^2 + exp(X%*%b2)*rnorm(n)
  list(X=X, Y=as.numeric(Y), d=2L, B_true=cbind(b1,b2))
}

gen_M6 <- function(n,p) {
  V1<-rnorm(n);V2<-rnorm(n);V3<-rnorm(n);W1<-rgamma(n,0.2,1);W2<-rgamma(n,0.2,1)
  X <- matrix(0,n,p)
  X[,1]<-W1; X[,2]<-V1+W2/2; X[,3]<- -V1+W2/2; X[,4]<-V2+V3; X[,5]<-V2-V3
  if(p>5) X[,6:p] <- matrix(rnorm(n*(p-5)),n,p-5)
  Y <- 1.5*(5+X[,1])*(2+X[,2]+X[,3]) + 0.5*rnorm(n)
  B <- matrix(0,p,2); B[1,1]<-1; B[2,2]<-1; B[3,2]<-1
  list(X=X, Y=as.numeric(Y), d=2L, B_true=B)
}

gen_M7 <- function(n,p,c=0) {
  base <- gen_M3(n,p); Rp <- Rp_reflector(p); Z <- base$X %*% Rp
  m <- floor(c*n)
  if(m>0) { idx <- sample.int(n,m); Z[idx,2:p] <- matrix(runif(m*(p-1),-20,20),m,p-1) }
  list(X=Z %*% t(Rp), Y=base$Y, d=1L, B_true=base$B_true)
}

gen_M8 <- function(n,p) {
  gam<-0.75; tau<-sqrt(1+2*gam+3*gam^2); Rp<-Rp_reflector(p)
  A<-rnorm(n); eta<-rnorm(n); Z <- matrix(rnorm(n*p),n,p)
  Z[,1]<-A; Z[,2]<-(1+gam*A^2)*eta/tau; X<-Z%*%t(Rp); Y<-A+0.5*rnorm(n)
  list(X=X, Y=as.numeric(Y), d=1L, B_true=matrix(Rp[,1],p,1))
}

# NOTE: gen_model dispatcher replaced by new_models.R block below



# ----- helper: 從有限混合分佈抽樣 -----
rmix <- function(n, comps, probs) {
  # comps: list of functions, each takes n and returns n samples
  # probs: mixing probabilities (sum to 1)
  k <- sample.int(length(probs), n, replace=TRUE, prob=probs)
  out <- numeric(n)
  for (j in seq_along(probs)) {
    idx <- which(k == j)
    if (length(idx)) out[idx] <- comps[[j]](length(idx))
  }
  out
}

# ============================================================================
# Model I (Xia 2007): Quadratic Location-Scale.  d = 2
#   X ~ N_p(0, Σ), Σ_ij = 0.5^|i-j|
#   Y = 2 (β1'X)^2 + exp(β2'X) ε
#   β1 = (1, 2, 0, ..., 0, 2)/3, β2 = (0, 0, 3, 4, 0, ..., 0)/5
# ============================================================================
gen_ModelI <- function(n, p) {
  Sigma <- outer(seq_len(p), seq_len(p), function(i, j) 0.5^abs(i - j))
  X <- MASS::mvrnorm(n, rep(0, p), Sigma)
  b1 <- rep(0, p); b1[c(1, 2, p)] <- c(1, 2, 2); b1 <- b1 / 3
  b2 <- rep(0, p); b2[c(3, 4)] <- c(3, 4); b2 <- b2 / 5
  Y <- 2 * (X %*% b1)^2 + exp(X %*% b2) * rnorm(n)
  list(X = X, Y = as.numeric(Y), d = 2L, B_true = cbind(b1, b2))
}

# ============================================================================
# Model L15 (Cook & Ni 2005): Skewed-predictor product.  d = 2
#   V1, V2, V3 ~ N(0,1); W1, W2 ~ Gamma(0.2, 1) i.i.d.
#   X1 = W1, X2 = V1 + W2/2, X3 = -V1 + W2/2
#   X4 = V2 + V3, X5 = V2 - V3
#   X6...Xp ~ N(0,1) independent
#   Y = 1.5 (5 + X1)(2 + X2 + X3) + 0.5 ε
#   central subspace = span(e1, e2+e3)
# ============================================================================
gen_ModelL15 <- function(n, p) {
  V1 <- rnorm(n); V2 <- rnorm(n); V3 <- rnorm(n)
  W1 <- rgamma(n, 0.2, 1); W2 <- rgamma(n, 0.2, 1)
  X <- matrix(0, n, p)
  X[, 1] <- W1
  X[, 2] <- V1 + W2 / 2
  X[, 3] <- -V1 + W2 / 2
  X[, 4] <- V2 + V3
  X[, 5] <- V2 - V3
  if (p > 5) X[, 6:p] <- matrix(rnorm(n * (p - 5)), n, p - 5)
  Y <- 1.5 * (5 + X[, 1]) * (2 + X[, 2] + X[, 3]) + 0.5 * rnorm(n)
  B <- matrix(0, p, 2); B[1, 1] <- 1; B[2, 2] <- 1; B[3, 2] <- 1
  list(X = X, Y = as.numeric(Y), d = 2L, B_true = B)
}

# ============================================================================
# Model VII: Fixed predictor contamination (5%).  d = 1
#   X ~ N_p(0, Σ), Σ_ij = 0.5^|i-j|
#   Y = (β1' X) ε,  β1 same as Model I
#   前 floor(0.05 n) 個觀測，X_{i,5}...X_{i,p} 換成 Uniform(40, 60)
# ============================================================================
gen_ModelVII <- function(n, p, contam = 0.05) {
  Sigma <- outer(seq_len(p), seq_len(p), function(i, j) 0.5^abs(i - j))
  X <- MASS::mvrnorm(n, rep(0, p), Sigma)
  b1 <- rep(0, p); b1[c(1, 2, p)] <- c(1, 2, 2); b1 <- b1 / 3
  Y <- as.numeric((X %*% b1) * rnorm(n))
  m <- floor(contam * n)
  if (m > 0 && p >= 5) {
    X[seq_len(m), 5:p] <- matrix(runif(m * (p - 4), 40, 60), m, p - 4)
  }
  list(X = X, Y = Y, d = 1L, B_true = matrix(b1, p, 1))
}

# Model C1 = Model VII with 10% contamination  (d = 1)
gen_ModelC1 <- function(n, p) gen_ModelVII(n, p, contam = 0.10)

# Model C2 = Model VII with 15% contamination  (d = 1)
gen_ModelC2 <- function(n, p) gen_ModelVII(n, p, contam = 0.15)

# ============================================================================
# Model II (Zhu et al. 2006): Linear mean + heteroscedastic error.  d = 2
#   X ~ N_p(0, Σ), Σ_ij = 0.5^|i-j|
#   β1 = (1,1,1,1,0,...,0)/2, β2 = (0,...,0,1,1,1,1)/2
#   Y = β1'X + (β2'X)^2 * 0.5 * ε
# ============================================================================
gen_ModelII <- function(n, p) {
  stopifnot(p >= 8)
  Sigma <- outer(seq_len(p), seq_len(p), function(i, j) 0.5^abs(i - j))
  X <- MASS::mvrnorm(n, rep(0, p), Sigma)
  b1 <- rep(0, p); b1[1:4] <- 1; b1 <- b1 / 2
  b2 <- rep(0, p); b2[(p - 3):p] <- 1; b2 <- b2 / 2
  Y <- as.numeric(X %*% b1 + (X %*% b2)^2 * 0.5 * rnorm(n))
  list(X = X, Y = Y, d = 2L, B_true = cbind(b1, b2))
}

# ============================================================================
# Model V (Li & Wang 2007): Non-symmetric oscillation.  d = 2
#   同 Model II 的 X
#   β1 = (1,1,1,0,...,0)/√3, β2 = (0,0,0,1,3,0,...,0)/√10
#   Y = 3 sin(β1'X / 4) + 0.2 (1 + (β2'X)^2) ε
# ============================================================================
gen_ModelV <- function(n, p) {
  stopifnot(p >= 5)
  Sigma <- outer(seq_len(p), seq_len(p), function(i, j) 0.5^abs(i - j))
  X <- MASS::mvrnorm(n, rep(0, p), Sigma)
  b1 <- rep(0, p); b1[1:3] <- 1; b1 <- b1 / sqrt(3)
  b2 <- rep(0, p); b2[4] <- 1; b2[5] <- 3; b2 <- b2 / sqrt(10)
  Y <- as.numeric(3 * sin((X %*% b1) / 4) + 0.2 * (1 + (X %*% b2)^2) * rnorm(n))
  list(X = X, Y = Y, d = 2L, B_true = cbind(b1, b2))
}

# ============================================================================
# Model TAQ3: Two-index quantile-width design.  d = 2
#   G ~ uniform on {1,...,5}, Y = G + ε_Y with ε_Y ~ Uniform(-0.2, 0.2)
#   A_j | G = h ~ mixture:
#     (1 - e_hj) Uniform[-a_hj, a_hj]  (interior body)
#     + (e_hj/2) [ Uniform[-3, -2.9] + Uniform[2.9, 3] ]  (tail masses)
#   Envelope 固定在 [-3, 3]，但 interior quantile widths 隨 slice 變
#   Z = R (A, W)^T,  X = Z (already standardized)
#   B_true = 前 d columns of R
#
#   Wu et al. (2020) 沒有給明確 e_hj, a_hj 常數；
#   這裡用 a 隨 slice 遞增 (0.5, 1.0, 1.5, 2.0, 2.5)，
#   e_hj = 0.10 (10% tail 質量在 |A|>2.9 的位置)
# ============================================================================
gen_ModelTAQ3 <- function(n, p) {
  stopifnot(p >= 2)
  Rp <- Rp_reflector(p)
  G <- sample.int(5L, n, replace = TRUE)
  Y <- G + runif(n, -0.2, 0.2)
  a_h <- c(0.5, 1.0, 1.5, 2.0, 2.5)     # interior half-width per slice
  e_h <- rep(0.10, 5)                    # tail mass per slice
  gen_Aj <- function() {
    a <- a_h[G]; e <- e_h[G]
    U <- runif(n)
    which_body <- U < (1 - e)
    which_pos <- U >= (1 - e) & U < (1 - e / 2)
    A <- numeric(n)
    A[which_body] <- runif(sum(which_body), -a[which_body], a[which_body])
    A[which_pos]  <- runif(sum(which_pos), 2.9, 3.0)
    A[!which_body & !which_pos] <- runif(sum(!which_body & !which_pos), -3.0, -2.9)
    A
  }
  A1 <- gen_Aj(); A2 <- gen_Aj()
  W <- matrix(rnorm(n * (p - 2)), n, p - 2)
  Zmat <- cbind(A1, A2, W)
  X <- Zmat %*% t(Rp)
  list(X = X, Y = Y, d = 2L, B_true = Rp[, 1:2, drop = FALSE])
}

# ============================================================================
# Model TAM2: Quartic endpoint-reach design.  d = 1
#   Y ~ Uniform[-1, 1]
#   常數: a = 0.5, e = 0.18, η = 0.05
#   b(v) = sqrt( (v - (1-e) a^2/3) / e - η^2/3 )
#   v(y) = 1 + 1.5 (y^4 - 1/5)
#   A | Y = y ~ 混合:
#     (1 - e) Uniform[-a, a]
#     + (e/2) Uniform[b(v(y))-η, b(v(y))+η]
#     + (e/2) Uniform[-b(v(y))-η, -b(v(y))+η]
#   中心 quantile 半寬固定，但邊界 reach 隨 y 變
# ============================================================================
gen_ModelTAM2 <- function(n, p) {
  Rp <- Rp_reflector(p)
  a_const <- 0.5; e_const <- 0.18; eta_const <- 0.05
  Y <- runif(n, -1, 1)
  v_of <- function(y) 1 + 1.5 * (y^4 - 1/5)
  b_of <- function(v) sqrt(pmax(0, (v - (1 - e_const) * a_const^2 / 3) / e_const - eta_const^2 / 3))
  v_y <- v_of(Y); b_y <- b_of(v_y)
  U <- runif(n)
  A <- numeric(n)
  in_body <- U < (1 - e_const)
  in_pos  <- U >= (1 - e_const) & U < (1 - e_const / 2)
  in_neg  <- !in_body & !in_pos
  A[in_body] <- runif(sum(in_body), -a_const, a_const)
  A[in_pos]  <- runif(sum(in_pos),  b_y[in_pos] - eta_const, b_y[in_pos] + eta_const)
  A[in_neg]  <- runif(sum(in_neg), -b_y[in_neg] - eta_const, -b_y[in_neg] + eta_const)
  W <- matrix(rnorm(n * (p - 1)), n, p - 1)
  Zmat <- cbind(A, W)
  X <- Zmat %*% t(Rp)
  list(X = X, Y = Y, d = 1L, B_true = matrix(Rp[, 1], p, 1))
}

# ============================================================================
# Model TAM3: Two-direction endpoint-reach design.  d = 2
#   A1 = 1.5 Y + E1,  E1 ~ Uniform[-sqrt(0.75), sqrt(0.75)]      (location)
#   A2 | Y 同 TAM2 但 v(y) = 1 + 1.5 (y^2 - 1/3)                (endpoint reach)
#   其他同 TAM2
# ============================================================================
gen_ModelTAM3 <- function(n, p) {
  stopifnot(p >= 2)
  Rp <- Rp_reflector(p)
  a_const <- 0.5; e_const <- 0.18; eta_const <- 0.05
  Y <- runif(n, -1, 1)
  # direction 1: location
  A1 <- 1.5 * Y + runif(n, -sqrt(0.75), sqrt(0.75))
  # direction 2: endpoint-reach with v(y) = 1 + 1.5(y^2 - 1/3)
  v_of <- function(y) 1 + 1.5 * (y^2 - 1/3)
  b_of <- function(v) sqrt(pmax(0, (v - (1 - e_const) * a_const^2 / 3) / e_const - eta_const^2 / 3))
  v_y <- v_of(Y); b_y <- b_of(v_y)
  U <- runif(n)
  A2 <- numeric(n)
  in_body <- U < (1 - e_const)
  in_pos  <- U >= (1 - e_const) & U < (1 - e_const / 2)
  in_neg  <- !in_body & !in_pos
  A2[in_body] <- runif(sum(in_body), -a_const, a_const)
  A2[in_pos]  <- runif(sum(in_pos),  b_y[in_pos] - eta_const, b_y[in_pos] + eta_const)
  A2[in_neg]  <- runif(sum(in_neg), -b_y[in_neg] - eta_const, -b_y[in_neg] + eta_const)
  W <- matrix(rnorm(n * (p - 2)), n, p - 2)
  Zmat <- cbind(A1, A2, W)
  X <- Zmat %*% t(Rp)
  list(X = X, Y = Y, d = 2L, B_true = Rp[, 1:2, drop = FALSE])
}

# ============================================================================
# 註冊到 gen_model dispatcher
# 需要更新原本的 gen_model 函式，加入新模型
# ============================================================================
gen_model <- function(name, n, p, c = 0) switch(name,
                                                # 原本的 Standalone models
                                                M0 = gen_M0(n, p),  M1 = gen_M1(n, p),  M2 = gen_M2(n, p),  M3 = gen_M3(n, p),
                                                M4 = gen_M4(n, p),  M5 = gen_M5(n, p),  M6 = gen_M6(n, p),
                                                M7 = gen_M7(n, p, c),  M8 = gen_M8(n, p),
                                                # 主論文 Group 1
                                                ModelI   = gen_ModelI(n, p),
                                                ModelL15 = gen_ModelL15(n, p),
                                                ModelVII = gen_ModelVII(n, p, contam = 0.05),
                                                # 主論文 Group 2
                                                TAQ3 = gen_ModelTAQ3(n, p),
                                                TAM2 = gen_ModelTAM2(n, p),
                                                TAM3 = gen_ModelTAM3(n, p),
                                                # 主論文 Group 3
                                                C1 = gen_ModelC1(n, p),
                                                C2 = gen_ModelC2(n, p),
                                                # 主論文 Group 4
                                                ModelII = gen_ModelII(n, p),
                                                ModelV  = gen_ModelV(n, p),
                                                stop("unknown model: ", name))

# ============================================================================
# 快速測試每個新模型
# ============================================================================
test_new_models <- function() {
  set.seed(1)
  for (nm in c("ModelI", "ModelL15", "ModelVII", "TAQ3", "TAM2", "TAM3",
               "C1", "C2", "ModelII", "ModelV")) {
    dat <- gen_model(nm, n = 200, p = 10)
    cat(sprintf("%-9s  n=%d  p=%d  d=%d  B_true dim=%s  Y range=[%.2f, %.2f]\n",
                nm, nrow(dat$X), ncol(dat$X), dat$d,
                paste(dim(dat$B_true), collapse="x"),
                min(dat$Y), max(dat$Y)))
  }
}


# ---- alpha per model ----
model_alpha <- function(name) {
  # 只有 M1 用 alpha=0.50，其他所有 model (含主論文新增) 用 alpha=0.20
  if (name == "M1") 0.50 else 0.20
}


########################################################################
#  PART 6 ── 模擬引擎（一次 replication + checkpoint runner）
########################################################################

METHODS_ALL       <- c("SIR","SAVE","pHd","DR","CUME","csMAVE",
                       "SIIE-Q(H)","SIIE-Q(HK)","DSIIE")
METHODS_FAST      <- setdiff(METHODS_ALL, "csMAVE")
METHODS_INTERVAL  <- c("SIIE-Q(H)","SIIE-Q(HK)","DSIIE")

one_rep <- function(model_id, design_id, rep, seed, n, p, H, K, alpha,
                    c_val=0, methods=METHODS_ALL) {
  set.seed(seed)
  dat <- gen_model(model_id, n, p, c=c_val)
  d <- if(is.null(dat$d)||dat$d==0L) 1L else dat$d
  fits <- fit_all(dat$X, dat$Y, d, H, K, alpha, methods)
  P_true <- if(is.null(dat$B_true)) NULL else tcrossprod(qr.Q(qr(dat$B_true)))
  do.call(rbind, lapply(methods, function(nm) {
    fit <- fits[[nm]]; fail <- isTRUE(fit$fail)||is.null(fit$projector)
    ra <- if(fail||is.null(P_true)) NA_real_ else r_ave(P_true, fit$projector)
    rm <- if(fail||is.null(dat$B_true)||NCOL(dat$B_true)<2) NA_real_ else r_min(dat$B_true, fit$projector)
    data.frame(model_id=model_id, design_id=design_id, replication=rep,
               data_seed=seed, method_id=nm, n=n, p=p, H=H, K=K,
               alpha=alpha, c=c_val, d=d,
               r_ave=ra, r_min=rm,
               fit_status=if(fail)"fail"else"success",
               stringsAsFactors=FALSE)
  }))
}

run_cell <- function(model, design_id, n, p, H, K, alpha, c_val=0, R,
                     methods, seed_base=0L, tag="") {
  fname <- sprintf("results/cell_%s_%s_n%d_p%d_H%d_K%d_a%03d_c%03d%s.rds",
                   design_id, model, n, p, H, K,
                   round(alpha*100), round(c_val*100),
                   if(nzchar(tag)) paste0("_",tag) else "")
  if (file.exists(fname)) { cat("SKIP:", basename(fname), "\n"); return(invisible()) }
  cat(sprintf("[%s] %s %s n=%d p=%d H=%d K=%d a=%.2f c=%.2f R=%d methods=%d\n",
              format(Sys.time(),"%H:%M:%S"), design_id, model, n, p, H, K,
              alpha, c_val, R, length(methods)))
  recs <- vector("list", R); t0 <- Sys.time()
  for (r in seq_len(R)) {
    recs[[r]] <- tryCatch(one_rep(model, design_id, r, seed_base+r, n, p, H, K,
                                  alpha, c_val, methods=methods),
                          error=function(e) NULL)
  }
  out <- do.call(rbind, recs); saveRDS(out, fname)
  cat(sprintf("  saved (%d rows) elapsed=%s\n", nrow(out), format(Sys.time()-t0)))
}

load_records <- function() {
  files <- list.files("results", pattern="^cell_.*\\.rds$", full.names=TRUE)
  if(!length(files)) stop("no cell files"); do.call(rbind, lapply(files, readRDS))
}

# alias: run_paper_models() 內部呼叫 run_cell2，與 run_cell 功能完全相同
run_cell2 <- run_cell


########################################################################
#  PART 7 ── 跑模擬
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  7-A) Table 3 + Table 5 + Table 6 的資料來源：SIM-PRIMARY
#       M1–M8, n=400, p=10, H=4, K=2
########################################################################

run_for_table3_5_6 <- function() {
  cat("\n========== Table 3 / 5 / 6: SIM-PRIMARY ==========\n")
  models <- c("M1","M2","M3","M4","M5","M6","M7","M8")
  # fast methods (R=50)
  for (i in seq_along(models)) {
    m <- models[i]; alpha <- model_alpha(m); cv <- if(m=="M7") 0.10 else 0
    run_cell(m, "SIM-PRIMARY", 400L, 10L, 4L, 2L, alpha, cv,
             R=R_PRIMARY, methods=METHODS_ALL, seed_base=10000L+i*1000L, tag="fast")
  }
  # csMAVE (R=15, 因為慢)
  for (i in seq_along(models)) {
    m <- models[i]; alpha <- model_alpha(m); cv <- if(m=="M7") 0.10 else 0
    run_cell(m, "SIM-PRIMARY", 400L, 10L, 4L, 2L, alpha, cv,
             R=R_CSMAVE, methods="csMAVE", seed_base=15000L+i*1000L, tag="csmave")
  }
}

########################################################################
#  7-B) Table 4 的資料來源：SIM-P50
#       M1–M4, n=800, p=50, H=4, K=2
########################################################################

run_for_table4 <- function() {
  cat("\n========== Table 4: SIM-P50 ==========\n")
  for (i in seq_along(c("M1","M2","M3","M4"))) {
    m <- c("M1","M2","M3","M4")[i]; alpha <- model_alpha(m)
    run_cell(m, "SIM-P50", 800L, 50L, 4L, 2L, alpha, 0,
             R=R_P50, methods=METHODS_ALL, seed_base=20000L+i*1000L)
  }
}

########################################################################
#  7-C) Table 7 的資料來源：M0 null permutation test
########################################################################

run_for_table7 <- function() {
  cat("\n========== Table 7: M0 null calibration ==========\n")
  frob_stat <- numeric(R_NULL); reject <- logical(R_NULL)
  for (r in seq_len(R_NULL)) {
    set.seed(90000+r); dat <- gen_M0(400L,10L)
    fit <- try(dsiie(dat$X, dat$Y, d=1L, H=4L, K=2L, alpha=0.20), silent=TRUE)
    if(inherits(fit,"try-error")) { frob_stat[r]<-NA; reject[r]<-NA; next }
    Tobs <- sqrt(sum(fit$M^2))
    Tperm <- numeric(B_NULL)
    for(b in seq_len(B_NULL)) {
      fb <- try(dsiie(dat$X, sample(dat$Y), d=1L, H=4L, K=2L, alpha=0.20), silent=TRUE)
      Tperm[b] <- if(inherits(fb,"try-error")) NA else sqrt(sum(fb$M^2))
    }
    frob_stat[r] <- Tobs
    pval <- (1+sum(Tperm>=Tobs,na.rm=TRUE))/(1+sum(!is.na(Tperm)))
    reject[r] <- pval <= 0.05
    if(r%%20==0) cat(sprintf("  M0 rep %d/%d  reject=%.3f\n", r, R_NULL, mean(reject[1:r],na.rm=TRUE)))
  }
  res <- list(frob_stat=frob_stat, reject=reject,
              rejection_rate=mean(reject,na.rm=TRUE),
              rejection_se=sqrt(mean(reject,na.rm=TRUE)*(1-mean(reject,na.rm=TRUE))/R_NULL))
  saveRDS(res, "results/M0_null.rds")
  cat(sprintf("  Rejection rate = %.3f\n", res$rejection_rate))
}

########################################################################
#  7-D) Table 8 的資料來源：runtime
########################################################################

run_for_table8 <- function() {
  cat("\n========== Table 8: Runtime (M2 + M4) ==========\n")
  cases <- list(c(400L,10L), c(800L,50L), c(2000L,50L))
  meths <- c("DSIIE","SIIE-Q(H)","SIIE-Q(HK)")
  gens  <- list(M2 = gen_M2, M4 = gen_M4)
  out <- list()
  for (mdl in names(gens)) {
    gen_fun <- gens[[mdl]]
    for(cs in cases) { n<-cs[1]; p<-cs[2]
    for(nm in meths) {
      dat <- gen_fun(n,p)  # warm up
      switch(nm, "DSIIE"=try(dsiie(dat$X,dat$Y,1L,4L,2L,0.20),silent=TRUE),
             "SIIE-Q(H)"=try(siie_q(dat$X,dat$Y,1L,4L,0.20),silent=TRUE),
             "SIIE-Q(HK)"=try(siie_q(dat$X,dat$Y,1L,8L,0.20),silent=TRUE))
      times <- numeric(R_RUNTIME)
      for(r in seq_len(R_RUNTIME)) {
        set.seed(99000+n+p+which(meths==nm)*7L+r+which(names(gens)==mdl)*131L)
        dat <- gen_fun(n,p)
        t0 <- proc.time()[[3]]
        switch(nm, "DSIIE"=try(dsiie(dat$X,dat$Y,1L,4L,2L,0.20),silent=TRUE),
               "SIIE-Q(H)"=try(siie_q(dat$X,dat$Y,1L,4L,0.20),silent=TRUE),
               "SIIE-Q(HK)"=try(siie_q(dat$X,dat$Y,1L,8L,0.20),silent=TRUE))
        times[r] <- proc.time()[[3]]-t0
      }
      key <- paste0(mdl,"_",n,"_",p,"_",nm)
      out[[key]] <- list(model=mdl,n=n,p=p,method=nm,
                         median=median(times),iqr=IQR(times))
      cat(sprintf("  %s n=%d p=%d %s median=%.4fs\n", mdl, n, p, nm, median(times)))
    }
    }
  }
  saveRDS(out, "results/runtime.rds")
}

########################################################################
#  7-E) Figure 1 的資料來源：SIM-N trajectories
#       M1–M4, n ∈ {200, 400, 800, 1600}
########################################################################

run_for_figure1 <- function() {
  cat("\n========== Figure 1: SIM-N trajectories ==========\n")
  for (i in seq_along(c("M1","M2","M3","M4"))) {
    m <- c("M1","M2","M3","M4")[i]; alpha <- model_alpha(m)
    for (j in seq_along(c(200L,400L,800L,1600L))) {
      n_i <- c(200L,400L,800L,1600L)[j]
      run_cell(m, "SIM-N", n_i, 10L, 4L, 2L, alpha, 0,
               R=R_TRAJ, methods=METHODS_ALL,
               seed_base=30000L+i*10000L+j*1000L)
    }
  }
}

########################################################################
#  7-F) Figure 2 的資料來源：四個 sensitivity study
########################################################################

run_for_figure2 <- function() {
  cat("\n========== Figure 2 Panel A: H sensitivity (M3) ==========\n")
  for (Hv in c(2L,3L,4L,5L))
    run_cell("M3","SIM-H",400L,10L,Hv,2L,0.20,0,R=R_SENS_H,
             methods=METHODS_INTERVAL, seed_base=40000L+Hv*1000L)
  
  cat("\n========== Figure 2 Panel B: K sensitivity (M1) ==========\n")
  for (Kv in c(1L,2L,3L,4L))
    run_cell("M1","SIM-K",400L,10L,4L,Kv,0.50,0,R=R_SENS_K,
             methods=METHODS_INTERVAL, seed_base=50000L+Kv*1000L)
  
  cat("\n========== Figure 2 Panel C: alpha sensitivity (M1, M2) ==========\n")
  for (m in c("M1","M2"))
    for (av in c(0.10,0.20,0.30,0.50))
      run_cell(m,"SIM-A",400L,10L,4L,2L,av,0,R=R_SENS_A,
               methods=METHODS_INTERVAL,
               seed_base=60000L+which(c("M1","M2")==m)*10000L+as.integer(av*100)*100L)
  
  cat("\n========== Figure 2 Panel D: contamination (M7) ==========\n")
  for (cv in c(0,0.05,0.10))
    run_cell("M7","SIM-C",400L,10L,4L,2L,0.20,cv,R=R_SENS_C,
             methods=METHODS_ALL, seed_base=70000L+as.integer(cv*100)*100L)
}

########################################################################
#  7-G) 附加：SIM-RANK（Table 4 的 rank-deficient stress test）
########################################################################

run_for_rank <- function() {
  cat("\n========== SIM-RANK ==========\n")
  for (i in seq_along(c("M1","M2","M3","M4"))) {
    m <- c("M1","M2","M3","M4")[i]; alpha <- model_alpha(m)
    run_cell(m,"SIM-RANK",400L,50L,4L,2L,alpha,0,R=R_RANK,
             methods=METHODS_ALL, seed_base=80000L+i*1000L)
  }
}


########################################################################
#  PART 8 ── 產生 LaTeX 表格與 PDF 圖
########################################################################

fmt <- function(x, digits=3) formatC(x, digits=digits, format="f")
mean_sd <- function(v) { v<-v[is.finite(v)]; list(mean=mean(v),sd=sd(v),n=length(v)) }

# ---------- Table 3 (tab:primary): §7.1 ----------
# 產生四個 sub-table (a)(b)(c)(d)，對應 Standalone Table 3 的 2x2 分類
# 使用者可以修改 TABLE3_GROUPS 自訂分組
TABLE3_GROUPS <- list(
  a = list(
    caption = "(a) Standalone base models (M1--M4)",
    models  = c("M1","M2","M3","M4")),
  b = list(
    caption = "(b) Standalone extended models (M5--M8)",
    models  = c("M5","M6","M7","M8")),
  c = list(
    caption = "(c) Paper literature/contamination models (Groups 1 \\& 3)",
    models  = c("ModelI","ModelL15","ModelVII","C1","C2")),
  d = list(
    caption = "(d) Paper theory-aligned/heteroscedastic (Groups 2 \\& 4)",
    models  = c("TAQ3","TAM2","TAM3","ModelII","ModelV"))
)

write_table3 <- function(recs) {
  cat("Writing Table 3 sub-tables (tab_primary_a/b/c/d.tex) ...\n")
  methods <- METHODS_ALL
  sub <- recs[recs$design_id=="SIM-PRIMARY",]
  
  # 為每個 group 產生獨立的 .tex body
  for (gid in names(TABLE3_GROUPS)) {
    g <- TABLE3_GROUPS[[gid]]
    models_g <- g$models
    f <- sprintf("tex/tab_primary_%s.tex", gid)
    cat(sprintf("%% Table 3 sub-table (%s): %s\n", gid, g$caption), file=f)
    for(m in models_g) {
      mean_cells <- sapply(methods, function(nm) {
        v <- sub$r_ave[sub$model_id==m & sub$method_id==nm]; v <- v[is.finite(v)]
        if (length(v)==0) "--" else fmt(mean(v))
      })
      sd_cells <- sapply(methods, function(nm) {
        v <- sub$r_ave[sub$model_id==m & sub$method_id==nm]; v <- v[is.finite(v)]
        if (length(v)==0) "--" else sprintf("(%s)", fmt(sd(v)))
      })
      cat(m, paste(mean_cells, collapse=" & "), "\\\\\n",
          sep=" & ", file=f, append=TRUE)
      cat("   ", paste(sd_cells, collapse=" & "), "\\\\\n",
          sep=" & ", file=f, append=TRUE)
    }
    cat(sprintf("  saved %s (%d models)\n", f, length(models_g)))
  }
  
  # 也產生完整的 tab_primary.tex (所有模型合在一起)，作為 backup
  all_models <- unique(unlist(lapply(TABLE3_GROUPS, function(g) g$models)))
  # 只保留有資料的
  all_models <- all_models[all_models %in% unique(sub$model_id)]
  f <- "tex/tab_primary.tex"
  cat("% Table 3: 全部模型 (可用作合併版)\n", file=f)
  for(m in all_models) {
    mean_cells <- sapply(methods, function(nm) {
      v <- sub$r_ave[sub$model_id==m & sub$method_id==nm]; v <- v[is.finite(v)]
      if (length(v)==0) "--" else fmt(mean(v))
    })
    sd_cells <- sapply(methods, function(nm) {
      v <- sub$r_ave[sub$model_id==m & sub$method_id==nm]; v <- v[is.finite(v)]
      if (length(v)==0) "--" else sprintf("(%s)", fmt(sd(v)))
    })
    cat(m, paste(mean_cells, collapse=" & "), "\\\\\n",
        sep=" & ", file=f, append=TRUE)
    cat("   ", paste(sd_cells, collapse=" & "), "\\\\\n",
        sep=" & ", file=f, append=TRUE)
  }
  cat(sprintf("  saved tex/tab_primary.tex (%d models total)\n", length(all_models)))
}

# ---------- Table 4 (tab:p50): §7.2 ----------
write_table4 <- function(recs) {
  cat("Writing Table 4 (tab:p50) ...\n")
  methods <- METHODS_ALL
  # 原本的 Standalone 模型 + 主論文 10 個新模型
  models <- c("M1","M2","M3","M4",
              "ModelI","ModelL15","ModelVII",
              "TAQ3","TAM2","TAM3",
              "C1","C2","ModelII","ModelV")
  M <- array(NA_real_, dim=c(length(models),9,2), dimnames=list(models,methods,c("mean","sd")))
  sub <- recs[recs$design_id=="SIM-P50",]
  for(m in models) for(nm in methods) {
    v <- sub$r_ave[sub$model_id==m & sub$method_id==nm]; s <- mean_sd(v)
    M[m,nm,"mean"]<-s$mean; M[m,nm,"sd"]<-s$sd
  }
  f <- "tex/tab_p50.tex"; cat("% Table 4: p=50 study n=800\n",file=f)
  for(m in models) {
    cat(m, paste(sapply(methods, function(nm)
      if(is.na(M[m,nm,"mean"])) "--" else fmt(M[m,nm,"mean"])), collapse=" & "),
      "\\\\\n", sep=" & ", file=f, append=TRUE)
    cat("   ", paste(sapply(methods, function(nm)
      if(is.na(M[m,nm,"sd"])) "--" else sprintf("(%s)",fmt(M[m,nm,"sd"]))), collapse=" & "),
      "\\\\\n", sep=" & ", file=f, append=TRUE)
  }
}

# ---------- Table 5 (tab:paired): §7.3 ----------
write_table5 <- function(recs) {
  cat("Writing Table 5 (tab:paired) ...\n")
  sub <- recs[recs$design_id=="SIM-PRIMARY" &
                recs$method_id %in% c("DSIIE","SIIE-Q(H)","SIIE-Q(HK)"),]
  f <- "tex/tab_paired.tex"; cat("% Table 5: Paired DSIIE contrasts\n",file=f)
  # 原本 M1-M8 + 主論文 10 個新模型
  paired_models <- c(paste0("M",1:8),
                     "ModelI","ModelL15","ModelVII",
                     "TAQ3","TAM2","TAM3",
                     "C1","C2","ModelII","ModelV")
  for(m in paired_models) {
    d1<-sub[sub$model_id==m & sub$method_id=="DSIIE",]
    d2<-sub[sub$model_id==m & sub$method_id=="SIIE-Q(H)",]
    d3<-sub[sub$model_id==m & sub$method_id=="SIIE-Q(HK)",]
    k <- intersect(intersect(d1$data_seed,d2$data_seed),d3$data_seed)
    if(!length(k)) next
    dC <- d1$r_ave[match(k,d1$data_seed)] - d2$r_ave[match(k,d2$data_seed)]
    dF <- d1$r_ave[match(k,d1$data_seed)] - d3$r_ave[match(k,d3$data_seed)]
    dC<-dC[is.finite(dC)]; dF<-dF[is.finite(dF)]; nc<-length(dC); nf<-length(dF)
    mc<-mean(dC); sc<-sd(dC); mf<-mean(dF); sf<-sd(dF)
    cat(sprintf("%s & %s\\;[%s,\\,%s] & %s\\;[%s,\\,%s] & %d \\\\\n",
                m, fmt(mc),fmt(mc-1.96*sc/sqrt(nc)),fmt(mc+1.96*sc/sqrt(nc)),
                fmt(mf),fmt(mf-1.96*sf/sqrt(nf)),fmt(mf+1.96*sf/sqrt(nf)),
                min(nc,nf)), file=f, append=TRUE)
  }
}

# ---------- Table 6 (tab:rmin): §7.4 ----------
write_table6 <- function(recs) {
  cat("Writing Table 6 (tab:rmin) ...\n")
  methods <- METHODS_ALL
  # 原本 Standalone 的 d=2 模型 + 主論文 6 個 d=2 模型
  models <- c("M4","M5","M6",
              "ModelI","ModelL15","TAQ3","TAM3","ModelII","ModelV")
  M <- array(NA_real_, dim=c(length(models),9,2), dimnames=list(models,methods,c("mean","sd")))
  sub <- recs[recs$design_id=="SIM-PRIMARY" & recs$model_id %in% models,]
  for(m in models) for(nm in methods) {
    v <- sub$r_min[sub$model_id==m & sub$method_id==nm]; s <- mean_sd(v)
    M[m,nm,"mean"]<-s$mean; M[m,nm,"sd"]<-s$sd
  }
  f <- "tex/tab_rmin.tex"; cat("% Table 6: r_min for M4-M6\n",file=f)
  for(m in models) {
    cat(m, paste(sapply(methods, function(nm)
      if(is.na(M[m,nm,"mean"])) "--" else fmt(M[m,nm,"mean"])), collapse=" & "),
      "\\\\\n", sep=" & ", file=f, append=TRUE)
    cat("   ", paste(sapply(methods, function(nm)
      if(is.na(M[m,nm,"sd"])) "--" else sprintf("(%s)",fmt(M[m,nm,"sd"]))), collapse=" & "),
      "\\\\\n", sep=" & ", file=f, append=TRUE)
  }
}

# ---------- Table 7 (tab:null): §7.5 ----------
write_table7 <- function() {
  cat("Writing Table 7 (tab:null) ...\n")
  n0 <- readRDS("results/M0_null.rds")
  rej<-n0$rejection_rate; se<-n0$rejection_se
  fro<-mean(n0$frob_stat,na.rm=TRUE); frose<-sd(n0$frob_stat,na.rm=TRUE)
  R_n<-length(n0$frob_stat)
  f<-"tex/tab_null.tex"; cat("% Table 7: Null calibration\n",file=f)
  cat(sprintf("DSIIE permutation rejection rate at 0.05 & %s & %s & [%s,\\,%s] & B=%d \\\\\n",
              fmt(rej),fmt(se),fmt(max(0,rej-1.96*se)),fmt(min(1,rej+1.96*se)),R_n), file=f,append=TRUE)
  cat(sprintf("Mean $\\|\\Mhat\\|_F$ under M0 & %s & %s & -- & R=%d \\\\\n",
              fmt(fro),fmt(frose/sqrt(R_n)),R_n), file=f,append=TRUE)
}

# ---------- Table 8 (tab:runtime): §7.5 ----------
write_table8 <- function() {
  cat("Writing Table 8 (tab:runtime) ...\n")
  rt <- readRDS("results/runtime.rds")
  f<-"tex/tab_runtime.tex"; cat("% Table 8: Runtime (M2 + M4)\n",file=f)
  # 相容舊格式 (無 model 前綴)
  has_prefix <- any(grepl("^M[0-9]+_", names(rt)))
  for (mdl in c("M2","M4")) {
    cat(sprintf("%% ---- Model %s ----\n", mdl), file=f, append=TRUE)
    for(n in c(400L,800L,2000L)) { p<-if(n==400L) 10L else 50L
    for(nm in c("DSIIE","SIIE-Q(H)","SIIE-Q(HK)")) {
      key <- if (has_prefix) paste0(mdl,"_",n,"_",p,"_",nm)
             else            paste0(n,"_",p,"_",nm)
      if(!key %in% names(rt)) next
      x<-rt[[key]]
      cat(sprintf("%s & $(%d,%d)$ & %s & %s & %s & -- \\\\\n",
                  mdl, n, p, nm, fmt(x$median), fmt(x$iqr)),
          file=f, append=TRUE)
    }
    }
    if (!has_prefix) break   # 舊檔只有 M2
  }
}

# ---------- Figure 1 (fig:trajectories): §7.6 ----------
write_figure1 <- function(recs) {
  cat("Writing Figure 1 (fig:trajectories) ...\n")
  sub <- recs[recs$design_id=="SIM-N",]; if(!nrow(sub)) return()
  ns <- c(200,400,800,1600)
  methods <- c("SIR","SAVE","pHd","DR","CUME","SIIE-Q(H)","SIIE-Q(HK)","DSIIE")
  cols <- c(SIR="black",SAVE="red",pHd="orange",DR="purple",
            CUME="brown",`SIIE-Q(H)`="darkgreen",`SIIE-Q(HK)`="steelblue",DSIIE="blue")
  ltys <- c(SIR=1,SAVE=2,pHd=3,DR=4,CUME=5,`SIIE-Q(H)`=1,`SIIE-Q(HK)`=2,DSIIE=1)
  lwds <- c(SIR=1,SAVE=1,pHd=1,DR=1,CUME=1,`SIIE-Q(H)`=1.4,`SIIE-Q(HK)`=1.4,DSIIE=2.2)
  pdf("figures/DSIIE_sample_size_trajectories.pdf", width=10, height=8)
  par(mfrow=c(2,2), mar=c(4,4,2.5,1))
  for(m in c("M1","M2","M3","M4")) {
    plot(NA,xlim=range(ns),ylim=c(0,1),log="x",xlab="n",
         ylab=expression(mean~r[ave]),main=m,xaxt="n")
    axis(1,at=ns,labels=ns); grid(col="grey85")
    for(nm in methods) {
      ys <- sapply(ns, function(nn) {
        v<-sub$r_ave[sub$model_id==m & sub$method_id==nm & sub$n==nn]
        mean(v[is.finite(v)])
      })
      lines(ns,ys,col=cols[nm],lty=ltys[nm],lwd=lwds[nm])
      points(ns,ys,col=cols[nm],pch=20,cex=0.8)
    }
    if(m=="M1") legend("bottomright",legend=methods,col=cols[methods],
                       lty=ltys[methods],lwd=lwds[methods],bty="n",cex=0.7,ncol=2)
  }
  dev.off()
}

# ---------- Figure 2 (fig:sensitivity): §7.6 ----------
write_figure2 <- function(recs) {
  cat("Writing Figure 2 (fig:sensitivity) ...\n")
  pdf("figures/DSIIE_sensitivity_panels.pdf", width=11, height=9)
  par(mfrow=c(2,2), mar=c(4,4,2.5,1))
  mi <- METHODS_INTERVAL
  cols <- c(DSIIE="blue",`SIIE-Q(H)`="darkgreen",`SIIE-Q(HK)`="steelblue")
  lwd <- c(DSIIE=2.2,`SIIE-Q(H)`=1.4,`SIIE-Q(HK)`=1.4)
  
  # Panel A: H
  sub <- recs[recs$design_id=="SIM-H",]
  plot(NA,xlim=c(2,5),ylim=c(0,1),xaxt="n",xlab="H",ylab=expression(mean~r[ave]),
       main="H sensitivity, M3"); axis(1,at=2:5); grid()
  for(nm in mi) {
    ys<-sapply(2:5, function(h) mean(sub$r_ave[sub$method_id==nm & sub$H==h],na.rm=TRUE))
    lines(2:5,ys,col=cols[nm],lwd=lwd[nm],lty=if(nm=="DSIIE")1 else 2)
    points(2:5,ys,col=cols[nm],pch=20)
  }
  legend("bottomright",legend=mi,col=cols[mi],lwd=lwd[mi],lty=c(1,2,2),bty="n")
  
  # Panel B: K
  sub <- recs[recs$design_id=="SIM-K",]
  plot(NA,xlim=c(1,4),ylim=c(0,1),xaxt="n",xlab="K",ylab=expression(mean~r[ave]),
       main="K sensitivity, M1"); axis(1,at=1:4); grid()
  for(nm in mi) {
    ys<-sapply(1:4, function(k) mean(sub$r_ave[sub$method_id==nm & sub$K==k],na.rm=TRUE))
    lines(1:4,ys,col=cols[nm],lwd=lwd[nm],lty=if(nm=="DSIIE")1 else 2)
    points(1:4,ys,col=cols[nm],pch=20)
  }
  
  # Panel C: alpha
  as_<-c(0.10,0.20,0.30,0.50); sub<-recs[recs$design_id=="SIM-A" & recs$model_id=="M1",]
  plot(NA,xlim=range(as_),ylim=c(0,1),xlab=expression(alpha),
       ylab=expression(mean~r[ave]),main="alpha sensitivity, M1"); grid()
  for(nm in mi) {
    ys<-sapply(as_, function(av) mean(sub$r_ave[sub$method_id==nm & abs(sub$alpha-av)<1e-6],na.rm=TRUE))
    lines(as_,ys,col=cols[nm],lwd=lwd[nm],lty=if(nm=="DSIIE")1 else 2)
    points(as_,ys,col=cols[nm],pch=20)
  }
  
  # Panel D: contamination
  cs_<-c(0,0.05,0.10); sub<-recs[recs$design_id=="SIM-C",]
  plot(NA,xlim=range(cs_),ylim=c(0,1),xlab="contamination c",
       ylab=expression(mean~r[ave]),main="Contamination, M7"); grid()
  ms<-c("DSIIE","SIIE-Q(H)","SIIE-Q(HK)","SIR","SAVE","DR")
  c2<-c(DSIIE="blue",`SIIE-Q(H)`="darkgreen",`SIIE-Q(HK)`="steelblue",
        SIR="black",SAVE="red",DR="purple")
  l2<-c(DSIIE=1,`SIIE-Q(H)`=2,`SIIE-Q(HK)`=2,SIR=3,SAVE=3,DR=3)
  w2<-c(DSIIE=2.2,`SIIE-Q(H)`=1.4,`SIIE-Q(HK)`=1.4,SIR=1,SAVE=1,DR=1)
  for(nm in ms) {
    ys<-sapply(cs_, function(cv) mean(sub$r_ave[sub$method_id==nm & abs(sub$c-cv)<1e-6],na.rm=TRUE))
    lines(cs_,ys,col=c2[nm],lwd=w2[nm],lty=l2[nm]); points(cs_,ys,col=c2[nm],pch=20)
  }
  legend("bottomleft",legend=ms,col=c2[ms],lwd=w2[ms],lty=l2[ms],bty="n",cex=0.7)
  dev.off()
}


########################################################################
#  PART 9 ── 主程式：依序跑全部、產出全部
########################################################################

main <- function(include_paper_models = TRUE, include_csmave_paper = FALSE) {
  t0 <- Sys.time()

  # 跑模擬（有 checkpoint，重跑會跳過已完成的 cell）
  run_for_table3_5_6()   # → Table 3, 5, 6 (M1-M8)
  run_for_table4()       # → Table 4 (M1-M4 at p=50)
  run_for_figure1()      # → Figure 1
  run_for_figure2()      # → Figure 2
  run_for_rank()         # → Table 4 附加
  run_for_table7()       # → Table 7
  run_for_table8()       # → Table 8

  # 主論文 10 個新模型 (ModelI/L15/VII, TAQ3/TAM2/TAM3, C1/C2, ModelII/V)
  # 需要時間較長，可以設 include_paper_models=FALSE 跳過
  if (isTRUE(include_paper_models)) {
    run_paper_models(R = R_PAPER, include_csmave = include_csmave_paper)
  }

  # 讀取所有結果
  recs <- load_records()
  cat(sprintf("\nLoaded %d records from %d design cells\n",
              nrow(recs), length(unique(recs$design_id))))

  # 產生 LaTeX 表格片段 → tex/ 資料夾
  write_table3(recs)
  write_table4(recs)
  write_table5(recs)
  write_table6(recs)
  write_table7()
  write_table8()
  if (isTRUE(include_paper_models)) write_paper_tables()

  # 產生 PDF 圖 → figures/ 資料夾
  write_figure1(recs)
  write_figure2(recs)

  cat(sprintf("\n===== ALL DONE  total=%s =====\n", format(Sys.time()-t0)))
  cat("輸出：\n")
  cat("  tex/tab_primary_[a|b|c|d].tex → 貼入 Table 3 (tab:primary)\n")
  cat("  tex/tab_p50.tex      → 貼入 Table 4 (tab:p50)\n")
  cat("  tex/tab_paired.tex   → 貼入 Table 5 (tab:paired)\n")
  cat("  tex/tab_rmin.tex     → 貼入 Table 6 (tab:rmin)\n")
  cat("  tex/tab_null.tex     → 貼入 Table 7 (tab:null)\n")
  cat("  tex/tab_runtime.tex  → 貼入 Table 8 (tab:runtime)\n")
  cat("  tex/tab_paper_G[1-4].tex → 主論文 4 Group 分表（若有跑 paper models）\n")
  cat("  figures/DSIIE_sample_size_trajectories.pdf → Figure 1\n")
  cat("  figures/DSIIE_sensitivity_panels.pdf       → Figure 2\n")
}

# ============================================================================
# 使用方式：
#   source("DSIIE_simulation_all.R")
#   main()                               # 跑完整版（含主論文 10 個新模型）
#   main(include_paper_models = FALSE)   # 只跑 M0-M8（快，約數小時）
#   main(include_csmave_paper = TRUE)    # 連 csMAVE 在新模型也跑（很慢）
# ============================================================================


########################################################################
#  PART 10 ── 主論文 4 Groups 模型的模擬 driver
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  主論文 (DSIIE_Chen_Wu_.pdf) 的 Group 1-4 全部 10 個模型
#  Group 1: n=200, p=10, H=10, K=2 (Model I, L15, VII)
#  Group 2: n=400, p=10, H=5,  K=2 (TAQ3, TAM2, TAM3)
#  Group 3: n=200, p=10, H=10, K=2 (C1, C2)
#  Group 4: n=200, p=10, H=10, K=2 (Model II, V)
########################################################################

R_PAPER <- 100L    # 主論文的 replications (原文各不同：Table 1 = 50, Table 2 = 200, Table 3 = 100)
# 這裡用統一 100，想跑滿可以個別加大

paper_group_config <- function(name) {
  # 統一用 Standalone Table 3 的規格：n=400, p=10, H=4, K=2, alpha=0.20
  # 這樣新模型可以直接放進 Table 3 的四個分類 (a) (b) (c) (d) 一起顯示
  switch(name,
         ModelI    = list(n=400L, p=10L, H=4L, K=2L, alpha=0.20, group=1L),
         ModelL15  = list(n=400L, p=10L, H=4L, K=2L, alpha=0.20, group=1L),
         ModelVII  = list(n=400L, p=10L, H=4L, K=2L, alpha=0.20, group=1L),
         TAQ3      = list(n=400L, p=10L, H=4L, K=2L, alpha=0.20, group=2L),
         TAM2      = list(n=400L, p=10L, H=4L, K=2L, alpha=0.20, group=2L),
         TAM3      = list(n=400L, p=10L, H=4L, K=2L, alpha=0.20, group=2L),
         C1        = list(n=400L, p=10L, H=4L, K=2L, alpha=0.20, group=3L),
         C2        = list(n=400L, p=10L, H=4L, K=2L, alpha=0.20, group=3L),
         ModelII   = list(n=400L, p=10L, H=4L, K=2L, alpha=0.20, group=4L),
         ModelV    = list(n=400L, p=10L, H=4L, K=2L, alpha=0.20, group=4L)
  )
}

run_paper_models <- function(R = R_PAPER, include_csmave = TRUE) {
  # 跑主論文 10 個模型，涵蓋 Standalone 的三個實驗規格：
  #   SIM-PRIMARY  (n=400, p=10)  → Table 3, Table 5, Table 6
  #   SIM-P50      (n=800, p=50)  → Table 4
  #   SIM-N        (n=200/400/800/1600, p=10) → Figure 1
  models_paper <- c("ModelI","ModelL15","ModelVII",
                    "TAQ3","TAM2","TAM3",
                    "C1","C2",
                    "ModelII","ModelV")
  R_csmave <- 15L
  R_p50 <- max(40L, R %/% 2L)    # P50 少一些 replications 節省時間
  R_traj <- max(40L, R %/% 2L)   # trajectory 也少一些
  
  for (i in seq_along(models_paper)) {
    m <- models_paper[i]
    cfg <- paper_group_config(m)
    
    # ------ SIM-PRIMARY (n=400, p=10) → Table 3 / 5 / 6 ------
    run_cell2(m, "SIM-PRIMARY",
              cfg$n, cfg$p, cfg$H, cfg$K, cfg$alpha, 0,
              R = R, methods = METHODS_FAST,
              seed_base = 100000L + i * 1000L, tag = "fast")
    if (include_csmave) {
      run_cell2(m, "SIM-PRIMARY",
                cfg$n, cfg$p, cfg$H, cfg$K, cfg$alpha, 0,
                R = R_csmave, methods = "csMAVE",
                seed_base = 150000L + i * 1000L, tag = "csmave")
    }
    
    # ------ SIM-P50 (n=800, p=50) → Table 4 ------
    run_cell2(m, "SIM-P50",
              800L, 50L, cfg$H, cfg$K, cfg$alpha, 0,
              R = R_p50, methods = METHODS_FAST,
              seed_base = 200000L + i * 1000L)
    if (include_csmave) {
      run_cell2(m, "SIM-P50",
                800L, 50L, cfg$H, cfg$K, cfg$alpha, 0,
                R = R_csmave, methods = "csMAVE",
                seed_base = 250000L + i * 1000L, tag = "csmave")
    }
    
    # ------ SIM-N trajectories → Figure 1 (n∈{200,400,800,1600}, p=10) ------
    for (j in seq_along(c(200L, 400L, 800L, 1600L))) {
      n_i <- c(200L, 400L, 800L, 1600L)[j]
      run_cell2(m, "SIM-N",
                n_i, 10L, cfg$H, cfg$K, cfg$alpha, 0,
                R = R_traj, methods = METHODS_FAST,
                seed_base = 300000L + i * 10000L + j * 1000L)
    }
  }
}

########################################################################
#  PART 11 ── 產生主論文 4 個 Group 的表格 (Table 1-4)
########################################################################

write_paper_tables <- function() {
  recs <- load_records()
  # 主論文 10 個模型是透過 run_paper_models() 寫入 SIM-PRIMARY 的，
  # 這裡依 group 把它們分別抓出來寫表
  paper_group_models <- list(
    `1` = c("ModelI", "ModelL15", "ModelVII"),
    `2` = c("TAQ3", "TAM2", "TAM3"),
    `3` = c("C1", "C2"),
    `4` = c("ModelII", "ModelV")
  )
  paper_recs <- recs[recs$design_id == "SIM-PRIMARY" &
                       recs$model_id %in% unlist(paper_group_models), ]
  if (!nrow(paper_recs)) { cat("No paper records found (SIM-PRIMARY)\n"); return() }
  methods <- METHODS_ALL

  # 逐個 Group 寫一張表
  for (g in 1:4) {
    models_g <- paper_group_models[[as.character(g)]]
    sub <- paper_recs[paper_recs$model_id %in% models_g, ]
    if (!nrow(sub)) next
    f <- sprintf("tex/tab_paper_G%d.tex", g)
    cat(sprintf("%% Paper Table for Group %d\n", g), file = f)
    for (m in models_g) {
      row_mean <- sapply(methods, function(nm) {
        v <- sub$r_ave[sub$model_id == m & sub$method_id == nm]
        v <- v[is.finite(v)]
        if (length(v) == 0) "--" else fmt(mean(v))
      })
      row_sd <- sapply(methods, function(nm) {
        v <- sub$r_ave[sub$model_id == m & sub$method_id == nm]
        v <- v[is.finite(v)]
        if (length(v) == 0) "--" else sprintf("(%s)", fmt(sd(v)))
      })
      cat(m, paste(row_mean, collapse=" & "), "\\\\\n", sep=" & ",
          file=f, append=TRUE)
      cat("   ", paste(row_sd, collapse=" & "), "\\\\\n", sep=" & ",
          file=f, append=TRUE)
    }
    cat("Wrote", f, "\n")
  }
}