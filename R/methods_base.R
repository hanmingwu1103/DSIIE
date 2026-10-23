# DSIIE research code: base definitions.
# Extracted verbatim (PART 1-6: utilities, DSIIE, SIIE, comparators, metrics, generators)
# from the authors' DSIIE_simulation_all.R (md5 ef880bd18eab2e6c66ffd2838a3851b1); directory creation removed.
# The dr::dr(method = 'dr') branch of fit_DR is not used: DR is fit_DR_LW in methods_corrected.R.

suppressPackageStartupMessages({ library(MASS); library(dr); library(MAVE) })



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
