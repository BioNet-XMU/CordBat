# Title: The CordBat algorithm
# Author: Fanjing Guo
# Date: 2022.12.28

library(igraph)
library(lava)
library(mixOmics)
library(MASS)
library(car)
library(DMwR2)

# ----------------------------------------------
# Community detection
# ----------------------------------------------
ComtyDet <- function(G, InputCOM, minCOMSize){
  # G: correlation coefficient matrix (including all metabolites)
  # InputCOM: a list contains all input communities, each set
  # is a vector (Metabolite Index)
  # minCOMSize: the maximum number of metabolites in a community, the
  # number of metabolites in each community should not be larger than this size
  nextCOM <- list()
  k <- 0
  
  InputCOM.size <- lengths(InputCOM)
  ii <- (InputCOM.size <= minCOMSize)
  if (sum(ii) != 0) {
    k <- sum(ii)
    nextCOM <- InputCOM[ii]
  }
  
  # get communities whose sizes are larger than minCOMSize
  InputCOM <- InputCOM[!ii]
  
  # Community detection
  for (i in c(1: length(InputCOM))) {
    InputG <- abs(G[InputCOM[[i]], InputCOM[[i]]])
    InputG <- InputG - diag(diag(InputG))
    Graph.InputG <- graph_from_adjacency_matrix(InputG, 
                                                weighted = TRUE, 
                                                mode = "undirected")
    COMTY <- cluster_louvain(Graph.InputG)
    cmty.id <- membership(COMTY)
    cmty.size <- as.vector(sizes(COMTY))
    NumComty <- length(cmty.size)
    
    for (n in c(1: NumComty)) {
      Node.Idx <- (cmty.id == n)
      Node.InitID <- InputCOM[[i]][Node.Idx]
      k <- k + 1
      nextCOM[[k]] <- Node.InitID
    }
  }
  
  return(nextCOM)
}


# --------------------------
# Get all communities
# --------------------------
getAllCom <- function(X){
  # X: the data input to get communities 
  # based on the features (column)
  
  N <- nrow(X) # the number of samples
  p <- ncol(X) # the number of features
  
  # get pearson correlation coefficient matrix
  G <- cor(X)
  G <- 1 / 2 * (G + t(G))
  
  initCOM <- list(c(1:p))
  lastCOM <- initCOM
  maxCOMSize <- min(round(N * 0.4), 50)
  cmty.sizemax <- p
  
  Nochange.times <- 0
  while(cmty.sizemax > maxCOMSize){
    com.num.old <- length(lastCOM)
    nextCOM <- ComtyDet(G, lastCOM, maxCOMSize)
    lastCOM <- nextCOM
    com.num.new <- length(lastCOM)
    ComtySize <- lengths(lastCOM)
    cmty.sizemax <- max(ComtySize)
    
    if (com.num.new == com.num.old) {
      Nochange.times <- Nochange.times + 1
    }
    else {
      Nochange.times <- 0
    }
    
    if (Nochange.times == 3) {
      break
    }
  }
  
  # a community including only one metabolite will be merged with another
  NumComty <- length(lastCOM)
  ComtySize <- lengths(lastCOM)
  # get corresponding community ID for each feature
  ComtyID <- rep(0, p)
  for (i in c(1: NumComty)) {
    cmtyi.metID <- lastCOM[[i]]
    ComtyID[cmtyi.metID] <- i
  }
  
  cmty.sizeis1 <- which(ComtySize == 1)
  Num.cmty.sizeis1 <- length(cmty.sizeis1)
  if (Num.cmty.sizeis1 != 0) {
    for (i in c(1: Num.cmty.sizeis1)){
      if (length(lastCOM[[cmty.sizeis1[i]]]) == 1) {
        S1.metID <- lastCOM[[cmty.sizeis1[i]]]
        allCorwithS1 <- G[, S1.metID]
        allCorwithS1[S1.metID] <- 0 # ignore correlation with itself
        Maxcor.metID <- which(allCorwithS1 == max(allCorwithS1))
        mergeCmty <- ComtyID[Maxcor.metID]
        mC.metIDs <- lastCOM[[mergeCmty]]
        mC.metIDs <- c(mC.metIDs, S1.metID)
        lastCOM[[mergeCmty]] <- mC.metIDs
        ComtyID[S1.metID] <- mergeCmty
      }
    }
    ComtySize <- lengths(lastCOM)
    cmty.sizeis1 <- which(ComtySize == 1)
    # delete
    lastCOM <- lastCOM[-cmty.sizeis1]
  }
  
  ComtySize <- lengths(lastCOM)
  cmty.sizeis2 <- which(ComtySize == 2)
  Num.cmty.sizeis2 <- length(cmty.sizeis2)
  if (Num.cmty.sizeis2 != 0) {
    for (i in c(1: Num.cmty.sizeis2)){
      if (length(lastCOM[[cmty.sizeis2[i]]]) == 2) {
        S2.metID <- lastCOM[[cmty.sizeis2[i]]]
        allCorwithS2 <- G[, S2.metID]
        allCorwithS2[S2.metID, ] <- 0 # ignore correlation with itself
        Maxcor.metID <- which(allCorwithS2 == max(allCorwithS2), 
                              arr.ind = T)[1]
        mergeCmty <- ComtyID[Maxcor.metID]
        mC.metIDs <- lastCOM[[mergeCmty]]
        mC.metIDs <- c(mC.metIDs, S2.metID)
        lastCOM[[mergeCmty]] <- mC.metIDs
      }
    }
    ComtySize <- lengths(lastCOM)
    cmty.sizeis2 <- which(ComtySize == 2)
    # delete
    lastCOM <- lastCOM[-cmty.sizeis2]
  }
  
  
  return(lastCOM)
}


# ----------------------------------------
# select a proper fold number for CV
# --------------------------------------
selfoldforCV <- function(N){
  foldstosel <- c(2:9)
  Num.quo <- N %/% foldstosel
  Num.rem <- N %% foldstosel
  
  # if no fold meets requirements, then change N
  selIdx <- (Num.rem == 0 & Num.quo >= 10)
  folds <- foldstosel[selIdx]
  folds.num <- length(folds)
  
  if (folds.num != 0) {
    d5 <- abs(rep(5, folds.num) - folds)
    fold <- folds[d5 == min(d5)]
    if (length(fold) > 1) {
      if (N > 300) {
        fold <- max(fold)
      }
      else {
        fold <- min(fold)
      }
    }
  }
  else {
    fold <- 1
  }
  
  return(fold)
}

# ------------------
# soft threshold
# ------------------
soft <- function(x, lambda){
  s <- sign(x) * max(abs(x) - lambda, 0)
  return(s)
}

# -------------------
# coordinate descent 
# -------------------
CDfgL <- function(V, beta_i, u, rho){
  p_1 <- ncol(V)
  
  # initialize
  beta.new <- rep(0, p_1)
  finished <- rep(FALSE, p_1)
  eps <- 1e-4
  times0 <- 0
  times1 <- 0
  
  while(TRUE){
    beta.old <- beta_i
    for (j in c(1: p_1)) {
      df <- V %*% beta_i - u
      x <- beta_i[j] - df[j] / V[j, j]
      beta_i[j] <- soft(x, rho / V[j, j])
    }
    beta.new <- beta_i
    
    zeroIdx <- (beta_i == 0)
    if (any(zeroIdx)) {
      if (all(beta.new[zeroIdx] == beta.old[zeroIdx])) {
        times0 <- times0 + 1
      }
      else{
        times0 <- 0
      }
      if (times0 >= 5) {
        finished[zeroIdx] <- TRUE
      }
    }
    
    if (any(!zeroIdx)) {
      if (all(abs((beta.new[!zeroIdx] - beta.old[!zeroIdx])
                  / beta.old[!zeroIdx]) < eps)) {
        times1 <- times1 + 1
      }
      else{
        times1 <- 0
      }
      if (times1 >= 3) {
        finished[!zeroIdx] <- TRUE
      }
    }
    
    
    if (all(finished)) {
      
      break
    }
  }
  return(beta.new)
}

# -------------------------------
# graphical lasso algorithm
# -------------------------------
graphicalLasso <- function(X, rho){
  N <- nrow(X)
  p <- ncol(X)
  # centered and scaling
  X <- scale(X, center = TRUE, scale = TRUE)
  # get covariance matrix
  S <- cov(X)
  
  # initialize
  Theta <- matrix(0, p, p)
  W <- S + rho * diag(1, p, p)
  B <- matrix(0, p-1, p)
  t <- 1e-5
  
  while (TRUE) {
    W.old <- W
    for (i in c(1: p)) {
      idx <- c(1: p)
      idx <- idx[-i]
      
      W_11 <- W[idx, idx]
      V <- W_11
      
      s_12 <- S[idx, i]
      u <- s_12
      
      B[, i] <- CDfgL(V, B[, i], u, rho)
      w_12 <- V %*% B[, i]
      W[idx, i] <- w_12
      W[i, idx] <- W[idx, i]
    }
    W.new <- W
    dW <- W.new - W.old
    S_ndiag <- c(S[upper.tri(S, diag = F)], 
                 S[lower.tri(S, diag = F)])
    if (mean(abs(as.vector(dW))) / mean(abs(S_ndiag)) < t) {
      break
    } 
  }
  
  for (i in c(1: p)) {
    idx <- c(1: p)
    idx <- idx[-i]
    Theta[i, i] <- 1 / (W[i, i] - t(W[idx, i]) %*% B[, i])
    Theta[idx, i] <- - B[, i] * Theta[i, i]
    Theta[i, idx] <-  Theta[idx, i]
  }
  
  c.mat <- list(Theta = Theta, W = W)
  return(c.mat)
}

# -------------------------------------------------------------
# Stability approach to regularization selection(StARS) 
# for High dimensional graphical models
# -------------------------------------------------------------
StARS <- function(X, b, M, print.detail = T){
  N <- nrow(X)
  p <- ncol(X)
  beta <- 0.05
  D_var <- 0
  Sel.rho <- 0.5
  
  # determine direction
  rho <- 0.5
  sum.fai <- matrix(0, p, p)
  for (i in c(1: M)) {
    subsampIdx <- sample(c(1: N), b, replace = F)
    Si <- X[subsampIdx, ]
    c.mat <- graphicalLasso(Si, rho)
    Theta <- c.mat$Theta
    # the element larger than 1e-6 means an edge exists
    Theta <- Theta * (abs(Theta) > 1e-6)
    fai <- 1 * (Theta != 0)
    sum.fai <- sum.fai + fai
  }
  thetb <- 1 / M * sum.fai
  ksib <- 2 * thetb * (matrix(1, p, p) - thetb)
  Db.5 <- sum(ksib[upper.tri(ksib, diag = FALSE)]) / (p * (p - 1) / 2)
  
  if (Db.5 > beta) {
    rhos <- seq(from = 0.9, to = 0.5, by = -0.1)
  }
  else {
    rho <- 0.1
    sum.fai <- matrix(0, p, p)
    for (i in c(1: M)) {
      subsampIdx <- sample(c(1: N), b)
      Si <- X[subsampIdx, ]
      c.mat <- graphicalLasso(Si, rho)
      Theta <- c.mat$Theta
      # the element larger than 1e-6 means an edge exists
      Theta <- Theta * (abs(Theta) > 1e-6)
      fai <- 1 * (Theta != 0)
      sum.fai <- sum.fai + fai
    }
    thetb <- 1 / M * sum.fai
    ksib <- 2 * thetb * (matrix(1, p, p) - thetb)
    Db.1 <- sum(ksib[upper.tri(ksib, diag = FALSE)]) / (p * (p - 1) / 2)
    
    if (Db.1 > beta) {
      rhos <- seq(from = 0.5, to = 0.1, by = -0.1)
    }
    else  {
      rho <- 0.05
      sum.fai <- matrix(0, p, p)
      for (i in c(1: M)) {
        subsampIdx <- sample(c(1: N), b)
        Si <- X[subsampIdx, ]
        c.mat <- graphicalLasso(Si, rho)
        Theta <- c.mat$Theta
        # the element larger than 1e-6 means an edge exists
        Theta <- Theta * (abs(Theta) > 1e-6)
        fai <- 1 * (Theta != 0)
        sum.fai <- sum.fai + fai
      }
      thetb <- 1 / M * sum.fai
      ksib <- 2 * thetb * (matrix(1, p, p) - thetb)
      Db.05 <- sum(ksib[upper.tri(ksib, diag = FALSE)]) / (p * (p - 1) / 2)
      
      if (Db.05 > beta) {
        rhos <- seq(from = 0.1, to = 0.05, by = -0.01)
      }
      else {
        rhos <- seq(from = 0.05, to = 0.01, by = -0.01)
      }
    }
  }
  
  # first selection
  for (r in c(1: length(rhos))) {
    rho <- rhos[r]
    sum.fai <- matrix(0, p, p)
    for (i in c(1: M)) {
      subsampIdx <- sample(c(1: N), b)
      Si <- X[subsampIdx, ]
      c.mat <- graphicalLasso(Si, rho)
      Theta <- c.mat$Theta
      # the element larger than 1e-6 means an edge exists
      Theta <- Theta * (abs(Theta) > 1e-6)
      fai <- 1 * (Theta != 0)
      sum.fai <- sum.fai + fai
    }
    
    thetb <- 1 / M * sum.fai
    ksib <- 2 * thetb * (matrix(1, p, p) - thetb)
    Db <- sum(ksib[upper.tri(ksib, diag = FALSE)]) / (p * (p - 1) / 2)
    
    if (Db <= beta) {
      if (Db > D_var) {
        D_var <- Db
      }
    }
    else {
      if (r - 1 != 0) {
        Sel.rho <- rhos[r - 1]
      }
      else {
        if (rhos[r] * 100 >= 10) {
          Sel.rho <- rhos[r] + 0.1
        }
        else {
          Sel.rho <- rhos[r] + 0.01
        }
      }
      
      break
    }
  }

  if (print.detail) {
    cat('Select rho =', Sel.rho, '\n')
  }
  return(c(Sel.rho, D_var))
}


# --------------------------------------------
# CV + BIC select rho
# --------------------------------------------
selrho.useCVBIC <- function(X, print.detail = T) {
  N <- nrow(X)
  fold <- selfoldforCV(N)
  CVset.size <- N / fold
  
  rhos <- seq(from = 0.1, to = 0.9, by = 0.1)
  CVerr1 <- matrix(0, fold, length(rhos))
  CVerr2 <- rep(0, length(rhos))
  
  for (r in c(1: length(rhos))) {
    rho <- rhos[r]
    if (fold != 1) {
      for (i in c(1: fold)) {
        start.index <- (i-1) * CVset.size + 1
        end.index <- i * CVset.size
        X.cv <- X[c(start.index: end.index), ]
        X.tr <- X
        X.tr <- X.tr[-c(start.index: end.index), ]
        
        c.mat <- graphicalLasso(X.tr, rho)
        Theta <- c.mat$Theta
        
        # compute error for CV set
        X.cv.sca <- scale(X.cv, center = TRUE, scale = TRUE)
        S.cv <- cov(X.cv.sca)
        
        k <- sum(Theta[upper.tri(Theta, diag = FALSE)] != 0)
        CVerr1[i, r] <- k * log(CVset.size) - CVset.size * (log(det(Theta)) - tr(S.cv %*% Theta))
      }
    }
    else {
      c.mat <- graphicalLasso(X, rho)
      Theta <- c.mat$Theta
      
      # compute error for CV set
      X.sca <- scale(X, center = TRUE, scale = TRUE)
      S <- cov(X.sca)
      
      k <- sum(Theta[upper.tri(Theta, diag = FALSE)] != 0)
      CVerr2[r] <- k * log(N) - 2 * (log(det(Theta)) - tr(S %*% Theta))
    }
    
  }
  
  if (fold != 1) {
    CVerr1 <- colMeans(CVerr1)
    MinCVerr <- min(CVerr1)
    rho.cv <- rhos[CVerr1 == MinCVerr]
  }
  else {
    MinCVerr <- min(CVerr2)
    rho.cv <- rhos[CVerr2 == MinCVerr]
  }
  
  
  
  if (print.detail) {
    cat('CVBIC: select rho =', rho.cv, '\n')
  }
  
  return(c(rho.cv, MinCVerr))
}

# -----------------------------------------
# update correction coefficients a and b
# -----------------------------------------
update.CorrectCoef <- function(X0.glist, X1.glist, Theta.list, 
                               a.i, b.i, penal.ksi, penal.gamma) {
  G <- length(X0.glist)
  p <- ncol(X0.glist[[1]])
  
  N0_gvec <- rep(0, G)
  N1_gvec <- rep(0, G)
  N_gvec <- rep(0, G)
  
  for (g in c(1: G)) {
    N0_gvec[g] <- nrow(X0.glist[[g]])
    N1_gvec[g] <- nrow(X1.glist[[g]])
    N_gvec[g] <- N0_gvec[g] + N1_gvec[g]
  }
  
  a.o <- a.i
  b.o <- b.i
  
  for (j in c(1: p)) {
    
    # update a
    tmp1.gvec <- rep(0, G)
    tmp2.gvec <- rep(0, G)
    for (g in c(1: G)) {
      A <- diag(a.o)
      B_gi <- rep(1, N1_gvec[g]) %*% t(b.o)
      
      X1.gi.cor <- X1.glist[[g]] %*% A + B_gi
      X.gi <- rbind(X0.glist[[g]], X1.gi.cor)
      
      X.gi.sca <- scale(X.gi, center = TRUE, scale = TRUE)
      X.gi.sca.attr <- attributes(X.gi.sca)
      Mu_g <- X.gi.sca.attr$`scaled:center`
      Sigma_g <- X.gi.sca.attr$`scaled:scale`
      
      Y_g <- X.gi.sca[(N0_gvec[g] + 1): N_gvec[g], ]
      
      Y_g[, j] <- rep(1, N1_gvec[g]) * b.o[j]
      Y_g[, j] <- Y_g[, j] - rep(1, N1_gvec[g]) * Mu_g[j]
      Y_g[, j] <- Y_g[, j] / (rep(1, N1_gvec[g]) * Sigma_g[j])
      Z_g <- Y_g
      
      tmp <- Z_g %*% Theta.list[[g]][, j]
      tmp1_g <- 2 / (N_gvec[g] * Sigma_g[j]) * sum(X1.glist[[g]][, j] * tmp)
      tmp1.gvec[g] <- tmp1_g
      tmp2_g <- 2 * Theta.list[[g]][j, j] / (N_gvec[g] * Sigma_g[j]^2) * sum(X1.glist[[g]][, j]^2)
      tmp2.gvec[g] <- tmp2_g
    }
    a.o[j] <- 1 + soft(- sum(tmp1.gvec) - sum(tmp2.gvec), penal.ksi) / sum(tmp2.gvec)
    
    # update b
    tmp3.gvec <- rep(0, G)
    tmp4.gvec <- rep(0, G)
    for (g in c(1: G)) {
      A <- diag(a.o)
      B_gi <- rep(1, N1_gvec[g]) %*% t(b.o)
      X1.gi.cor <- X1.glist[[g]] %*% A + B_gi
      X.gi <- rbind(X0.glist[[g]], X1.gi.cor)
      
      X.gi.sca <- scale(X.gi, center = TRUE, scale = TRUE)
      X.gi.sca.attr <- attributes(X.gi.sca)
      Mu_g <- X.gi.sca.attr$`scaled:center`
      Sigma_g <- X.gi.sca.attr$`scaled:scale`
      Y_g <- X.gi.sca[(N0_gvec[g] + 1): N_gvec[g], ]
      
      Y_g[, j] <- X1.glist[[g]][, j] * a.o[j]
      Y_g[, j] <- Y_g[, j] - rep(1, N1_gvec[g]) * Mu_g[j]
      Y_g[, j] <- Y_g[, j] / (rep(1, N1_gvec[g]) * Sigma_g[j])
      Z_g <- Y_g
      
      tmp <- Z_g %*% Theta.list[[g]][, j]
      tmp3_g <- 2 / (N_gvec[g] * Sigma_g[j]) * sum(tmp)
      tmp3.gvec[g] <- tmp3_g
      tmp4_g <- 2 * N1_gvec[g] * Theta.list[[g]][j, j] / (N_gvec[g] * Sigma_g[j]^2)
      tmp4.gvec[g] <- tmp4_g
    }
    b.o[j] <- soft(- sum(tmp3.gvec), penal.gamma) / sum(tmp4.gvec)
    
  }
  
  coef.out <- list(coef.a = a.o, 
                   coef.b = b.o)
  return(coef.out)
}



# ------------------------------------
# use GGM for batch effect correction
# ------------------------------------
BEgLasso <- function(X0.glist, X1.glist, penal.rho, penal.ksi,
                     penal.gamma, eps) {
  G <- length(X0.glist)
  p <- ncol(X0.glist[[1]])
  
  N0_gvec <- rep(0, G)
  N1_gvec <- rep(0, G)
  
  X0.m.gmat <- matrix(0, p, G)
  X1.m.gmat <- matrix(0, p, G)
  for (g in c(1: G)) {
    N0_gi <- nrow(X0.glist[[g]])
    N0_gvec[g] <- N0_gi
    
    X0.m.gmat[, g] <- colMeans(X0.glist[[g]])
    
    N1_gi <- nrow(X1.glist[[g]])
    N1_gvec[g] <- N1_gi
    
    X1.m.gmat[, g] <- colMeans(X1.glist[[g]])
  }
  
  # Initialize
  Theta.list <- list()
  B.list <- list()
  coef.as <- matrix(0, p, G)
  for (g in c(1: G)) {
    Theta.list[[g]] <- matrix(0, p, p)
    B.list[[g]] <- matrix(0, p-1, p)
    coef.as[, g] <- X0.m.gmat[, g] / X1.m.gmat[, g]
  }
  coef.a <- rowMeans(coef.as)
  coef.b <- rep(0, p)
  
  coef.A <- diag(coef.a)
  X1.cor.glist <- list()
  for (g in c(1: G)) {
    coef.B.gi <- rep(1, N1_gvec[g]) %*% t(coef.b)
    X1.gi.cor <- X1.glist[[g]] %*% coef.A + coef.B.gi
    X1.cor.glist[[g]] <- X1.gi.cor
  }
  
  X.glist <- list()
  W.list <- list()
  for (g in c(1: G)) {
    X.gi <- rbind(X0.glist[[g]], X1.cor.glist[[g]])
    X.gi.sca <- scale(X.gi, center = TRUE, scale = TRUE)
    S0_gi <- cov(X.gi.sca)
    W_i <- S0_gi + penal.rho * diag(1, p)
    W.list[[g]] <- W_i
  }

  t <- eps
  
  # iterations
  times0.gvec <- rep(0, G)
  times1.gvec <- rep(0, G)
  finished.gmat <- matrix(F, (p - 1) * p / 2, G)
  
  while(TRUE) {
    W_old.list <- W.list
    coef.A <- diag(coef.a)
    S.list <- list()
    for (g in c(1: G)) {
      coef.B.gi <- rep(1, N1_gvec[g]) %*% t(coef.b)
      X1.gi.cor <- X1.glist[[g]] %*% coef.A + coef.B.gi
      
      X.gi <- rbind(X0.glist[[g]], X1.gi.cor)
      X.gi.sca <- scale(X.gi, center = TRUE, scale = TRUE)
      S_i <- cov(X.gi.sca)
      S.list[[g]] <- S_i
    }
    
    for (j in c(1: p)) {
      idx <- c(1: p)
      idx <- idx[-j]
      
      for (g in c(1: G)) {
        Wi_11 <- W.list[[g]][idx, idx]
        V_i <- Wi_11
        
        si_12 <- S.list[[g]][idx, j]
        u_i <- si_12
        
        B.list[[g]][, j] <- CDfgL(V_i, B.list[[g]][, j], u_i, penal.rho)
        wi_12 <- V_i %*% B.list[[g]][, j]
        W.list[[g]][idx, j] <- wi_12
        W.list[[g]][j, idx] <- W.list[[g]][idx, j]
      }
    }
    
    for (j in c(1: p)) {
      idx <- c(1: p)
      idx <- idx[-j]
      
      for (g in c(1: G)) {
        Theta.list[[g]][j, j] <- (W.list[[g]][j, j] - t(W.list[[g]][idx, j]) %*% B.list[[g]][, j])
        Theta.list[[g]][idx, j] <- - B.list[[g]][, j] * Theta.list[[g]][j, j]
        Theta.list[[g]][j, idx] <-  Theta.list[[g]][idx, j]
      }
    }
    
    W_new.list <- W.list
    W_new.ndiag.list <- list()
    W_old.ndiag.list <- list()
    zeroIdx.list <- list()
    for (g in c(1: G)) {
      W_new.ndiag.list[[g]] <- W_new.list[[g]][upper.tri(W_new.list[[g]], diag = F)]
      W_old.ndiag.list[[g]] <- W_old.list[[g]][upper.tri(W_new.list[[g]], diag = F)]
      zeroIdx_g <- (W_old.ndiag.list[[g]] < 1e-5)
      zeroIdx.list[[g]] <- zeroIdx_g
    }
    
    for (g in c(1: G)) {
      zeroIdx_g <- zeroIdx.list[[g]]
      if (any(zeroIdx_g)) {
        if (all(W_new.ndiag.list[[g]][zeroIdx_g] < 1e-5)) {
          times0.gvec[[g]] <- times0.gvec[[g]] + 1
        }
        else {
          times0.gvec[[g]] <- 0
        }
        if (times0.gvec[[g]] >= 5) {
          finished.gmat[zeroIdx_g, g] <- TRUE
        }
      }
      
      if (any(!zeroIdx_g)) {
        dW_g <- max(abs((W_new.ndiag.list[[g]][!zeroIdx_g] - W_old.ndiag.list[[g]][!zeroIdx_g]) 
                        / W_old.ndiag.list[[g]][!zeroIdx_g])) 
        
        if (dW_g < t) {
          times1.gvec[[g]] <- times1.gvec[[g]] + 1
        }
        else {
          times1.gvec[[g]] <- 0
        }
        
        if (times1.gvec[[g]] >= 3) {
          finished.gmat[!zeroIdx_g, g] <- TRUE
        }
      }
    }
    
    if (all(as.vector(finished.gmat))) {
      break
    }
    else {
      coef.update <- update.CorrectCoef(X0.glist, X1.glist, Theta.list, 
                                        coef.a, coef.b, penal.ksi, penal.gamma)
      coef.a <- coef.update$coef.a
      coef.b <- coef.update$coef.b
    }
  }
  
  # replace a negative coefficient a with a small positive value
  neg.Idx <- (coef.a < 0)
  coef.a[neg.Idx] <- 0.05
  
  coef.A <- diag(coef.a)
  
  for (g in c(1: G)) {
    coef.B.gi <- rep(1, N1_gvec[g]) %*% t(coef.b)
    X1.gi.cor <- X1.glist[[g]] %*% coef.A + coef.B.gi
    X1.cor.glist[[g]] <- X1.gi.cor
  }
  
  para.out <- list(Theta = Theta.list,  
                   X1.cor = X1.cor.glist,
                   coef.a = coef.a, 
                   coef.b = coef.b)
  
  return(para.out)
}



# ----------------------------------------------------------------
# find best parameters (ksi and gamma) using Mahalanobis distance
# ----------------------------------------------------------------
findBestPara <- function(X0.glist, X1.glist, penal.rho, eps) {
  G <- length(X0.glist)
  p <- ncol(X0.glist[[1]])
  
  N0_gvec <- rep(0, G)
  N1_gvec <- rep(0, G)
  N_gvec <- rep(0, G)
  for (g in c(1: G)) {
    N0_gvec[g] <- nrow(X0.glist[[g]])
    N1_gvec[g] <- nrow(X1.glist[[g]])
    N_gvec[g] <- N0_gvec[g] + N1_gvec[g]
  }
  
  Sel.ksi <- 0
  Sel.gamma <- 0
  
  ksis <- c(1, 0.5, 0.3, 0.1)
  gammas <- c(1, 0.5, 0.3, 0.1)
  
  MinAvedist <- Inf
  
  for (k in c(1: length(ksis))) {
    for (g in c(1: length(gammas))) {
      ksi.k <- ksis[k]
      gamma.g <- gammas[g]
      
      Allpara <- BEgLasso(X0.glist, X1.glist, penal.rho, ksi.k, gamma.g, eps)
      X1.cor.glist <- Allpara$X1.cor
      Theta.list <- Allpara$Theta
      
      r <- 0.5
      
      ebic.gvec <- rep(0, G)
      for (i in c(1: G)) {
        X.gi <- rbind(X0.glist[[i]], X1.cor.glist[[i]])
        
        # get empirical covariance matrix
        X.gi.sca <- scale(X.gi, center = T, scale = T)
        S_i <- cov(X.gi.sca)
        Theta_i <- Theta.list[[i]]
        E.num.gi <- sum(Theta_i[upper.tri(Theta_i, diag = FALSE)] != 0)
        
        # EBIC
        ebic.gvec[i] <- - N_gvec[i] * (log(det(Theta_i)) - tr(S_i %*% Theta_i)) + 
          E.num.gi * log(N_gvec[i]) + 4 * E.num.gi * r * log(p)
      }
      
      ebic <- sum(ebic.gvec)
      
      totdist <- ebic
      
      
      if (totdist < MinAvedist) {
        MinAvedist <- totdist
        Sel.ksi <- ksi.k
        Sel.gamma <- gamma.g
        
      }
    }
  }
  
  
  penterm <- list(penal.ksi = Sel.ksi,
                  penal.gamma = Sel.gamma,
                  MinAvedist = MinAvedist)
  
  return(penterm)
}

# -------------------------
# Delete outliers in data 
# --------------------------
DelOutlier <- function(X) {
  pca.dat <- pca(X, ncomp = 3, center = TRUE, scale = TRUE)
  pca.dat.varX <- pca.dat$variates$X
  delsampIdx <- c()
  for (i in c(1: 3)) {
    pc.i <- pca.dat.varX[, i]
    pc.i.m <- mean(pc.i)
    pc.i.sd <- sd(pc.i)
    pc.i.min <- pc.i.m - 3 * pc.i.sd
    pc.i.max <- pc.i.m + 3 * pc.i.sd
    delsampIdx <- c(delsampIdx, which(pc.i < pc.i.min | pc.i > pc.i.max))
    
  }
  delsampIdx <- unique(delsampIdx)
  
  if (length(delsampIdx) == 0) {
    X.out <- X
  }
  else {
    X.out <- X[-delsampIdx, ]
  }
  
  Del.result <- list(delsampIdx = delsampIdx, 
                     X.out = X.out)
  return(Del.result)
}

# --------------------------
# Impute outliers in data
# -------------------------
ImputeOutlier <- function(X) {
  p <- ncol(X)
  X.out <- X
  for (i in c(1: p)) {
    dat.i <- X[, p]
    dat.i.m <- mean(dat.i)
    dat.i.sd <- sd(dat.i)
    dat.i.max <- dat.i.m + 3 * dat.i.sd
    dat.i.min <- dat.i.m - 3 * dat.i.sd
    dat.i[dat.i < dat.i.min | dat.i > dat.i.max] <- NA
    
    X.out[, p] <- dat.i
  }
  
  na.num <- sum(is.na(X.out))
  if (na.num != 0) {
    X.out <- as.data.frame(X.out)
    X.out <- knnImputation(X.out)
  }
  
  X.out <- as.matrix(X.out)
  
  return(X.out)
}

# -------------------
# get reference batch
# ---------------------
getRefbat <- function(Data, batch) {
  # cumulative RSD
  batch <- as.factor(batch)
  batch.num <- nlevels(batch)
  batch.lev <- levels(batch)
  
  rsd <- seq(0.01, 1, 0.01)
  rsd.cumfreq <- data.frame(batch = factor(rep(batch.lev, each = length(rsd)), 
                                           levels = c(1: batch.num)), 
                            rsd = rep(rsd, batch.num), 
                            cumfreq = rep(0, length(rsd) * batch.num))
  
  p <- ncol(Data)
  for (i in c(1: batch.num)) {
    dat.i <- Data[batch == batch.lev[i], ]
    m.val <- apply(dat.i, MARGIN = 2, FUN = mean)
    sd.val <- apply(dat.i, MARGIN = 2, FUN = sd)
    rsd.val <- sd.val / m.val
    cumfreq.i <- rep(0, length(rsd))
    
    for (k in c(1: length(rsd))) {
      cumfreq.i[k] <- sum(rsd.val <= rsd[k]) / p
    }
    rsd.cumfreq$cumfreq[rsd.cumfreq$batch == batch.lev[i]] <- cumfreq.i
  }
  
  # compute the area under the curve
  AUCFC_rsd <- rep(0, batch.num) # area under the cumulative frequency curve of RSD for each batch
  delta_x <- 0.01
  for (i in c(1: batch.num)) {
    cumfreq.bati <- rsd.cumfreq[rsd.cumfreq$batch == batch.lev[i], ]$cumfreq
    n <- length(cumfreq.bati)
    
    y.1 <- cumfreq.bati[1]
    y.n <- cumfreq.bati[n]
    y.is <- cumfreq.bati[-c(1, n)]
    
    AUCFC_rsd[i] <- delta_x * (0.5 * (y.1 + y.n) + sum(y.is))
  }
  
  batch.rank <- as.numeric(batch.lev[order(AUCFC_rsd, decreasing = T)])
  AUCFC_rsd.sort <- sort(AUCFC_rsd, decreasing = T)
  
  # plot
  linesize <- 1
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", 
                 "#F0E442", "#0072B2", "#D55E00", "#999999", 
                 "#33CC99", "#990000", "#CC99FF", "#FF00FF", 
                 '#5500FF', '#66FF66', '#0088A8', '#A500CC')
  
  lab.txt <- paste0('Area under the curve: \n(decreasing order)\n')
  for (i in c(1: batch.num)) {
    lab.txt <- paste0(lab.txt, 
                      'Batch ', batch.rank[i], ': ', round(AUCFC_rsd.sort[i], 4), '\n')
  }
  
  ggplot(data = rsd.cumfreq, aes(x = rsd, y = cumfreq)) +
    geom_line(aes(color = batch), linewidth = linesize) +
    geom_vline(xintercept = 0.2, linetype = 2) +
    geom_vline(xintercept = 0.3, linetype = 2) + 
    geom_text(aes(x = 0.85, y = 0.02 * (batch.num + 2)), 
              label = lab.txt, color = "black", size = 5) + 
    xlab('RSD') + ylab('Cumulative Frequency') + xlim(0, 1) + ylim(0, 1) +
    scale_color_manual(values = cbPalette) + scale_linetype_discrete() + theme_bw() +
    theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12, face = 'bold'),
          legend.title = element_text(size = 15), legend.text = element_text(size = 12),
          plot.title = element_text(size = 18, hjust = 0.5)) +
    labs(color = 'Batch', linetype = NA, title = 'RSD')
}

# ----------------------------------------
# Batch effect correction (BEC) - CordBat
# ----------------------------------------
CordBat <- function(X, batch, group = NULL, grouping = F, ref.batch, 
                    eps = 1e-5, print.detail = T) {
  # X is the log-transformed data matrix of all batches
  # batch is the batch IDs for samples, a vector
  # group is the group IDs for samples, a vector
  # ref.batch is the ID of the reference batch
  if (is.null(group)) {
    containQC <- NA
    group <- rep(1, nrow(X))
  }
  else {
    if (sum(group == 'QC') != 0) {
      containQC <- T
      
      X.init <- X
      batch.init <- batch
      group.init <- group
      
      X_QC <- X[group == 'QC', ]
      X <- X[group != 'QC', ]
      
      QC_batch <- batch[group == 'QC']
      batch <- batch[group != 'QC']
      group <- group[group != 'QC']
    }
    else {
      containQC <- F
    }
    
    if (!grouping) {
      group <- rep(1, nrow(X))
    }
  }
  
  X <- as.matrix(X)
  p <- ncol(X)
  
  batch.f <- as.factor(batch)
  batch.num <- nlevels(batch.f)
  batch.level <- levels(batch.f)
  
  # delete outliers
  batch.new <- batch
  group.new <- group
  X.delout <- X
  delsampIdx <- c()
  for (i in c(1: batch.num)) {
    X.bati <- X[batch == batch.level[i], ]
    bati.idx <- which(batch == batch.level[i])
    X.bati <- DelOutlier(X.bati)
    delsamp.bati <- X.bati$delsampIdx
    dat.bati <- X.bati$X.out
    dat.bati <- ImputeOutlier(dat.bati)
    
    if (length(delsamp.bati) != 0) {
      bati.delinitIdx <- bati.idx[delsamp.bati]
      X.delout[bati.idx[-delsamp.bati], ] <- dat.bati
      delsampIdx <- c(delsampIdx, bati.delinitIdx)
    }
    else {
      X.delout[bati.idx, ] <- dat.bati
    }
  }
  
  if (length(delsampIdx) != 0) {
    X.nodel <- X.delout
    X.delout <- X.delout[-delsampIdx, ]
    batch.new <- batch.new[-delsampIdx]
    group.new <- group.new[-delsampIdx]
  }
  else {
    X.nodel <- X
  }
  
  # transform to factor
  batch.new <- as.factor(batch.new)
  batch.num <- nlevels(batch.new)
  batch.level <- levels(batch.new)
  
  group.new <- as.factor(group.new)
  group.num <- nlevels(group.new)
  grp.level <- levels(group.new)
  
  # initialize output
  Theta.list <- list()
  for (i in c(1: group.num)) {
    Theta.list[[i]] <- matrix(0, p, p)
  }
  a <- rep(1, p)
  b <- rep(0, p)
  
  para <- list(Theta = Theta.list,  
               coef.a = a,
               coef.b = b)

  Xcor.para <- list()
  for (i in c(1: batch.num)) {
    Xcor.para[[i]] <- para
  }
  
  # all data
  X.cor <- matrix(0, nrow(X), ncol(X))
  X.cor[batch == ref.batch, ] <- X[batch == ref.batch, ]
  
  # data delete outlier
  X.cor.1 <- matrix(0, nrow(X.delout), ncol(X.delout))
  X.cor.1[batch.new == ref.batch, ] <- X.delout[batch.new == ref.batch, ]
  
  # if data contain QCs
  if ((!is.na(containQC)) & containQC) {
    X.cor.withQC <- matrix(0, nrow(X.init), ncol(X.init))
    X.cor.withQC[batch.init == ref.batch, ] <- X.init[batch.init == ref.batch, ]
  }
  
  #### community detection ------------------------
  Xb0.mat <- X.delout[batch.new == ref.batch, ]
  
  COM <- getAllCom(Xb0.mat)
  cat("Community detection: ", length(COM), "communities", "\n",
      "Size: ", lengths(COM), "\n")
  
  # batch effect correction for communities
  for (i in c(1: length(COM))) {
    metID <- COM[[i]]
    Xb0.COMi <- Xb0.mat[, metID]
    Nb0.COMi <- nrow(Xb0.COMi)
    
    Xb0.COMi.glist <- list()
    Nb0.COMi.gvec <- rep(0, group.num)
    for (g in c(1: group.num)) {
      Xb0.gi.COMi <- X.delout[batch.new == ref.batch & group.new == grp.level[g], metID]
      Xb0.COMi.glist[[g]] <- Xb0.gi.COMi
      Nb0.gi.COMi <- nrow(Xb0.gi.COMi)
      Nb0.COMi.gvec[g] <- Nb0.gi.COMi
    }
    
    if (length(metID) > 5) {
      rhos <- rep(0, group.num)
      for (g in c(1: group.num)) {
        rho_g <- StARS(Xb0.COMi.glist[[g]], round(0.7 * Nb0.COMi.gvec[g]), 100, 
                       print.detail = print.detail)
        rho_g <- rho_g[1]
        rhos[g] <- rho_g
      }
      
      rho <- mean(rhos)
    }
    else {
      rhos <- rep(0, group.num)
      for (g in c(1: group.num)) {
        rho_g <- selrho.useCVBIC(Xb0.COMi.glist[[g]], print.detail = print.detail)
        rho_g <- rho_g[1]
        rhos[g] <- rho_g
      }
      
      rho <- mean(rhos)
    }
    
    if (print.detail) {
      cat('Set rho = ', rho, '\n')
    }
    
    for (k in c(1: batch.num)) {
      if (batch.level[k] != ref.batch) {
        Xb1.Batk.COMi.glist <- list()
        for (g in c(1: group.num)) {
          Xb1.gi.Batk <- X.delout[batch.new == batch.level[k] & group.new == grp.level[g], ]
          Xb1.gi.Batk.COMi <- Xb1.gi.Batk[, metID]
          Xb1.Batk.COMi.glist[[g]] <- Xb1.gi.Batk.COMi
        }
        
        penterm <- findBestPara(Xb0.COMi.glist, Xb1.Batk.COMi.glist, rho, eps)
        ksi <- penterm$penal.ksi
        gama <- penterm$penal.gamma
        
        if (print.detail) {
          cat('Batch ', k, ' correction begin......')
        }
        
        para.out <- BEgLasso(Xb0.COMi.glist, Xb1.Batk.COMi.glist, rho, ksi, gama, eps)
        X1.cor.glist <- para.out$X1.cor
        Theta.list <- para.out$Theta
        coef.a <- para.out$coef.a
        coef.b <- para.out$coef.b
        
        if (print.detail) {
          cat('finshed', '\n')
        }
        
        for (g in c(1: group.num)) {
          X.cor.1[batch.new == batch.level[k] & group.new == grp.level[g], metID] <- X1.cor.glist[[g]]
          Xcor.para[[k]]$Theta[[g]][metID, metID] <- Theta.list[[g]]
        }
        Xcor.para[[k]]$coef.a[metID] <- coef.a
        Xcor.para[[k]]$coef.b[metID] <- coef.b
      
        Xb1.nodel <- X.nodel[batch == batch.level[k], ]
        N1 <- nrow(Xb1.nodel)
        coef.A <- diag(coef.a)
        coef.B <- rep(1, N1) %*% t(coef.b)
        X.cor[batch == batch.level[k], metID] <- Xb1.nodel[, metID] %*% coef.A + coef.B
        
        if ((!is.na(containQC)) & containQC) {
          X.cor.withQC[batch.init == batch.level[k] & group.init != 'QC', metID] <- 
            Xb1.nodel[, metID] %*% coef.A + coef.B
          
          QC.batk <- X_QC[QC_batch == batch.level[k], ]
          Nqc.batk <- nrow(QC.batk)
          coef.A <- diag(coef.a)
          coef.B <- rep(1, Nqc.batk) %*% t(coef.b)
          X.cor.withQC[batch.init == batch.level[k] & group.init == 'QC', metID] <- 
            QC.batk[, metID] %*% coef.A + coef.B
        }
        else {
          X.cor.withQC <- NULL
        }
      }
      
    }
    
    cat('Finish correction of community ', i, '\n')
  }
  
  Xcor <- list(batch.level = batch.level,
               delsampIdx = delsampIdx,
               batch.new = batch.new, 
               group.new = group.new, 
               X.delout = X.delout, 
               X.cor = X.cor, 
               X.cor.1 = X.cor.1, 
               X.cor.withQC = X.cor.withQC, 
               Xcor.para = Xcor.para)
  
  return(Xcor)
}

