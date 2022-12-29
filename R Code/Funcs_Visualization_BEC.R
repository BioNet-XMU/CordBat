# Title: Functions of Visualization
# Author: FanJing Guo
# Date: 2022.12.28

library(ggplot2)
library(mixOmics)
library(gridExtra)
library(igraph)
library(networkD3)
library(RBGL)
library(Rgraphviz)
library(networkD3)
library(ggsci)

#---------------------------------------------------------------------
# Principal component analysis (PCA)
#---------------------------------------------------------------------
PCA_scatter <- function(data, batch, group = NULL, expl.var, xlim, ylim, 
                         batch.legend.title = 'Batch', 
                         group.legend.title = 'Group', 
                         legend.pos = 'right', 
                         title = NULL, title.cex = 3, 
                         fon = 'sans') {
  data <- as.data.frame(data)
  batch <- as.factor(batch)
  
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", "#D55E00", "#999999", 
                 "#33CC99", "#990000", "#CC99FF", "#FF00FF", '#5500FF', '#66FF66', '#0088A8', '#A500CC')
  
  
  pMain <- ggplot(data = data, aes(x = data[, 1], y = data[, 2])) + 
    geom_point(aes(color = batch), shape = 16, size = 5) + 
    xlab(paste0('PC1: ', round(as.numeric(expl.var[1])*100, 1), '%')) + 
    ylab(paste0('PC2: ', round(as.numeric(expl.var[2])*100, 1), '%')) + 
    scale_color_manual(values = cbPalette) + 
    theme_bw() + 
    theme(legend.position = legend.pos, 
          axis.title.x = element_text(size = 25, family = fon, face = 'bold'), 
          axis.title.y = element_text(size = 25, family = fon, face = 'bold'), 
          axis.text.x = element_text(size = 18, family = fon, face = "bold"), 
          axis.text.y = element_text(size = 18, family = fon, face = "bold"), 
          panel.border = element_rect(size = 2)) + 
    xlim(xlim[1], xlim[2]) + ylim(ylim[1], ylim[2])
  
  if (legend.pos != 'none') {
    pMain + theme(legend.title = element_text(size = 20, family = fon, face = 'bold'), 
                  legend.text = element_text(size = 18, family = fon, face = 'bold')) + 
      labs(color = batch.legend.title, shape = group.legend.title, 
           title = title)
  }
  else {
    pMain + labs(title = title)
  }
  
}

PCA_scatter1 <- function(data, batch = NULL, group, expl.var, xlim, ylim, 
                         batch.legend.title = 'Batch', 
                         group.legend.title = 'Group', 
                         legend.pos = 'right', 
                         title = NULL, title.cex = 3, 
                         fon = 'sans') {
  data <- as.data.frame(data)
  group <- as.factor(group)
  grp.lev <- levels(group)
  
  dat_grp1 <- data[group == grp.lev[1], ]
  dat_grp2 <- data[group == grp.lev[2], ]
  grp1 <- group[group == grp.lev[1]]
  grp2 <- group[group == grp.lev[2]]
  
  cbPalette <- c("red", "blue", "green", "yellow", "gray", "black")
  
  
  pMain <- ggplot() + geom_point(data = dat_grp1, aes(x = dat_grp1[, 1], y = dat_grp1[, 2], fill = grp1), 
                                 shape = 21, color = 'black', size = 5, stroke = 1.5) +
    geom_point(data = dat_grp2, aes(x = dat_grp2[, 1], y = dat_grp2[, 2], fill = grp2), 
               shape = 21, color = 'black', size = 5, stroke = 1.5) + 
    xlab(paste0('PC1: ', round(as.numeric(expl.var[1])*100, 1), '%')) + 
    ylab(paste0('PC2: ', round(as.numeric(expl.var[2])*100, 1), '%')) + 
    # scale_fill_manual(values = cbPalette) + 
    scale_fill_manual(values = cbPalette) + 
    theme_bw() + 
    theme(legend.position = legend.pos, 
          axis.title.x = element_text(size = 25, family = fon, face = 'bold'), 
          axis.title.y = element_text(size = 25, family = fon, face = 'bold'), 
          axis.text.x = element_text(size = 18, family = fon, face = "bold"), 
          axis.text.y = element_text(size = 18, family = fon, face = "bold"), 
          panel.border = element_rect(size = 2), 
          plot.title = element_text(size = 25, family = fon, face = 'bold', hjust = 0.5)) + 
    xlim(xlim[1], xlim[2]) + ylim(ylim[1], ylim[2])
  
  if (legend.pos != 'none') {
    pMain + theme(legend.title = element_text(size = 20, family = fon, face = 'bold'), 
                  legend.text = element_text(size = 18, family = fon, face = 'bold')) + 
      labs(fill = group.legend.title, shape = batch.legend.title, 
           title = title)
  }
  else {
    pMain + labs(title = title)
  }
  
}

PCA_scatter_withQC <- function(data, batch, group, QC_label = 'QC', 
                                expl.var, xlim, ylim, 
                                batch.legend.title = 'Batch', 
                                group.legend.title = 'Group', 
                                legend.pos = 'right', 
                                title = NULL, title.cex = 3, 
                                fon = 'sans') {
  data <- as.data.frame(data)
  batch <- as.factor(batch)
  group <- as.factor(group)
  sbjdata <- data[group != QC_label, ]
  sbjbatch <- batch[group != QC_label]
  sbjgroup <- as.factor(rep('Sample', nrow(sbjdata)))
  
  QCdata <- data[group == QC_label, ]
  QCbatch <- batch[group == QC_label]
  QCgroup <- as.factor(rep('QC', nrow(QCdata)))
  
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", "#D55E00", "#999999", 
                 "#33CC99", "#990000", "#CC99FF", "#FF00FF", '#5500FF', '#66FF66', '#0088A8', '#A500CC')
  shapePool <- c(17, 16)
  
  
  pMain <- ggplot() + geom_point(data = sbjdata, aes(x = sbjdata[, 1], y = sbjdata[, 2], 
                                                     color = sbjbatch, shape = sbjgroup), 
                                 size = 5) + 
    geom_point(data = QCdata, aes(x = QCdata[, 1], y = QCdata[, 2], shape = QCgroup), 
               size = 5, color = 'black') + 
    xlab(paste0('PC1: ', round(as.numeric(expl.var[1])*100, 1), '%')) + 
    ylab(paste0('PC2: ', round(as.numeric(expl.var[2])*100, 1), '%')) + 
    scale_color_manual(values = cbPalette) + 
    scale_shape_manual(values = shapePool) + theme_bw() + 
    theme(legend.position = legend.pos, 
          axis.title.x = element_text(size = 30, family = fon, face = 'bold'), 
          axis.title.y = element_text(size = 30, family = fon, face = 'bold'), 
          axis.text.x = element_text(size = 25, family = fon, face = "bold"), 
          axis.text.y = element_text(size = 25, family = fon, face = "bold"), 
          panel.border = element_rect(linewidth = 3), 
          plot.title = element_text(hjust = 0.5, size = rel(title.cex), 
                                    family = fon, face = 'bold')) + 
    xlim(xlim[1], xlim[2]) + ylim(ylim[1], ylim[2])
  
  if (legend.pos != 'none') {
    pMain + theme(legend.title = element_text(size = 30, family = fon, face = 'bold'), 
                  legend.text = element_text(size = 25, family = fon, face = 'bold')) + 
      labs(color = batch.legend.title, shape = group.legend.title, 
           title = title)
  }
  else {
    pMain + labs(title = title)
  }
}

# ------------------------
# DE metabolites analysis
# -------------------------
DEana <- function(X, batch, top.rank = 1) {
  N <- nrow(X)
  p <- ncol(X)
  
  batch <- as.factor(batch)
  batch.num <- nlevels(batch)
  batch.level <- levels(batch)
  DEmetIdx <- c(1: p)
  P_val <- matrix(0, batch.num - 1, p)
  for (i in c(2: batch.num)) {
    
    for (j in 1: p) {
      # F-test
      vartt <- var.test(X[batch == batch.level[1], j], 
                        X[batch == batch.level[i], j])
      varPval <- vartt$p.value
      alpha <- 0.05
      if (varPval < alpha){
        varEqual <- F
      } else{
        varEqual <- T
      }
      stutt <- t.test(X[batch == batch.level[1], j], 
                      X[batch == batch.level[i], j], 
                      var.equal = varEqual)
      P_val[i - 1, j] <- stutt$p.value
    }
    
    if (length(intersect(DEmetIdx, which(P_val < 0.01))) == 0) {
      DEmetIdx <- intersect(DEmetIdx, which(P_val[i -1, ] < 0.05))
    }
    else {
      DEmetIdx <- intersect(DEmetIdx, which(P_val[i - 1, ] < 0.01))
    }
  }
  
  P_val <- apply(P_val, MARGIN = 2, FUN = max)
  DE.P_val <- P_val[DEmetIdx]
  sort.DE.P_val <- order(DE.P_val)
  DEmetIdx <- DEmetIdx[sort.DE.P_val]
  sel.Met <- DEmetIdx[top.rank]
  
  return(sel.Met)
}

DEmet_plot_merge <- function(dat.bef, dat.cor, batch, InjOrd, fon = 'sans') {
  N <- nrow(dat.bef)
  met.names <- colnames(dat.bef)
  met1.max <- max(c(dat.bef[, 1], dat.cor[, 1]))
  met1.min <- min(c(dat.bef[, 1], dat.cor[, 1]))
  met2.max <- max(c(dat.bef[, 2], dat.cor[, 2]))
  met2.min <- min(c(dat.bef[, 2], dat.cor[, 2]))
  
  dat.bef <- dat.bef[order(InjOrd), ]
  dat.cor <- dat.cor[order(InjOrd), ]
  batch <- batch[order(InjOrd)]
  InjOrd <- sort(InjOrd)
  
  dat.bef <- c(dat.bef[, 1], dat.bef[, 2])
  dat.cor <- c(dat.cor[, 1], dat.cor[, 2])
  dat <- c(dat.bef, dat.cor)
  dat.type <- factor(c(rep('Raw', 2 * N), rep('CordBat', 2 * N)), 
                     levels = c('Raw', 'CordBat'))
  metabo <- factor(c(rep(met.names[1], N), rep(met.names[2], N), 
                     rep(met.names[1], N), rep(met.names[2], N)), 
                   levels = c(met.names[1], met.names[2]))
  dat.bat <- factor(rep(batch, 4), levels = c(1: max(batch)))
  dat.ord <- rep(InjOrd, 4)
  
  dat <- data.frame(dat.type = dat.type, 
                    metabo = metabo, 
                    dat.ord = dat.ord, 
                    dat.bat = dat.bat, 
                    met.dat = dat)
  
  
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", "#D55E00", "#999999", 
                 "#33CC99", "#990000", "#CC99FF", "#FF00FF", '#5500FF', '#66FF66', '#0088A8', '#A500CC')
  
  # plot
  blank_data <- data.frame(metabo = c(rep(met.names, each = 2)), 
                           x = 0, 
                           y = c(met1.min, met1.max, met2.min, met2.max))
  
  ggplot(data = dat, aes(x = dat.ord, y = met.dat)) + 
    geom_point(aes(color = dat.bat), size = 3) + 
    scale_color_manual(values = cbPalette) + 
    geom_blank(data = blank_data, aes(x = x, y = y)) + 
    facet_grid(metabo ~ dat.type, scales = 'free_y') + 
    theme_bw() + 
    theme(strip.text = element_text(size = 15, colour = 'black', family = fon, face = 'bold'), 
          strip.background = element_rect(fill = 'transparent', size = 2), 
          strip.placement = 'outside', 
          axis.title.x = element_text(size = 25, family = fon, face = 'bold'), 
          axis.title.y = element_text(size = 25, family = fon, face = 'bold'), 
          axis.text.x = element_text(size = 18, family = fon, face = "bold"), 
          axis.text.y = element_text(size = 18, family = fon, face = "bold"), 
          panel.border = element_rect(fill = NA, colour = 'black'),
          axis.line = element_line(colour = 'black', size = 1), 
          legend.position = 'right', 
          legend.title = element_text(size = 20, family = fon, face = 'bold'), 
          legend.text = element_text(size = 18, family = fon, face = 'bold')) + 
    xlab('Sample Index') + ylab('Metabolite Intensity') + 
    labs(color = 'Batch', title = NULL)
  
}


# --------------------------------
# Data structure preserving
# --------------------------------
getdiffPccMat <- function(X1, X2) {
  N.1 <- nrow(X1)
  p <- ncol(X1)
  X1.pcc <- cor(X1)
  
  N.2 <- nrow(X2)
  p <- ncol(X2)
  X2.pcc <- cor(X2)
  
  X1.pcc.ndiag <- X1.pcc - diag(diag(X1.pcc))
  X2.pcc.ndiag <- X2.pcc - diag(diag(X2.pcc))
  
  # Fisher-z transformation
  X1.fz <- 1 / 2 * log((matrix(1, p, p) + X1.pcc.ndiag) / 
                         (matrix(1, p, p) - X1.pcc.ndiag))
  X2.fz <- 1 / 2 * log((matrix(1, p, p) + X2.pcc.ndiag) / 
                         (matrix(1, p, p) - X2.pcc.ndiag))
  
  # Z test
  SE <- sqrt(1 / (N.1 - 3) + 1 / (N.2 - 3))
  dfz <- abs(X1.fz - X2.fz)
  Z <- dfz / SE
  
  # transform Z to adjacent matrix A
  pcc.result <- list(X1.pcc.ndiag = X1.pcc.ndiag, 
                     X2.pcc.ndiag = X2.pcc.ndiag,
                     diffMat = Z)
  return(pcc.result)
}

getstdPcc <- function(X.intraCor, X.bat, Notdiff.ratio) {
  p <- ncol(X.intraCor)
  X.bat <- as.factor(X.bat)
  batch.num <- nlevels(X.bat)
  # 1. Compute PCC for each batch
  PCC.list <- list()
  for (i in c(1: batch.num)) {
    Dat.bati <- X.intraCor[X.bat == i, ]
    PCC.bati <- cor(Dat.bati)
    PCC.list[[i]] <- PCC.bati
  }
  
  # 2. Compute Z between batches
  diff.result.num <- batch.num * (batch.num - 1) / 2
  Z.list <- list()
  P.list <- list()
  difflist2bat <- matrix(0, diff.result.num, 2)
  count <- 0
  for (i in c(1: (batch.num - 1))) {
    Dat.bati <- X.intraCor[X.bat == i, ]
    for (j in c((i + 1): batch.num)) {
      Dat.batj <- X.intraCor[X.bat == j, ]
      diffPCC <- getdiffPccMat(X1 = Dat.bati, 
                               X2 = Dat.batj)
      Z.bati_batj <- diffPCC$diffMat
      
      count <- count + 1
      Z.list[[count]] <- Z.bati_batj
      # Z to P
      P.list[[count]] <- 2 * pnorm(Z.bati_batj, lower.tail = F)
      difflist2bat[count, ] <- c(i, j)
    }
  }
  
  # adjust P values
  P.mat <- matrix(0, diff.result.num, p * (p - 1) / 2)
  count <- 0
  for (i in c(1: (p - 1))) {
    for (j in c((i + 1): p)) {
      count <- count + 1
      for (k in c(1: length(P.list))) {
        P.k <- P.list[[k]]
        P.mat[k, count] <- P.k[i, j]
      }
    }
  }
  adjust.P.mat <- matrix(p.adjust(P.mat, method = "BH"), 
                         diff.result.num, p * (p - 1) / 2)
  
  # 3. determine the retaining PCCs
  std.PCC <- matrix(NA, p, p)
  count <- 0
  for (i in c(1: (p - 1))) {
    for (j in c((i + 1): p)) {
      count <- count + 1
      Notdiff.num <- sum(adjust.P.mat[, count] >= 0.05)
      Notdiff.listid <- which(adjust.P.mat[, count] >= 0.05) 

      if (Notdiff.num >= round(diff.result.num * Notdiff.ratio)) {
        Not.diff.bat <- difflist2bat[Notdiff.listid, ]
        Not.diff.bat <- unique(as.vector(Not.diff.bat))
        
        std.PCC[i, j] <- 0
        for (n in c(1: length(Not.diff.bat))) {
          pcc.val <- PCC.list[[Not.diff.bat[n]]][i, j]
          std.PCC[i, j] <- std.PCC[i, j] + pcc.val
        }
        std.PCC[i, j] <- std.PCC[i, j] / length(Not.diff.bat)
        std.PCC[j, i] <- std.PCC[i, j]
      }
    }
  }
  
  retain.PCC.num <- length(std.PCC[!is.na(std.PCC)]) / 2
  
  std.PCC.result <- list(retain.PCC.num = retain.PCC.num, 
                         std.PCC.Mat = std.PCC)
  
  return(std.PCC.result)
}


getPdiff2Std <- function(std.pcc, X.dat, conf_alpha = 0.05) {
  N <- nrow(X.dat)
  p <- ncol(X.dat)
  X.pcc <- cor(X.dat)
  
  tot.pcc.num <- length(std.pcc[!is.na(std.pcc)]) / 2
  std.pcc1 <- std.pcc
  std.pcc1[is.na(std.pcc)] <- 0
  X.pcc[is.na(std.pcc)] <- 0
  
  # Fisher-z transformation
  std.fz <- 1 / 2 * log((matrix(1, p, p) + std.pcc1) / (matrix(1, p, p) - std.pcc1))
  X.fz <- 1 / 2 * log((matrix(1, p, p) + X.pcc) / (matrix(1, p, p) - X.pcc))
  
  # Z test
  SE <- sqrt(2 / (N - 3))
  dfz <- abs(std.fz - X.fz)
  Z <- dfz / SE
  # Z to P
  P <- 2 * pnorm(Z, lower.tail = F)
  P <- matrix(p.adjust(P, method = "BH"),
              nrow(P), ncol(P))
  
  diff.pcc.num <- sum(P[!is.na(std.pcc)] < conf_alpha) / 2
  diff.pcc.ratio <- diff.pcc.num / tot.pcc.num
  
  diff2std.result <- list(diff.pcc.num = diff.pcc.num, 
                          diff.ratio = diff.pcc.ratio, 
                          P.Mat = P)
  return(diff2std.result)
}

getLargestComp <- function(P.mat, conf_alpha = 0.05) {
  # transform Z to adjacent matrix A
  Adiff <- 1 * (P.mat < conf_alpha)
  nlogP.mat <- -log10(P.mat)
  Pdiff <- nlogP.mat * (P.mat < conf_alpha)
  
  G <- graph_from_adjacency_matrix(Adiff, weighted = NULL, mode = "undirected")
  G <- as_graphnel(G)
  connectComp <- connectedComp(G)
  connectComp.size <- lengths(connectComp)
  c.Comp.nodesid <- c()
  for (i in c(1: length(connectComp.size))) {
    if (connectComp.size[i] != 1) {
      Compi.nodes <- connectComp[[i]]
      nodes.id <- as.numeric(Compi.nodes)
      c.Comp.nodesid <- c(c.Comp.nodesid, nodes.id)
    }
  }
  P.connectComp <- Pdiff[c.Comp.nodesid, c.Comp.nodesid]
  
  if (max(connectComp.size) != 1) {
    largestCompIdx <- which(connectComp.size == max(connectComp.size))
    if (length(largestCompIdx) > 1) {
      largestComp.nodes <- connectComp[[largestCompIdx[1]]]
    }
    else {
      largestComp.nodes <- connectComp[[largestCompIdx]]
    }
    
    node.id <- as.numeric(largestComp.nodes)
    
    P.largestComp <- Pdiff[node.id, node.id]
  }
  else {
    node.id <- NULL
    P.largestComp <- NULL
  }
  
  LargestComp <- list(P.connectComp = P.connectComp, 
                      LCnodes = node.id, 
                      P.largestComp = P.largestComp)
  return(LargestComp)
  
}

P_to_dataframe <- function(P) {
  p <- ncol(P)
  
  nodes <- data.frame(node = c(1: p) - 1, 
                      group = rep(1, p))
  node1 <- c()
  node2 <- c()
  w <- c()
  for (i in c(1: (p - 1))) {
    nodesId <- which(P[i, ] != 0)
    nodesId <- nodesId[nodesId > i]
    nodeNum <- length(nodesId)
    if (nodeNum != 0) {
      node1 <- c(node1, rep(i - 1, nodeNum))
      node2 <- c(node2, nodesId - 1)
      w <- c(w, P[i, nodesId])
    }
  }
  
  edges <- data.frame(node1 = node1, 
                      node2 = node2, 
                      weight = w)
  
  networkDat <- list(nodes = nodes, 
                     edges = edges)
  
  return(networkDat)
}


Graph_plot_forP <- function(P) {
  NetworkDat <- P_to_dataframe(P)
  Nodes <- NetworkDat$nodes
  Edges <- NetworkDat$edges
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # plot
  forceNetwork(Links = Edges, Nodes = Nodes, Source = "node1",
               Target = "node2", Value = "weight", NodeID = "node", 
               Group = "group", opacity = 1, zoom = TRUE)
}


# ---------------------------
# add noise to QC samples
# ---------------------------
addnoise_to_QC <- function(dat, SNR, dolog = T) {
  N <- nrow(dat)
  p <- ncol(dat)
  
  if (dolog) {
    dat.log <- log2(dat)
  }
  else {
    dat.log <- dat
  }
  
  Noise <- matrix(0, N, p)
  for (i in c(1: p)) {
    Noise[, i] <- rnorm(N, mean = 0, sd = 1)
  }
  
  dat.m <- apply(dat.log, MARGIN = 2, FUN = mean)
  noi.sd <- dat.m * 2^(- SNR / 20)
  Noise <- Noise %*% diag(noi.sd)
  
  dat.addnoise <- dat.log + Noise
  
  addnoise.result <- list(dat = dat.addnoise, 
                          addnoise = Noise)
  return(addnoise.result)
}


# --------------------------------
# Biological effect preserving
# ----------------------------------
# F-test and t-test
getttestResult <- function(dat_g1, dat_g2) {
  p <- ncol(dat_g1)
  
  dat.tt.tval <- rep(0, p)
  dat.tt.Pval <- rep(0, p)
  for (i in 1: p) {
    # F test
    vartt <- var.test(dat_g1[, i], dat_g2[, i])
    varPval <- vartt$p.value
    alpha <- 0.05
    if (varPval < alpha){
      varEqual <- F
    } else{
      varEqual <- T
    }
    stutt <- t.test(dat_g1[, i], dat_g2[, i], var.equal = varEqual)
    dat.tt.tval[i] <- stutt$statistic
    dat.tt.Pval[i] <- stutt$p.value
  }
  
  tt.result <- list(tt.tval = dat.tt.tval, 
                    tt.Pval = dat.tt.Pval)
  return(tt.result)
}

# get overlapped metabolites vs the ground truth
getOLperc <- function(GTmetID, Pval) {
  Pval.sortID <- order(Pval)
  
  top.num <- seq(2, 120, 2)
  OLmetNum <- rep(0, length(top.num))
  OLperc <- rep(0, length(top.num))
  for (i in c(1: length(top.num))) {
    DEmet.Pval.topID <- Pval.sortID[c(1: top.num[i])]
    OLmetId <- intersect(GTmetID, DEmet.Pval.topID)
    OLmetNum[i] <- length(OLmetId)
    OLperc[i] <- length(OLmetId) / length(DEmet.Pval.topID)
  }
  
  OLresult <- list(OLmetNum = OLmetNum,  
                   OLperc = OLperc)
  return(OLresult)
}

# SVM with ensemble
SVMensemble <- function(dat.tr, tr.grp, dat.te, Tn = 1000) {
  # group is encoded with 0 and 1
  # grp.lev[1] = 0, grp.lev[2] = 1
  tr.grp <- as.factor(tr.grp)
  grp.lev <- levels(tr.grp)
  
  grp1.trnum <- sum(tr.grp == grp.lev[1])
  grp2.trnum <- sum(tr.grp == grp.lev[2])
  
  if (grp1.trnum < grp2.trnum) {
    minor.grp <- as.numeric(grp.lev[1])
    major.grp <- as.numeric(grp.lev[2])
  }
  else {
    minor.grp <- as.numeric(grp.lev[2])
    major.grp <- as.numeric(grp.lev[1])
  }
  
  dat.tr.minor <- dat.tr[tr.grp == minor.grp, ]
  dat.tr.major <- dat.tr[tr.grp == major.grp, ]
  N.minor <- nrow(dat.tr.minor)
  N.major <- nrow(dat.tr.major)
  
  pred.all <- matrix(0, nrow(dat.te), Tn)
  for (i in c(1: Tn)) {
    major.subsetIdx <- sample(c(1: N.major), N.minor)
    dat.tr.major.sub <- dat.tr.major[major.subsetIdx, ]
    dat.tr.i <- rbind(dat.tr.major.sub, dat.tr.minor)
    # probability
    dat.tr.i.grp <- as.factor(c(rep(major.grp, N.minor), rep(minor.grp, N.minor)))
    model.i <- svm(dat.tr.i, dat.tr.i.grp, probability = T)
    pred.i <- predict(model.i, dat.te, probability = T)
    pred.prob.i <- attr(pred.i, "probabilities")
    colIdx <- which(as.numeric(colnames(pred.prob.i)) == 1)
    pred.prob.i <- pred.prob.i[, colIdx]
    pred.all[, i] <- pred.prob.i
    
    if (i %% 1000 == 0) {
      cat(i, '\n')
    }
    
  }
  
  pred <- apply(pred.all, MARGIN = 1, FUN = median)
  
  return(pred)
}

# CV for SVM ensemble
CVforSVMensemble <- function(dat, grp, fold = 7) {
  dat.g1 <- dat[grp == 1, ]
  dat.g2 <- dat[grp == 2, ]
  N.g1 <- nrow(dat.g1)
  N.g2 <- nrow(dat.g2)
  
  # disturb sample order
  dat.g1 <- dat.g1[sample(c(1: N.g1), N.g1), ]
  dat.g2 <- dat.g2[sample(c(1: N.g2), N.g2), ]
  
  CVsize1_6.g1 <- round(N.g1 / fold)
  CVsize1_6.g2 <- round(N.g2 / fold)
  
  clas.pred <- c()
  clas.GT <- c()
  cv.err <- rep(0, fold)
  accuracy <- rep(0, fold)
  for (i in c(1: fold)) {
    start.index1 <- (i-1) * CVsize1_6.g1 + 1
    if (i != fold) {
      end.index1 <- i * CVsize1_6.g1
    }
    else {
      end.index1 <- N.g1
    }
    dat.cv.g1 <- dat.g1[c(start.index1: end.index1), ]
    
    start.index2 <- (i-1) * CVsize1_6.g2 + 1
    if (i != fold) {
      end.index2 <- i * CVsize1_6.g2
    }
    else {
      end.index2 <- N.g2
    }
    dat.cv.g2 <- dat.g2[c(start.index2: end.index2), ]
    
    dat.cv <- rbind(dat.cv.g1, dat.cv.g2)
    cv.grp <- c(rep(1, nrow(dat.cv.g1)), rep(0, nrow(dat.cv.g2)))
    
    dat.tr.g1 <- dat.g1[-c(start.index1: end.index1), ]
    dat.tr.g2 <- dat.g2[-c(start.index2: end.index2), ]
    dat.tr <- rbind(dat.tr.g1, dat.tr.g2)
    tr.grp <- c(rep(1, nrow(dat.tr.g1)), rep(0, nrow(dat.tr.g2)))
    tr.grp <- as.factor(tr.grp)
    
    clas.pred.i <- SVMensemble(dat.tr, tr.grp, dat.cv)
    clas.pred <- c(clas.pred, clas.pred.i)
    clas.GT <- c(clas.GT, cv.grp)
    
    # MSE
    cv.err[i] <- 1 / (nrow(dat.cv)) * sum((clas.pred.i - as.numeric(cv.grp))^2)
    accuracy[i] <- 1 - cv.err[i]
  }
  
  clas.roc <- roc(as.numeric(clas.GT), as.numeric(clas.pred), quiet = T)
  auc.val <- as.numeric(clas.roc$auc)
  roc.spec <- clas.roc$specificities
  roc.sens <- clas.roc$sensitivities

  clas1.GT.ID <- which(clas.GT == 1)
  clas1.ID <- which(clas.pred > 0.5)
  
  clasROC.result <- list(clas1.GT.ID = clas1.GT.ID, 
                         clas1.ID = clas1.ID, 
                         cv.err = cv.err, 
                         accuracy = accuracy, 
                         auc.val = auc.val, 
                         spec = roc.spec, 
                         sens = roc.sens)
  
  return(clasROC.result)
}