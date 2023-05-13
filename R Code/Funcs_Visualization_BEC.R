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
          axis.text.x = element_text(size = 22, family = fon, face = "bold", color = 'black'), 
          axis.text.y = element_text(size = 22, family = fon, face = "bold", color = 'black'), 
          axis.ticks.x.bottom = element_line(linewidth = 1, lineend = 10, color = 'black'), 
          axis.ticks.y.left = element_line(linewidth = 1, lineend = 10, color = 'black'), 
          panel.border = element_rect(linewidth = 2), 
          panel.grid = element_blank()) + 
    xlim(xlim[1], xlim[2]) + ylim(ylim[1], ylim[2])
  
  if (legend.pos != 'none') {
    pMain + theme(legend.title = element_text(size = 25, family = fon, face = 'bold'), 
                  legend.text = element_text(size = 20, family = fon, face = 'bold')) + 
      labs(color = batch.legend.title, shape = group.legend.title, 
           title = title)
  }
  else {
    pMain + labs(title = title)
  }
  
}

PCA_scatter1 <- function(data, group, expl.var, xlim, ylim, 
                         legend.pos = c(0.92, 0.92), 
                         title = NULL, title.cex = 3, 
                         fon = 'sans') {
  dat <- as.data.frame(data)
  group <- as.factor(group)
  grp.lev <- levels(group)
  
  cbPalette <- c("red", "blue", "green", "yellow", "gray", "black")
  
  pMain <- ggplot() + geom_point(data = dat, aes(x = dat[, 1], y = dat[, 2], 
                                                 fill = group), 
                                 shape = 21, color = 'black', size = 5, stroke = 1.5) + 
    xlab(paste0('PC1: ', round(as.numeric(expl.var[1])*100, 1), '%')) + 
    ylab(paste0('PC2: ', round(as.numeric(expl.var[2])*100, 1), '%')) + 
    scale_fill_manual(values = cbPalette) + 
    theme_bw() + 
    theme(legend.position = legend.pos, 
          axis.title.x = element_text(size = 25, family = fon, face = 'bold'), 
          axis.title.y = element_text(size = 25, family = fon, face = 'bold'), 
          axis.text.x = element_text(size = 22, family = fon, face = "bold", color = 'black'), 
          axis.text.y = element_text(size = 22, family = fon, face = "bold", color = 'black'), 
          axis.ticks = element_line(color = 'black', linewidth = 1), 
          axis.ticks.length = unit(2, 'mm'), 
          panel.border = element_rect(linewidth = 2), 
          panel.grid = element_blank(), 
          plot.title = element_text(size = 25, family = fon, face = 'bold', hjust = 0.5)) + 
    xlim(xlim[1], xlim[2]) + ylim(ylim[1], ylim[2])
  
  if (legend.pos[1] != 'none') {
    pMain + theme(legend.background = element_rect(color = 'black', linewidth = 1), 
                  legend.text = element_text(size = 18, family = fon, face = 'bold')) + 
      labs(fill = NULL, title = title)
  }
  else {
    pMain + labs(title = title)
  }
}

PCA_scatter_withQC <- function(data, batch, group, QC_label = 'QC', 
                                expl.var, xlim, ylim, 
                                batch.legend.title = 'Combined \nBatch', 
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
          axis.text.x = element_text(size = 25, family = fon, face = "bold", color = 'black'), 
          axis.text.y = element_text(size = 25, family = fon, face = "bold", color = 'black'), 
          axis.ticks = element_line(color = 'black', linewidth = 1, lineend = 20), 
          panel.border = element_rect(linewidth = 3), 
          panel.grid = element_blank(), 
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
  dat.type <- factor(c(rep('Uncorrect', 2 * N), rep('CordBat', 2 * N)), 
                     levels = c('Uncorrect', 'CordBat'))
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
    theme(strip.text = element_text(size = 25, colour = 'black', family = fon, face = 'bold'), 
          strip.background = element_blank(), 
          strip.placement = 'outside', 
          axis.title.x = element_text(size = 25, family = fon, face = 'bold'), 
          axis.title.y = element_text(size = 25, family = fon, face = 'bold', vjust = 2), 
          axis.text.x = element_text(size = 20, family = fon, face = "bold", color = 'black'), 
          axis.text.y = element_blank(), 
          axis.ticks.x = element_line(color = 'black', linewidth = 1, lineend = 10), 
          axis.ticks.y = element_blank(), 
          panel.border = element_rect(fill = NA, colour = 'black', linewidth = 2),
          panel.grid = element_blank(), 
          legend.position = 'none', 
          plot.margin = margin(2, 2, 5, 10)) + 
    xlab('Sample Index') + ylab('Intensity (a.u.)') + 
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
  retain.PCC.perc <- retain.PCC.num / (p * (p - 1) / 2)
  
  std.PCC.result <- list(retain.PCC.num = retain.PCC.num, 
                         retain.PCC.perc = retain.PCC.perc, 
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
addnoise_to_QC <- function(dat, SNR, dolog = F) {
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

setArrowPos <- function(x0, y0, pos, ang, len = NA) {
  
  if (is.na(len)) {
    len <- 9
  }
  ex <- 1.2
  ey <- 1.2
  dx <- len * sin(ang)
  dy <- len * cos(ang)
  
  arrow.x <- switch(pos, 
                    top = c(x0, x0), 
                    bottom = c(x0, x0), 
                    left = c(x0 - ex - len, x0 - ex), 
                    right = c(x0 + ex + len, x0 + ex), 
                    topleft = c(x0 - ex - dx, x0 - ex), 
                    topright = c(x0 + ex + dx, x0 + ex), 
                    bottomleft = c(x0 - ex - dx, x0 - ex), 
                    bottomright = c(x0 + ex + dx, x0 + ex))
  
  arrow.y <- switch(pos, 
                    top = c(y0 + ey + len, y0 + ey), 
                    bottom = c(y0 - ey - len, y0 - ey), 
                    left = c(y0, y0), 
                    right = c(y0, y0), 
                    topleft = c(y0 + ey + dy, y0 + ey), 
                    topright = c(y0 + ey + dy, y0 + ey), 
                    bottomleft = c(y0 - ey - dy, y0 - ey), 
                    bottomright = c(y0 - ey - dy, y0 - ey))
  
  arrow.pos <- data.frame(x = arrow.x, 
                          y = arrow.y)
  
  return(arrow.pos)
}

library(dplyr)
plotCumRSD_m <- function(cumPerc.list, 
                         tar.method.name = 'CordBat', 
                         rsd.gap = 0.001, 
                         mark.rsd = 20, 
                         arrow.pos = c('topleft', 'topright', 'bottomright'), 
                         arrow.angle = c(pi/4, pi/20, pi/3), 
                         arrow.len = c(9, 9, 9)) {
  
  CdB.cumPerc <- cumPerc.list[[tar.method.name]]
  SR.cumPerc <- cumPerc.list$SERRF
  SR03.cumPerc <- cumPerc.list$SERRF_03
  SR045.cumPerc <- cumPerc.list$SERRF_045
  SR06.cumPerc <- cumPerc.list$SERRF_06
  
  # area
  rsd.num <- 1 / rsd.gap
  area.CdB <- (CdB.cumPerc[1] + CdB.cumPerc[rsd.num] + 
                 2 * sum(CdB.cumPerc[c(2: (rsd.num - 1))])) * rsd.gap / 2
  
  area.SR <- (SR.cumPerc[1] + SR.cumPerc[rsd.num] + 
                2 * sum(SR.cumPerc[c(2: (rsd.num - 1))])) * rsd.gap / 2
  
  area.SR03 <- (SR03.cumPerc[1] + SR03.cumPerc[rsd.num] + 
                2 * sum(SR03.cumPerc[c(2: (rsd.num - 1))])) * rsd.gap / 2
  
  area.SR045 <- (SR045.cumPerc[1] + SR045.cumPerc[rsd.num] + 
                2 * sum(SR045.cumPerc[c(2: (rsd.num - 1))])) * rsd.gap / 2
  
  area.SR06 <- (SR06.cumPerc[1] + SR06.cumPerc[rsd.num] + 
                2 * sum(SR06.cumPerc[c(2: (rsd.num - 1))])) * rsd.gap / 2
  
  labelCdB <- paste0(tar.method.name, ', Area = ', round(area.CdB, 2))
  labelSR <- paste0('SERRF, Area = ', round(area.SR, 2))
  labelSR03 <- paste0('SERRF(0.030), Area = ', round(area.SR03, 2))
  labelSR045 <- paste0('SERRF(0.045), Area = ', round(area.SR045, 2))
  labelSR06 <- paste0('SERRF(0.060), Area = ', round(area.SR06, 2))
  
  dat.type <- factor(c(rep(labelCdB, rsd.num), 
                       rep(labelSR, rsd.num), 
                       rep(labelSR03, rsd.num), 
                       rep(labelSR045, rsd.num), 
                       rep(labelSR06, rsd.num)), 
                     levels = c(labelCdB, labelSR, 
                                labelSR03, labelSR045, 
                                labelSR06))
  
  all.cumPerc <- c(CdB.cumPerc, SR.cumPerc, 
                   SR03.cumPerc, SR045.cumPerc, 
                   SR06.cumPerc)
  
  rsd.val <- seq(rsd.gap, 1, rsd.gap)
  dat.cumPerc <- data.frame(rsd.val = rep(rsd.val * 100, 5), 
                            cumPerc = all.cumPerc * 100)
  
  # plot
  fon <- 'sans'
  idx20 <- which(dat.cumPerc$rsd.val == mark.rsd)
  p.idx <- idx20
  cpers20 <- dat.cumPerc[p.idx, ]
  set_col <- c('red', 'green', 'magenta', 'blue', 'black')
  set_linetype <- c(2, 1, 1, 2, 1)
  
  a1.pos <- setArrowPos(cpers20[1, 1], cpers20[1, 2], pos = arrow.pos[1], 
                        ang = arrow.angle[1], len = arrow.len[1])
  a2.pos <- setArrowPos(cpers20[2, 1], cpers20[2, 2], pos = arrow.pos[2], 
                        ang = arrow.angle[2], len = arrow.len[2])
  a3.pos <- setArrowPos(cpers20[4, 1], cpers20[4, 2], pos = arrow.pos[3], 
                        ang = arrow.angle[3], len = arrow.len[3])
  
  a.pos.list <- list(a1 = a1.pos, 
                     a2 = a2.pos, 
                     a3 = a3.pos)
  text.pos.list <- list()
  dy <- 3.5
  dx <- 6
  for (i in c(1: length(a.pos.list))) {
    pos.x.i <- a.pos.list[[i]]$x
    pos.y.i <- a.pos.list[[i]]$y
    
    if (length(grep('top', arrow.pos[i])) != 0) {
      text.pos.list[[i]] <- data.frame(x = pos.x.i[1], 
                                       y = pos.y.i[1] + dy)
    }
    else if (length(grep('bottom', arrow.pos[i])) != 0) {
      text.pos.list[[i]] <- data.frame(x = pos.x.i[1], 
                                       y = pos.y.i[1] - dy)
    }
    else if (length(grep('left', arrow.pos[i])) != 0) {
      text.pos.list[[i]] <- data.frame(x = pos.x.i[1] - dx, 
                                       y = pos.y.i[1])
    }
    else {
      text.pos.list[[i]] <- data.frame(x = pos.x.i[1] + dx, 
                                       y = pos.y.i[1])
    }
  }

  ggplot() +
    geom_line(data = dat.cumPerc, 
              aes(x = rsd.val, y = cumPerc, color = dat.type, linetype = dat.type), 
              linewidth = 1.3) + 
    geom_vline(xintercept = mark.rsd, linetype = 2, linewidth = 1, color = '#00FFFFFF') + 
    geom_segment(data = a1.pos, aes(x = x[1], xend = x[2], y = y[1], yend = y[2]), 
                 arrow = arrow(angle = 25, type = 'open', length = unit(3, 'mm')), 
                 linewidth = 1.2, lineend = 'round', linejoin = 'mitre', color = 'red') + 
    annotate(geom = 'text', x = text.pos.list[[1]]$x, y = text.pos.list[[1]]$y, 
             label = paste0(round(cpers20[1, 2], 1), '%'), 
             size = 6, color = 'red', family = 'sans', fontface = 'bold') + 
    geom_segment(data = a2.pos, aes(x = x[1], xend = x[2], y = y[1], yend = y[2]), 
                 arrow = arrow(angle = 25, type = 'open', length = unit(3, 'mm')), 
                 linewidth = 1.2, lineend = 'round', linejoin = 'mitre', color = 'green') + 
    annotate(geom = 'text', x = text.pos.list[[2]]$x, y = text.pos.list[[2]]$y, 
             label = paste0(round(cpers20[2, 2], 1), '%'), 
             size = 6, color = 'green', family = 'sans', fontface = 'bold') + 
    geom_segment(data = a3.pos, aes(x = x[1], xend = x[2], y = y[1], yend = y[2]), 
                 arrow = arrow(angle = 25, type = 'open', length = unit(3, 'mm')), 
                 linewidth = 1.2, lineend = 'round', linejoin = 'mitre', color = 'blue') + 
    annotate(geom = 'text', x = text.pos.list[[3]]$x, y = text.pos.list[[3]]$y, 
             label = paste0(round(cpers20[4, 2], 1), '%'), 
             size = 6, color = 'blue', family = 'sans', fontface = 'bold') + 
    scale_color_manual(values = set_col) + 
    scale_linetype_manual(values = set_linetype) + 
    scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) + 
    theme_bw() + guides(linetype = "none") + 
    theme(legend.position = c(0.7, 0.17), 
          axis.title = element_text(size = 25, family = fon, face = 'bold'), 
          axis.text = element_text(size = 20, family = fon, face = 'bold', color = 'black'),
          legend.box.spacing = unit(1, 'cm'), 
          legend.text = element_text(size = 18, family = fon, 
                                     margin = margin(b = 5)),
          panel.border = element_rect(linewidth = 2.5), 
          panel.grid = element_blank()) + 
    xlab('RSD (%)') + 
    ylab('Percentage of metabolite (%)') + ylim(-2, 102) +
    labs(color = NULL, title = NULL)
  
}

# --------------------------------
# Biological effect preserving
# ----------------------------------
# F-test and t-test
getttestResult <- function(dat_g1, dat_g2, BH.adj = T) {
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
  
  if (BH.adj) {
    dat.tt.Pval <- p.adjust(dat.tt.Pval, method = 'BH')
  } 
  
  tt.result <- list(tt.tval = dat.tt.tval, 
                    tt.Pval = dat.tt.Pval)
  return(tt.result)
}

# get overlapped metabolites vs the ground truth
getOLperc <- function(GT.ID, sort.Idx, gap = 1, max.topRank = NA) {
  
  GT.num <- length(GT.ID)
  if (is.na(max.topRank)) {
    max.topRank <- GT.num
  }
  
  top.num <- seq(gap, max.topRank, gap)
  OL.Num <- rep(0, length(top.num))
  OLperc <- rep(0, length(top.num))
  for (i in c(1: length(top.num))) {
    DE.topID <- sort.Idx[c(1: top.num[i])]
    OLId <- intersect(GT.ID, DE.topID)
    OL.Num[i] <- length(OLId)
    OLperc[i] <- length(OLId) / length(DE.topID)
  }
  
  OLresult <- list(topRank = top.num, 
                   OL.Num = OL.Num,  
                   OLperc = OLperc)
  return(OLresult)
}

# score plot for PLS-DA model
PLSDA_scatter <- function(data, group, xlim, ylim, 
                          legend.pos = 'none', 
                          title = NULL, title.cex = 3, fon = 'sans') {
  dat <- as.data.frame(data)
  group <- as.factor(group)
  grp.lev <- levels(group)
  
  set_col <- c("blue", "red", "green", "yellow", "gray", "black")
  set_shape <- c(16, 17, 18, 15, 19, 20)
  
  x.len <- max(dat[, 1]) - min(dat[, 1])
  y.len <- max(dat[, 2]) - min(dat[, 2])
  dx <- x.len / 8
  dy <- y.len / 18
  
  pMain <- ggplot() + geom_point(data = dat, aes(x = dat[, 1], y = dat[, 2], 
                                                 color = group, shape = group), 
                                 size = 5) + 
    geom_vline(xintercept = 0, linetype = 2) + 
    geom_hline(yintercept = 0, linetype = 2) + 
    stat_ellipse(data = dat[group == grp.lev[1], ], 
                 aes(x = dat[group == grp.lev[1], 1], y = dat[group == grp.lev[1], 2]), 
                 color = 'black', fill = NA, type = 'norm', geom = 'polygon', 
                 level = 0.95, linetype = 2) + 
    annotate(geom = 'text', x = max(dat[, 1]) + 1 - dx, y = min(dat[, 2]) - 1 + dy, 
             label = title, size = 8, color = 'black', family = 'sans', fontface = 'bold') + 
    scale_color_manual(values = set_col) + 
    scale_shape_manual(values = set_shape) + 
    theme_bw() + 
    theme(legend.position = legend.pos, 
          axis.title.x = element_text(size = 25, family = fon, face = 'bold'), 
          axis.title.y = element_text(size = 25, family = fon, face = 'bold'), 
          axis.text.x = element_text(size = 22, family = fon, face = "bold", color = 'black'), 
          axis.text.y = element_text(size = 22, family = fon, face = "bold", color = 'black'), 
          axis.ticks = element_line(color = 'black', linewidth = 1), 
          axis.ticks.length = unit(2, 'mm'), 
          panel.border = element_rect(linewidth = 2), 
          panel.grid = element_blank(), 
          plot.title = element_text(size = 25, family = fon, face = 'bold', hjust = 0.5)) + 
    xlim(xlim[1], xlim[2]) + ylim(ylim[1], ylim[2]) + 
    xlab(NULL) + ylab(NULL)
  
  if (legend.pos[1] != 'none') {
    pMain + theme(legend.title = element_text(size = 18, family = fon, face = 'bold'), 
                  legend.text = element_text(size = 18, family = fon, face = 'bold')) + 
      labs(color = NULL, shape = NULL, title = NULL)
  }
  else {
    pMain + labs(title = NULL)
  }
}

# SVM with ensemble: for unbalanced group
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
  
  pred.res <- list(pred.all = pred.all, 
                   pred.med = pred)
  return(pred.res)
}

# CV for SVM ensemble
CVforSVMensemble <- function(dat, grp, fold = 7, Tn = 100) {
  grp.f <- as.factor(grp)
  grp.lev <- levels(grp.f)
  
  # scale
  dat <- scale(dat, center = T, scale = T)
  
  dat.g1 <- dat[grp == grp.lev[1], ]
  dat.g2 <- dat[grp == grp.lev[2], ]
  
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
    cv.grp <- c(rep(0, nrow(dat.cv.g1)), rep(1, nrow(dat.cv.g2)))
    
    dat.tr.g1 <- dat.g1[-c(start.index1: end.index1), ]
    dat.tr.g2 <- dat.g2[-c(start.index2: end.index2), ]
    dat.tr <- rbind(dat.tr.g1, dat.tr.g2)
    tr.grp <- c(rep(0, nrow(dat.tr.g1)), rep(1, nrow(dat.tr.g2)))
    tr.grp <- as.factor(tr.grp)
    
    Ntr.g1 <- nrow(dat.tr.g1)
    Ntr.g2 <- nrow(dat.tr.g2)
    
    # dat.tr.sca <- scale(dat.tr, center = T, scale = T)
    # dat.cv.sca <- scale(dat.cv, center = T, scale = T)
    
    if (max(Ntr.g1, Ntr.g2) / min(Ntr.g1, Ntr.g2) >= 1.5) {
      clas.pred.i <- SVMensemble(dat.tr, tr.grp, dat.cv, Tn = Tn)
      # pred.all.i <- clas.pred.i$pred.all
      clas.pred.i <- clas.pred.i$pred
    }
    else {
      clas.pred.i <- matrix(0, nrow(dat.cv), Tn)
      for (k in c(1: Tn)) {
        model.i <- svm(dat.tr, tr.grp, probability = T)
        pred.i <- predict(model.i, dat.cv, probability = T)
        pred.prob.i <- attr(pred.i, "probabilities")
        colIdx <- which(as.numeric(colnames(pred.prob.i)) == 1)
        clas.pred.i[, k] <- pred.prob.i[, colIdx]
      }
      clas.pred.i <- apply(clas.pred.i, MARGIN = 1, FUN = median)
    }
    
    clas.pred <- c(clas.pred, clas.pred.i)
    clas.GT <- c(clas.GT, cv.grp)
    
    # MSE
    cv.err[i] <- 1 / (nrow(dat.cv)) * sum((clas.pred.i - as.numeric(cv.grp))^2)
    accuracy[i] <- 1 - cv.err[i]
    
    cat(i, '\n')
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

