# Title: All Evaluation
# Author: Fanjing Guo
# Date: 2023.3.6

library(mixOmics)
library(car)
library(webshot)
library(statTarget)
library(MetNormalizer)
library(ropls)
library(pROC)
library(e1071)
library(ggbreak)

# -------------------------------------
# 1. Removal of Inter-batch variation
# -------------------------------------
RMbatVariation <- function(dat.bef, dat.bef.intraCor = NULL, dat.cor, 
                           batch, group, order, QC_label = 'QC', 
                           Fval.cutoff = 100, feature.list = NULL, 
                           save.path0) {
  bat.f <- as.factor(batch)
  batch.num <- nlevels(bat.f)
  
  # 1.1 PCA
  # Before correction
  pca.bef <- pca(dat.bef, ncomp = 2, center = TRUE, scale = TRUE)
  pca.bef.varX <- pca.bef$variates$X
  pca.bef.expl <- pca.bef$prop_expl_var$X
  
  p1 <- PCA_scatter(data = pca.bef.varX, batch = batch, 
              expl.var = pca.bef.expl,
              xlim = c(min(pca.bef.varX[, 1]) - 1, max(pca.bef.varX[, 1]) + 1), 
              ylim = c(min(pca.bef.varX[, 2]) - 1, max(pca.bef.varX[, 2]) + 1), 
              legend.pos = 'none')
  print(p1)
  ggsave(paste0(save.path0, 'PCA_bef_', gsub('-', '', Sys.Date()), '.emf'), plot = p1, 
         width = 7.5, height = 7, units = 'in', dpi = 300)
  
  # maha dist in PC space
  bef.md <- rep(0, batch.num * (batch.num - 1))
  l <- 0
  for (i in c(1: batch.num)) {
    PCdat.bati <- pca.bef.varX[batch == i, ]
    PCdat.bati.m <- colMeans(PCdat.bati)
    PCdat.bati.s <- cov(PCdat.bati)
    for (j in c(1: batch.num)) {
      if (j != i) {
        l <- l + 1
        PCdat.batj <- pca.bef.varX[batch == j, ]
        bef.md[l] <- mean(mahalanobis(PCdat.batj, PCdat.bati.m, PCdat.bati.s))
      }
    }
  }
  bef.md.mean <- mean(bef.md)
  
  # After correction
  pca.cor <- pca(dat.cor, ncomp = 2, center = TRUE, scale = TRUE)
  pca.cor.varX <- pca.cor$variates$X
  pca.cor.expl <- pca.cor$prop_expl_var$X
  
  p2 <- PCA_scatter(data = pca.cor.varX, batch = batch, 
              expl.var = pca.cor.expl,
              xlim = c(min(pca.cor.varX[, 1]) - 1, max(pca.cor.varX[, 1]) + 1), 
              ylim = c(min(pca.cor.varX[, 2]) - 1, max(pca.cor.varX[, 2]) + 1), 
              legend.pos = 'right')
  
  print(p2)
  ggsave(paste0(save.path0, 'PCA_cor_', gsub('-', '', Sys.Date()), '.emf'), plot = p2, 
         width = 8.5, height = 7, units = 'in', dpi = 300)
  
  # maha dist in PC space
  cor.md <- rep(0, batch.num * (batch.num - 1))
  l <- 0
  for (i in c(1: batch.num)) {
    PCdat.bati <- pca.cor.varX[batch == i, ]
    PCdat.bati.m <- colMeans(PCdat.bati)
    PCdat.bati.s <- cov(PCdat.bati)
    for (j in c(1: batch.num)) {
      if (j != i) {
        l <- l + 1
        PCdat.batj <- pca.cor.varX[batch == j, ]
        cor.md[l] <- mean(mahalanobis(PCdat.batj, PCdat.bati.m, PCdat.bati.s))
      }
    }
  }
  cor.md.mean <- mean(cor.md)
  
  # 1.2 ANOVA
  if (sum(group == QC_label) != 0) {
    dat.bef.sbj <- dat.bef[group != QC_label, ]
    dat.cor.sbj <- dat.cor[group != QC_label, ]
    group.sbj <- group[group != QC_label]
    batch.sbj <- batch[group != QC_label]
  }
  else {
    dat.bef.sbj <- dat.bef
    dat.cor.sbj <- dat.cor
    group.sbj <- group
    batch.sbj <- batch
  }
  
  p <- ncol(dat.bef.sbj)
  
  grp.sbj.f <- as.factor(group.sbj)
  bat.sbj.f <- as.factor(batch.sbj)
  
  bef.Sum_sq <- data.frame(BE = rep(0, p))
  bef.Fval <- data.frame(BE = rep(0, p))
  
  cor.Sum_sq <- data.frame(BE = rep(0, p))
  cor.Fval <- data.frame(BE = rep(0, p))
  
  for (i in c(1: p)) {
    # fit model before correction
    bef.mod.i <- lm(dat.bef.sbj[, i] ~ bat.sbj.f * grp.sbj.f)
    bef.aov.i <- Anova(bef.mod.i, type = 3)
    bef.Sum_sq[i, ] <- bef.aov.i$`Sum Sq`[2]
    bef.Fval[i, ] <- bef.aov.i$`F value`[2]
    
    # fit model after CordBat correction
    cor.mod.i <- lm(dat.cor.sbj[, i] ~ bat.sbj.f * grp.sbj.f)
    cor.aov.i <- Anova(cor.mod.i, type = 3)
    cor.Sum_sq[i, ] <- cor.aov.i$`Sum Sq`[2]
    cor.Fval[i, ] <- cor.aov.i$`F value`[2]
    
    # cat(i, '\n')
  }
  
  bef.BEcut <- bef.Fval$BE
  cor.BEcut <- cor.Fval$BE
  bef.BEmet.num <- sum(bef.BEcut > 1.7)
  cor.BEmet.num <- sum(cor.BEcut > 1.7)
  cor.BEcut <- cor.BEcut[order(bef.BEcut)]
  bef.BEcut <- sort(bef.BEcut)
  
  cutoff <- Fval.cutoff
  bef.BEcut[bef.BEcut > cutoff] <- cutoff
  cor.BEcut[cor.BEcut > cutoff] <- cutoff
  
  dat.BE <- data.frame(dat.name = factor(c(rep('Uncorrect', p), 
                                           rep('CordBat', p)), 
                                         levels = c('Uncorrect', 'CordBat')), 
                       met.id = rep(c(1: p), 2),  
                       Fval = c(bef.BEcut, cor.BEcut))
  
  gaps <- seq(10, 100, 10)
  abs_5 <- abs(p %/% gaps - 5)
  set_gap <- gaps[which(abs_5 == min(abs_5))[1]]
  
  p3 <- ggplot(data = dat.BE, aes(x = dat.name, y = met.id, fill = Fval)) + 
    geom_raster() + theme_bw() + 
    scale_fill_gradientn(colors = c('#39489f', '#39bbec', '#f9ed36', '#f38466', '#b81f25')) + 
    scale_y_continuous(breaks = c(p, seq((p - set_gap + 1), ((p %% set_gap) + 1), -set_gap)), 
                       labels = c(1, seq(set_gap, p, set_gap))) +
    scale_x_discrete(expand = expansion(mult = c(0, 0)), position = 'top') + 
    theme(axis.text.x.top = element_text(size = 25, family = 'sans', face = 'bold', 
                                         vjust = -3, color = 'black'), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 25, family = 'sans', face = 'bold', vjust = 2), 
          axis.text.y = element_text(size = 22, family = 'sans', face = 'bold', color = 'black'),
          axis.ticks = element_blank(),
          legend.position = 'right', legend.justification = c(1.3, 0.95), 
          legend.title = element_text(size = 20, family = 'sans', face = 'bold'), 
          legend.text = element_text(size = 15, family = 'sans', face = 'bold'), 
          legend.box.spacing = unit(0.3, 'cm'),
          legend.background = element_rect(fill = 'transparent'), 
          panel.border = element_blank(), 
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 25, 
                                    family = 'sans', face = 'bold')) + 
    labs(y = 'Index of sorted metabolite', fill = 'F-value', title = NULL)
  
  print(p3)
  ggsave(paste0(save.path0, 'ANOVA_', gsub('-', '', Sys.Date()), '.emf'), plot = p3, 
         width = 6.5, height = 7, units = 'in', dpi = 300)
  
  
  # 1.3 inter-batch DE metabolites 
  if (!is.null(dat.bef.intraCor)) {
    dat.bef <- dat.bef.intraCor
  }
  
  if (is.null(feature.list)) {
    feature.list <- c(1: p)
  }
  DEmetIdx <- DEana(dat.bef, batch, top.rank = c(1, 2))
  
  met.bef <- dat.bef[, DEmetIdx]
  colnames(met.bef) <- feature.list[DEmetIdx]
  
  met.cor <- dat.cor[, DEmetIdx]
  colnames(met.cor) <- feature.list[DEmetIdx]
  
  p4 <- DEmet_plot_merge(dat.bef = met.bef, 
                   dat.cor = met.cor, 
                   batch = batch, 
                   InjOrd = order)
  
  print(p4)
  ggsave(paste0(save.path0, 'BatDEmet_', gsub('-', '', Sys.Date()), '.emf'), plot = p4, 
         width = 7, height = 7, units = 'in', dpi = 300)
  
  rmBE.res <- data.frame(maha.bef = bef.md.mean,
                         maha.cor = cor.md.mean,
                         bef.BEmetNum = bef.BEmet.num,
                         cor.BEmetNum = cor.BEmet.num)
  write.table(rmBE.res, paste0(save.path0, 'result_', gsub('-', '', Sys.Date()), '.txt'))
  
}


# -----------------------------------
# 2. Preservation of data structure
# -----------------------------------
DelOut_Imp <- function(dat, batch, group = NULL) {
  
  bat.f <- as.factor(batch)
  batch.lev <- levels(bat.f)
  batch.num <- nlevels(bat.f)
  
  # delete outliers
  batch.new <- batch
  group.new <- group
  dat.delout <- dat
  delsampIdx <- c()
  
  for (i in c(1: batch.num)) {
    dat.bati <- dat[batch == batch.lev[i], ]
    bati.idx <- which(batch == batch.lev[i])
    dat.bati <- DelOutlier(dat.bati)
    delsamp.bati <- dat.bati$delsampIdx
    dat.bati <- dat.bati$X.out
    dat.bati <- ImputeOutlier(dat.bati)
    
    if (length(delsamp.bati) != 0) {
      bati.delinitIdx <- bati.idx[delsamp.bati]
      dat.delout[bati.idx[-delsamp.bati], ] <- dat.bati
      delsampIdx <- c(delsampIdx, bati.delinitIdx)
    }
    else {
      dat.delout[bati.idx, ] <- dat.bati
    }
  }
  
  if (length(delsampIdx) != 0) {
    dat.delout <- dat.delout[-delsampIdx, ]
    batch.new <- batch.new[-delsampIdx]
    group.new <- group.new[-delsampIdx]
  }
  
  delout.res <- list(dat.delout = dat.delout, 
                     batch.new = batch.new, 
                     group.new = group.new)
  return(delout.res)
}

PresDatStruct <- function(dat.bef, dat.bef.delout, dat.cor, batch.delout, 
                          NotdifBat.ratio = 0.8, delays = c(0.2, 0.5), save.path0) {
  # dat.bef.delout: the data before batch effect correction, in which the outliers 
  # in each batch are deleted to obtain more stable PCC standard set (since the PCC 
  # is sensitive to outliers).
  
  stdPcc <- getstdPcc(dat.bef.delout, batch.delout, Notdiff.ratio = NotdifBat.ratio)
  stdPcc.num <- stdPcc$retain.PCC.num
  stdPcc.perc <- stdPcc$retain.PCC.perc
  stdPcc <- stdPcc$std.PCC.Mat
  
  colnames(dat.bef) <- c(1: ncol(dat.bef))
  colnames(dat.cor) <- c(1: ncol(dat.cor))
  
  confs <- c(0.1, 0.05, 0.01)
  PresDS.res <- data.frame(stdPcc.num = rep(stdPcc.num, length(confs)), 
                           stdPcc.perc = rep(stdPcc.perc, length(confs)), 
                           conf_level = confs, 
                           bef.dPccNum = rep(0, length(confs)), 
                           bef.dPccPerc = rep(0, length(confs)), 
                           cor.dPccNum = rep(0, length(confs)), 
                           cor.dPccPerc = rep(0, length(confs)))
  
  for (i in c(1: length(confs))) {
    conf <- confs[i]
    
    # before correction
    P.bef <- getPdiff2Std(stdPcc, dat.bef, conf_alpha = conf)
    bef.dPCC.num <- P.bef$diff.pcc.num
    PresDS.res$bef.dPccNum[i] <- bef.dPCC.num
    bef.dPCC.perc <- P.bef$diff.ratio
    PresDS.res$bef.dPccPerc[i] <- bef.dPCC.perc
    P.bef <- P.bef$P.Mat
    
    bef_vs_std.largestComp <- getLargestComp(P.bef, conf_alpha = conf)
    bef_vs_std.Pcc <- bef_vs_std.largestComp$P.connectComp
    # p1 <- Graph_plot_forP(bef_vs_std.Pcc)
    # 
    file.name <- paste0(save.path0, 'bef_dPCCnet_', gsub('[.]', '_', conf))
    # saveNetwork(p1, paste0(file.name, '.html'))
    # webshot::webshot(paste0(file.name, '.html'), file = paste0(file.name, '_', 
    #                                                            gsub('-', '', Sys.Date()), '.png'),
    #                  delay = delays[1], vwidth = 500, vheight = 500)
    unlink(paste0(file.name, '.html'))
    unlink(paste0(file.name, '_files'), recursive = T)
    
    # print(p1)
    
    # after correction
    P.cor <- getPdiff2Std(stdPcc, dat.cor, conf_alpha = conf)
    cor.dPCC.num <- P.cor$diff.pcc.num
    PresDS.res$cor.dPccNum[i] <- cor.dPCC.num
    cor.dPCC.perc <- P.cor$diff.ratio
    PresDS.res$cor.dPccPerc[i] <- cor.dPCC.perc
    P.cor <- P.cor$P.Mat
    
    if (cor.dPCC.num != 0) {
      cor_vs_std.largestComp <- getLargestComp(P.cor, conf_alpha = conf)
      cor_vs_std.Pcc <- cor_vs_std.largestComp$P.connectComp
      p2 <- Graph_plot_forP(cor_vs_std.Pcc)

      file.name <- paste0(save.path0, 'cor_dPCCnet_', gsub('[.]', '_', conf))
      saveNetwork(p2, paste0(file.name, '.html'))
      webshot::webshot(paste0(file.name, '.html'), file = paste0(file.name, '_',
                                                                 gsub('-', '', Sys.Date()), '.png'),
                       delay = delays[2], vwidth = 350, vheight = 350)
      unlink(paste0(file.name, '.html'))
      unlink(paste0(file.name, '_files'), recursive = T)
      
      # print(p2)
    }
  }
  
  write.table(PresDS.res, file = paste0(save.path0, 'Data structure pres_', 
                                        gsub('-', '', Sys.Date()), '.txt'))
}

# --------------------------------------
# 3. Comparison with QC-based method
# --------------------------------------
############## QC-RLSC BEC #####################
QCRLSC_cor <- function(dat, samp.info, feature.list, QC_label = 'QC') {
  # generate pheno and data file
  samp.name <- samp.info$SampName
  samp.order <- samp.info$RunOrder
  batch <- samp.info$Batch
  group <- samp.info$Group
  
  # standardize sample name
  qc.num <- sum(group == QC_label)
  if (length(grep('QC', samp.name)) != qc.num) {
    if (length(grep('qc', samp.name, ignore.case = T)) == 0) {
      samp.name <- paste0(samp.name, '#QC')
    }
    else {
      samp.name <- gsub('qc', 'QC', samp.name, ignore.case = T)
    }
  }
  
  # standardize group label
  group[group == QC_label] <- NA
  grp.f <- as.factor(group[!is.na(group)])
  grp.num <- nlevels(grp.f)
  grp.lev <- levels(grp.f)
  std.grp <- group
  for (i in c(1: grp.num)) {
    std.grp[group == grp.lev[i]] <- i
  }
  group <- std.grp
  
  samPeno <- data.frame(sample = samp.name, 
                        batch = batch, 
                        class = group, 
                        order = samp.order)
  write.csv(samPeno, 'pheno.csv', row.names = F)
  samPeno <- 'pheno.csv'
  
  samFile <- data.frame(name = feature.list, 
                        t(dat))
  colnames(samFile)[-1] <- samp.name
  write.csv(samFile, 'file.csv', row.names = F)
  samFile <- 'file.csv'
  
  shiftCor(samPeno, samFile, MLmethod = "QCRLSC", QCspan = 0.75, 
           imputeM = "KNN", coCV = 100, plot = F)
  
  read.path <- paste0(getwd(), '/statTarget/shiftCor/After_shiftCor/shift_all_cor.csv')
  QRcor_Dat <- read.csv(read.path, header = TRUE)
  QRcor_Dat <- as.matrix(QRcor_Dat[, -c(1, 2)])
  
  unlink(paste0(getwd(), '/statTarget'), recursive = T)
  unlink(paste0(getwd(), '/pheno.csv'))
  unlink(paste0(getwd(), '/file.csv'))
  
  return(QRcor_Dat)
}

################ SERRF format data generating ########################
GenSRformatdat <- function(dat, samp.info, feature.list, QC_label = 'QC', 
                           wrfile.path) {
  # generate data file with SERRF format
  samp.name <- samp.info$SampName
  samp.order <- samp.info$RunOrder
  batch <- samp.info$Batch
  group <- samp.info$Group
  
  # standardize batch iD
  bat.f <- as.factor(batch)
  bat.num <- nlevels(bat.f)
  bat.lev <- levels(bat.f)
  for (i in c(1: bat.num)) {
    batch[batch == bat.lev[i]] <- LETTERS[i]
  }
  
  # standardize group label
  group[group != QC_label] <- 'sample'
  group[group == QC_label] <- 'qc'
  
  samp_dat <- data.frame(batch = batch, 
                         sampleType = group, 
                         time = samp.order, 
                         label = samp.name, 
                         dat)
  colnames(samp_dat)[-c(1: 4)] <- feature.list
  samp_dat <- t(samp_dat)
  samp_dat <- as.data.frame(samp_dat)

  p <- length(feature.list)
  feat.no <- c(1: p)
  all_dat <- data.frame(c1 = c(rep('', 3), 'No', feat.no), 
                        c2 = rownames(samp_dat), 
                        samp_dat)  
  
  # write to the excel file
  write.table(all_dat, wrfile.path, row.names = F, col.names = F, sep = ',')
  
}

###################### SVR BEC ##############################
SVR_cor <- function(dat, samp.info, mets.name = NULL, mzs.list = NULL, 
                    rts.list = NULL, QC_label = 'QC') {
  
  # generate pheno and data file
  samp.name <- samp.info$SampName
  samp.order <- samp.info$RunOrder
  group <- samp.info$Group
  
  # standardize group labels
  group[group != QC_label] <- 'Subject'
  group[group == QC_label] <- 'QC'
  
  SampInfo <- data.frame(sample.name = samp.name, 
                         injection.order = samp.order, 
                         class = group)
  write.csv(SampInfo, 'sample.info.csv', row.names = F)
  SampInfo <- 'sample.info.csv'
  
  p <- ncol(dat)
  if (is.null(mets.name)) {
    mets.name <- rep('', p)
    for (i in c(1: p)) {
      mets.name[i] <- paste0('Unknown-', i)
    }
  }
  if (is.null(mzs.list)) {
    mzs.list <- rep(0, p)
  }
  else {
    mzs.list <- as.numeric(gsub('[^0-9,.]', '', mzs.list))
  }
  if (is.null(rts.list)) {
    rts.list <- rep(0, p)
  }
  
  SampDat <- data.frame(name = mets.name, 
                        mz = mzs.list, 
                        rt = rts.list, 
                        t(dat))
  colnames(SampDat)[-c(1: 3)] <- samp.name
  write.csv(SampDat, 'data.csv', row.names = F)
  SampDat <- 'data.csv'
  
  metNor(ms1.data.name = SampDat, sample.info.name = SampInfo)
  
  unlink(SampInfo)
  unlink(SampDat)
  
  read.path <- paste0(getwd(), '/svr_normalization_result/data_svr_normalization.csv')
  SVRcor.res <- read.csv(read.path, header = T)
  SVRcor1 <- as.matrix(t(SVRcor.res[, -c(1: 6)]))
  N <- nrow(SVRcor1)
  p <- ncol(SVRcor1)
  N.sbj <- sum(group == 'Subject')
  SVRcor_Dat <- matrix(0, N, p)
  SVRcor_Dat[group == 'Subject', ] <- SVRcor1[c(1: N.sbj), ]
  SVRcor_Dat[group == 'QC', ] <- SVRcor1[c((N.sbj + 1): N), ]
  
  unlink(paste0(getwd(), '/svr_normalization_result'), recursive = T)
  
  return(SVRcor_Dat)
}

################## Evaluation ###################
# PCA considering whether the QCs are pooled samples
PCA_1 <- function(dat, group, QC_label = 'QC', pooledQC = T) {
  if (pooledQC) {
    pca.list <- pca(dat, ncomp = 2, center = T, scale = T)
    pca.varX <- pca.list$variates$X
    pca.expl <- pca.list$prop_expl_var$X
  }
  else {
    dat.scale <- scale(dat, center = T, scale = T)
    
    dat.sbj <- dat.scale[group != QC_label, ]
    dat.qc <- dat.scale[group == QC_label, ]
    
    pca.sbj <- pca(dat.sbj, ncomp = 2, center = T, scale = T)
    pca.expl <- pca.sbj$prop_expl_var$X
    pca.sbj.varx <- pca.sbj$variates$X
    pca.sbj.load <- pca.sbj$loadings$X
    pca.qc.varx <- dat.qc %*% pca.sbj.load
    pca.varX <- matrix(0, nrow(dat), 2)
    pca.varX[group != QC_label, ] <- pca.sbj.varx
    pca.varX[group == QC_label, ] <- pca.qc.varx
  }
  
  pca.res <- list(varX = pca.varX, 
                  expl = pca.expl)
  
  return(pca.res)
}

# get mahalanobis distance between batches
MahaDist_bat <- function(pca.varX, batch) {
  
  bat.f <- as.factor(batch)
  batch.num <- nlevels(bat.f)
  batch.lev <- levels(bat.f)
  maha.dist <- c()
  for (i in c(1: batch.num)) {
    PCdat.bati <- pca.varX[batch == batch.lev[i], ]
    PCdat.bati.m <- colMeans(PCdat.bati)
    PCdat.bati.s <- cov(PCdat.bati)
    for (j in c(1: batch.num)) {
      if (j != i) {
        PCdat.batj <- pca.varX[batch == batch.lev[j], ]
        maha.dist <- c(maha.dist, mahalanobis(PCdat.batj, PCdat.bati.m, PCdat.bati.s))
      }
    }
  }
  
  return(maha.dist)
}

# add noise to data with ratio of metabolite intensity
addnoise_ratio <- function(dat, noise.ratio, dolog = F) {
  N <- nrow(dat)
  p <- ncol(dat)
  
  if (dolog) {
    dat.log <- log2(dat)
  }
  else {
    dat.log <- dat
  }
  
  # generate noise
  dat.m <- apply(dat.log, MARGIN = 2, FUN = mean)
  Noise <- matrix(0, N, p)
  for (i in c(1: p)) {
    Noise[, i] <- rnorm(N, mean = 0, sd = noise.ratio * dat.m[i])
  }
  
  dat.addnoise <- dat.log + Noise
  
  addnoise.result <- list(dat = dat.addnoise, 
                          addnoise = Noise)
  return(addnoise.result)
}

# get cumulative percentages of metabolites with RSDs
getCumPerc_RSD <- function(dat, inv.log = T, 
                           rsd.val = seq(0.001, 1, 0.001)) {
  
  if (inv.log) {
    dat <- 2^dat
  }
  
  met.m <- apply(dat, MARGIN = 2, FUN = mean)
  met.sd <- apply(dat, MARGIN = 2, FUN = sd)
  met.rsd <- met.sd / met.m
  
  cumPerc <- rep(0, length(rsd.val))
  for (i in c(1: length(rsd.val))) {
    cumPerc[i] <- sum(met.rsd <= rsd.val[i]) / length(met.rsd)
  }
  
  return(cumPerc)
}

GenDatNoise_SRformat <- function(dat, dat.name, noise.ratio = c(0.03, 0.045, 0.06), samp.info, 
                                 feature.list, QC_label = 'QC', dolog = F){
  dir.name <- paste0('SRformat_file_addNoise_', dat.name)
  unlink(dir.name, recursive = T)
  dir.create(dir.name)
  
  group <- samp.info$Group
  dat.qc <- dat[group == QC_label, ]
  dat.addNoise <- dat
  for (i in c(1: length(noise.ratio))) {
    
    addNoise <- addnoise_ratio(dat.qc, noise.ratio[i], dolog = dolog)
    NoiDat.i <- addNoise$dat
    dat.addNoise[group == QC_label, ] <- NoiDat.i
    
    GenSRformatdat(dat.addNoise, samp.info, feature.list, QC_label = QC_label, 
                   wrfile.path = paste0(getwd(), '/', dir.name, 
                                        '/SRformat_Dat_addNoise_', i, '.csv'))
  }
  
  cat('Files in directory: ', dir.name, '\n')
}

# main function - Compare with QC-based methods (RLSC, SVR, SERRF)
CmpWithQCmethod <- function(dat.list, datNoi.SRcor.list, batch, group, 
                            target.method.name = 'CordBat', 
                            QC_label = 'QC', pooledQC = T, 
                            PCA.print = T, MahaDist.print = T, 
                            CumRSD.sbj.print = T, CumRSD.qc.print = T, 
                            mark.rsd.sbj = 20, mark.rsd.qc = 20, 
                            set_arrow.sbj = NULL, set_arrow.qc = NULL, 
                            save.path0) {
  
  dat.name <- names(dat.list)
  dat.num <- length(dat.list)
  
  # Uncor_Dat <- dat.list$Uncorrect # [[1]]
  # QRcor_Dat <- dat.list$RLSC # [[2]]
  # SVRcor_Dat <- dat.list$SVR # [[3]]
  # SRcor_Dat <- dat.list$SERRF # [[4]]
  # CdBcor_Dat <- dat.list$CordBat # [[5]]
  
  # 3.1 PCA with QCs
  if (any(c(PCA.print, MahaDist.print))) {
    pca.varX.list <- list()
    for (i in c(1: dat.num)) {
      dat.i <- dat.list[[i]]
      pca.i.list <- PCA_1(dat.i, group = group, QC_label = QC_label, pooledQC = pooledQC)
      pca.i.varX <- pca.i.list$varX
      pca.varX.list[[i]] <- pca.i.varX
      pca.i.expl <- pca.i.list$expl
      
      if (PCA.print) {
        if (i != 3) {
          legend.pos <- 'none'
          set_width <- 7.5
        }
        else {
          legend.pos <- 'right'
          set_width <- 8.5
        }
        
        p1 <- PCA_scatter_withQC(data = pca.i.varX, batch = batch, QC_label = QC_label, 
                                 group = group, expl.var = pca.i.expl,
                                 xlim = c(min(pca.i.varX[, 1]) - 1, max(pca.i.varX[, 1]) + 1), 
                                 ylim = c(min(pca.i.varX[, 2]) - 1, max(pca.i.varX[, 2]) + 1), 
                                 legend.pos = legend.pos, title = dat.name[i])
        print(p1)
        ggsave(paste0(save.path0, 'PCA_', dat.name[i], '_', gsub('-', '', Sys.Date()), '.emf'), 
               plot = p1, width = set_width, height = 7.8, units = 'in', dpi = 300)
      }
    }
  }
  
  # 3.2 Mahalanobis distance between batches for study samples
  if (MahaDist.print) {
    batch.sbj <- batch[group != QC_label]
    mds.med <- rep(0, dat.num)
    maha_dist.list <- list()
    for (i in c(1: dat.num)) {
      pca.varX.i <- pca.varX.list[[i]]
      pca.sbj.varX.i <- pca.varX.i[group != QC_label, ]
      maha_dist.list[[i]] <- MahaDist_bat(pca.sbj.varX.i, batch.sbj)
      mds.med[i] <- median(maha_dist.list[[i]])
    }
    
    # plot
    dattype.vec <- unlist(lapply(dat.name, function(x) rep(x, length(maha_dist.list[[1]]))))
    dattype <- factor(dattype.vec, levels = dat.name)
    
    dat.md <- data.frame(dattype = dattype, 
                         md = unlist(maha_dist.list))
    
    max.mds_75 <- max(sapply(maha_dist.list, FUN = function(x) quantile(x, 0.75)))
    min.mds_75 <- min(sapply(maha_dist.list, FUN = function(x) quantile(x, 0.75)))
    
    if (max.mds_75 / min.mds_75 > 30) {
      f <- 1
    }
    else {
      f <- 2.3
    }
    
    fon <- 'sans'
    set_col <- c('#F39B7FFF', '#3C5488FF', '#00A087FF', 
                 '#4DBBD5FF', '#E64B35FF')
    p2 <- ggplot(data = dat.md, aes(x = dattype, y = md)) + 
      geom_boxplot(aes(fill = dattype), color = 'black', size = 1, outlier.color = NA) + 
      scale_fill_manual(values = set_col) + theme_bw() + 
      xlab(NULL) + ylab(NULL) + 
      theme(axis.text.x = element_text(size = 25, family = fon, face = 'bold', 
                                       hjust = 0.7, vjust = 0.85, angle = 30, color = 'black'), 
            axis.text.y = element_text(size = 25, family = fon, face = 'bold', 
                                       color = 'black'), 
            axis.ticks = element_line(color = 'black', linewidth = 1, lineend = 20), 
            legend.position = 'none', 
            plot.title = element_text(hjust = 0.5, size = 30, family = fon, face = 'bold'), 
            panel.border = element_rect(linewidth = 3), 
            panel.grid = element_blank()) + 
      coord_cartesian(ylim = c(0, f * max.mds_75)) + 
      labs(title = 'Mahalanobis distance')
    
    print(p2)
    ggsave(paste0(save.path0, 'MahaDist_all_', gsub('-', '', Sys.Date()), '.emf'), 
           plot = p2, width = 7.5, height = 8.5, units = 'in', dpi = 300)
    
    Md.res <- data.frame(dat = dat.name,
                         mahadist = mds.med)
    write.table(Md.res, paste0(save.path0, 'MD result_', gsub('-', '', Sys.Date()), '.txt'))
  }
  

  # 3.3 Add noise to QC samples (CordBat vs SERRF)
  if (any(c(CumRSD.sbj.print, CumRSD.qc.print))) {
    CdBcor_Dat <- dat.list[[target.method.name]]
    SRcor_Dat <- dat.list$SERRF
    SRcor_NoiDat_03 <- SRcor_NoiDat.list[[1]]
    SRcor_NoiDat_045 <- SRcor_NoiDat.list[[2]]
    SRcor_NoiDat_06 <- SRcor_NoiDat.list[[3]]
    
    # get cumulative percentages with RSDs
    rsd.val <- seq(0.001, 1, 0.001)
    
    if (CumRSD.sbj.print) {
      # a. study samples
      CdB.sbj.cumPerc <- getCumPerc_RSD(CdBcor_Dat[group != QC_label, ], 
                                        rsd.val = rsd.val)
      SR.sbj.cumPerc <- getCumPerc_RSD(SRcor_Dat[group != QC_label, ], 
                                       rsd.val = rsd.val)
      SR03.sbj.cumPerc <- getCumPerc_RSD(SRcor_NoiDat_03[group != QC_label, ], 
                                         rsd.val = rsd.val)
      SR045.sbj.cumPerc <- getCumPerc_RSD(SRcor_NoiDat_045[group != QC_label, ], 
                                          rsd.val = rsd.val)
      SR06.sbj.cumPerc <- getCumPerc_RSD(SRcor_NoiDat_06[group != QC_label, ], 
                                         rsd.val = rsd.val)
      
      sbj.cumPerc.list <- list(CordBat = CdB.sbj.cumPerc, 
                               SERRF = SR.sbj.cumPerc, 
                               SERRF_03 = SR03.sbj.cumPerc, 
                               SERRF_045 = SR045.sbj.cumPerc, 
                               SERRF_06 = SR06.sbj.cumPerc)
      
      if (is.null(set_arrow.sbj)) {
        arrow.pos <- c('left', 'topleft', 'bottomright')
        arrow.ang <- c(0, pi/20, pi/2.4)
        arrow.len <- c(7, 9, 8)
      }
      else {
        arrow.pos <- set_arrow.sbj$arrow.pos
        arrow.ang <- set_arrow.sbj$arrow.ang
        arrow.len <- set_arrow.sbj$arrow.len
      }
      
      p3 <- plotCumRSD_m(sbj.cumPerc.list, arrow.pos = arrow.pos, 
                         mark.rsd = mark.rsd.sbj, 
                         arrow.angle = arrow.ang, arrow.len = arrow.len)
      print(p3)
      ggsave(paste0(save.path0, 'CumRSD_sbj_', gsub('-', '', Sys.Date()), '.emf'), 
             plot = p3, width = 7.5, height = 7, units = 'in', dpi = 300)
    }
    
    if (CumRSD.qc.print) {
      # b. QC samples
      CdB.qc.cumPerc <- getCumPerc_RSD(CdBcor_Dat[group == QC_label, ], 
                                       rsd.val = rsd.val)
      SR.qc.cumPerc <- getCumPerc_RSD(SRcor_Dat[group == QC_label, ], 
                                      rsd.val = rsd.val)
      SR03.qc.cumPerc <- getCumPerc_RSD(SRcor_NoiDat_03[group == QC_label, ], 
                                        rsd.val = rsd.val)
      SR045.qc.cumPerc <- getCumPerc_RSD(SRcor_NoiDat_045[group == QC_label, ], 
                                         rsd.val = rsd.val)
      SR06.qc.cumPerc <- getCumPerc_RSD(SRcor_NoiDat_06[group == QC_label, ], 
                                        rsd.val = rsd.val)
      
      qc.cumPerc.list <- list(CordBat = CdB.qc.cumPerc, 
                              SERRF = SR.qc.cumPerc, 
                              SERRF_03 = SR03.qc.cumPerc, 
                              SERRF_045 = SR045.qc.cumPerc, 
                              SERRF_06 = SR06.qc.cumPerc)
      
      if (is.null(set_arrow.qc)) {
        arrow.pos <- c('bottomleft', 'left', 'bottomright')
        arrow.ang <- c(pi/2.5, 0, pi/2.4)
        arrow.len <- c(10, 8, 10)
      }
      else {
        arrow.pos <- set_arrow.qc$arrow.pos
        arrow.ang <- set_arrow.qc$arrow.ang
        arrow.len <- set_arrow.qc$arrow.len
      }
      
      p4 <- plotCumRSD_m(qc.cumPerc.list, arrow.pos = arrow.pos, 
                         mark.rsd = mark.rsd.qc, 
                         arrow.angle = arrow.ang, arrow.len = arrow.len)
      print(p4)
      ggsave(paste0(save.path0, 'CumRSD_qc_', gsub('-', '', Sys.Date()), '.emf'), 
             plot = p4, width = 7.5, height = 7, units = 'in', dpi = 300)
    }
  }
}

# --------------------------------------
# 4. Preservation of biological effects
# --------------------------------------
# get P-value for unbalanced groups
getPval_unbalGrp <- function(dat.grp1, dat.grp2, Tn = 100, BH.adj = T) {
  N.1 <- nrow(dat.grp1)
  N.2 <- nrow(dat.grp2)
  p <- ncol(dat.grp1)
  
  if (N.1 != N.2) {
    ttest.pval.mat <- matrix(0, Tn, p)
    if (N.1 < N.2) {
      dat.minorGrp <- dat.grp1
      dat.majorGrp <- dat.grp2
    }
    else {
      dat.minorGrp <- dat.grp2
      dat.majorGrp <- dat.grp1
    }
    
    N.minor <- min(N.1, N.2)
    N.major <- max(N.1, N.2)
    
    for (i in c(1: Tn)) {
      dat.sampMajor <- dat.majorGrp[sample(c(1: N.major), N.minor), ]
      ttest.result <- getttestResult(dat.sampMajor, dat.minorGrp, BH.adj = BH.adj)
      ttest.pval.mat[i, ] <- ttest.result$tt.Pval
    }
    
    ttest.pval <- apply(ttest.pval.mat, MARGIN = 2, FUN = median)
  }
  else {
    ttest.result <- getttestResult(dat.grp1, dat.grp2, BH.adj = BH.adj)
    ttest.pval <- ttest.result$tt.Pval
  }
  
  return(ttest.pval)
}

# Determine potential biomarkers in each batch
PotenBiomarker <- function(dat.sbj, batch.sbj, group.sbj, 
                           sel.type = 1, 
                           alp = 0.1, BH.adj = T, 
                           toprank = NA, min.batch.num = 2) {
  
  p <- ncol(dat.sbj)
  
  bat.f <- as.factor(batch.sbj)
  batch.num <- nlevels(bat.f)
  batch.lev <- levels(bat.f)
  
  grp.f <- as.factor(group.sbj)
  group.num <- nlevels(grp.f)
  group.lev <- levels(grp.f)
  if (group.num > 2) {
    cat('Use only ', group.lev[1], ' and ', group.lev[2], ' group for t-test.')
  }
  
  DEmet_grp.VIP <- matrix(0, p, batch.num)
  DEmet_grp.Pval <- matrix(0, p, batch.num)
  
  for (i in c(1: batch.num)) {
    dat.i <- dat.sbj[batch.sbj == batch.lev[i], ]
    grp.i <- group.sbj[batch.sbj == batch.lev[i]]
    plsda.i <- opls(dat.i, grp.i, orthoI = 0, 
                    predI = 3, fig.pdfC = 'none', info.txtC = 'none')
    bati.VIP <- getVipVn(plsda.i)
    DEmet_grp.VIP[, i] <- bati.VIP
    
    N1.i <- sum(grp.i == group.lev[1])
    N2.i <- sum(grp.i == group.lev[2])
    
    if (max(N1.i, N2.i) / min(N1.i, N2.i) > 1.5) {
      ttest.pval <- getPval_unbalGrp(dat.i[grp.i == group.lev[1], ], 
                                     dat.i[grp.i == group.lev[2], ], 
                                     BH.adj = BH.adj)
    }
    else {
      ttest.result <- getttestResult(dat.i[grp.i == group.lev[1], ],
                                     dat.i[grp.i == group.lev[2], ], 
                                     BH.adj = BH.adj)
      ttest.pval <- ttest.result$tt.Pval
    }

    DEmet_grp.Pval[, i] <- ttest.pval
  }
  
  if (is.na(toprank)) {
    toprank <- round(0.2 * p)
  }
  
  DEmet_eachBat <- matrix(0, p, batch.num)
  
  if (sel.type == 1) {
    # firstly select features p-value < 0.1, secondly select VIP top rank
    for (i in c(1: batch.num)) {
      vi <- DEmet_grp.VIP[, i]
      pi <- DEmet_grp.Pval[, i]
      
      pi_1.num <- sum(pi < alp)
      pi_1.idx <- which(pi < alp)
      vi_1 <- vi[pi_1.idx]
      
      mi <- rep(0, p)
      if (pi_1.num > toprank) {
        vi.sortidx <- pi_1.idx[order(vi_1, decreasing = T)]
        vi.topidx <- vi.sortidx[c(1: toprank)]
        mi[vi.topidx] <- 1
      }
      else {
        mi[pi < alp] <- 1
      }
      
      DEmet_eachBat[, i] <- mi
    }
  }
  
  if (sel.type == 2) {
    # firstly select features VIP > 1, secondly select p-value top rank
    for (i in c(1: batch.num)) {
      vi <- DEmet_grp.VIP[, i]
      pi <- DEmet_grp.Pval[, i]
      
      vi_1.num <- sum(vi > 1)
      vi_1.idx <- which(vi > 1)
      pi_1 <- pi[vi_1.idx]
      
      mi <- rep(0, p)
      if (vi_1.num > toprank) {
        pi.sortidx <- vi_1.idx[order(pi_1, decreasing = F)]
        pi.topidx <- pi.sortidx[c(1: toprank)]
        mi[pi.topidx] <- 1
      }
      else {
        mi[vi > 1] <- 1
      }
      
      DEmet_eachBat[, i] <- mi
    }
  }

  
  DEmet_lt2Bat.Idx <- which(rowSums(DEmet_eachBat) >= min.batch.num)
  
  return(DEmet_lt2Bat.Idx)
}

# Identifying biomarkers in whole data
identBiomarker <- function(dat.sbj, group.sbj, potenBiomarks = NULL, 
                           rankMet.method = 'both', 
                           rankWeight = c(0.8, 0.2), 
                           BH.adj = T, max.topRank = NULL) {
  
  p <- ncol(dat.sbj)

  grp.f <- as.factor(group.sbj)
  group.num <- nlevels(grp.f)
  group.lev <- levels(grp.f)
  if (group.num > 2) {
    cat('Use only ', group.lev[1], ' and ', group.lev[2], ' group for t-test.')
  }
  
  # VIP
  dat.plsda <- opls(dat.sbj, group.sbj, orthoI = 0, 
                    predI = 3, fig.pdfC = 'none', info.txtC = 'none')
  met.VIP <- getVipVn(dat.plsda)
  
  # t-test
  N.1 <- sum(group.sbj == group.lev[1])
  N.2 <- sum(group.sbj == group.lev[2])
  if (max(N.1, N.2) / min(N.1, N.2) >= 1.5) {
    ttest.pval <- getPval_unbalGrp(dat.sbj[group.sbj == group.lev[1], ], 
                                   dat.sbj[group.sbj == group.lev[2], ], 
                                   BH.adj = BH.adj)
  }
  else {
    ttest.result <- getttestResult(dat.sbj[group.sbj == group.lev[1], ],
                                   dat.sbj[group.sbj == group.lev[2], ], 
                                   BH.adj = BH.adj)
    ttest.pval <- ttest.result$tt.Pval
  }
  
  # get overlapped feature number vs potential biomarkers
  if (rankMet.method == 'VIP') {
    sort.MetIdx <- order(met.VIP, decreasing = T)
  }
  
  if (rankMet.method == 'ttest') {
    sort.MetIdx <- order(ttest.pval, decreasing = F)
  }
  
  if (rankMet.method == 'both') {
    VIP.sort.MetIdx <- order(met.VIP, decreasing = T)
    VIP.metRank <- order(VIP.sort.MetIdx)
    
    Pval.sort.MetIdx <- order(ttest.pval, decreasing = F)
    Pval.metRank <- order(Pval.sort.MetIdx)
    
    met.mRank <- rankWeight[1] * VIP.metRank + rankWeight[2] * Pval.metRank
    sort.MetIdx <- order(met.mRank)
  }

  if (!is.null(max.topRank)) {
    if (max.topRank > p) {
      max.topRank <- p
    }
  }

  if (!is.null(potenBiomarks)) {
    OL.DEmet.res <- getOLperc(potenBiomarks, sort.MetIdx, max.topRank = max.topRank)
  }
  else {
    OL.DEmet.res <- NULL
  }
  
  identBioM.res <- list(sort.Met = sort.MetIdx, 
                        OL.DEmet = OL.DEmet.res)
  
  return(identBioM.res)
}

PotenCorrBioM <- function(dat.sbj, batch.sbj, group.sbj, 
                          BH.adj = T, toprank = NA, min.batch.num = 2) {
  p <- ncol(dat.sbj)
  
  bat.f <- as.factor(batch.sbj)
  batch.num <- nlevels(bat.f)
  batch.lev <- levels(bat.f)
  
  grp.f <- as.factor(group.sbj)
  group.num <- nlevels(grp.f)
  group.lev <- levels(grp.f)
  if (group.num > 2) {
    cat('Use only ', group.lev[1], ' and ', group.lev[2], ' group.')
  }
  
  pcc.num <- p * (p - 1) / 2
  if (is.na(toprank)) {
    toprank <- round(0.1 * pcc.num)
  }
  
  DEpcc_eachBat <- matrix(0, pcc.num, batch.num)
  for (i in c(1: batch.num)) {
    dat.i <- dat.sbj[batch.sbj == batch.lev[i], ]
    grp.i <- group.sbj[batch.sbj == batch.lev[i]]
    dat.i.g1 <- dat.i[grp.i == group.lev[1], ]
    dat.i.g2 <- dat.i[grp.i == group.lev[2], ]
    
    dat.delout <- DelOutlier(dat.i.g1)
    dat.i.g1.delout <- dat.delout$X.out
    dat.i.g1.delout <- ImputeOutlier(dat.i.g1.delout)
    
    dat.delout <- DelOutlier(dat.i.g2)
    dat.i.g2.delout <- dat.delout$X.out
    dat.i.g2.delout <- ImputeOutlier(dat.i.g2.delout)
    
    diff.res <- getdiffPccMat(dat.i.g1.delout, dat.i.g2.delout)
    diff.Z <- diff.res$diffMat
    diff.P <- 2 * pnorm(diff.Z, lower.tail = F)
    
    if (BH.adj) {
      diff.P <- matrix(p.adjust(diff.P, method = "BH"), nrow(diff.P), ncol(diff.P))
    }
    
    diff.Pu <- diff.P[upper.tri(diff.P)]
    sort.PccIdx <- order(diff.Pu)
    Pcc.topIdx <- sort.PccIdx[c(1: toprank)]
    
    DEpcc_eachBat[Pcc.topIdx, i] <- 1
  }
  
  DEpcc_lt2Bat.Idx <- which(rowSums(DEpcc_eachBat) >= min.batch.num)
  
  return(DEpcc_lt2Bat.Idx)
}

identCorrBioM <- function(dat.sbj, group.sbj, PotenCorrBioMarks, 
                          BH.adj = T, max.topRank) {
  p <- ncol(dat.sbj)
  
  grp.f <- as.factor(group.sbj)
  group.num <- nlevels(grp.f)
  group.lev <- levels(grp.f)
  if (group.num > 2) {
    cat('Use only ', group.lev[1], ' and ', group.lev[2], ' group.')
  }
  
  dat.g1 <- dat.sbj[group.sbj == group.lev[1], ]
  dat.g2 <- dat.sbj[group.sbj == group.lev[2], ]
  
  dat.delout <- DelOutlier(dat.g1)
  dat.g1.delout <- dat.delout$X.out
  dat.g1.delout <- ImputeOutlier(dat.g1.delout)
  
  dat.delout <- DelOutlier(dat.g2)
  dat.g2.delout <- dat.delout$X.out
  dat.g2.delout <- ImputeOutlier(dat.g2.delout)
  
  diff.res <- getdiffPccMat(dat.g1.delout, dat.g2.delout)
  diff.Z <- diff.res$diffMat
  diff.P <- 2 * pnorm(diff.Z, lower.tail = F)
  
  if (BH.adj) {
    diff.P <- matrix(p.adjust(diff.P, method = "BH"), nrow(diff.P), ncol(diff.P))
  }
  
  diff.Pu <- diff.P[upper.tri(diff.P)]
  sort.PccIdx <- order(diff.Pu)
  
  pcc.num <- p * (p - 1) / 2
  
  if (max.topRank > pcc.num) {
    max.topRank <- pcc.num
  }
  
  if (!is.null(PotenCorrBioMarks)) {
    OL.DEpcc.res <- getOLperc(PotenCorrBioMarks, sort.PccIdx, 
                              max.topRank = max.topRank)
  }
  else {
    OL.DEpcc.res <- NULL
  }
  
  identBioM.res <- list(sort.Pcc = sort.PccIdx, 
                        OL.DEpcc = OL.DEpcc.res)
  
  return(identBioM.res)

}

PotenCorrBioM_1 <- function(dat.sbj, batch.sbj, group.sbj, 
                            BH.adj = T, sel.type = 1, alp = 0.1, 
                            VIP.thresh = 1, 
                            toprank = NA, min.batch.num) {
  p <- ncol(dat.sbj)
  
  bat.f <- as.factor(batch.sbj)
  batch.num <- nlevels(bat.f)
  batch.lev <- levels(bat.f)
  
  grp.f <- as.factor(group.sbj)
  group.num <- nlevels(grp.f)
  group.lev <- levels(grp.f)
  if (group.num > 2) {
    cat('Use only ', group.lev[1], ' and ', group.lev[2], ' group.')
  }
  
  pcc.num <- p * (p - 1) / 2
  DEpcc_grp.VIP <- matrix(0, pcc.num, batch.num)
  DEpcc_grp.Pval <- matrix(0, pcc.num, batch.num)
  
  cat('Get potential correlation biomarkers in batch: \n')
  for (i in c(1: batch.num)) {
    dat.i <- dat.sbj[batch.sbj == batch.lev[i], ]
    grp.i <- group.sbj[batch.sbj == batch.lev[i]]
    N.i <- nrow(dat.i)
    
    samp.disturb.i <- matrix(0, N.i, pcc.num)
    ref.samp.i <- dat.i[grp.i == group.lev[1], ]
    ref.pcc.i <- cor(ref.samp.i)
    
    for (k in c(1: N.i)) {
      dpcc.k <- cor(rbind(ref.samp.i, dat.i[k, ])) - ref.pcc.i
      samp.disturb.i[k, ] <- dpcc.k[upper.tri(dpcc.k)]
    }
    
    samp.disturb.i <- samp.disturb.i + 2 # make all values positive
    
    if (is.numeric(sel.type)) {
      # VIP
      plsda.i <- opls(samp.disturb.i, grp.i, orthoI = 0, 
                      predI = 3, fig.pdfC = 'none', info.txtC = 'none')
      bati.VIP <- getVipVn(plsda.i)
      DEpcc_grp.VIP[, i] <- bati.VIP
      
      # t-test
      samp.disturb.i.g1 <- samp.disturb.i[grp.i == group.lev[1], ]
      samp.disturb.i.g2 <- samp.disturb.i[grp.i == group.lev[2], ]
      N1.i <- nrow(samp.disturb.i.g1)
      N2.i <- nrow(samp.disturb.i.g2)
      
      if (max(N1.i, N2.i) / min(N1.i, N2.i) > 1.5) {
        ttest.pval <- getPval_unbalGrp(samp.disturb.i.g1, samp.disturb.i.g2, 
                                       BH.adj = BH.adj)
      }
      else {
        ttest.result <- getttestResult(samp.disturb.i.g1, samp.disturb.i.g2, 
                                       BH.adj = BH.adj)
        ttest.pval <- ttest.result$tt.Pval
      }
      
      DEpcc_grp.Pval[, i] <- ttest.pval
    }
    else if (sel.type == 'VIP') {
      # VIP
      plsda.i <- opls(samp.disturb.i, grp.i, orthoI = 0, 
                      predI = 3, fig.pdfC = 'none', info.txtC = 'none')
      bati.VIP <- getVipVn(plsda.i)
      DEpcc_grp.VIP[, i] <- bati.VIP
    }
    else if (sel.type == 'ttest') {
      # t-test
      samp.disturb.i.g1 <- samp.disturb.i[grp.i == group.lev[1], ]
      samp.disturb.i.g2 <- samp.disturb.i[grp.i == group.lev[2], ]
      N1.i <- nrow(samp.disturb.i.g1)
      N2.i <- nrow(samp.disturb.i.g2)
      
      if (max(N1.i, N2.i) / min(N1.i, N2.i) > 1.5) {
        ttest.pval <- getPval_unbalGrp(samp.disturb.i.g1, samp.disturb.i.g2, 
                                       BH.adj = BH.adj)
      }
      else {
        ttest.result <- getttestResult(samp.disturb.i.g1, samp.disturb.i.g2, 
                                       BH.adj = BH.adj)
        ttest.pval <- ttest.result$tt.Pval
      }
      
      DEpcc_grp.Pval[, i] <- ttest.pval
    }
    else {
      stop('Find no method to determine potential biomarkers!')
    }
    
    cat(i, '\n')
  }
  
  DEpcc_eachBat <- matrix(0, pcc.num, batch.num)
  
  if (is.na(toprank)) {
    toprank <- round(0.05 * pcc.num)
  }
  
  if (sel.type == 'VIP') {
    for (i in c(1: batch.num)) {
      mi <- rep(0, pcc.num)
      
      vi <- DEpcc_grp.VIP[, i]
      vi.sortidx <- order(vi, decreasing = T)
      vi.topidx <- vi.sortidx[c(1: toprank)]
      mi[vi.topidx] <- 1
      DEpcc_eachBat[, i] <- mi
    }
  }
  
  if (sel.type == 'ttest') {
    for (i in c(1: batch.num)) {
      mi <- rep(0, pcc.num)
      
      pi <- DEpcc_grp.Pval[, i]
      pi.sortidx <- order(pi, decreasing = F)
      pi.topidx <- pi.sortidx[c(1: toprank)]
      mi[pi.topidx] <- 1
      DEpcc_eachBat[, i] <- mi
    }
  }
  
  if (sel.type == 1) {
    # firstly select features p-value < 0.1, secondly select VIP top rank
    for (i in c(1: batch.num)) {
      vi <- DEpcc_grp.VIP[, i]
      pi <- DEpcc_grp.Pval[, i]
      
      pi_1.num <- sum(pi < alp)
      pi_1.idx <- which(pi < alp)
      vi_1 <- vi[pi_1.idx]
      
      mi <- rep(0, pcc.num)
      if (pi_1.num > toprank) {
        vi.sortidx <- pi_1.idx[order(vi_1, decreasing = T)]
        vi.topidx <- vi.sortidx[c(1: toprank)]
        mi[vi.topidx] <- 1
      }
      else {
        mi[pi < alp] <- 1
      }
      
      DEpcc_eachBat[, i] <- mi
    }
  }
  
  if (sel.type == 2) {
    # firstly select features VIP > 1, secondly select p-value top rank
    for (i in c(1: batch.num)) {
      vi <- DEpcc_grp.VIP[, i]
      pi <- DEpcc_grp.Pval[, i]
      
      vi_1.num <- sum(vi > VIP.thresh)
      vi_1.idx <- which(vi > VIP.thresh)
      pi_1 <- pi[vi_1.idx]
      
      mi <- rep(0, pcc.num)
      if (vi_1.num > toprank) {
        pi.sortidx <- vi_1.idx[order(pi_1, decreasing = F)]
        pi.topidx <- pi.sortidx[c(1: toprank)]
        mi[pi.topidx] <- 1
      }
      else {
        mi[vi > VIP.thresh] <- 1
      }
      
      DEpcc_eachBat[, i] <- mi
    }
  }
  
  DEpcc_lt2Bat.Idx <- which(rowSums(DEpcc_eachBat) >= min.batch.num)
  
  return(DEpcc_lt2Bat.Idx)
}

identCorrBioM_1 <- function(dat.sbj, group.sbj, PotenCorrBioMarks, 
                            rankPcc.method = 'both', 
                            rankWeight = c(0.5, 0.5), 
                            BH.adj = T, max.topRank) {
  N <- nrow(dat.sbj)
  p <- ncol(dat.sbj)
  
  grp.f <- as.factor(group.sbj)
  group.num <- nlevels(grp.f)
  group.lev <- levels(grp.f)
  if (group.num > 2) {
    cat('Use only ', group.lev[1], ' and ', group.lev[2], ' group.')
  }
  
  pcc.num <- p * (p - 1) / 2
  
  samp.disturb <- matrix(0, N, pcc.num)
  ref.samp <- dat.sbj[group.sbj == group.lev[1], ]
  ref.pcc <- cor(ref.samp)
  
  for (k in c(1: N)) {
    dpcc.k <- cor(rbind(ref.samp, dat.sbj[k, ])) - ref.pcc
    samp.disturb[k, ] <- dpcc.k[upper.tri(dpcc.k)]
  }
  
  samp.disturb <- samp.disturb + 2
  
  samp.disturb.g1 <- samp.disturb[group.sbj == group.lev[1], ]
  samp.disturb.g2 <- samp.disturb[group.sbj == group.lev[2], ]
  N.1 <- nrow(samp.disturb.g1)
  N.2 <- nrow(samp.disturb.g2)
  
  if (rankPcc.method == 'VIP') {
    sd.plsda <- opls(samp.disturb, group.sbj, orthoI = 0, 
                    predI = 3, fig.pdfC = 'none', info.txtC = 'none')
    sd.VIP <- getVipVn(sd.plsda)
    sort.PccIdx <- order(sd.VIP, decreasing = T)
  }
  
  if (rankPcc.method == 'ttest') {
    if (max(N.1, N.2) / min(N.1, N.2) >= 1.5) {
      ttest.pval <- getPval_unbalGrp(samp.disturb.g1, samp.disturb.g2, 
                                     BH.adj = BH.adj)
    }
    else {
      ttest.result <- getttestResult(samp.disturb.g1, samp.disturb.g2, 
                                     BH.adj = BH.adj)
      ttest.pval <- ttest.result$tt.Pval
    }
    sort.PccIdx <- order(ttest.pval)
  }
  
  if (rankPcc.method == 'both') {
    # VIP
    sd.plsda <- opls(samp.disturb, group.sbj, orthoI = 0, 
                     predI = 3, fig.pdfC = 'none', info.txtC = 'none')
    sd.VIP <- getVipVn(sd.plsda)
    VIP.sort.PccIdx <- order(sd.VIP, decreasing = T)
    VIP.pccRank <- order(VIP.sort.PccIdx)
    
    # t-test
    if (max(N.1, N.2) / min(N.1, N.2) >= 1.5) {
      ttest.pval <- getPval_unbalGrp(samp.disturb.g1, samp.disturb.g2, 
                                     BH.adj = BH.adj)
    }
    else {
      ttest.result <- getttestResult(samp.disturb.g1, samp.disturb.g2, 
                                     BH.adj = BH.adj)
      ttest.pval <- ttest.result$tt.Pval
    }
    Pval.sort.PccIdx <- order(ttest.pval, decreasing = F)
    Pval.pccRank <- order(Pval.sort.PccIdx)
    
    # sel.idx <- which(ttest.pval < alp)
    # sort.PccIdx <- sel.idx[order(sd.VIP[sel.idx], decreasing = T)]
    pcc.mRank <- rankWeight[1] * VIP.pccRank + rankWeight[2] * Pval.pccRank
    sort.PccIdx <- order(pcc.mRank)
  }
  
  if (max.topRank > length(sort.PccIdx)) {
    max.topRank <- length(sort.PccIdx)
  }
  
  if (!is.null(PotenCorrBioMarks)) {
    OL.DEpcc.res <- getOLperc(PotenCorrBioMarks, sort.PccIdx, 
                              max.topRank = max.topRank)
  }
  else {
    OL.DEpcc.res <- NULL
  }
  
  identBioM.res <- list(sort.Pcc = sort.PccIdx, 
                        OL.DEpcc = OL.DEpcc.res)
  
  return(identBioM.res)
  
}

set_axisBreak <- function(minVal = 0, maxVal) {
  dVal <- maxVal - minVal
  candi.br <- c(1, 2, 5)
  while (dVal / candi.br[3] > 10) {
    candi.br <- 10 * candi.br
  }
  
  br.num <- round(dVal / candi.br)
  br.num_5 <- abs(br.num - 5)
  sel.br <- candi.br[which(br.num_5 == min(br.num_5))[1]]
  
  return(sel.br)
}

# Main function: Preservation of biological effects
PresBioEffect <- function(dat.sbj.list, dat.intraCor = NULL, batch.sbj, group.sbj, 
                          tar.method.name = 'CordBat', 
                          PCA.print = T, OLDEmet.print = T, OLDEpcc.print = T, 
                          PLSDA.print = T, ClasAcc.print = T, ClasROC.print = T, 
                          selBioM.type = c(1, 1), BH.adj = c(F, T), alp = c(0.1, 0.1), 
                          r1 = c(0.2, 0.05), r2 = c(2, 0.5), rank.method = c('both', 'VIP'), 
                          min.batch.num = c(2, 2), useDEforClas = T, DE.ratio = 0.3, 
                          save.path0) {
  
  dat.num <- length(dat.sbj.list)
  dat.name <- names(dat.sbj.list)
  
  # 4.1 PCA (group)
  if (PCA.print) {
    for (i in c(1: dat.num)) {
      dat.i <- dat.sbj.list[[i]]
      pca.i <- pca(dat.i, ncomp = 2, center = T, scale = T)
      pca.i.expl <- pca.i$prop_expl_var$X
      pca.i.varX <- pca.i$variates$X
      
      p1 <- PCA_scatter1(data = pca.i.varX, group = group.sbj,
                         expl.var = pca.i.expl,
                         xlim = c(min(pca.i.varX[, 1]) - 1, max(pca.i.varX[, 1]) + 1),
                         ylim = c(min(pca.i.varX[, 2]) - 1, max(pca.i.varX[, 2]) + 1),
                         title = dat.name[i])
      print(p1)
      ggsave(paste0(save.path0, 'PCA_grp_', dat.name[i], '_', gsub('-', '', Sys.Date()), '.emf'),
             plot = p1, width = 7.5, height = 7.6, units = 'in', dpi = 300)
    }
  }
  
  # 4.2 identifying DE metabolites between groups
  if (OLDEmet.print) {
    p <- ncol(dat.sbj.list[[1]])
    
    if (is.null(dat.intraCor)) {
      dat.bef <- dat.sbj.list[['Uncorrect']]
    }
    else {
      dat.bef <- dat.intraCor
    }
    DEmet.GT <- PotenBiomarker(dat.bef, batch.sbj, group.sbj,
                               sel.type = selBioM.type[1], BH.adj = BH.adj[1], 
                               toprank = round(r1[1] * p), alp = alp[1], 
                               min.batch.num = min.batch.num[1])
    
    OLmetNum.list <- list()
    if (selBioM.type[1] == 1) {
      rW <- c(0.8, 0.2)
    }
    
    if (selBioM.type[1] == 2) {
      rW <- c(0.65, 0.35)
    }
    
    for (i in c(1: dat.num)) {
      dat.i <- dat.sbj.list[[i]]
      OLmet.res <- identBiomarker(dat.i, group.sbj, potenBiomarks = DEmet.GT,
                                  rankMet.method = rank.method[1],
                                  rankWeight = rW, BH.adj = BH.adj[1],
                                  max.topRank = round(r2[1] * length(DEmet.GT)))
      OLmetNum.list[[i]] <- OLmet.res$OL.DEmet$OL.Num
    }
    topRank <- OLmet.res$OL.DEmet$topRank
    
    # plot
    dattype.vec <- unlist(lapply(dat.name, function(x) rep(x, length(topRank))))
    dattype <- factor(dattype.vec, levels = dat.name)
    
    dat.OLmetNum <- data.frame(dattype = dattype,
                               top.n = rep(topRank, dat.num),
                               OLnum = unlist(OLmetNum.list))
    
    max.olnum <- max(unlist(OLmetNum.list))
    
    gap_x <- set_axisBreak(maxVal = max(topRank))
    br <- seq(0, round(max(topRank)/gap_x) * gap_x, gap_x)
    
    gap_y <- set_axisBreak(maxVal = max.olnum)
    
    fon <- 'sans'
    all_col <- c('black', 'blue', 'yellow', 'magenta', 'cyan', 'red')
    set_col <- c(all_col[c(1: (dat.num - 1))], 'red')
    p2 <- ggplot(data = dat.OLmetNum, aes(x = top.n, y = OLnum, group = dattype, color = dattype)) +
      geom_line(linewidth = 1.5) + 
      scale_color_manual(values = set_col) +
      scale_x_continuous(breaks = br,
                         limits = c(0, max(topRank))) +
      scale_y_continuous(breaks = seq(0, round(max.olnum/gap_y) * gap_y, gap_y),
                         limits = c(0, max.olnum + 1)) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 22, family = fon, face = 'bold', color = 'black'),
            axis.text.y = element_text(size = 22, family = fon, face = 'bold', color = 'black'),
            axis.title.x = element_text(size = 25, family = fon, face = 'bold'),
            axis.title.y = element_text(size = 25, family = fon, face = 'bold'),
            axis.ticks = element_line(color = 'black', linewidth = 1),
            axis.ticks.length = unit(2, 'mm'),
            panel.border = element_rect(fill = NA, linewidth = 2),
            panel.grid = element_blank(),
            legend.position = c(0.16, 0.8),
            legend.text = element_text(size = 18, family = fon, face = 'bold', margin = margin(b = 3))) +
      xlab('Top N') + ylab('# overlapped metabolites') +
      labs(group = NULL, color = NULL, title = NULL)
    
    print(p2)
    ggsave(paste0(save.path0, 'OLDEmet_all_', gsub('-', '', Sys.Date()), '.emf'),
           plot = p2, width = 7.5, height = 7, units = 'in', dpi = 300)
  }


  # 4.3 identifying DE correlations between groups
  if (OLDEpcc.print) {
    p <- ncol(dat.sbj.list[[1]])
    pcc.num <- p * (p - 1) / 2
    
    if (is.null(dat.intraCor)) {
      dat.bef <- dat.sbj.list[['Uncorrect']]
    }
    else {
      dat.bef <- dat.intraCor
    }
    DEpcc.GT <- PotenCorrBioM_1(dat.bef, batch.sbj, group.sbj, BH.adj = BH.adj[2],
                                sel.type = selBioM.type[2], alp = alp[2], 
                                toprank = round(r1[2] * pcc.num), 
                                min.batch.num = min.batch.num[2])
    
    if (selBioM.type[2] == 1) {
      rW <- c(0.65, 0.35)
      
    }
    
    if (selBioM.type[2] == 2) {
      rW <- c(0.65, 0.35)
    }
    
    OLpccNum.list <- list()
    
    cat('Get overlapped biomarkers in data: \n')
    for (i in c(1: dat.num)) {
      dat.i <- dat.sbj.list[[i]]
      OLpcc.res <- identCorrBioM_1(dat.i, group.sbj, PotenCorrBioMarks = DEpcc.GT,
                                   rankPcc.method = rank.method[2], rankWeight = rW,
                                   BH.adj = BH.adj[2],
                                   max.topRank = round(r2[2] * length(DEpcc.GT)))
      OLpccNum.list[[i]] <- OLpcc.res$OL.DEpcc$OL.Num
      
      cat(dat.name[i], '\n')
    }
    topRank <- OLpcc.res$OL.DEpcc$topRank
    
    # plot
    dattype.vec <- unlist(lapply(dat.name, function(x) rep(x, length(topRank))))
    dattype <- factor(dattype.vec, levels = dat.name)
    
    dat.OLpccNum <- data.frame(dattype = dattype,
                               top.n = rep(topRank, dat.num),
                               OLnum = unlist(OLpccNum.list))
    
    max.olnum <- max(unlist(OLpccNum.list))
    
    gap_x <- set_axisBreak(maxVal = max(topRank))
    br <- seq(0, round(max(topRank)/gap_x) * gap_x, gap_x)
    
    gap_y <- set_axisBreak(maxVal = max.olnum)
    
    fon <- 'sans'
    all_col <- c('black', 'blue', 'yellow', 'magenta', 'cyan', 'red')
    set_col <- c(all_col[c(1: (dat.num - 1))], 'red')
    p3 <- ggplot(data = dat.OLpccNum, aes(x = top.n, y = OLnum, group = dattype, color = dattype)) +
      geom_line(linewidth = 1.5) +
      scale_color_manual(values = set_col) + 
      scale_x_continuous(breaks = br,
                         limits = c(0, max(topRank))) +
      scale_y_continuous(breaks = seq(0, round(max.olnum/gap_y) * gap_y, gap_y),
                         limits = c(0, max.olnum + 1)) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 22, family = fon, face = 'bold', color = 'black'),
            axis.text.y = element_text(size = 22, family = fon, face = 'bold', color = 'black'),
            axis.title.x = element_text(size = 25, family = fon, face = 'bold'),
            axis.title.y = element_text(size = 25, family = fon, face = 'bold'),
            axis.ticks = element_line(color = 'black', linewidth = 1),
            axis.ticks.length = unit(2, 'mm'),
            panel.border = element_rect(fill = NA, linewidth = 2),
            panel.grid = element_blank(),
            legend.position = c(0.16, 0.8),
            legend.text = element_text(size = 18, family = fon, face = 'bold', margin = margin(b = 3))) +
      xlab('Top N') + ylab('# overlapped correlations') +
      labs(group = NULL, color = NULL, title = NULL)
    
    print(p3)
    ggsave(paste0(save.path0, 'OLDEpcc_all_', gsub('-', '', Sys.Date()), '.emf'),
           plot = p3, width = 7.5, height = 7, units = 'in', dpi = 300)
  }

  # 4.4 PLSDA model
  if (PLSDA.print) {
    p <- ncol(dat.sbj.list[[1]])
    
    grp.f <- as.factor(group.sbj)
    group.num <- nlevels(grp.f)
    group.lev <- levels(grp.f)
    
    if (group.num > 2) {
      cat('Use only ', group.lev[1], ' and ', group.lev[2], ' group.\n')
    }
    
    R2Y.vals <- rep(0, dat.num)
    Q2.vals <- rep(0, dat.num)
    
    # sel.num <- round(0.2 * p)
    for (i in c(1: dat.num)) {
      dat.i <- dat.sbj.list[[i]]
      
      # sortMet.i <- identBiomarker(dat.i, group.sbj)
      # selMetIdx.i <- sortMet.i$sort.Met[c(1: sel.num)]
      # dat.i <- dat.i[, selMetIdx.i]
      
      # PLS-DA model
      plsda.i <- opls(dat.i, group.sbj, orthoI = 0,
                      predI = 2, fig.pdfC = 'none', info.txtC = 'none')
      pls.score.i <- plsda.i@scoreMN
      
      p4 <- PLSDA_scatter(pls.score.i, group = group.sbj,
                          xlim = c(min(pls.score.i[, 1]) - 1, max(pls.score.i[, 1]) + 1),
                          ylim = c(min(pls.score.i[, 2]) - 1, max(pls.score.i[, 2]) + 1),
                          title = dat.name[i])
      
      print(p4)
      ggsave(paste0(save.path0, 'PLSDA_', dat.name[i], '_', gsub('-', '', Sys.Date()), '.emf'),
             plot = p4, width = 7.3, height = 7, units = 'in', dpi = 300)
      
      R2Y.vals[i] <- plsda.i@summaryDF$`R2Y(cum)`
      Q2.vals[i] <- plsda.i@summaryDF$`Q2(cum)`
      
      cat(R2Y.vals[i], '  ', Q2.vals[i], '\n')
      
    }
    
    RQ.df <- data.frame(dat = dat.name,
                        R2Y = R2Y.vals,
                        Q2 = Q2.vals)
    
    write.table(RQ.df, paste0(save.path0, 'PLSDA_result_', gsub('-', '', Sys.Date()), '.txt'))
  }

  
  # 4.5 Classification accuracy
  if (any(c(ClasAcc.print, ClasROC.print))) {
    AUC.list <- list()
    ACC.list <- list()
    ROC.list <- list()
    
    p <- ncol(dat.sbj.list[[1]])
    sel.num <- round(DE.ratio * p)
    cat('Construct SVM classification model for data: \n')
    for (i in c(1: dat.num)) {
      dat.i <- dat.sbj.list[[i]]
      
      if (useDEforClas) {
        sortMet.i <- identBiomarker(dat.i, group.sbj)
        selMetIdx.i <- sortMet.i$sort.Met[c(1: sel.num)]
        
        dat.i <- dat.i[, selMetIdx.i]
      }
      
      cvSVMens.roc.i <- CVforSVMensemble(dat.i, group.sbj)
      
      roc.i <- data.frame(spec = rev(1 - cvSVMens.roc.i$spec), 
                          sens = rev(cvSVMens.roc.i$sens))
      ROC.list[[i]] <- roc.i
      
      auc.i <- cvSVMens.roc.i$auc.val
      AUC.list[[i]] <- auc.i
      
      acc.i <- cvSVMens.roc.i$accuracy
      ACC.list[[i]] <- acc.i
      
      cat(dat.name[i], '\n')
    }
    
    names(ROC.list) <- dat.name
    names(AUC.list) <- dat.name
    names(ACC.list) <- dat.name
    
    if (ClasAcc.print) {
      # plot
      # Prediction accuracy of the SVM models
      mean.ACCs <- sapply(ACC.list, FUN = mean)
      sd.ACCs <- sapply(ACC.list, FUN = sd)
      acc.val <- data.frame(m = mean.ACCs, 
                            sd = sd.ACCs)
      rownames(acc.val) <- dat.name
      
      dat.lev <- dat.name[order(mean.ACCs, decreasing = T)]
      dat.acc <- data.frame(dat.type = factor(rownames(acc.val), levels = dat.lev),
                            acc.mean = acc.val$m, 
                            acc.sd = acc.val$sd)
      
      fon <- 'sans'
      max.acc <- max(c(dat.acc$acc.mean - dat.acc$acc.sd, dat.acc$acc.mean + dat.acc$acc.sd))
      min.acc <- min(c(dat.acc$acc.mean - dat.acc$acc.sd, dat.acc$acc.mean + dat.acc$acc.sd))
      
      dacc <- max.acc - min.acc
      if (dacc * 100 < 10) {
        gap <- 0.01
      } else if (dacc * 10 < 5) {
        gap <- 0.02
      } else {
        gap <- 0.1
      }
      
      p5 <- ggplot(dat.acc, aes(fill = dat.type, y = dat.type, x = acc.mean)) +
        geom_bar(position = position_dodge(preserve = 'single'), stat = 'identity',
                 width = 0.7, colour = NA, linewidth = 1.5) + 
        geom_errorbar(aes(xmin = acc.mean - acc.sd, xmax = acc.mean + acc.sd),
                      position = position_dodge2(preserve = 'single', padding = 0.5),
                      linewidth = 1.5, width = 0.3) +
        scale_y_discrete(limits = rev, expand = expansion(0.16)) +
        scale_x_break(c(0.05, round(min.acc - 0.0001, 4)), 
                      ticklabels = seq(round(min.acc - 0.001, 2), round(max.acc + 0.001, 2), gap),
                      scales = 10) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, round(max.acc + 0.002, 3)),
                           breaks = c(0, seq(round(min.acc - 0.001, 2), round(max.acc + 0.001, 2), gap))) +
        scale_fill_npg() + 
        theme_bw() +
        theme(axis.text.x.top = element_blank(), 
              axis.text.x.bottom = element_text(size = 20, family = fon, face = 'bold', 
                                                color = 'black'), 
              axis.text.y = element_text(size = 22, family = fon, face = 'bold', 
                                         color = 'black'), 
              axis.title.y = element_blank(), 
              axis.title.x = element_text(size = 25, family = fon, face = 'bold', 
                                          hjust = 0.65), 
              axis.ticks.x.top = element_blank(), 
              axis.ticks.x.bottom = element_line(linewidth = 1, color = 'black'), 
              axis.ticks.y = element_line(linewidth = 1, color = 'black'), 
              axis.ticks.length.y = unit(2, 'mm'), 
              axis.ticks.length.x.bottom = unit(2, 'mm'), 
              panel.border = element_blank(), 
              panel.grid = element_blank(), 
              axis.line.x.bottom = element_line(linewidth = 1.5), 
              axis.line.x.top = element_blank(), 
              axis.line.y = element_line(linewidth = 2), 
              legend.position = 'none') + 
        xlab('Accuracy') + 
        labs(fill = 'Data', title = NULL) 
      
      print(p5)
      ## scale_x_break会导致图无法保存为emf文件，因此需要手动保存
      # ggsave(paste0(save.path0, 'ClasAcc_', gsub('-', '', Sys.Date()), '.pdf'),
      #        plot = p5, width = 8, height = 7, units = 'in', dpi = 300)
    }
    
    if (ClasROC.print) {
      # ROC curves for the SVM models
      TopAcc.dat <- dat.lev[c(1: 3)]
      top.tar <- grep(tar.method.name, TopAcc.dat)
      if (length(top.tar) == 0) {
        topdat.name <- c('Uncorrect', TopAcc.dat[c(1, 2)], tar.method.name)
      }
      else {
        topdat.name <- c('Uncorrect', TopAcc.dat[-top.tar], tar.method.name)
      }
      Top.ROC.list <- ROC.list[topdat.name]
      Top.AUCs <- unlist(AUC.list[topdat.name])
      
      topdat.label <- rep('', length(topdat.name))
      for (i in c(1: length(topdat.name))) {
        topdat.label[i] <- paste0(topdat.name[i], ', AUC = ', round(Top.AUCs[i], 2))
      }
      
      dat.roc <- data.frame(dat.type = factor(rep(topdat.label, 
                                                  each = length(Top.ROC.list[[1]][, 1])), 
                                              levels = topdat.label), 
                            spec = unlist(lapply(Top.ROC.list, FUN = function(x) x$spec)),
                            sens = unlist(lapply(Top.ROC.list, FUN = function(x) x$sens)))
      
      
      
      fon <- 'sans'
      set_col <- c('black', 'green', 'blue', 'red')
      p6 <- ggplot(data = dat.roc, aes(x = spec, y = sens, group = dat.type, color = dat.type)) +
        geom_line(linewidth = 1.5) + 
        scale_color_manual(values = set_col) + 
        # scale_x_continuous(breaks = seq(0, 1, 0.25)) + 
        # scale_y_continuous(breaks = seq(0, 1, 0.25)) + 
        xlab('1 - specificity') + ylab('sensitivity') + 
        xlim(0, 1) + ylim(0, 1) +
        theme_bw() +
        theme(axis.title = element_text(size = 25, family = fon, face = 'bold'),
              axis.text = element_text(size = 22, family = fon, face = 'bold', color = 'black'), 
              axis.ticks = element_line(linewidth = 1, color = 'black'), 
              axis.ticks.length = unit(2, 'mm'), 
              legend.position = c(0.65, 0.15), 
              legend.title = element_blank(),
              legend.text = element_text(size = 22, family = fon, face = 'bold', 
                                         margin = margin(b = 5)),
              plot.title = element_text(size = 25, family = fon, face = 'bold', hjust = 0.5),
              panel.border = element_rect(linewidth = 2), 
              panel.grid = element_blank()) +
        labs(color = NULL, title = NULL)
      
      print(p6)
      ggsave(paste0(save.path0, 'ROC_clas_', gsub('-', '', Sys.Date()), '.emf'),
             plot = p6, width = 7.5, height = 7, units = 'in', dpi = 300)
    }
  }
}











