# Title: Data IV BEC and evaluation
# Author: Fanjing Guo
# Date: 2023.2.28

# --------------------- Data pre-processing and batch effect correction ------------------- #
# Untargeted urine LC-HRMS metabolomics data
# X <- read.csv("raw6_filled_filtered_mzmine.csv", header = TRUE)
X <- read.csv("raw10_filled_filtered_mzmine.csv", header = TRUE)
mzs <- X$row.m.z
X_Dat <- as.matrix(X[, -c(1: 3)])
X_Dat <- X_Dat[, -312]
N <- ncol(X_Dat)
p <- nrow(X_Dat)

# delete features with same m/z value
NA.num <- rowSums(X_Dat == 0)
feat.order <- order(NA.num)
X_Dat <- X_Dat[feat.order, ]
mzs <- mzs[feat.order]
del.featIdx <- which(duplicated(round(mzs, 1)))  
X_Dat.del <- X_Dat[-del.featIdx, ]
mzs <- mzs[-del.featIdx]
N <- ncol(X_Dat.del)
p <- nrow(X_Dat.del)

# delete missing value and impute
NA.num <- rowSums(X_Dat.del == 0)
del.featIdx <- which(NA.num > 0.15 * N)
X_Dat.del1 <- X_Dat.del[-del.featIdx, ]
N <- ncol(X_Dat.del1)
p <- nrow(X_Dat.del1)
mzs.del <- mzs[-del.featIdx]
mzs.list <- paste0('m/z ', round(as.numeric(mzs.del), 3))

# impute missing value
X_Dat.del1[X_Dat.del1 == 0] <- NA

library(DMwR2)
X_Dat.del1 <- as.data.frame(t(X_Dat.del1))
X_Dat.del.imp <- knnImputation(X_Dat.del1)

X_info <- read.csv('Data IV_info.csv', header = T)
samp.name <- rownames(X_Dat.del.imp)
samp.name <- sapply(strsplit(samp.name, '[.]'), '[[', 3)
samp.order <- rep(0, N)
for (i in c(1: N)) {
  samp.order[i] <- X_info$RunOrder[X_info$SampName == samp.name[i]]
}
X_Dat.del.imp <- X_Dat.del.imp[order(samp.order), ]
samp.name <- samp.name[order(samp.order)]
samp.order <- sort(samp.order)
rownames(X_Dat.del.imp) <- samp.name
X_Dat.del.imp <- as.matrix(X_Dat.del.imp) 

# get sample information
group <- X_info$Group
batch <- X_info$Batch

# log-transformed
X_Dat.log.proc <- log2(X_Dat.del.imp)

# # pca
# library(mixOmics)
# pca.befintra <- pca(X_Dat.log.proc, ncomp = 2, center = T, scale = T)
# pca.befintra.varX <- pca.befintra$variates$X
# pca.befintra.expl <- pca.befintra$prop_expl_var$X
# 
# pc1.max <- max(pca.befintra.varX[, 1])
# pc1.min <- min(pca.befintra.varX[, 1])
# pc2.max <- max(pca.befintra.varX[, 2])
# pc2.min <- min(pca.befintra.varX[, 2])
# PCA_scatter_withQC(data = pca.befintra.varX, batch = batch.mer, group = group, 
#                    QC_label = 'QC', expl.var = pca.befintra.expl, 
#                    xlim = c(pc1.min, pc1.max), ylim = c(pc2.min, pc2.max))
# 
# dat.sbj <- X_Dat.log.proc[group != 'QC', ]
# pca.befintra <- pca(dat.sbj, ncomp = 2, center = T, scale = T)
# pca.befintra.varX <- pca.befintra$variates$X
# pca.befintra.expl <- pca.befintra$prop_expl_var$X
# 
# pc1.max <- max(pca.befintra.varX[, 1])
# pc1.min <- min(pca.befintra.varX[, 1])
# pc2.max <- max(pca.befintra.varX[, 2])
# pc2.min <- min(pca.befintra.varX[, 2])
# PCA_scatter1(data = pca.befintra.varX, group = group[group != 'QC'], 
#              expl.var = pca.befintra.expl, xlim = c(pc1.min, pc1.max), 
#              ylim = c(pc2.min, pc2.max))

# get subject samples
batch.sbj <- batch[group != 'QC']
group.sbj <- group[group != 'QC']
samp.name.sbj <- samp.name[group != 'QC']
samp.order.sbj <- samp.order[group != 'QC']
Sbj_Dat.log.proc <- X_Dat.log.proc[group != 'QC', ]

# intra-batch correction (LOESS - Only subject samples)
bat.f <- as.factor(batch.sbj)
batch.num <- nlevels(bat.f)
batch.lev <- levels(bat.f)
Sbj_Dat.intraCor <- Sbj_Dat.log.proc
for (i in c(1: batch.num)) {
  dat.i <- Sbj_Dat.log.proc[bat.f == batch.lev[i], ]
  ord.i <- samp.order.sbj[bat.f == batch.lev[i]]
  N.i <- nrow(dat.i)
  dat.i.cor <- dat.i
  for (j in c(1: p)) {
    dat <- data.frame(x = ord.i,
                      y = dat.i[, j])
    fitted.result <- loess(y ~ x, data = dat, span = 0.75)
    fit.metj <- fitted.result$fitted
    y.cor <- dat.i[, j] + median(dat.i[, j]) - fit.metj
    dat.i.cor[, j] <- y.cor
  }
  Sbj_Dat.intraCor[bat.f == batch.lev[i], ] <- dat.i.cor
  cat(i, '\n')
}

Sbj_Dat.intraCor <- as.matrix(Sbj_Dat.intraCor)
colnames(Sbj_Dat.intraCor) <- c(1: ncol(Sbj_Dat.intraCor))

# # PQN normalization - intra-batch
# X_Dat.intraCor1 <- X_Dat.intraCor
# for (i in c(1: batch.num)) {
#   dat.bat.i.CG <- X_Dat.intraCor[batch == batch.lev[i] & group == 'CG' , ]
#   N.i <- nrow(dat.bat.i.CG)
#   ref.samp <- apply(dat.bat.i.CG, MARGIN = 2, FUN = median)
#   ratio.by_ref <- matrix(ref.samp, nrow = N.i, ncol = p, byrow = T) / dat.bat.i.CG
#   norm.coef <- apply(ratio.by_ref, MARGIN = 1, FUN = median)
#   dat.bat.i.norm <- diag(norm.coef) %*% dat.bat.i.CG
#   X_Dat.intraCor1[batch == batch.lev[i] & group == 'CG', ] <- dat.bat.i.norm
# 
#   dat.bat.i.TG <- X_Dat.intraCor[batch == batch.lev[i] & group == 'TG' , ]
#   N.i <- nrow(dat.bat.i.TG)
#   ref.samp <- apply(dat.bat.i.TG, MARGIN = 2, FUN = median)
#   ratio.by_ref <- matrix(ref.samp, nrow = N.i, ncol = p, byrow = T) / dat.bat.i.TG
#   norm.coef <- apply(ratio.by_ref, MARGIN = 1, FUN = median)
#   dat.bat.i.norm <- diag(norm.coef) %*% dat.bat.i.TG
#   X_Dat.intraCor1[batch == batch.lev[i] & group == 'TG', ] <- dat.bat.i.norm
# 
#   dat.bat.i.QC <- X_Dat.intraCor[batch == batch.lev[i] & group == 'QC' , ]
#   N.i <- nrow(dat.bat.i.QC)
#   ref.samp <- apply(dat.bat.i.QC, MARGIN = 2, FUN = median)
#   ratio.by_ref <- matrix(ref.samp, nrow = N.i, ncol = p, byrow = T) / dat.bat.i.QC
#   norm.coef <- apply(ratio.by_ref, MARGIN = 1, FUN = median)
#   dat.bat.i.norm <- diag(norm.coef) %*% dat.bat.i.QC
#   X_Dat.intraCor1[batch == batch.lev[i] & group == 'QC', ] <- dat.bat.i.norm
# }

# pca
library(mixOmics)
pca.intra <- pca(X_Dat.intraCor, ncomp = 2, center = T, scale = T)
# pca.intra <- pca(X_Dat.intraCor1, ncomp = 2, center = T, scale = T)
pca.intra.varX <- pca.intra$variates$X
pca.intra.expl <- pca.intra$prop_expl_var$X

pc1.max <- max(pca.intra.varX[, 1])
pc1.min <- min(pca.intra.varX[, 1])
pc2.max <- max(pca.intra.varX[, 2])
pc2.min <- min(pca.intra.varX[, 2])
PCA_scatter_withQC(data = pca.intra.varX, batch = batch, group = group, 
                   QC_label = 'QC', expl.var = pca.intra.expl, 
                   xlim = c(pc1.min, pc1.max), ylim = c(pc2.min, pc2.max))

dat.sbj <- X_Dat.intraCor[group != 'QC', ]
pca.intra <- pca(dat.sbj, ncomp = 2, center = T, scale = T)
pca.intra.varX <- pca.intra$variates$X
pca.intra.expl <- pca.intra$prop_expl_var$X

pc1.max <- max(pca.intra.varX[, 1])
pc1.min <- min(pca.intra.varX[, 1])
pc2.max <- max(pca.intra.varX[, 2])
pc2.min <- min(pca.intra.varX[, 2])
PCA_scatter1(data = pca.intra.varX, group = group[group != 'QC'], 
             expl.var = pca.intra.expl, xlim = c(pc1.min, pc1.max), 
             ylim = c(pc2.min, pc2.max))

# Batch effect correction
# CorbBat BEC
getRefbat(X_Dat.intraCor, batch.mer)
getRefbat(Sbj_Dat.intraCor, batch.sbj)
CdBDat <- CordBat(X_Dat.intraCor, batch, ref.batch = 9, eps = 1e-4, 
                  print.detail = F)

CdBDat <- CordBat(Sbj_Dat.intraCor, batch.sbj, ref.batch = 9, eps = 1e-4, 
                  print.detail = F)
# CdBDat <- CordBat(X_Dat.intraCor1, batch, ref.batch = 9, eps = 1e-4, 
#                   print.detail = F)
CdBcor_Data_IV <- CdBDat$X.cor
Sbj_Dat.delout <- CdBDat$X.delout
batch.sbj.delout <- CdBDat$batch.new
# CdBcor_Data_IV <- CdBDat$X.cor.withQC
X_Dat.delout <- CdBDat$X.delout
batch.delout <- CdBDat$batch.new
write.csv(CdBcor_Data_IV, 'CdBcor_Data_IV_6.csv', row.names = F)
write.csv(CdBcor_Data_IV, 'CdBcor_Data_IV_6_LES_PQN.csv', row.names = F)
write.csv(CdBcor_Data_IV, 'CdBcor_Data_IV_6_LES_PQN_ref9.csv', row.names = F)
write.csv(CdBcor_Data_IV, 'CdBcor_Data_IV_6_LES_ref9.csv', row.names = F)
write.csv(CdBcor_Data_IV, 'CdBcor_Data_IV_6_LES_ref9_new.csv', row.names = F)
write.csv(CdBcor_Data_IV, 'CdBcor_sbjData_IV_6_LES_ref9.csv', row.names = F)

CdBcor_Data_IV <- read.csv('CdBcor_Data_IV_6.csv', header = T)
CdBcor_Data_IV <- read.csv('CdBcor_Data_IV_6_LES_PQN.csv', header = T)
CdBcor_Data_IV <- read.csv('CdBcor_Data_IV_6_LES_PQN_ref9.csv', header = T)
CdBcor_Data_IV <- read.csv('CdBcor_Data_IV_6_LES_ref9.csv', header = T)
CdBcor_Data_IV <- read.csv('CdBcor_Data_IV_6_LES_ref9_new.csv', header = T)

# ------------------------------- Evaluation ----------------------------------------- #
savefile.path0 <- 'F:/Graduate/Graduate Studies/Papers/1-GGM for Batch Effect Correction/New Results/Data IV/'

# ------------------------
# removal of batch effect
# ------------------------
RMbatVariation(X_Dat.intraCor, CdBcor_Data_IV, batch, group, samp.order, 
               QC_label = 'QC', Fval.cutoff = 20, feature.list = mzs.list, 
               save.path0 = paste0(savefile.path0, 'Eval 1/'))

RMbatVariation(X_Dat.intraCor, CdBcor_Data_IV, batch.mer, group, samp.order, 
               QC_label = 'QC', Fval.cutoff = 20, feature.list = mzs.list, 
               save.path0 = paste0(savefile.path0, 'Eval 1/'))

RMbatVariation(Sbj_Dat.intraCor, CdBcor_Data_IV, batch.sbj, group.sbj, samp.order.sbj, 
               QC_label = 'QC', Fval.cutoff = 20, feature.list = mzs.list, 
               save.path0 = paste0(savefile.path0, 'Eval 1/without QC/'))

RMbatVariation(X_Dat.log.proc, X_Dat.intraCor.mer, CdBcor_Data_IV, 
               batch.mer, group, samp.order, QC_label = 'QC', 
               Fval.cutoff = 20, feature.list = mzs.list, 
               save.path0 = paste0(savefile.path0, 'Eval 1/'))

# --------------------------------
# Preservation of data structure
# --------------------------------
PresDatStruct(X_Dat.intraCor, X_Dat.delout, CdBcor_Data_IV, batch.delout, 
              NotdifBat.ratio = 0.85, delays = c(0.2, 1), 
              save.path0 = paste0(savefile.path0, 'Eval 2/'))

PresDatStruct(Sbj_Dat.intraCor, Sbj_Dat.delout, CdBcor_Data_IV, batch.sbj.delout, 
              NotdifBat.ratio = 0.8, delays = c(0.2, 1), 
              save.path0 = paste0(savefile.path0, 'Eval 2/without QC/'))

# ------------------------------
# Compared with QC-based method
# ------------------------------
# merge adjacent 2 batches
# -------- batch 1-2: new batch 1 ------------
# -------- batch 3-4: new batch 2 ------------
# -------- batch 5-6: new batch 3 -----------
# -------- batch 7-8: new batch 4 ----------
# -------- batch 9-10: new batch 5 ----------
# -------- batch 11-13: new batch 6 ----------
batch.mer <- rep(0, length(batch))
batch.mer[batch == 1 | batch == 2] <- 1
batch.mer[batch == 3 | batch == 4] <- 2
batch.mer[batch == 5 | batch == 6] <- 3
batch.mer[batch == 7 | batch == 8] <- 4
batch.mer[batch == 9 | batch == 10] <- 5
batch.mer[batch == 11 | batch == 12 | batch == 13] <- 6

X_info.mer <- X_info
X_info.mer$Batch <- batch.mer

# remove the batch effect between the merged batches - LOESS & SERRF
p <- ncol(X_Dat.log.proc)
bat.f <- as.factor(batch.mer)
batch.num <- nlevels(bat.f)
batch.lev <- levels(bat.f)
X_Dat.intraCor.mer1 <- X_Dat.log.proc

# LOESS
for (i in c(1: batch.num)) {
  dat.i <- X_Dat.log.proc[bat.f == batch.lev[i], ]
  ord.i <- samp.order[bat.f == batch.lev[i]]
  N.i <- nrow(dat.i)
  dat.i.cor <- dat.i
  for (j in c(1: p)) {
    dat <- data.frame(x = ord.i,
                      y = dat.i[, j])
    fitted.result <- loess(y ~ x, data = dat, span = 0.75)
    fit.metj <- fitted.result$fitted
    y.cor <- dat.i[, j] + median(dat.i[, j]) - fit.metj
    dat.i.cor[, j] <- y.cor
  }
  X_Dat.intraCor.mer1[bat.f == batch.lev[i], ] <- dat.i.cor
  cat(i, '\n')
}

X_Dat.intraCor.mer1 <- as.matrix(X_Dat.intraCor.mer1)
colnames(X_Dat.intraCor.mer1) <- c(1: ncol(X_Dat.intraCor.mer1))

# SERRF
X_Dat.intraCor.mer <- X_Dat.intraCor.mer1
SRformat_filePath0 <- paste0(getwd(), '/SRformat_file_Data IV_eachBat/')
for (i in c(1: batch.num)) {
  dat.bati <- X_Dat.intraCor.mer1[batch.mer == batch.lev[i], ]
  bati.info <- X_info.mer[batch.mer == batch.lev[i], ]
  GenSRformatdat(dat.bati, bati.info, mzs.list, 
                 wrfile.path = paste0(SRformat_filePath0, 'SRformatDat_batch', 
                                      batch.lev[i], '.csv'))
}

source('app.R')
runApp()

for (i in c(1: batch.num)) {
  zip.i.path <- paste0(SRformat_filePath0, 'SERRF Result', i, '.zip')
  dat.bati.cor <- read.csv(unzip(zip.i.path, 'normalized by - SERRF.csv'), 
                           header = T)
  dat.bati.cor <- as.matrix(t(dat.bati.cor[, -1]))
  X_Dat.intraCor.mer[batch.mer == i, ] <- dat.bati.cor
}
write.csv(X_Dat.intraCor.mer, 'Data IV_SRintraCor.csv', row.names = F)
X_Dat.intraCor.mer <- read.csv('Data IV_SRintraCor.csv', header = T)
X_Dat.intraCor.mer <- as.matrix(X_Dat.intraCor.mer)

write.csv(X_Dat.intraCor.mer, 'Data IV_LES&SRintraCor.csv', row.names = F)
X_Dat.intraCor.mer <- read.csv('Data IV_LES&SRintraCor.csv', header = T)
X_Dat.intraCor.mer <- as.matrix(X_Dat.intraCor.mer)



# Batch effect correction
# CordBat
getRefbat(X_Dat.intraCor.mer, batch.mer)
# CdBDat <- CordBat(X_Dat.intraCor.mer, batch.mer, ref.batch = 1, eps = 1e-4, 
#                   print.detail = F)
# CdBDat <- CordBat(X_Dat.intraCor.mer1, batch.mer, ref.batch = 1, eps = 1e-4, 
#                   print.detail = F)
CdBDat <- CordBat(X_Dat.intraCor.mer, batch.mer, ref.batch = 5, eps = 1e-4, 
                  print.detail = F)
CdBcor_Data_IV <- CdBDat$X.cor

# write.csv(CdBcor_Data_IV, 'CdBcor_Data_IV_6_SR_ref1_mer.csv', row.names = F)
# CdBcor_Data_IV <- read.csv('CdBcor_Data_IV_6_SR_ref1_mer.csv', header = T)
# 
# write.csv(CdBcor_Data_IV, 'CdBcor_Data_IV_6_SRLES_ref1_mer.csv', row.names = F)
# CdBcor_Data_IV <- read.csv('CdBcor_Data_IV_6_SRLES_ref1_mer.csv', header = T)

write.csv(CdBcor_Data_IV, 'CdBcor_Data_IV_6_LESSR_ref5_mer.csv', row.names = F)
CdBcor_Data_IV <- as.matrix(read.csv('CdBcor_Data_IV_6_LESSR_ref5_mer.csv', header = T))

# QC-RLSC
library(statTarget)
QRcor_Data_IV <- QCRLSC_cor(X_Dat.log.proc, samp.info = X_info.mer, 
                            feature.list = mzs.list)
# QRcor_Data_IV <- QCRLSC_cor(X_Dat.intraCor.mer, samp.info = X_info.mer, 
#                             feature.list = mzs.list)

# SERRF
GenSRformatdat(X_Dat.log.proc, X_info.mer, mzs.list, 
               wrfile.path = 'SRformat_Data_IV.csv')
# GenSRformatdat(X_Dat.intraCor.mer1, X_info.mer, mzs.list, 
#                wrfile.path = 'SRformat_Data_IV.csv')
source('app.R')
runApp()
SRcor_Data_IV <- read.csv(unzip(paste0(getwd(), '/SERRF Result - Data IV.zip'), 
                                'normalized by - SERRF.csv'), header = T)
SRcor_Data_IV <- as.matrix(t(SRcor_Data_IV[, -1]))
# unlink('SERRF Result.zip')

# SVR
library(MetNormalizer)
SVRcor_Data_IV <- SVR_cor(X_Dat.log.proc, samp.info = X_info.mer, 
                          mzs.list = mzs.list)
# SVRcor_Data_IV <- SVR_cor(X_Dat.intraCor.mer, samp.info = X_info.mer, 
#                           mzs.list = mzs.list)

# add noise to QC
GenDatNoise_SRformat(X_Dat.log.proc, dat.name = 'Data IV', 
                     samp.info = X_info.mer, feature.list = mzs.list)
runApp() # download result file in the same directory

SRformat_filePath0 <- paste0(getwd(), '/SRformat_file_addNoise_Data IV/')
SRcor_NoiDat.list <- list()
for (i in c(1: 3)) {
  zip.i.path <- paste0(SRformat_filePath0, 'SERRF Result', i, '.zip')
  dat.cor.i <- read.csv(unzip(zip.i.path, 'normalized by - SERRF.csv'), 
                        header = T)
  SRcor_NoiDat.list[[i]] <- as.matrix(t(dat.cor.i[, -1]))
}
unlink('normalized by - none.csv')
unlink('normalized by - SERRF.csv')
unlink('QC-RSDs.csv')
unlink('Bar Plot and PCA plot.png')


dat.list <- list(Uncorrect = X_Dat.log.proc, 
                 RLSC = QRcor_Data_IV, 
                 SVR = SVRcor_Data_IV, 
                 SERRF = SRcor_Data_IV, 
                 CordBat = CdBcor_Data_IV)

# CmpWithQCmethod(dat.list, SRcor_NoiDat.list, batch = batch.mer, group = group, 
#                 save.path0 = paste0(savefile.path0, 'Eval 3/'))

CmpWithQCmethod(dat.list, SRcor_NoiDat.list, batch = batch.mer, group = group, 
                PCA.print = T, MahaDist.print = F, 
                CumRSD.sbj.print = F, CumRSD.qc.print = F,
                save.path0 = paste0(savefile.path0, 'Eval 3 - 1/'))


# -----------------------------
# Preserving biological effect
# -----------------------------
Sbj_Dat.log.proc <- X_Dat.log.proc[group != 'QC', ]
Sbj_Dat.intraCor.mer <- X_Dat.intraCor.mer[group != 'QC', ]
Sbj_Dat.intraCor.mer1 <- X_Dat.intraCor.mer1[group != 'QC', ]
batch.sbj.mer <- batch.mer[group != 'QC']
batch.sbj <- batch[group != 'QC']
group.sbj <- group[group != 'QC']

# CordBat
# getRefbat(Sbj_Dat.intraCor.mer1, batch.sbj.mer)
CdBDat <- CordBat(Sbj_Dat.intraCor.mer, batch.sbj.mer, group = group.sbj, 
                  grouping = T, ref.batch = 5, eps = 1e-4, print.detail = F)
CdBcor_Data_IV.sbj <- CdBDat$X.cor
write.csv(CdBcor_Data_IV.sbj, 'CdBcor_Data_IV_6_sbj_LESSR_ref5_mer.csv', row.names = F)
CdBcor_Data_IV.sbj <- as.matrix(read.csv('CdBcor_Data_IV_6_sbj_LESSR_ref5_mer.csv', header = T))

# CdBDat <- CordBat(Sbj_Dat.intraCor.mer1, batch.sbj.mer, group = group.sbj, 
#                   grouping = T, ref.batch = 5, eps = 1e-4, print.detail = F)
# CdBcor_Data_IV.sbj <- CdBDat$X.cor
# write.csv(CdBcor_Data_IV.sbj, 'CdBcor_Data_IV_6_sbj_LES_ref5_mer.csv', row.names = F)
# CdBcor_Data_IV.sbj <- as.matrix(read.csv('CdBcor_Data_IV_6_sbj_LES_ref5_mer.csv', header = T))

# ComBat BEC
library(sva)
sbj.mod <- model.matrix(~as.factor(group.sbj))
CBcor_Data_IV <- ComBat(dat = t(Sbj_Dat.intraCor.mer), batch = batch.sbj.mer, 
                        mod = sbj.mod, par.prior = T)
CBcor_Data_IV <- t(CBcor_Data_IV)
colnames(CBcor_Data_IV) <- c(1: ncol(CBcor_Data_IV))

# FAbatch BEC
library(bapred)
FADat <- fabatch(Sbj_Dat.log.proc, as.factor(group.sbj), as.factor(batch.sbj.mer))
FAcor_Data_IV <- FADat$xadj
colnames(FAcor_Data_IV) <- c(1: ncol(FAcor_Data_IV))

# dat.list1 <- list(Uncorrect = Sbj_Dat.intraCor.mer, 
#                   ComBat = CBcor_Data_IV, 
#                   FAbatch = FAcor_Data_IV, 
#                   RLSC = QRcor_Data_IV[group != 'QC', ], 
#                   SERRF = SRcor_Data_IV[group != 'QC', ], 
#                   CordBat = CdBcor_Data_IV[group != 'QC', ])

QRcor_Data_IV.sbj <- QRcor_Data_IV[group != 'QC', ]
SRcor_Data_IV.sbj <- SRcor_Data_IV[group != 'QC', ]

dat.list1 <- list(Uncorrect = Sbj_Dat.log.proc, 
                  ComBat = CBcor_Data_IV, 
                  FAbatch = FAcor_Data_IV, 
                  RLSC = QRcor_Data_IV.sbj, 
                  SERRF = SRcor_Data_IV.sbj, 
                  CordBat = CdBcor_Data_IV.sbj)

PresBioEffect(dat.list1, dat.intraCor = Sbj_Dat.intraCor.mer, 
              batch = batch.sbj.mer, group = group.sbj, 
              selBioM.type = c(1, 1), BH.adj = c(F, F), min.batch.num = c(2, 2), 
              save.path0 = paste0(savefile.path0, 'Eval 4/'))




# -----------------------------------
# Principal component analysis (PCA)
# -----------------------------------
library(mixOmics)
# Before correction
pca.bef <- pca(X_Dat.intraCor, ncomp = 2, center = TRUE, scale = TRUE)
pca.bef.varX <- pca.bef$variates$X
pca.bef.expl <- pca.bef$prop_expl_var$X

PCA_scatter(data = pca.bef.varX, batch = batch, 
            expl.var = pca.bef.expl,
            xlim = c(min(pca.bef.varX[, 1]) - 1, max(pca.bef.varX[, 1]) + 1), 
            ylim = c(min(pca.bef.varX[, 2]) - 1, max(pca.bef.varX[, 2]) + 1), 
            legend.pos = 'right')

# maha dist in PC space
bef.md <- rep(0, batch.num * (batch.num - 1))
l <- 0
for (i in c(1: batch.num)) {
  PCdat.bati <- pca.bef.varX[Sbj_Bat == i, ]
  PCdat.bati.m <- colMeans(PCdat.bati)
  PCdat.bati.s <- cov(PCdat.bati)
  for (j in c(1: batch.num)) {
    if (j != i) {
      l <- l + 1
      PCdat.batj <- pca.bef.varX[Sbj_Bat == j, ]
      bef.md[l] <- mean(mahalanobis(PCdat.batj, PCdat.bati.m, PCdat.bati.s))
    }
  }
}
bef.md.mean <- mean(bef.md)

# CordBat correction
pca.cor <- pca(CdBcor_Data_IV, ncomp = 2, center = TRUE, scale = TRUE)
pca.cor.varX <- pca.cor$variates$X
pca.cor.expl <- pca.cor$prop_expl_var$X

PCA_scatter(data = pca.cor.varX, batch = batch, 
            expl.var = pca.cor.expl,
            xlim = c(min(pca.cor.varX[, 1]) - 1, max(pca.cor.varX[, 1]) + 1), 
            ylim = c(min(pca.cor.varX[, 2]) - 1, max(pca.cor.varX[, 2]) + 1), 
            legend.pos = 'right')

PCA_scatter_withQC(data = pca.cor.varX, batch = batch.mer, group = group, 
            QC_label = 'QC', expl.var = pca.cor.expl,
            xlim = c(min(pca.cor.varX[, 1]) - 1, max(pca.cor.varX[, 1]) + 1), 
            ylim = c(min(pca.cor.varX[, 2]) - 1, max(pca.cor.varX[, 2]) + 1))

# maha dist in PC space
cor.md <- rep(0, batch.num * (batch.num - 1))
l <- 0
for (i in c(1: batch.num)) {
  PCdat.bati <- pca.cor.varX[Sbj_Bat == i, ]
  PCdat.bati.m <- colMeans(PCdat.bati)
  PCdat.bati.s <- cov(PCdat.bati)
  for (j in c(1: batch.num)) {
    if (j != i) {
      l <- l + 1
      PCdat.batj <- pca.cor.varX[Sbj_Bat == j, ]
      cor.md[l] <- mean(mahalanobis(PCdat.batj, PCdat.bati.m, PCdat.bati.s))
    }
  }
}
cor.md.mean <- mean(cor.md)

# ---------------------------------
# Two-way (batch & group) ANOVA 
# ----------------------------------
library(car)
Sbj.Bat <- as.factor(Sbj_Bat)
Sbj.Grp <- as.factor(Sbj_Grp)

bef.Sum_sq <- data.frame(BE = rep(0, p))
bef.Fval <- data.frame(BE = rep(0, p))

cor.Sum_sq <- data.frame(BE = rep(0, p))
cor.Fval <- data.frame(BE = rep(0, p))

for (i in c(1: p)) {
  # fit model before correction
  bef.mod.i <- lm(Sbj_Dat.intraCor[, i] ~ Sbj.Bat * Sbj.Grp)
  bef.aov.i <- Anova(bef.mod.i, type = 3)
  bef.Sum_sq[i, ] <- bef.aov.i$`Sum Sq`[2]
  bef.Fval[i, ] <- bef.aov.i$`F value`[2]
  
  # fit model after CordBat correction
  cor.mod.i <- lm(CdBcor_Data_I[, i] ~ Sbj.Bat * Sbj.Grp)
  cor.aov.i <- Anova(cor.mod.i, type = 3)
  cor.Sum_sq[i, ] <- cor.aov.i$`Sum Sq`[2]
  cor.Fval[i, ] <- cor.aov.i$`F value`[2]
  
  cat(i, '\n')
}

bef.BEcut <- bef.Fval$BE
cor.BEcut <- cor.Fval$BE
cor.BEcut <- cor.BEcut[order(bef.BEcut)]
bef.BEcut <- sort(bef.BEcut)

cutoff <- 100
bef.BEcut[bef.BEcut > cutoff] <- cutoff
cor.BEcut[cor.BEcut > cutoff] <- cutoff

dat.BE <- data.frame(dat.name = factor(c(rep('Uncorrect', p), 
                                         rep('CordBat', p)), 
                                       levels = c('Uncorrect', 'CordBat')), 
                     met.id = rep(c(1: p), 2),  
                     Fval = c(bef.BEcut, cor.BEcut))

ggplot(data = dat.BE, aes(x = dat.name, y = met.id, fill = Fval)) + 
  geom_raster() + theme_bw() + 
  scale_fill_gradientn(colors = c('#39489f', '#39bbec', '#f9ed36', '#f38466', '#b81f25')) + 
  scale_y_continuous(breaks = c(53, seq(44, 4, -10)), labels = c(1, seq(10, 50, 10))) +
  scale_x_discrete(expand = c(0, 0)) + 
  theme(axis.text.x = element_text(size = 18, family = 'sans', face = 'bold', vjust = 5), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 22, family = 'sans', face = 'bold', vjust = 2), 
        axis.text.y = element_text(size = 18, family = 'sans', face = 'bold'),
        axis.ticks = element_blank(),
        legend.position = 'top',
        legend.title = element_text(size = 18, family = 'sans', face = 'bold'), 
        legend.text = element_text(size = 12, family = 'sans', face = 'bold'), 
        legend.box.spacing = unit(-0.5, 'cm'), 
        legend.background = element_rect(fill = 'transparent'), 
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 25, 
                                  family = 'sans', face = 'bold')) + 
  labs(y = 'Sorted Index of Metabolites', fill = 'F value', title = NULL)


# ---------------------
# DE metabolites
# ---------------------
DEmetIdx <- DEana(Sbj_Dat.intraCor, Sbj_Bat, top.rank = 1)

DEmetIdx2 <- c(37, 47)
dat.bef <- Sbj_Dat.intraCor[, DEmetIdx2]
colnames(dat.bef) <- colnames(Sbj_Dat.log.proc)[DEmetIdx2]
dat.cor <- CdBcor_Data_I[, DEmetIdx2]
colnames(dat.cor) <- colnames(Sbj_Dat.log.proc)[DEmetIdx2]
DEmet_plot_merge(dat.bef = dat.bef, 
                 dat.cor = dat.cor, 
                 batch = Sbj_Bat, InjOrd = Sbj_info$order)


# -----------------
# Data structures
# ------------------
stdPcc <- getstdPcc(Sbj_Dat.delout, Sbj_Bat, Notdiff.ratio = 0.8)
stdPcc <- stdPcc$std.PCC.Mat

Dat.bef <- Sbj_Dat.intraCor
colnames(Dat.bef) <- c(1: ncol(Dat.bef))
Dat.cor <- CdBcor_Data_I
colnames(Dat.cor) <- c(1: ncol(Dat.cor))

conf <- 0.05
# conf <- 0.01
# conf <- 0.1

P.bef <- getPdiff2Std(stdPcc, Dat.bef, conf_alpha = conf)
P.bef <- P.bef$P.Mat

P.QF <- getPdiff2Std(stdPcc, Dat.cor, conf_alpha = conf)
P.QF <- P.QF$P.Mat

bef_vs_std.largestComp <- getLargestComp(P.bef, conf_alpha = conf)
bef_vs_std.Pcc <- bef_vs_std.largestComp$P.connectComp
Graph_plot_forP(bef_vs_std.Pcc)

QF_vs_std.largestComp <- getLargestComp(P.QF, conf_alpha = conf)
QF_vs_std.Pcc <- QF_vs_std.largestComp$P.connectComp
Graph_plot_forP(QF_vs_std.Pcc)



#################################################################
#########################################################################
############################################################################
# Untargeted urine LC-HRMS metabolomics data












