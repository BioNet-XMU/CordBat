# Title: Data III BEC and evaluation
# Author: Fanjing Guo
# Date: 2022.12.28

# --------------------- Data pre-processing and batch effect correction ------------------- #
# UHPLC-QTOF/MS metabolomcis dataset from an adenocarcinoma study
X <- read.csv("Data III_raw.csv", header = TRUE)
X_info <- X[, c(1: 5)]
X_Dat <- as.matrix(X[, -c(1: 5)])
N <- nrow(X_Dat)
p <- ncol(X_Dat)

# log2-transformed (preprocessing in Matlab)
X_Dat.log <- log2(X_Dat)
batch <- X_info$batch
group <- X_info$group
samp.order <- X_info$injection.order
X_Dat.log.proc <- X_Dat.log
met.list <- colnames(X_Dat.log.proc)
N <- nrow(X_Dat.log.proc)
p <- ncol(X_Dat.log.proc)

colnames(X_info) <- c('SampName', 'RunOrder', 'Class', 'Batch', 'Group')

# # intra-batch correction (Loess line using subject samples)
# Sbj_info <- X_info[X_Grp != 3, ]
# Sbj_Dat.log.proc <- X_Dat.log.proc[X_Grp != 3, ]
# p <- ncol(Sbj_Dat.log.proc)
# Sbj_Bat <- X_Bat[X_Grp != 3]
# Sbj_Bat.f <- as.factor(Sbj_Bat)
# Sbj_Grp <- X_Grp[X_Grp != 3]
# batch.num <- nlevels(Sbj_Bat.f)
# batch.lev <- levels(Sbj_Bat.f)
# Sbj_Dat.intraCor <- Sbj_Dat.log.proc
# for (i in c(1: batch.num)) {
#   dat.i <- Sbj_Dat.log.proc[Sbj_Bat.f == batch.lev[i], ]
#   ord.i <- Sbj_info$injection.order[Sbj_Bat.f == batch.lev[i]]
#   
#   dat.i.cor <- dat.i
#   for (j in c(1: p)) {
#     dat <- data.frame(x = ord.i, 
#                       y = dat.i[, j])
#     fitted.result <- loess(y ~ x, data = dat, span = 0.75)
#     fit.metj <- fitted.result$fitted
#     y.cor <- dat.i[, j] + median(dat.i[, j]) - fit.metj
#     
#     dat.i.cor[, j] <- y.cor
#   }
#   Sbj_Dat.intraCor[Sbj_Bat.f == batch.lev[i], ] <- dat.i.cor
#   cat(i, '\n')
# }
# 
# Sbj_Dat.intraCor <- as.matrix(Sbj_Dat.intraCor)
# colnames(Sbj_Dat.intraCor) <- c(1: ncol(Sbj_Dat.intraCor))

# intra-batch correction 1 - LOESS
bat.f <- as.factor(batch)
batch.num <- nlevels(bat.f)
batch.lev <- levels(bat.f)
X_Dat.intraCor1 <- X_Dat.log.proc
for (i in c(1: batch.num)) {
  dat.i <- X_Dat.log.proc[batch == batch.lev[i], ]
  ord.i <- samp.order[batch == batch.lev[i]]
  
  dat.i.cor <- dat.i
  for (j in c(1: p)) {
    dat <- data.frame(x = ord.i, 
                      y = dat.i[, j])
    fitted.result <- loess(y ~ x, data = dat, span = 0.75)
    fit.metj <- fitted.result$fitted
    y.cor <- dat.i[, j] + median(dat.i[, j]) - fit.metj
    
    dat.i.cor[, j] <- y.cor
  }
  X_Dat.intraCor1[batch == batch.lev[i], ] <- dat.i.cor
  cat(i, '\n')
}

X_Dat.intraCor1 <- as.matrix(X_Dat.intraCor1)
colnames(X_Dat.intraCor1) <- c(1: ncol(X_Dat.intraCor1))

# SERRF
X_Dat.intraCor <- X_Dat.intraCor1
SRformat_filePath0 <- paste0(getwd(), '/SRformat_file_Data III_eachBat/')
for (i in c(1: batch.num)) {
  dat.bati <- X_Dat.intraCor1[batch == batch.lev[i], ]
  bati.info <- X_info[batch == batch.lev[i], ]
  GenSRformatdat(dat.bati, bati.info, met.list, QC_label = 3, 
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
  X_Dat.intraCor[batch == i, ] <- dat.bati.cor
}

write.csv(X_Dat.intraCor, 'Data_III_LES&SRintraCor_ref2.csv', row.names = F)

# Batch effect correction
Sbj_Dat.intraVor <- X_Dat.intraCor[group != 3, ]
Sbj_Dat.log.proc <- X_Dat.log.proc[group != 3, ]
batch.sbj <- batch[group != 3]
group.sbj <- group[group != 3]

# CordBat BEC
# with QC
getRefbat(X_Dat.intraCor, batch)
CdBDat <- CordBat(X_Dat.intraCor, batch, ref.batch = 2, eps = 1e-3, print.detail = F)
CdBcor_Data_III <- CdBDat$X.cor
CdBcor_Data_III.sbj <- CdBcor_Data_III[group != 3, ]

write.csv(CdBcor_Data_III, 'CdBcor_Data_III_LESSR_ref2.csv', row.names = F)

# # without QC
# getRefbat(Sbj_Dat.intraCor, batch)
# CdBDat <- CordBat(Sbj_Dat.intraCor, batch.sbj, group = group.sbj, 
#                   grouping = T, ref.batch = 1, eps = 1e-3, print.detail = F)
# CdBcor_Data_III.sbj <- CdBDat$X.cor

# CdBcor_Data_III.sbj <- as.matrix(read.csv('CdBcor_Data_III_ref1.csv', header = T))

# ComBat BEC
library(sva)
X_Mod <- model.matrix(~as.factor(group.sbj))
Sbj_Dat.intraCor <- t(Sbj_Dat.intraCor)
CBcor_Data_III <- ComBat(dat = Sbj_Dat.intraCor, batch = batch.sbj, 
                         mod = X_Mod, par.prior = T)
CBcor_Data_III <- t(CBcor_Data_III)
colnames(CBcor_Data_III) <- c(1: ncol(CBcor_Data_III))
Sbj_Dat.intraCor <- t(Sbj_Dat.intraCor)

# FAbatch BEC
library(bapred)
FADat <- fabatch(Sbj_Dat.log.proc, as.factor(batch.sbj), as.factor(group.sbj))
FAcor_Data_III <- FADat$xadj
colnames(FAcor_Data_III) <- c(1: ncol(FAcor_Data_III))


# QC-based methods
# QC-RLSC BEC
library(statTarget)
QRcor_Data_III <- QCRLSC_cor(X_Dat.log.proc, samp.info = X_info, 
                             QC_label = 3, 
                             feature.list = met.list)

# SVR
library(MetNormalizer)
SVRcor_Data_III <- SVR_cor(X_Dat.log.proc, samp.info = X_info, 
                           QC_label = 3, 
                           mets.name = met.list)

# SERRF BEC
GenSRformatdat(X_Dat.log.proc, X_info, met.list, QC_label = 3, 
               wrfile.path = 'SRformat_Data_III.csv')
source('app.R')
runApp()
SRcor_Data_III <- read.csv(unzip(paste0(getwd(), '/SERRF Result.zip'), 
                                'normalized by - SERRF.csv'), header = T)
SRcor_Data_III <- as.matrix(t(SRcor_Data_III[, -1]))
unlink('SERRF Result.zip')

SRcor_Data_III <- read.csv("SERRFDat_real_1_ordered3.csv", header = TRUE)
SRcor_Data_III <- as.matrix(t(SRcor_Data_III[, -1]))
SRcor_Data_III <- log2(SRcor_Data_III)
colnames(SRcor_Data_III) <- c(1: ncol(SRcor_Data_III))
rownames(SRcor_Data_III) <- c(1: nrow(SRcor_Data_III))


# add noise to QC
GenDatNoise_SRformat(X_Dat.log.proc, dat.name = 'Data III', QC_label = 3, 
                     samp.info = X_info, feature.list = met.list)
runApp() # download result file in the same directory

SRformat_filePath0 <- paste0(getwd(), '/SRformat_file_addNoise_Data III/')
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

# ----------------------------------------- Evaluation ------------------------------------- #

savefile.path0 <- 'F:/Graduate/Graduate Studies/Papers/1-GGM for Batch Effect Correction/New Results/Data III/'

# --------------------------------
# Preservation of data structure
# --------------------------------
delout.res <- DelOut_Imp(X_Dat.intraCor, batch)
X_Dat.delout <- delout.res$dat.delout
batch.delout <- delout.res$batch.new
PresDatStruct(X_Dat.intraCor, X_Dat.delout, CdBcor_Data_III, batch.delout, 
              NotdifBat.ratio = 0.9, delays = c(0.2, 0.5), 
              save.path0 = paste0(savefile.path0, 'Eval 2/'))

# ------------------------------
# Compared with QC-based method
# ------------------------------
dat.list <- list(Uncorrect = X_Dat.log.proc, 
                 RLSC = QRcor_Data_III, 
                 SVR = SVRcor_Data_III, 
                 SERRF = SRcor_Data_III, 
                 CordBat = CdBcor_Data_III)

set_arrow1 <- data.frame(arrow.pos = c('left', 'left', 'topleft'), 
                         arrow.ang = c(0, 0, pi/2.4), 
                         arrow.len = c(7, 9, 15))

set_arrow2 <- data.frame(arrow.pos = c('bottomleft', 'left', 'right'), 
                         arrow.ang = c(pi/2.5, 0, 0), 
                         arrow.len = c(10, 8, 8))

CmpWithQCmethod(dat.list, SRcor_NoiDat.list, batch = batch, group = group, 
                QC_label = 3, pooledQC = F, 
                PCA.print = F, MahaDist.print = T, 
                CumRSD.sbj.print = F, CumRSD.qc.print = F, 
                mark.rsd.sbj = 30, mark.rsd.qc = 20, 
                set_arrow.sbj = set_arrow1, set_arrow.qc = set_arrow2, 
                save.path0 = paste0(savefile.path0, 'Eval 3/'))

# -----------------------------------
# Principal component analysis (PCA)
# -----------------------------------
library(mixOmics)

Sbj_Grp <- X_Grp[X_Grp != 3]
Sbj_Grp[Sbj_Grp == 1] <- 'CRC'
Sbj_Grp[Sbj_Grp == 2] <- 'CE'
Sbj_Grp.f <- factor(Sbj_Grp, levels = c('CRC', 'CE'))

# Before BEC
pca.bef <- pca(Sbj_Dat.log.proc, ncomp = 2, center = T, scale = T)
pca.bef.expl <- pca.bef$prop_expl_var$X
pca.bef.varX <- pca.bef$variates$X

PCA_scatter1(data = pca.bef.varX, group = Sbj_Grp.f, 
             expl.var = pca.bef.expl,
             xlim = c(min(pca.bef.varX[, 1]) - 1, max(pca.bef.varX[, 1]) + 1),
             ylim = c(min(pca.bef.varX[, 2]) - 1, max(pca.bef.varX[, 2]) + 1),
             title = NULL, legend.pos = 'none')

# CordBat correction
pca.CdBcor <- pca(CdBcor_Data_III, ncomp = 2, center = T, scale = T)
pca.CdBcor.expl <- pca.CdBcor$prop_expl_var$X
pca.CdBcor.varX <- pca.CdBcor$variates$X

PCA_scatter1(data = pca.CdBcor.varX, group = Sbj_Grp.f, 
             expl.var = pca.CdBcor.expl,
             xlim = c(min(pca.CdBcor.varX[, 1]) - 1, max(pca.CdBcor.varX[, 1]) + 1), 
             ylim = c(min(pca.CdBcor.varX[, 2]) - 1, max(pca.CdBcor.varX[, 2]) + 1),
             title = NULL)

# ComBat correction
pca.CBcor <- pca(CBcor_Data_III, ncomp = 2, center = T, scale = T)
pca.CBcor.expl <- pca.CBcor$prop_expl_var$X
pca.CBcor.varX <- pca.CBcor$variates$X

PCA_scatter1(data = pca.CBcor.varX, group = Sbj_Grp.f, 
             expl.var = pca.CBcor.expl,
             xlim = c(min(pca.CBcor.varX[, 1]) - 1, max(pca.CBcor.varX[, 1]) + 1), 
             ylim = c(min(pca.CBcor.varX[, 2]) - 1, max(pca.CBcor.varX[, 2]) + 1),
             title = NULL, legend.pos = 'none')

# FAbatch correction
pca.FAcor <- pca(FAcor_Data_III, ncomp = 2, center = T, scale = T)
pca.FAcor.expl <- pca.FAcor$prop_expl_var$X
pca.FAcor.varX <- pca.FAcor$variates$X

PCA_scatter1(data = pca.FAcor.varX, group = Sbj_Grp.f, 
             expl.var = pca.FAcor.expl,
             xlim = c(min(pca.FAcor.varX[, 1]) - 1, max(pca.FAcor.varX[, 1]) + 1), 
             ylim = c(min(pca.FAcor.varX[, 2]) - 1, max(pca.FAcor.varX[, 2]) + 1),
             title = NULL)

# SERRF correction
pca.SRcor <- pca(SRcor_Data_III[X_Grp != 3, ], ncomp = 2, center = T, scale = T)
pca.SRcor.expl <- pca.SRcor$prop_expl_var$X
pca.SRcor.varX <- pca.SRcor$variates$X

PCA_scatter1(data = pca.SRcor.varX, group = Sbj_Grp.f, 
             expl.var = pca.cor.expl,
             xlim = c(min(pca.SRcor.varX[, 1]) - 1, max(pca.SRcor.varX[, 1]) + 1), 
             ylim = c(min(pca.SRcor.varX[, 2]) - 1, max(pca.SRcor.varX[, 2]) + 1),
             title = NULL, legend.pos = 'none')

# QC-RLSC correction
pca.QRcor <- pca(QRcor_Data_III, ncomp = 2, center = T, scale = T)
pca.QRcor.expl <- pca.QRcor$prop_expl_var$X
pca.QRcor.varX <- pca.QRcor$variates$X

PCA_scatter1(data = pca.QRcor.varX, group = Sbj_Grp.f, 
             expl.var = pca.QRcor.expl,
             xlim = c(min(pca.QRcor.varX[, 1]) - 1, max(pca.QRcor.varX[, 1]) + 1), 
             ylim = c(min(pca.QRcor.varX[, 2]) - 1, max(pca.QRcor.varX[, 2]) + 1),
             title = NULL)


# --------------------
# Biological effect
# ---------------------
library(ropls)
library(pROC)

Sbj_Bat <- X_Bat[X_Grp != 3]
Sbj_Grp <- X_Grp[X_Grp != 3]

# Identification of potential biomarkers in each batch
batch.num_1 <- 3
DEmet_grp.VIP <- matrix(0, p, batch.num_1)
DEmet_grp.VIP.sortID <- matrix(0, p, batch.num_1)
DEmet_grp.Pval <- matrix(0, p, batch.num_1)
for (i in c(1: batch.num_1)) {
  Sbj_Dat.i <- Sbj_Dat.intraCor[Sbj_Bat == i, ]
  Sbj_Grp.i <- Sbj_Grp[Sbj_Bat == i]
  plsda.i <- opls(Sbj_Dat.i, Sbj_Grp.i, orthoI = 0, 
                  predI = 3, fig.pdfC = 'none', info.txtC = 'none')
  bati.VIP <- getVipVn(plsda.i)
  DEmet_grp.VIP[, i] <- bati.VIP
  DEmet_grp.VIP.sortID[, i] <- order(DEmet_grp.VIP[, i], decreasing = T)
  
  ttest.result <- getttestResult(Sbj_Dat.i[Sbj_Grp.i == 1, ], 
                                 Sbj_Dat.i[Sbj_Grp.i == 2, ])
  bef.ttest.pval <- ttest.result$tt.Pval
  DEmet_grp.Pval[, i] <- bef.ttest.pval
}
ave.Pval <- apply(DEmet_grp.Pval, MARGIN = 1, FUN = mean)

DEmet_pID1 <- which(DEmet_grp.Pval[, 1] < 0.1)
DEmet1 <- DEmet_pID1[order(DEmet_grp.VIP[DEmet_pID1, 1], decreasing = T)]
DEmet1 <- DEmet1[c(1: 200)]

DEmet_pID2 <- which(DEmet_grp.Pval[, 2] < 0.1)
DEmet2 <- DEmet_pID2[order(DEmet_grp.VIP[DEmet_pID2, 1], decreasing = T)]
DEmet2 <- DEmet2[c(1: 200)]

DEmet_pID3 <- which(DEmet_grp.Pval[, 3] < 0.1)
DEmet3 <- DEmet_pID3[order(DEmet_grp.VIP[DEmet_pID3, 1], decreasing = T)]
DEmet3 <- DEmet3[c(1: 200)]

DEmet1_3 <- unique(c(intersect(DEmet1, DEmet2), 
                     intersect(DEmet2, DEmet3), 
                     intersect(DEmet1, DEmet3)))
metId.sortScore <- DEmet1_3


# Identify potential biomarkers in whole data

# Before BEC
Sbj_Dat.bef <- X_Dat[X_Grp != 3, ]
plsda.bef <- opls(Sbj_Dat.bef, Sbj_Grp, orthoI = 0, 
                  predI = 3, fig.pdfC = 'none', info.txtC = 'none')
bef.VIP <- getVipVn(plsda.bef)
bef.VIP.sortID <- order(bef.VIP, decreasing = T)

ttest.result <- getttestResult(Sbj_Dat.bef[Sbj_Grp == 1, ], 
                               Sbj_Dat.bef[Sbj_Grp == 2, ])
bef.ttest.pval <- ttest.result$tt.Pval

bef.OLres <- getOLperc(metId.sortScore, bef.ttest.pval)
bef.OLnum <- bef.OLres$OLmetNum
bef.OLperc <- bef.OLres$OLperc

# CordBat correction
plsda.CdB <- opls(CdBcor_Data_III, Sbj_Grp, orthoI = 0, 
                 predI = 3, fig.pdfC = 'none', info.txtC = 'none')
CdB.VIP <- getVipVn(plsda.CdB)
CdB.VIP.sortID <- order(CdB.VIP, decreasing = T)

ttest.result <- getttestResult(CdB.Sbj_Dat[Sbj_Grp == 1, ], 
                               CdB.Sbj_Dat[Sbj_Grp == 2, ])
CdB.ttest.pval <- ttest.result$tt.Pval

CdB.OLres <- getOLperc(metId.sortScore, CdB.ttest.pval)
CdB.OLnum <- CdB.OLres$OLmetNum
CdB.OLperc <- CdB.OLres$OLperc

# ComBat correction
plsda.CB <- opls(CBcor_Data_III, Sbj_Grp, orthoI = 0, 
                 predI = 3, fig.pdfC = 'none', info.txtC = 'none')
CB.VIP <- getVipVn(plsda.CB)
CB.VIP.sortID <- order(CB.VIP, decreasing = T)

ttest.result <- getttestResult(CBcor_Data_III[Sbj_Grp == 1, ], 
                               CBcor_Data_III[Sbj_Grp == 2, ])
CB.ttest.pval <- ttest.result$tt.Pval

CB.OLres <- getOLperc(metId.sortScore, CB.ttest.pval)
CB.OLnum <- CB.OLres$OLmetNum
CB.OLperc <- CB.OLres$OLperc


# FAbatch correction
plsda.FA <- opls(FAcor_Data_III, Sbj_Grp, orthoI = 0, 
                 predI = 3, fig.pdfC = 'none', info.txtC = 'none')
FA.VIP <- getVipVn(plsda.FA)
FA.VIP.sortID <- order(FA.VIP, decreasing = T)

ttest.result <- getttestResult(FAcor_Data_III[Sbj_Grp == 1, ], 
                               FAcor_Data_III[Sbj_Grp == 2, ])
FA.ttest.pval <- ttest.result$tt.Pval

FA.OLres <- getOLperc(metId.sortScore, FA.ttest.pval)
FA.OLnum <- FA.OLres$OLmetNum
FA.OLperc <- FA.OLres$OLperc

# SERRF correction
SR.Sbj_Dat <- SRcor_Data_III[X_Grp != 3, ]
plsda.SR <- opls(SR.Sbj_Dat, Sbj_Grp, orthoI = 0, 
                 predI = 3, fig.pdfC = 'none', info.txtC = 'none')
SR.VIP <- getVipVn(plsda.SR)
SR.VIP.sortID <- order(SR.VIP, decreasing = T)

ttest.result <- getttestResult(SR.Sbj_Dat[Sbj_Grp == 1, ], 
                               SR.Sbj_Dat[Sbj_Grp == 2, ])
SR.ttest.pval <- ttest.result$tt.Pval

SR.OLres <- getOLperc(metId.sortScore, SR.ttest.pval)
SR.OLnum <- SR.OLres$OLmetNum
SR.OLperc <- SR.OLres$OLperc

# QC-RLSC correction
plsda.QR <- opls(QRcor_Data_III, Sbj_Grp, orthoI = 0, 
                 predI = 3, fig.pdfC = 'none', info.txtC = 'none')
QR.VIP <- rep(0, p)
QR.VIP[retain.met.id] <- getVipVn(plsda.QR)
QR.VIP.sortID <- order(QR.VIP, decreasing = T)

ttest.result <- getttestResult(QR.Sbj_Dat[Sbj_Grp == 1, ], 
                               QR.Sbj_Dat[Sbj_Grp == 2, ])
QR.ttest.pval <- rep(1, p)
QR.ttest.pval[retain.met.id] <- ttest.result$tt.Pval

QR.OLres <- getOLperc(metId.sortScore, QR.ttest.pval)
QR.OLnum <- QR.OLres$OLmetNum
QR.OLperc <- QR.OLres$OLperc


# plot 
# Number of shared metabolites between top n metabolites 
# in p-value list and potential biomarkers for each dataset
dattype <- factor(c(rep('Uncorrect', 60), 
                    rep('ComBat', 60), 
                    rep('FAbatch', 60), 
                    rep('SERRF', 60), 
                    rep('QC-RLSC', 60), 
                    rep('CordBat', 60)), 
                  levels = c('Uncorrect', 'ComBat', 'FAbatch', 'SERRF', 
                             'QC-RLSC', 'CordBat'))

dat <- data.frame(dattype = dattype, 
                  top.n = rep(seq(2, 120, 2), 6), 
                  OLnum = c(bef.OLnum, CB.OLnum, FA.OLnum, 
                            SR.OLnum, QR.OLnum, CdB.OLnum))

max.olnum <- max(c(bef.OLnum, CB.OLnum, FA.OLnum, 
                   SR.OLnum, QR.OLnum, CdB.OLnum))

fon <- 'sans'
set_col <- c('black', 'blue', 'yellow', 'cyan', 'magenta', 'red')
ggplot(data = dat, aes(x = top.n, y = OLnum, group = dattype, color = dattype)) + 
  geom_line(size = 1.5) + scale_color_manual(values = set_col) + 
  scale_x_continuous(breaks = seq(0, 120, 20), labels = seq(0, 120, 20), 
                     limits = c(0, (120 + 1))) + 
  scale_y_continuous(breaks = seq(0, 90, 10), labels = seq(0, 90, 10), 
                     limits = c(0, (max.olnum + 1))) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 18, family = fon, face = 'bold'), 
        axis.text.y = element_text(size = 18, family = fon, face = 'bold'), 
        axis.title.x = element_text(size = 25, family = fon, face = 'bold'), 
        axis.title.y = element_text(size = 25, family = fon, face = 'bold'), 
        panel.border = element_rect(fill = NA, size = 2), 
        legend.position = c(0.17, 0.75), 
        legend.text = element_text(size = 18, family = fon, face = 'bold')) + 
  xlab('Top n') + ylab('Number of overlapped metabolites') + 
  labs(group = NULL, color = NULL, title = NULL)


# Classification accuracy
# CV for SVM ensemble
library(e1071)
library(pROC)
Sbj_Dat.bef <- X_Dat.log.proc[X_Grp != 3, ]
Sbj_Grp <- X_Grp[X_Grp != 3]

bef.cvSVMens.roc <- CVforSVMensemble(Sbj_Dat.bef, Sbj_Grp)
CdB.cvSVMens.roc <- CVforSVMensemble(CdBcor_Data_III, Sbj_Grp)
CB.cvSVMens.roc <- CVforSVMensemble(CBcor_Data_III, Sbj_Grp)
FA.cvSVMens.roc <- CVforSVMensemble(FAcor_Data_III, Sbj_Grp)
SR.cvSVMens.roc <- CVforSVMensemble(SRcor_Data_III[X_Grp != 3, ], Sbj_Grp)
QR.cvSVMens.roc <- CVforSVMensemble(QRcor_Data_III, Sbj_Grp)


bef.roc <- data.frame(spec = rev(1 - bef.cvSVMens.roc$spec), 
                      sens = rev(bef.cvSVMens.roc$sens))
dat.bef <- rep('Uncorrect', nrow(bef.roc))
bef.auc <- bef.cvSVMens.roc$auc.val
bef.CVerr <- bef.cvSVMens.roc$cv.err
bef.acc <- bef.cvSVMens.roc$accuracy

CdB.roc <- data.frame(spec = rev(1 - CdB.cvSVMens.roc$spec), 
                     sens = rev(CdB.cvSVMens.roc$sens))
dat.CdB <- rep('CordBat', nrow(CdB.roc))
CdB.auc <- CdB.cvSVMens.roc$auc.val
CdB.CVerr <- CdB.cvSVMens.roc$cv.err
CdB.acc <- CdB.cvSVMens.roc$accuracy

CB.roc <- data.frame(spec = rev(1 - CB.cvSVMens.roc$spec), 
                     sens = rev(CB.cvSVMens.roc$sens))
dat.CB <- rep('ComBat', nrow(CB.roc))
CB.auc <- CB.cvSVMens.roc$auc.val
CB.CVerr <- CB.cvSVMens.roc$cv.err
CB.acc <- CB.cvSVMens.roc$accuracy

FA.roc <- data.frame(spec = rev(1 - FA.cvSVMens.roc$spec), 
                     sens = rev(FA.cvSVMens.roc$sens))
dat.FA <- rep('FAbatch', nrow(FA.roc))
FA.auc <- FA.cvSVMens.roc$auc.val
FA.CVerr <- FA.cvSVMens.roc$cv.err
FA.acc <- FA.cvSVMens.roc$accuracy

SR.roc <- data.frame(spec = rev(1 - SR.cvSVMens.roc$spec), 
                     sens = rev(SR.cvSVMens.roc$sens))
dat.SR <- rep('SERRF', nrow(SR.roc))
SR.auc <- SR.cvSVMens.roc$auc.val
SR.CVerr <- SR.cvSVMens.roc$cv.err
SR.acc <- SR.cvSVMens.roc$accuracy

QR.roc <- data.frame(spec = rev(1 - QR.cvSVMens.roc$spec), 
                     sens = rev(QR.cvSVMens.roc$sens))
dat.QR <- rep('QC-RLSC', nrow(QR.roc))
QR.auc <- QR.cvSVMens.roc$auc.val
QR.CVerr <- QR.cvSVMens.roc$cv.err
QR.acc <- QR.cvSVMens.roc$accuracy


# plot
# Prediction accuracies of the SVM models
acc.val <- data.frame(m = c(mean(bef.acc), mean(CdB.acc), mean(CB.acc),
                            mean(FA.acc), mean(SR.acc), mean(QR.acc)),
                      sd = c(sd(bef.acc), sd(CdB.acc), sd(CB.acc),
                             sd(FA.acc), sd(SR.acc), sd(QR.acc)))
rownames(acc.val) <- c('Uncorrect', 'CordBat', 'ComBat', 'FAbatch', 'SERRF', 'QC-RLSC')

library(ggbreak)
dat.acc <- data.frame(dat.type = factor(rownames(acc.val), 
                                        levels = c('CordBat', 'ComBat', 'SERRF', 'QC-RLSC', 'Uncorrect', 'FAbatch')),
                      acc.mean = acc.val$m, 
                      acc.sd = acc.val$sd)


fon <- 'sans'
max.acc <- max(c(dat.acc$acc.mean - dat.acc$acc.sd, dat.acc$acc.mean + dat.acc$acc.sd))
min.acc <- min(c(dat.acc$acc.mean - dat.acc$acc.sd, dat.acc$acc.mean + dat.acc$acc.sd))
ggplot(dat.acc, aes(fill = dat.type, y = dat.type, x = acc.mean)) +
  geom_bar(position = position_dodge(preserve = 'single'), stat = 'identity',
           width = 0.7, colour = 'black', size = 1.5) +
  geom_errorbar(aes(xmin = acc.mean - acc.sd, xmax = acc.mean + acc.sd),
                position = position_dodge2(preserve = 'single', padding = 0.5),
                size = 1.5, width = 0.3) +
  geom_vline(xintercept = mean(QF.acc), linetype = 2, size = 1.2) + 
  scale_y_discrete(limits = rev, expand = expansion(0.16)) +
  scale_x_break(c(0.05, round(min.acc - 0.0001, 4)), 
                ticklabels = seq(round(min.acc - 0.001, 2), round(max.acc + 0.001, 2), 0.01),
                scales = 10) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, round(max.acc + 0.002, 3)),
                     breaks = c(0, seq(round(min.acc - 0.001, 2), round(max.acc + 0.001, 2), 0.01))) +
  scale_fill_npg() +
  theme_bw() +
  theme(axis.text = element_text(size = 18, family = fon, face = 'bold'),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 25, family = fon, face = 'bold'),
        panel.border = element_blank(),
        axis.line = element_line(size = 2),
        legend.position = 'right',
        legend.direction = 'vertical',
        legend.title = element_text(size = 20, family = fon, face = 'bold'),
        legend.text = element_text(size = 18, family = fon, face = 'bold')) +
  xlab('Prediction Accuracy') + 
  labs(fill = 'Data', title = NULL)


# ROC curves for the SVM models
dat.ty <- factor(c(dat.bef, dat.CdB, dat.CB, dat.SR), 
                 levels = c('Uncorrect', 'ComBat', 'SERRF', 'CordBat'))
dat.roc <- data.frame(spec = c(bef.roc$spec, CdB.roc$spec, CB.roc$spec,
                               SR.roc$spec),
                      sens = c(bef.roc$sens, CdB.roc$sens, CB.roc$sens,
                               SR.roc$sens))

fon <- 'sans'
set_col <- c('black', 'green', 'blue', 'red')
ggplot(data = dat.roc, aes(x = spec, y = sens, group = dat.ty, color = dat.ty)) +
  geom_line(size = 1.5) +
  geom_text(aes(x = 0.76, y = 0.15),
            label = paste0("AUC: \n","Uncorrect: ", round(raw.auc, 4), "\n",
                           "ComBat: ", round(CB.auc, 4), "\n",
                           "SERRF: ", round(SR.auc, 4), "\n",
                           "CordBat: ", round(QF.auc, 4), "\n"),
            color = "black", size = 9) +
  xlab('1 - specificity') + ylab('sensitivity') + xlim(0, 1) + ylim(0, 1) +
  scale_color_manual(values = set_col) + theme_bw() +
  theme(axis.title = element_text(size = 25, family = fon, face = 'bold'),
        axis.text = element_text(size = 25, family = fon, face = 'bold'),
        legend.title = element_text(size = 25, family = fon, face = 'bold'),
        legend.text = element_text(size = 25, family = fon, face = 'bold'),
        plot.title = element_text(size = 25, family = fon, face = 'bold', hjust = 0.5),
        panel.border = element_rect(size = 2)) +
  labs(color = 'Method', title = 'ROC for CV SVM classification')



