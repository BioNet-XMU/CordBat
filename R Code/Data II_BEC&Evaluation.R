# Title: Data II BEC and evaluation
# Author: Fanjing Guo
# Date: 2022.12.28

# --------------------- Data pre-processing and batch effect correction ------------------- #
# LC-TQMS metabolomics dataset from a trauma study
X <- read.csv("Data II_Raw.csv", header = TRUE)
X_info <- X[, c(1: 4)]
X_Dat <- as.matrix(X[, -c(1: 4)])
N <- nrow(X_Dat)
p <- ncol(X_Dat)

# pre-processing
X_Dat[X_Dat == 'N/A'] <- NA
X_Dat <- apply(X_Dat, MARGIN = 2, FUN = as.numeric)
# delete metabolites with >= 30% missing values
na.num <- colSums(is.na(X_Dat))
delmetIdx <- which(na.num / N >= 0.5)
if (length(delmetIdx) != 0) {
  X_Dat.proc <- X_Dat[, -delmetIdx]
}
p <- ncol(X_Dat.proc)

# impute missing values using knn
library(DMwR2)
X_Dat.proc <- as.data.frame(X_Dat.proc)
X_Dat.proc <- knnImputation(X_Dat.proc)

# detect and delete outliers
library(mixOmics)
X_Dat.log.proc <- log2(X_Dat.proc)
pca.dat <- pca(X_Dat.log.proc, ncomp = 2, center = TRUE, scale = TRUE)
pca.dat.varX <- pca.dat$variates$X
delsampIdx <- c()
for (i in c(1: 2)) {
  pc.i <- pca.dat.varX[, i]
  pc.i.m <- mean(pc.i)
  pc.i.sd <- sd(pc.i)
  pc.i.min <- pc.i.m - 3 * pc.i.sd
  pc.i.max <- pc.i.m + 3 * pc.i.sd
  delsampIdx <- c(delsampIdx, which(pc.i < pc.i.min | pc.i > pc.i.max))
}
delsampIdx <- unique(delsampIdx)

X_Dat.proc <- X_Dat.proc[-delsampIdx, ]
X_Dat.log.proc <- X_Dat.log.proc[-delsampIdx, ]
X_info <- X_info[-delsampIdx, ]

X_Dat.proc <- as.matrix(X_Dat.proc)
X_Dat.log.proc <- as.matrix(X_Dat.log.proc)


X_Bat <- X_info$Batch
X_Grp <- X_info$Group

N <- nrow(X_Dat.log.proc)
p <- ncol(X_Dat.log.proc)
X_Dat.log.proc <- as.matrix(X_Dat.log.proc)

# merge adjacent 4 batches
# -------- batch 1-4: new batch 1 ------------
# -------- batch 5-8: new batch 2 ------------
# -------- batch 9-12: new batch 3 -----------
# -------- batch 13-16: new batch 4 ----------
X_Bat.new <- rep(0, length(X_Bat))
X_Bat.new[X_Bat == 1 | X_Bat == 2 | X_Bat == 3 | X_Bat == 4] <- 1
X_Bat.new[X_Bat == 5 | X_Bat == 6 | X_Bat == 7 | X_Bat == 8] <- 2
X_Bat.new[X_Bat == 9 | X_Bat == 10 | X_Bat == 11 | X_Bat == 12] <- 3
X_Bat.new[X_Bat == 13 | X_Bat == 14 | X_Bat == 15 | X_Bat == 16] <- 4

X_info.new <- X_info
X_info.new$Batch <- X_Bat.new

# use SERRF for intra-batch correction
X_Dat.intraCor <- read.csv("Data_II_intraCor_withSERRF.csv", header = TRUE)
X_Dat.intraCor <- as.matrix(X_Dat.intraCor[, -c(1: 4)])

# Batch effect correction
# CordBat BEC
CdBDat <- CordBat(X_Dat.intraCor, X_Bat.new, group = NULL, ref.batch = 1, eps = 1e-5)
CdBcor_Data_II <- CdBDat$X.cor

# QC-RLSC BEC
library(statTarget)
samPeno <- 'Data II_pheno.csv'
samFile <- 'Data II_file.csv'
shiftCor(samPeno, samFile, MLmethod = "QCRLSC", QCspan = 0.75, 
         imputeM = "KNN", coCV = 100)

QRcor_Data_II <- read.csv("Data_II_QRcor.csv", header = TRUE)
QRcor_Data_II <- as.matrix(QRcor_Data_II[, -c(1, 2)])

# SERRF BEC
# through the web tool 
SRcor_Data_II <- read.csv("Data_II_SRcor.csv", header = TRUE)
SRcor_Data_II <- as.matrix(t(SRcor_Data_II.1[, -1]))

# SVR BEC
library(MetNormalizer)
metNor(ms1.data.name = 'data.csv', 
       sample.info.name = 'sample.info.csv')

SVRcor <- read.csv("Data_II_SVRcor.csv", header = TRUE)
SVRcor <- as.matrix(t(SVRcor[, -c(1:6)]))

SVRcor_Data_II <- X_Dat.log.proc
N.sbj <- sum(X_Grp != 'qc')
SVRcor_Data_II[X_Grp != 'qc', ] <- SVRcor[c(1: N.sbj), ]
SVRcor_Data_II[X_Grp == 'qc', ] <- SVRcor.1[(N.sbj + 1): N, ]


# --------------------------------------- Evaluation --------------------------------------- #
# -----------------------------------
# Principal component analysis (PCA)
# -----------------------------------
library(mixOmics)

# Before BEC
X_Dat.scale.bef <- scale(X_Dat.intraCor, center = T, scale = T)

Sbj_Dat <- X_Dat.scale.bef[X_Grp != 'qc', ]
QC_Dat <- X_Dat.scale.bef[X_Grp == 'qc', ]

pca.sbj.bef <- pca(Sbj_Dat, ncomp = 2, center = T, scale = T)
pca.bef.expl <- pca.sbj.bef$prop_expl_var$X
pca.sbj.bef.varx <- pca.sbj.bef$variates$X
pca.sbj.bef.load <- pca.sbj.bef$loadings$X
pca.qc.bef.varx <- QC_Dat %*% pca.sbj.bef.load
pca.bef.varX <- matrix(0, nrow(X_Dat.scale.bef), 2)
pca.bef.varX[X_Grp != 'qc', ] <- pca.sbj.bef.varx
pca.bef.varX[X_Grp == 'qc', ] <- pca.qc.bef.varx

PCA_scatter_withQC(data = pca.bef.varX, batch = X_Bat.new, QC_label = 'qc', 
                    group = X_Grp, expl.var = pca.bef.expl,
                    xlim = c(min(pca.bef.varX[, 1]) - 1, max(pca.bef.varX[, 1]) + 1), 
                    ylim = c(min(pca.bef.varX[, 2]) - 1, max(pca.bef.varX[, 2]) + 1), 
                    legend.pos = 'none')

# After CordBat correction
X_Dat.scale.cor <- scale(CdBcor_Data_II, center = T, scale = T)

Sbj_Dat <- X_Dat.scale.cor[X_Grp != 'qc', ]
QC_Dat <- X_Dat.scale.cor[X_Grp == 'qc', ]

pca.sbj.CdBcor <- pca(Sbj_Dat, ncomp = 2, center = F, scale = F)
pca.CdBcor.expl <- pca.sbj.CdBcor$prop_expl_var$X
pca.sbj.CdBcor.varx <- pca.sbj.CdBcor$variates$X
pca.sbj.CdBcor.load <- pca.sbj.CdBcor$loadings$X
pca.qc.CdBcor.varx <- QC_Dat %*% pca.sbj.CdBcor.load
pca.CdBcor.varX <- matrix(0, nrow(X_Dat.scale.cor), 2)
pca.CdBcor.varX[X_Grp != 'qc', ] <- pca.sbj.CdBcor.varx
pca.CdBcor.varX[X_Grp == 'qc', ] <- pca.qc.CdBcor.varx

PCA_scatter_withQC(data = pca.CdBcor.varX, batch = X_Bat.new, QC_label = 'qc', 
                    group = X_Grp, expl.var = pca.CdBcor.expl,
                    xlim = c(min(pca.CdBcor.varX[, 1]) - 1, max(pca.CdBcor.varX[, 1]) + 1), 
                    ylim = c(min(pca.CdBcor.varX[, 2]) - 1, max(pca.CdBcor.varX[, 2]) + 1))


# After QC-RLSC correction
X_Dat.scale.cor <- scale(QRcor_Data_II, center = T, scale = T)

Sbj_Dat <- X_Dat.scale.cor[X_Grp != 'qc', ]
QC_Dat <- X_Dat.scale.cor[X_Grp == 'qc', ]

pca.sbj.QRcor <- pca(Sbj_Dat, ncomp = 2, center = F, scale = F)
pca.QRcor.expl <- pca.sbj.QRcor$prop_expl_var$X
pca.sbj.QRcor.varx <- pca.sbj.QRcor$variates$X
pca.sbj.QRcor.load <- pca.sbj.QRcor$loadings$X
pca.qc.QRcor.varx <- QC_Dat %*% pca.sbj.QRcor.load
pca.QRcor.varX <- matrix(0, nrow(X_Dat.scale.cor), 2)
pca.QRcor.varX[X_Grp != 'qc', ] <- pca.sbj.QRcor.varx
pca.QRcor.varX[X_Grp == 'qc', ] <- pca.qc.QRcor.varx

PCA_scatter_withQC(data = pca.QRcor.varX, batch = X_Bat.new, QC_label = 'qc', 
                    group = X_Grp, expl.var = pca.QRcor.expl,
                    xlim = c(min(pca.QRcor.varX[, 1]) - 1, max(pca.QRcor.varX[, 1]) + 1), 
                    ylim = c(min(pca.QRcor.varX[, 2]) - 1, max(pca.QRcor.varX[, 2]) + 1), 
                    legend.pos = 'right')


# After SERRF correction
X_Dat.scale.cor <- scale(SRcor_Data_II, center = T, scale = T)

Sbj_Dat <- X_Dat.scale.cor[X_Grp != 'qc', ]
QC_Dat <- X_Dat.scale.cor[X_Grp == 'qc', ]

pca.sbj.SRcor <- pca(Sbj_Dat, ncomp = 2, center = F, scale = F)
pca.SRcor.expl <- pca.sbj.SRcor$prop_expl_var$X
pca.sbj.SRcor.varx <- pca.sbj.SRcor$variates$X
pca.sbj.SRcor.load <- pca.sbj.SRcor$loadings$X
pca.qc.SRcor.varx <- QC_Dat %*% pca.sbj.SRcor.load
pca.SRcor.varX <- matrix(0, nrow(X_Dat.scale.cor), 2)
pca.SRcor.varX[X_Grp != 'qc', ] <- pca.sbj.SRcor.varx
pca.SRcor.varX[X_Grp == 'qc', ] <- pca.qc.SRcor.varx

PCA_scatter_withQC(data = pca.SRcor.varX, batch = X_Bat.new, QC_label = 'qc', 
                    group = X_Grp, expl.var = pca.SRcor.expl,
                    xlim = c(min(pca.SRcor.varX[, 1]) - 1, max(pca.SRcor.varX[, 1]) + 1), 
                    ylim = c(min(pca.SRcor.varX[, 2]) - 1, max(pca.SRcor.varX[, 2]) + 1), 
                    legend.pos = 'right')


# After SVR correction
X_Dat.scale.cor <- scale(SVRcor_Data_II, center = T, scale = T)

Sbj_Dat <- X_Dat.scale.cor[X_Grp != 'qc', ]
QC_Dat <- X_Dat.scale.cor[X_Grp == 'qc', ]

pca.sbj.SVRcor <- pca(Sbj_Dat, ncomp = 2, center = F, scale = F)
pca.SVRcor.expl <- pca.sbj.SVRcor$prop_expl_var$X
pca.sbj.SVRcor.varx <- pca.sbj.SVRcor$variates$X
pca.sbj.SVRcor.load <- pca.sbj.SVRcor$loadings$X
pca.qc.SVRcor.varx <- QC_Dat %*% pca.sbj.SVRcor.load
pca.SVRcor.varX <- matrix(0, nrow(X_Dat.scale.cor), 2)
pca.SVRcor.varX[X_Grp != 'qc', ] <- pca.sbj.SVRcor.varx
pca.SVRcor.varX[X_Grp == 'qc', ] <- pca.qc.SVRcor.varx

PCA_scatter_withQC(data = pca.SVRcor.varX, batch = X_Bat.new, QC_label = 'qc', 
                    group = X_Grp, expl.var = pca.SVRcor.expl,
                    xlim = c(min(pca.SVRcor.varX[, 1]) - 1, max(pca.SVRcor.varX[, 1]) + 1), 
                    ylim = c(min(pca.SVRcor.varX[, 2]) - 1, max(pca.SVRcor.varX[, 2]) + 1), 
                    legend.pos = 'right')


## Mahalanobis distance in PC space
# subject samples

# Before BEC
bef.sbj.md <- c()
for (i in c(1: batch.num)) {
  PCdat.bati <- pca.bef.varX[X_Grp != 'qc' & X_Bat.new == i, ]
  PCdat.bati.m <- colMeans(PCdat.bati)
  PCdat.bati.s <- cov(PCdat.bati)
  for (j in c(1: batch.num)) {
    if (j != i) {
      PCdat.batj <- pca.bef.varX[X_Grp != 'qc' & X_Bat.new == j, ]
      bef.sbj.md <- c(bef.sbj.md, mahalanobis(PCdat.batj, PCdat.bati.m, PCdat.bati.s))
    }
  }
}
bef.sbj.md.med <- median(bef.sbj.md)

# CordBat correction
CdBcor.sbj.md <- c()
for (i in c(1: batch.num)) {
  PCdat.bati <- pca.CdBcor.varX[X_Grp != 'qc' & X_Bat.new == i, ]
  PCdat.bati.m <- colMeans(PCdat.bati)
  PCdat.bati.s <- cov(PCdat.bati)
  for (j in c(1: batch.num)) {
    if (j != i) {
      PCdat.batj <- pca.CdBcor.varX[X_Grp != 'qc' & X_Bat.new == j, ]
      CdBcor.sbj.md <- c(CdBcor.sbj.md, mahalanobis(PCdat.batj, PCdat.bati.m, PCdat.bati.s))
    }
  }
}
CdBcor.sbj.md.med <- median(CdBcor.sbj.md)

# QC-RLSC correction
QRcor.sbj.md <- c()
for (i in c(1: batch.num)) {
  PCdat.bati <- pca.QRcor.varX[X_Grp != 'qc' & X_Bat.new == i, ]
  PCdat.bati.m <- colMeans(PCdat.bati)
  PCdat.bati.s <- cov(PCdat.bati)
  for (j in c(1: batch.num)) {
    if (j != i) {
      PCdat.batj <- pca.QRcor.varX[X_Grp != 'qc' & X_Bat.new == j, ]
      QRcor.sbj.md <- c(QRcor.sbj.md, mahalanobis(PCdat.batj, PCdat.bati.m, PCdat.bati.s))
    }
  }
}
QRcor.sbj.md.med <- median(QRcor.sbj.md)

# SERRF correction
SRcor.sbj.md <- c()
for (i in c(1: batch.num)) {
  PCdat.bati <- pca.SRcor.varX[X_Grp != 'qc' & X_Bat.new == i, ]
  PCdat.bati.m <- colMeans(PCdat.bati)
  PCdat.bati.s <- cov(PCdat.bati)
  for (j in c(1: batch.num)) {
    if (j != i) {
      
      PCdat.batj <- pca.SRcor.varX[X_Grp != 'qc' & X_Bat.new == j, ]
      SRcor.sbj.md <- c(SRcor.sbj.md, mahalanobis(PCdat.batj, PCdat.bati.m, PCdat.bati.s))
    }
  }
}
SRcor.sbj.md.med <- median(SRcor.sbj.md)

# SVR correction
SVRcor.sbj.md <- c()
for (i in c(1: batch.num)) {
  PCdat.bati <- pca.SVRcor.varX[X_Grp != 'qc' & X_Bat.new == i, ]
  PCdat.bati.m <- colMeans(PCdat.bati)
  PCdat.bati.s <- cov(PCdat.bati)
  for (j in c(1: batch.num)) {
    if (j != i) {
      PCdat.batj <- pca.SVRcor.varX[X_Grp != 'qc' & X_Bat.new == j, ]
      SVRcor.sbj.md <- c(SVRcor.sbj.md, mahalanobis(PCdat.batj, PCdat.bati.m, PCdat.bati.s))
    }
  }
}
SVRcor.sbj.md.med <- median(SVRcor.sbj.md)

# plot
dattype <- factor(c(rep('Uncorrected', length(bef.sbj.md)), 
                    rep('QC-RLSC', length(QRcor.sbj.md)), 
                    rep('SVR', length(SVRcor.sbj.md)), 
                    rep('SERRF', length(SRcor.sbj.md)), 
                    rep('CordBat', length(CdBcor.sbj.md))), 
                  levels = c('Uncorrected', 'QC-RLSC', 'SVR','SERRF', 'CordBat'))

dat.md <- data.frame(dattype = dattype, 
                     md = c(bef.sbj.md, QRcor.sbj.md, SVRcor.sbj.md, 
                            SRcor.sbj.md, CdBcor.sbj.md))

fon <- 'sans'
set_col <- c('#F39B7FFF', '#3C5488FF', '#00A087FF', 
             '#4DBBD5FF', '#E64B35FF')
ggplot(data = dat.md, aes(x = dattype, y = md)) + 
  geom_boxplot(aes(fill = dattype), color = 'black', size = 1, outlier.color = NA) + 
  scale_fill_manual(values = set_col) + theme_bw() + 
  xlab(NULL) + ylab('Mahalanobis distance') + 
  theme(axis.text.y = element_text(size = 25, family = fon, face = 'bold', 
                                   vjust = 0.5, hjust = 0.5), 
        axis.title.y = element_text(size = 30, family = fon, face = 'bold'),
        axis.text.x = element_blank(), 
        legend.position = c(0.8, 0.8), 
        legend.background = element_rect(color = 'black', linewidth = 1), 
        legend.title = element_text(size = 25, family = fon, face = 'bold'), 
        legend.text = element_text(size = 25, family = fon, face = 'bold'), 
        panel.border = element_rect(linewidth = 3)) + 
  coord_cartesian(ylim = c(0, 11)) + 
  labs(fill = 'Method')


# ------------------
# cumulative RSD
# ------------------
rsd.val <- seq(0.001, 1, 0.001)

# CordBat correction
CdBdat <- 2^CdBcor_Data_II
CdB.Sbj <- CdBdat[X_Grp != 'qc', ]
CdB.QC <- CdBdat[X_Grp == 'qc', ]
Sbj.met.m <- apply(CdB.Sbj, MARGIN = 2, FUN = mean)
Sbj.met.sd <- apply(CdB.Sbj, MARGIN = 2, FUN = sd)
Sbj.met.rsd <- Sbj.met.sd / Sbj.met.m
CdB.Sbj.rsd.med <- median(Sbj.met.rsd)
CdB.Sbj.cumfreq <- rep(0, length(rsd.val))
for (i in c(1: length(rsd.val))) {
  CdB.Sbj.cumfreq[i] <- sum(Sbj.met.rsd <= rsd.val[i]) / length(Sbj.met.rsd)
}

QC.met.m <- apply(CdB.QC, MARGIN = 2, FUN = mean)
QC.met.sd <- apply(CdB.QC, MARGIN = 2, FUN = sd)
QC.met.rsd <- QC.met.sd / QC.met.m
CdB.QC.rsd.med <- median(QC.met.rsd)
CdB.QC.cumfreq <- rep(0, length(rsd.val))
for (i in c(1: length(rsd.val))) {
  CdB.QC.cumfreq[i] <- sum(QC.met.rsd <= rsd.val[i]) / length(QC.met.rsd)
}


# SERRF correction
SRdat <- 2^SRcor_Data_II
SR.Sbj <- SRdat[X_Grp != 'qc', ] 
SR.QC <- SRdat[X_Grp == 'qc', ] 
Sbj.met.m <- apply(SR.Sbj, MARGIN = 2, FUN = mean)
Sbj.met.sd <- apply(SR.Sbj, MARGIN = 2, FUN = sd)
Sbj.met.rsd <- Sbj.met.sd / Sbj.met.m
SR.Sbj.rsd.med <- median(Sbj.met.rsd)
SR.Sbj.cumfreq <- rep(0, length(rsd.val))
for (i in c(1: length(rsd.val))) {
  SR.Sbj.cumfreq[i] <- sum(Sbj.met.rsd <= rsd.val[i]) / length(Sbj.met.rsd)
}

QC.met.m <- apply(SR.QC, MARGIN = 2, FUN = mean)
QC.met.sd <- apply(SR.QC, MARGIN = 2, FUN = sd)
QC.met.rsd <- QC.met.sd / QC.met.m
SR.QC.rsd.med <- median(QC.met.rsd)
SR.QC.cumfreq <- rep(0, length(rsd.val))
for (i in c(1: length(rsd.val))) {
  SR.QC.cumfreq[i] <- sum(QC.met.rsd <= rsd.val[i]) / length(QC.met.rsd)
}

# -------------------
# add noise to QC
# -------------------
QCdat <- X_Dat.intraCor[X_Grp == 'qc', ]
SNR.vals <- seq(100, 80, -10)

for (i in c(1: length(SNR.vals))) {
  QCdat <- X_Dat.intraCor[X_Grp == 'qc', ]
  QCdat.addnoise <- addnoise_to_QC(QCdat, SNR = SNR.vals[i], dolog = F)
  QCdat <- QCdat.addnoise$dat
  
  X_Dat.addnoise <- X_Dat.intraCor
  X_Dat.addnoise[X_Grp == 'qc', ] <- QCdat
  
  dat.addnoise <- data.frame(X_info.new, X_Dat.addnoise)
  write.csv(dat.addnoise, file = paste0('Data_II_addnoisetoQC_SNR', SNR.vals[i], '.csv'))
}

# add noise SNR = 80
snr <- 80
SRdat_SNR <- read.csv(paste0("Data_II_SRcor_SNR", snr, ".csv"), header = TRUE)
SRdat_SNR <- as.matrix(t(SRdat_SNR[, -1]))

SRdat <- 2^SRdat_SNR
SRn.Sbj <- SRdat[X_Grp != 'qc', ] 
SRn.QC <- SRdat[X_Grp == 'qc', ] 
Sbj.met.m <- apply(SRn.Sbj, MARGIN = 2, FUN = mean)
Sbj.met.sd <- apply(SRn.Sbj, MARGIN = 2, FUN = sd)
Sbj.met.rsd <- Sbj.met.sd / Sbj.met.m
SRn.Sbj.rsd.med <- median(Sbj.met.rsd)
SRn.Sbj.cumfreq <- rep(0, length(rsd.val))
for (i in c(1: length(rsd.val))) {
  SRn.Sbj.cumfreq[i] <- sum(Sbj.met.rsd <= rsd.val[i]) / length(Sbj.met.rsd)
}

QC.met.m <- apply(SRn.QC, MARGIN = 2, FUN = mean)
QC.met.sd <- apply(SRn.QC, MARGIN = 2, FUN = sd)
QC.met.rsd <- QC.met.sd / QC.met.m
SRn.QC.rsd.med <- median(QC.met.rsd)
SRn.QC.cumfreq <- rep(0, length(rsd.val))
for (i in c(1: length(rsd.val))) {
  SRn.QC.cumfreq[i] <- sum(QC.met.rsd <= rsd.val[i]) / length(QC.met.rsd)
}

SRn80.Sbj.cumfreq <- SRn.Sbj.cumfreq
SRn80.QC.cumfreq <- SRn.QC.cumfreq

# add noise SNR = 90
snr <- 90
SRdat_SNR <- read.csv(paste0("Data_II_SRcor_SNR", snr, ".csv"), header = TRUE)
SRdat_SNR <- as.matrix(t(SRdat_SNR[, -1]))

SRdat <- 2^SRdat_SNR
SRn.Sbj <- SRdat[X_Grp != 'qc', ] 
SRn.QC <- SRdat[X_Grp == 'qc', ] 
Sbj.met.m <- apply(SRn.Sbj, MARGIN = 2, FUN = mean)
Sbj.met.sd <- apply(SRn.Sbj, MARGIN = 2, FUN = sd)
Sbj.met.rsd <- Sbj.met.sd / Sbj.met.m
SRn.Sbj.rsd.med <- median(Sbj.met.rsd)
SRn.Sbj.cumfreq <- rep(0, length(rsd.val))
for (i in c(1: length(rsd.val))) {
  SRn.Sbj.cumfreq[i] <- sum(Sbj.met.rsd <= rsd.val[i]) / length(Sbj.met.rsd)
}

QC.met.m <- apply(SRn.QC, MARGIN = 2, FUN = mean)
QC.met.sd <- apply(SRn.QC, MARGIN = 2, FUN = sd)
QC.met.rsd <- QC.met.sd / QC.met.m
SRn.QC.rsd.med <- median(QC.met.rsd)
SRn.QC.cumfreq <- rep(0, length(rsd.val))
for (i in c(1: length(rsd.val))) {
  SRn.QC.cumfreq[i] <- sum(QC.met.rsd <= rsd.val[i]) / length(QC.met.rsd)
}

SRn90.Sbj.cumfreq <- SRn.Sbj.cumfreq
SRn90.QC.cumfreq <- SRn.QC.cumfreq


# add noise SNR = 100
snr <- 100
SRdat_SNR <- read.csv(paste0("Data_II_SRcor_SNR", snr, ".csv"), header = TRUE)
SRdat_SNR <- as.matrix(t(SRdat_SNR[, -1]))

SRdat <- 2^SRdat_SNR
SRn.Sbj <- SRdat[X_Grp != 'qc', ] 
SRn.QC <- SRdat[X_Grp == 'qc', ] 
Sbj.met.m <- apply(SRn.Sbj, MARGIN = 2, FUN = mean)
Sbj.met.sd <- apply(SRn.Sbj, MARGIN = 2, FUN = sd)
Sbj.met.rsd <- Sbj.met.sd / Sbj.met.m
SRn.Sbj.rsd.med <- median(Sbj.met.rsd)
SRn.Sbj.cumfreq <- rep(0, length(rsd.val))
for (i in c(1: length(rsd.val))) {
  SRn.Sbj.cumfreq[i] <- sum(Sbj.met.rsd <= rsd.val[i]) / length(Sbj.met.rsd)
}

QC.met.m <- apply(SRn.QC, MARGIN = 2, FUN = mean)
QC.met.sd <- apply(SRn.QC, MARGIN = 2, FUN = sd)
QC.met.rsd <- QC.met.sd / QC.met.m
SRn.QC.rsd.med <- median(QC.met.rsd)
SRn.QC.cumfreq <- rep(0, length(rsd.val))
for (i in c(1: length(rsd.val))) {
  SRn.QC.cumfreq[i] <- sum(QC.met.rsd <= rsd.val[i]) / length(QC.met.rsd)
}

SRn100.Sbj.cumfreq <- SRn.Sbj.cumfreq
SRn100.QC.cumfreq <- SRn.QC.cumfreq

# plot

# subject samples
all.sbj.cumfreq <- c(SRn100.Sbj.cumfreq, SRn90.Sbj.cumfreq, 
                     SRn80.Sbj.cumfreq, SR.Sbj.cumfreq, 
                     CdB.Sbj.cumfreq)

dat.sbj <- data.frame(rsd.val = rep(rsd.val * 100, 5), 
                      cumfreq = all.sbj.cumfreq * 100)

# area
area.CdB.sbj <- (CdB.Sbj.cumfreq[1] + CdB.Sbj.cumfreq[1000] + 
                  2 * sum(CdB.Sbj.cumfreq[c(2: 999)])) * 0.001 / 2

area.SR.sbj <- (SR.Sbj.cumfreq[1] + SR.Sbj.cumfreq[1000] + 
                  2 * sum(SR.Sbj.cumfreq[c(2: 999)])) * 0.001 / 2

area.SRn80.sbj <- (SRn80.Sbj.cumfreq[1] + SRn80.Sbj.cumfreq[1000] + 
                     2 * sum(SRn80.Sbj.cumfreq[c(2: 999)])) * 0.001 / 2

area.SRn90.sbj <- (SRn90.Sbj.cumfreq[1] + SRn90.Sbj.cumfreq[1000] + 
                     2 * sum(SRn90.Sbj.cumfreq[c(2: 999)])) * 0.001 / 2

area.SRn100.sbj <- (SRn100.Sbj.cumfreq[1] + SRn100.Sbj.cumfreq[1000] + 
                      2 * sum(SRn100.Sbj.cumfreq[c(2: 999)])) * 0.001 / 2

label80 <- paste0('SERRF SNR-80, Area = ', round(area.SRn80.sbj, 2))
label90 <- paste0('SERRF SNR-90, Area = ', round(area.SRn90.sbj, 2))
label100 <- paste0('SERRF SNR-100, Area = ', round(area.SRn100.sbj, 2))
labelSR <- paste0('SERRF, Area = ', round(area.SR.sbj, 2))
labelCdB <- paste0('CordBat, Area = ', round(area.CdB.sbj, 2))

dat.type <- factor(c(rep(label100, length(rsd.val)), 
                     rep(label90, length(rsd.val)), 
                     rep(label80, length(rsd.val)), 
                     rep(labelSR, length(rsd.val)), 
                     rep(labelCdB, length(rsd.val))), 
                   levels = c(label80, label90, 
                              label100, labelSR, 
                              labelCdB))

fon <- 'sans'
idx20 <- which(dat.sbj$rsd.val == 20)
p.idx <- idx20
cbPalette <- c("#CC79A7", "blue", "#E69F00", "green", "red")
set_alpha <- c(0.3, 1, 0.3, rep(1, 2))
set_shape <- c(NA, 24, NA, rep(24, 2))
set_linetype <- c(1, 4, 1, 1, 2)

library(dplyr)
ggplot() +
  geom_line(data = dat.sbj, 
            aes(x = rsd.val, y = cumfreq, color = dat.type, alpha = dat.type, 
                linetype = dat.type), 
            linewidth = linesize) + 
  geom_vline(xintercept = 20, linetype = 2, size = 1) +
  geom_point(data = dat.sbj[p.idx, ], 
             aes(x = rsd.val, 
                 y = cumfreq, 
                 fill = dat.type[p.idx], 
                 alpha = dat.type[p.idx], 
                 shape = dat.type[p.idx]), 
             size = 5, color = 'black', stroke = 1) + 
  annotate(geom = 'segment', x = 13, y = 28, xend = 18.5, yend = 24, 
           arrow = arrow(length = unit(0.5, 'cm')), linewidth = 1) + 
  annotate(geom = 'text', x = 11, y = 30, 
           label = paste0(round(QF1.Sbj.cumfreq[rsd.val == 0.2] * 100, 2), '%'), 
           size = 5, color = 'red') + 
  annotate(geom = 'segment', x = 25, y = 34, xend = 21.5, yend = 24, 
           arrow = arrow(length = unit(0.5, 'cm')), linewidth = 1) + 
  annotate(geom = 'text', x = 26, y = 36, 
           label = paste0(round(SR.Sbj.cumfreq[rsd.val == 0.2] * 100, 2), '%'), 
           size = 5, color = '#008B45FF') + 
  annotate(geom = 'segment', x = 30, y = 15, xend = 22, yend = 14, 
           arrow = arrow(length = unit(0.5, 'cm')), linewidth = 1) + 
  annotate(geom = 'text', x = 37, y = 15, 
           label = paste0(round(SRn90.Sbj.cumfreq[rsd.val == 0.2] * 100, 2), '%'), 
           size = 5, color = 'blue') + 
  scale_color_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette) + 
  scale_alpha_manual(values = set_alpha) + 
  scale_shape_manual(values = set_shape) + 
  scale_linetype_manual(values = set_linetype) + 
  scale_x_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) + 
  theme_bw() + guides(fill = "none", alpha = "none", shape = "none", linetype = "none") + 
  theme(legend.position = c(0.27, 0.83), 
        axis.title = element_text(size = 25, family = fon, face = 'bold'), 
        axis.text = element_text(size = 18, family = fon, face = 'bold'),
        legend.background = element_rect(color = 'black', linewidth = 1), 
        legend.box.spacing = unit(1, 'cm'), 
        legend.text = element_text(size = 15, family = fon),
        panel.border = element_rect(linewidth = 2.5)) + 
  xlab('RSD (%)') + 
  ylab('Cumulative Percentages (%)') + ylim(0, 100) +
  labs(color = NULL, title = NULL)


# QC samples
all.qc.cumfreq <- c(SRn100.QC.cumfreq, SRn90.QC.cumfreq, 
                    SRn80.QC.cumfreq, SR.QC.cumfreq, 
                    CdB.QC.cumfreq)

dat.qc <- data.frame(rsd.val = rep(rsd.val * 100, 5), 
                     cumfreq = all.qc.cumfreq * 100)

# area
area.CdB.qc <- (CdB.QC.cumfreq[1] + CdB.QC.cumfreq[1000] + 
                 2 * sum(CdB.QC.cumfreq[c(2: 999)])) * 0.001 / 2

area.SR.qc <- (SR.QC.cumfreq[1] + SR.QC.cumfreq[1000] + 
                 2 * sum(SR.QC.cumfreq[c(2: 999)])) * 0.001 / 2

area.SRn80.qc <- (SRn80.QC.cumfreq[1] + SRn80.QC.cumfreq[1000] + 
                    2 * sum(SRn80.QC.cumfreq[c(2: 999)])) * 0.001 / 2

area.SRn90.qc <- (SRn90.QC.cumfreq[1] + SRn90.QC.cumfreq[1000] + 
                    2 * sum(SRn90.QC.cumfreq[c(2: 999)])) * 0.001 / 2

area.SRn100.qc <- (SRn100.QC.cumfreq[1] + SRn100.QC.cumfreq[1000] + 
                     2 * sum(SRn100.QC.cumfreq[c(2: 999)])) * 0.001 / 2

label80 <- paste0('SERRF SNR-80, Area = ', round(area.SRn80.qc, 2))
label90 <- paste0('SERRF SNR-90, Area = ', round(area.SRn90.qc, 2))
label100 <- paste0('SERRF SNR-100, Area = ', round(area.SRn100.qc, 2))
labelSR <- paste0('SERRF, Area = ', round(area.SR.qc, 2))
labelCdB <- paste0('CordBat, Area = ', round(area.CdB.qc, 2))

dat.type <- factor(c(rep(label100, length(rsd.val)), 
                     rep(label90, length(rsd.val)), 
                     rep(label80, length(rsd.val)), 
                     rep(labelSR, length(rsd.val)), 
                     rep(labelCdB, length(rsd.val))), 
                   levels = c(label80, label90, 
                              label100, labelSR, 
                              labelCdB))

fon <- 'sans'
idx20 <- which(dat.qc$rsd.val == 20)
p.idx <- idx20
cbPalette <- c("#CC79A7", "blue", "#E69F00", "green", "red")
set_alpha <- c(0.3, 1, 0.3, rep(1, 2))
set_shape <- c(NA, 24, NA, rep(24, 2))
set_linetype <- c(1, 4, 1, 1, 2)

ggplot() +
  geom_line(data = dat.qc, 
            aes(x = rsd.val, y = cumfreq, color = dat.type, alpha = dat.type, 
                linetype = dat.type), 
            linewidth = linesize) + 
  geom_vline(xintercept = 20, linetype = 2, size = 1) +
  geom_point(data = dat.qc[p.idx, ], 
             aes(x = rsd.val, 
                 y = cumfreq, 
                 fill = dat.type[p.idx], 
                 alpha = dat.type[p.idx], 
                 shape = dat.type[p.idx]), 
             size = 5, color = 'black', stroke = 1.5) + 
  annotate(geom = 'segment', x = 13, y = 85, xend = 17.5, yend = 88, 
           arrow = arrow(length = unit(0.5, 'cm')), linewidth = 1) + 
  annotate(geom = 'text', x = 11, y = 83,
           label = paste0(round(QF1.QC.cumfreq[rsd.val == 0.2] * 100, 2), '%'),
           size = 5, color = 'red') +
  annotate(geom = 'segment', x = 12, y = 100, xend = 18.5, yend = 100,
           arrow = arrow(length = unit(0.5, 'cm')), linewidth = 1) +
  annotate(geom = 'text', x = 6, y = 100,
           label = paste0(round(SR.QC.cumfreq[rsd.val == 0.2] * 100, 2), '%'),
           size = 5, color = '#008B45FF') +
  annotate(geom = 'segment', x = 27, y = 24, xend = 21.5, yend = 19,
           arrow = arrow(length = unit(0.5, 'cm')), linewidth = 1) +
  annotate(geom = 'text', x = 29, y = 26,
           label = paste0(round(SRn90.QC.cumfreq[rsd.val == 0.2] * 100, 2), '%'),
           size = 5, color = 'blue') +
  scale_color_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette) + 
  scale_alpha_manual(values = set_alpha) + 
  scale_shape_manual(values = set_shape) + 
  scale_linetype_manual(values = set_linetype) + 
  scale_x_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) + 
  theme_bw() + guides(fill = "none", alpha = 'none', shape = 'none', linetype = "none") + 
  theme(axis.title = element_text(size = 25, family = fon, face = 'bold'), 
        axis.text = element_text(size = 18, family = fon, face = 'bold'), 
        legend.position = c(0.73, 0.18), 
        legend.background = element_rect(linewidth = 1, color = 'black'), 
        legend.text = element_text(size = 15, family = fon),
        panel.border = element_rect(linewidth = 2)) + 
  xlab('RSD (%)') + ylab('Cumulative Percentages (%)') + ylim(0, 100) +
  labs(color = NULL, title = NULL)


