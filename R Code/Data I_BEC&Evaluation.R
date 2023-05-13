# Title: Data I BEC and evaluation
# Author: Fanjing Guo
# Date: 2022.12.28

# --------------------- Data pre-processing and batch effect correction ------------------- #
# LC-MS/MS metabolomics dataset from an atherosclerosis study
X <- read.csv("Data I_Raw_filtered.csv", header = TRUE)
X_info <- X[, c(1: 4)]
X_Dat <- X[, -c(1: 4)]

# add group information
X_clinic <- read.csv("Data I_clinical.csv", header = TRUE)
patID <- X_clinic$Pat.ID
samp.patID <- gsub('[^0-9]', '', X_info$sample)

# HTN
HTN.label <- rep(NA, nrow(X_info))
for (i in c(1: length(patID))) {
  HTN.label[which(samp.patID == patID[i])] <- X_clinic$HTN[i]
}
HTN.isna.Idx <- which(is.na(HTN.label))
QC.Idx <- which(X_info$type == 'QC')
HTN.label[intersect(HTN.isna.Idx, QC.Idx)] <- 'QC'
delsampIdx <- which(is.na(HTN.label))
HTN.label <- HTN.label[-delsampIdx]

# Gender
Gen.label <- rep(NA, nrow(X_info))
for (i in c(1: length(patID))) {
  Gen.label[which(samp.patID == patID[i])] <- X_clinic$Gender[i]
}
Gen.isna.Idx <- which(is.na(Gen.label))
QC.Idx <- which(X_info$type == 'QC')
Gen.label[intersect(Gen.isna.Idx, QC.Idx)] <- 'QC'
delsampIdx <- which(is.na(Gen.label))
Gen.label <-Gen.label[-delsampIdx]

# delete samples without group information
X_Dat <- X_Dat[-delsampIdx, ]
X_info <- X_info[-delsampIdx, ]

# add group information
# grp.label <- HTN.label
grp.label <- Gen.label
X_info <- cbind(X_info, grp.label)

# pre-processing
N <- nrow(X_Dat)

# delete metabolites with >= 50% missing values
Sbj_Dat.proc <- X_Dat[X_info$type == 'S', ]
Sbj_info <- X_info[X_info$type == 'S', ]

na.num <- colSums(is.na(Sbj_Dat.proc))
delmetIdx <- which(na.num / N >= 0.5)
if (length(delmetIdx) != 0) {
  Sbj_Dat.proc <- Sbj_Dat.proc[, -delmetIdx]
}
p <- ncol(Sbj_Dat.proc)

# impute missing values using knn
library(DMwR2)
Sbj_Dat.proc <- as.data.frame(Sbj_Dat.proc)
Sbj_Dat.proc <- knnImputation(Sbj_Dat.proc)

# detect and delete outliers
library(mixOmics)
Sbj_Dat.log.proc <- as.matrix(log2(Sbj_Dat.proc))
pca.dat <- pca(Sbj_Dat.log.proc, ncomp = 2, center = TRUE, scale = TRUE)
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

Sbj_Dat.log.proc <- Sbj_Dat.log.proc[-delsampIdx, ]
Sbj_info <- Sbj_info[-delsampIdx, ]

Sbj_Bat <- Sbj_info$batch
Sbj_Grp <- Sbj_info$HTN

N <- nrow(Sbj_Dat.log.proc)
p <- ncol(Sbj_Dat.log.proc)

# intra-batch correction (LOESS - Only subject samples)
Sbj_Bat.f <- as.factor(Sbj_Bat)
batch.num <- nlevels(Sbj_Bat.f)
batch.lev <- levels(Sbj_Bat.f)
Sbj_Dat.intraCor <- Sbj_Dat.log.proc
for (i in c(1: batch.num)) {
  dat.i <- Sbj_Dat.log.proc[Sbj_Bat.f == batch.lev[i], ]
  ord.i <- Sbj_info$order[Sbj_Bat.f == batch.lev[i]]
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
  Sbj_Dat.intraCor[Sbj_Bat.f == batch.lev[i], ] <- dat.i.cor
  cat(i, '\n')
}

Sbj_Dat.intraCor <- as.matrix(Sbj_Dat.intraCor)
colnames(Sbj_Dat.intraCor) <- c(1: ncol(Sbj_Dat.intraCor))

# Batch effect correction
# CorbBat BEC
CdBDat <- CordBat(Sbj_Dat.intraCor, Sbj_Bat, Sbj_Grp, ref.batch = ref.batch, eps = 1e-3)
CdBcor_Data_I <- CdBDat$X.cor
Sbj_Dat.delout <- CdBDat$X.delout


# ------------------------------- Evaluation ----------------------------------------- #
# -----------------------------------
# Principal component analysis (PCA)
# -----------------------------------
library(mixOmics)
# Before correction
pca.bef <- pca(Sbj_Dat.log.proc, ncomp = 2, center = TRUE, scale = TRUE)
pca.bef.varX <- pca.bef$variates$X
pca.bef.expl <- pca.bef$prop_expl_var$X

PCA_scatter(data = pca.bef.varX, batch = Sbj_Bat, 
            expl.var = pca.bef.expl,
            xlim = c(min(pca.bef.varX[, 1]) - 1, max(pca.bef.varX[, 1]) + 1), 
            ylim = c(min(pca.bef.varX[, 2]) - 1, max(pca.bef.varX[, 2]) + 1), 
            legend.pos = 'none')

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
pca.cor <- pca(QFcor_Data_I, ncomp = 2, center = TRUE, scale = TRUE)
pca.cor.varX <- pca.cor$variates$X
pca.cor.expl <- pca.cor$prop_expl_var$X

PCA_scatter(data = pca.cor.varX, batch = Sbj_Bat, 
             expl.var = pca.cor.expl,
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


###########################################################################################
# Other evaluation
savefile.path0 <- 'F:/Graduate/Graduate Studies/Papers/1-GGM for Batch Effect Correction/New Results/Data I/'

# ------------------------------
# Compared with QC-based method
# ------------------------------
# delete metabolites with >= 50% missing values
na.num <- colSums(is.na(X_Dat))
delmetIdx <- which(na.num / N >= 0.5)
X_Dat.proc <- X_Dat
if (length(delmetIdx) != 0) {
  X_Dat.proc <- X_Dat[, -delmetIdx]
}
X_Dat.proc <- as.matrix(X_Dat.proc)
p <- ncol(X_Dat.proc)

# impute missing values using knn
library(DMwR2)
X_Dat.proc <- as.data.frame(X_Dat.proc)
X_Dat.proc <- knnImputation(X_Dat.proc)

# log2 transformed
X_Dat.log.proc <- as.matrix(log2(X_Dat.proc))

# detect and delete outliers
library(mixOmics)
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

X_Dat.log.proc <- X_Dat.log.proc[-delsampIdx, ]
X_info <- X_info[-delsampIdx, ]

N <- nrow(X_Dat.log.proc)
p <- ncol(X_Dat.log.proc)

# delete the last 3 sample (not QCs)
# - the sample in the end of a batch should be a QC sample
X_Dat.log.proc <- X_Dat.log.proc[c(1: (N - 3)), ]
X_info <- X_info[c(1: (N - 3)), ]

N <- nrow(X_Dat.log.proc)
p <- ncol(X_Dat.log.proc)

colnames(X_info) <- c('SampName', 'RunOrder', 'SampType', 'Batch', 'Group')
batch <- X_info$Batch
group <-  X_info$Group
samp.order <- X_info$RunOrder
samp.name_type = data.frame(SampName = X_info$SampName, 
                            SampType = X_info$SampType)
samp.name <- apply(samp.name_type, MARGIN = 1, 
                   FUN = function(x) paste0(x[2], '-', x[1]))
samp.name.qc <- samp.name[group == 'QC'][1]
samp.name[group == 'QC'] <- paste0(samp.name.qc, '-', c(1: sum(group == 'QC')))
X_info$SampName <- samp.name

met.list <- colnames(X_Dat.log.proc)

# LOESS - intra-batch correction 1
p <- ncol(X_Dat.log.proc)
bat.f <- as.factor(batch)
batch.num <- nlevels(bat.f)
batch.lev <- levels(bat.f)
X_Dat.intraCor1 <- X_Dat.log.proc

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
  X_Dat.intraCor1[bat.f == batch.lev[i], ] <- dat.i.cor
  cat(i, '\n')
}

X_Dat.intraCor1 <- as.matrix(X_Dat.intraCor1)

# SERRF
X_Dat.intraCor <- X_Dat.intraCor1
SRformat_filePath0 <- paste0(getwd(), '/SRformat_file_Data I_eachBat/')
for (i in c(1: batch.num)) {
  dat.bati <- X_Dat.intraCor1[batch == batch.lev[i], ]
  bati.info <- X_info[batch == batch.lev[i], ]
  GenSRformatdat(dat.bati, bati.info, met.list, 
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

write.csv(X_Dat.intraCor, 'Data I_LES&SRintraCor.csv', row.names = F)
X_Dat.intraCor <- read.csv('Data I_LES&SRintraCor.csv', header = T)
X_Dat.intraCor <- as.matrix(X_Dat.intraCor)

# Batch effect correction
# CordBat
getRefbat(X_Dat.intraCor, batch)
CdBDat <- CordBat(X_Dat.intraCor, batch, ref.batch = 3, eps = 1e-3, 
                  print.detail = F)
CdBcor_Data_I <- CdBDat$X.cor

write.csv(CdBcor_Data_I, 'CdBcor_Data_I_LESSR_ref3.csv', row.names = F)
CdBcor_Data_I <- as.matrix(read.csv('CdBcor_Data_I_LESSR_ref3.csv', header = T))

# QC-based method
# QC-RLSC
QRcor_Data_I <- QCRLSC_cor(X_Dat.log.proc, samp.info = X_info, 
                           feature.list = met.list)

# SERRF
GenSRformatdat(X_Dat.log.proc, X_info, met.list, 
               wrfile.path = 'SRformat_Data_I.csv')

source('app.R')
runApp()

SRcor_Data_I <- read.csv(unzip(paste0(getwd(), '/SERRF Result_Data I.zip'), 
                                'normalized by - SERRF.csv'), header = T)
SRcor_Data_I <- as.matrix(t(SRcor_Data_I[, -1]))
# unlink('SERRF Result.zip')

# SVR
SVRcor_Data_I <- SVR_cor(X_Dat.log.proc, samp.info = X_info, 
                          mets.name = met.list)


# add noise to QC
GenDatNoise_SRformat(X_Dat.log.proc, dat.name = 'Data I', 
                     samp.info = X_info, feature.list = met.list)
runApp() # download result file in the same directory

SRformat_filePath0 <- paste0(getwd(), '/SRformat_file_addNoise_Data I/')
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

# Evaluation
dat.list <- list(Uncorrect = X_Dat.log.proc, 
                 RLSC = QRcor_Data_I, 
                 SVR = SVRcor_Data_I, 
                 SERRF = SRcor_Data_I, 
                 CordBat = CdBcor_Data_I)

set_arrow1 <- data.frame(arrow.pos = c('left', 'topleft', 'topleft'), 
                         arrow.ang = c(0, pi/20, pi/2.4), 
                         arrow.len = c(7, 9, 13))

set_arrow2 <- data.frame(arrow.pos = c('bottomleft', 'left', 'right'), 
                         arrow.ang = c(pi/2.5, 0, 0), 
                         arrow.len = c(10, 8, 8))

CmpWithQCmethod(dat.list, SRcor_NoiDat.list, batch = batch, group = group, 
                mark.rsd.sbj = 30, 
                set_arrow.sbj = set_arrow1, set_arrow.qc = set_arrow2, 
                save.path0 = paste0(savefile.path0, 'Eval 3/'))


# -----------------------------
# Preserving biological effect
# -----------------------------
Sbj_Dat.log.proc <- X_Dat.log.proc[group != 'QC', ]
Sbj_Dat.intraCor <- X_Dat.intraCor[group != 'QC', ]
batch.sbj <- batch[group != 'QC']
group.sbj <- group[group != 'QC']

# CordBat
getRefbat(Sbj_Dat.intraCor, batch.sbj)
CdBDat <- CordBat(Sbj_Dat.intraCor, batch.sbj, group = group.sbj, 
                  grouping = T, ref.batch = 3, eps = 1e-3, print.detail = F)
CdBcor_Data_I.sbj <- CdBDat$X.cor

write.csv(CdBcor_Data_I.sbj, 'CdBcor_Data_I_sbj_LESSR_ref3.csv', row.names = F)
CdBcor_Data_I.sbj <- as.matrix(read.csv('CdBcor_Data_I_sbj_LESSR_ref3.csv', header = T))

# Other data-driven method
# ComBat
library(sva)
sbj.mod <- model.matrix(~as.factor(group.sbj))
CBcor_Data_I <- ComBat(dat = t(Sbj_Dat.intraCor), batch = batch.sbj, 
                        mod = sbj.mod, par.prior = T)
CBcor_Data_I <- t(CBcor_Data_I)
colnames(CBcor_Data_I) <- c(1: ncol(CBcor_Data_I))

# FAbatch
library(bapred)
FADat <- fabatch(Sbj_Dat.log.proc, as.factor(group.sbj), as.factor(batch.sbj))
FAcor_Data_I <- FADat$xadj
colnames(FAcor_Data_I) <- c(1: ncol(FAcor_Data_I))

# QC-based method
QRcor_Data_I.sbj <- QRcor_Data_I[group != 'QC', ]
SRcor_Data_I.sbj <- SRcor_Data_I[group != 'QC', ]

dat.list1 <- list(Uncorrect = Sbj_Dat.log.proc, 
                  ComBat = CBcor_Data_I, 
                  RLSC = QRcor_Data_I.sbj, 
                  SERRF = SRcor_Data_I.sbj, 
                  CordBat = CdBcor_Data_I.sbj)

PresBioEffect(dat.list1, dat.intraCor = Sbj_Dat.intraCor, 
              batch = batch.sbj, group = group.sbj, 
              PCA.print = F, OLDEmet.print = F, OLDEpcc.print = F, 
              PLSDA.print = F, 
              selBioM.type = c(1, 1), BH.adj = c(F, F), 
              r1 = c(0.2, 0.1), r2 = c(2, 1), 
              rank.method = c('both', 'both'), min.batch.num = c(3, 3), 
              useDEforClas = T, DE.ratio = 0.5, 
              save.path0 = paste0(savefile.path0, 'Eval 4/'))

