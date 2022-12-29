# Title: Data I BEC and evaluation
# Author: Fanjing Guo
# Date: 2022.12.28

# --------------------- Data pre-processing and batch effect correction ------------------- #
# LC-MS/MS metabolomics dataset from an atherosclerosis study
X <- read.csv("Data I_Raw_filtered.csv", header = TRUE)
X_info <- X[, c(1: 4)]
X_Dat <- X[, -c(1: 4)]

# add HTN group information
HTN.label <- rep(NA, nrow(X_info))
X_info <- cbind(X_info, HTN.label)
colnames(X_info)[5] <- 'HTN' 
X_clinic <- read.csv("Data I_clinical.csv", header = TRUE)
patID <- X_clinic$Pat.ID
for (i in c(1: length(patID))) {
  X_info$HTN[which(X_info$sample == patID[i])] <- X_clinic$HTN[i]
}
delsampIdx <- c()
HTN.isna.Idx <- which(is.na(X_info$HTN))
for (i in c(1: length(HTN.isna.Idx))) {
  Samp.type <- X_info$type[HTN.isna.Idx[i]]
  if (Samp.type == 'S') {
    delsampIdx <- c(delsampIdx, HTN.isna.Idx[i])
  }
  else if (Samp.type == 'QC') {
    X_info$HTN[HTN.isna.Idx[i]] <- 'QC'
  }
  else {
    samp <- X_info$sample[HTN.isna.Idx[i]]
    repsamp <- strsplit(samp, '*', fixed = T)[[1]][1]
    repsamp.id <- which(X_info$sample == repsamp)
    if (is.na(X_info$HTN[repsamp.id])) {
      delsampIdx <- c(delsampIdx, HTN.isna.Idx[i])
    }
    else {
      X_info$HTN[HTN.isna.Idx[i]] <- X_info$HTN[repsamp.id]
    }
  }
}
X_Dat <- X_Dat[-delsampIdx, ]
X_info <- X_info[-delsampIdx, ]

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



