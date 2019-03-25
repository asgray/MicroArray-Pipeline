# Anthony Gray
# 410.671
# Final Project

#---------------------------------------------------
# METHODS
#---------------------------------------------------
# setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library.setup <- function(){
  source("http://www.bioconductor.org/biocLite.R")
  biocLite("GEOquery")
  library(Biobase)
  library(GEOquery)
  library(limma)
  library(gplots)
  library(raster)
  library(impute)
  library(MASS)
  old.par <- par(mar = c(0, 0, 0, 0))
}

download.gse28521 <- function(){
  # GSE28521
  # download and store data
  gse28521 <- getGEO('GSE28521', destdir='.', GSEMatrix = F)
  gsmplatforms <- lapply(GSMList(gse28521),function(x) {Meta(x)$platform_id})
  gsmnorms <- lapply(GSMList(gse28521),function(x) {Meta(x)$data_processing[1]})
  #length(unique(unlist(gsmplatforms))) == 1
  #[1] TRUE
  # All platforms are the same, == GPL6883
  #length(unique(unlist(gsmnorms))) == 1
  #[1] TRUE
  # all data is normalized
    # Log2 transformation and Quantile normalization uzing the R Lumi package.
  
  gpl6883 <- getGEO('GPL6883', destdir = '.')
  gsmlist <- GSMList(gse28521)
}

# data frame generator ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gen.dat <- function(){
  # extract expression values from GSE, assemble into dataframe
  dat <- c()
  anote <- c()
  labs  <- c()
  # gene ids
  rowlabs <- Table(gsmlist[[1]])[,1]
  for (i in c(1:79)) {
    # values for each sample
    vals <- as.numeric(Table(gsmlist[[i]])[,2])
    # name of each sample
    labs[i] <- gsmlist[[i]]@header$geo_accession
    # add vals as columns
    dat <- cbind(dat, vals)
    
    # retrieves autistic/control and tissue type
    desc <- gsmlist[[i]]@header$characteristics_ch1
    # A or C
    aorc <- toupper(substring(desc[1], 17 ,17))
    # brain region
    reg <- substring(desc[2], 24, 24)
    anote[i] <- paste(aorc, reg, sep = '_')
  }
  # concatenate sample name and tpe code. A = autistic, C = control
  # T = temporal cortex, C = cerebellum, F = frontal cortex
  colnames(dat) <- paste(labs, anote, sep = '_')
  rownames(dat) <- rowlabs
  
  if (! full.run){
    # no correlation at all between cerebullum and other tissues, cerebullum removed
    return(dat[,22:79])
  } else {
    return(dat)
  }
}

# identifying outliers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
outlier.tests <- function(dat){
  # heatmap
  dat.cor <- cor(dat, use = 'pairwise.complete.obs')
  layout(matrix(c(1,1,1,1,1,1,1,1,2,2), 5, 2, byrow = TRUE))
  par(oma = c(5,7,1,1))
  cx <- rev(colorpanel(25,"red","black","blue"))
  leg <- seq(min(dat.cor, na.rm = T), max(dat.cor, na.rm = T), length = 10)
  image(dat.cor, main = "GSE28521 Pearson Correlations", axes = F, col = cx)
  axis(1, at = seq(0, 1, length = ncol(dat.cor)),
       label = dimnames(dat.cor)[[2]], cex.axis = 0.9, las = 2)
  axis(2, at = seq(0, 1, length = ncol(dat.cor)),
       label = dimnames(dat.cor)[[2]], cex.axis = 0.9, las = 2)
  image(as.matrix(leg), col = cx, axes = F)
  tmp <- round(leg, 2)
  axis(1, at = seq(0, 1, length = length(leg)), labels = tmp, cex.axis = 1)
  
  Sys.sleep(1)
  
  # dendrogram
  par(mfrow = c(1,1))
  dd <- dist(scale(dat.cor), method = 'euclidean')
  hc <- hclust(dd, method = 'ward.D2')
  plot(hc, main = 'GSE28521 Cluster Dendrogram')
  
  # CV plot
  cvs <- apply(dat, 2, cv, na.rm = T)
  mns <- apply(dat, 2, mean, na.rm = T)
  plot(x = mns, y = cvs, main = 'CV vs. Mean (GSE28521)',
       xlab = 'Mean Values',
       ylab = 'Coeffecient of Variation', pch = 25, col = 'green',
       bg = 'purple')
  text(mns, cvs, labels = names(mns), pos = 4)
  
  # avg cor plot
  dfmeans <- apply(dat.cor, 1, mean)
  plot(dfmeans, main = 'GSE28521 Mean Correlations',
       ylab = 'Mean Correlation', xlab = '', pch = 24,
       col = 'blue', bg = 'red', axes = F)
  axis(1, at = c(1:ncol(dat)), labels = names(dfmeans), las = 2)
  axis(2)
  
  # boxplot
  boxplot(dfmeans, main = 'GSE28521 Means')
}

# ANOVA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
aov.all.genes4 <- function(x,s1,s2,s3,s4) {
  x1 <- as.numeric(x[s1])
  x2 <- as.numeric(x[s2])
  x3 <- as.numeric(x[s3])
  x4 <- as.numeric(x[s4])
  fac <- c(rep('A',length(x1)), rep('B', length(x2)), rep('C', length(x3)), rep('D', length(x4)))
  a.dat <- data.frame(as.factor(fac),c(x1,x2,x3,x4))
  names(a.dat) <- c('factor', 'express')
  p.out <- summary(aov(express~factor, a.dat))[[1]][1,5]
  #p.out <- summary(aov(express~factor, a.dat))[[1]][1,4]	# use to get F-statistic
  return(p.out)
}

aov.all.genes6 <- function(x,s1,s2,s3,s4, s5, s6) {
  x1 <- as.numeric(x[s1])
  x2 <- as.numeric(x[s2])
  x3 <- as.numeric(x[s3])
  x4 <- as.numeric(x[s4])
  x5 <- as.numeric(x[s5])
  x6 <- as.numeric(x[s6])
  fac <- c(rep('A',length(x1)), rep('B', length(x2)), rep('C', length(x3)), rep('D', length(x4)), rep('E', length(x5)), rep('F', length(x6)))
  a.dat <- data.frame(as.factor(fac),c(x1,x2,x3,x4, x5, x6))
  names(a.dat) <- c('factor', 'express')
  p.out <- summary(aov(express~factor, a.dat))[[1]][1,5]
  #p.out <- summary(aov(express~factor, a.dat))[[1]][1,4]	# use to get F-statistic
  return(p.out)
}

run.anova <- function(dat){
  if (! full.run){
    # number vectors for each category
    # w/out outliers
    AF <- c(1:15)
    CF <- c(16:31)
    AT <- c(32:41)
    CT <- c(42:54)
    
    
    # w/ outliers
    # AF <- c(1:16)
    # CF <- c(17:32)
    # AT <- c(33:45)
    # CT <- c(46:58)
    # contains p-values between 4 groups
    aov.run <- apply(dat, 1, aov.all.genes4, s1 = AF, s2 = CF, s3 = AT, s4 = CT)
    hist(aov.run, xlab = 'p-values', main = 'ANOVA of Temporal and Frontal Expression (GSE28521)')
  } else {
    # w/out outliers
    AC <- c(1:8)
    CC <- c(9:18)
    AF <- c(19:33)
    CF <- c(34:49)
    AT <- c(50:59)
    CT <- c(60:72)
  
    # w/ outliers
    # AC <- c(1:10)
    # CC <- c(11:21)
    # AF <- c(22:37)
    # CF <- c(38:53)
    # AT <- c(54:66)
    # CT <- c(67:79)
    aov.run <- apply(dat, 1, aov.all.genes6, s1 = AC, s2 = CC, s3 = AF, s4 = CF, s5 = AT, s6 = CT)
    hist(aov.run, xlab = 'p-values', main = 'ANOVA of Temporal and Frontal Expression (GSE28521)')
  }
  return(aov.run)
}

# PCA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pca.with.scree <- function(sigs) {
  dat.pca <- prcomp(t(sigs))
  dat.pca.var <- round(dat.pca$sdev^2 / sum(dat.pca$sdev^2)*100,2)
  plot(c(1:length(dat.pca.var)), dat.pca.var, type="b",
       xlab="# of Principal Components", ylab="% variance", pch=21,
       col=1, bg=3, cex=1.5,
       main = 'Scree Plot for Significant GSE28521 Genes')
  return(dat.pca)
}

pca.plots <- function(dat.pca){
  par(mfrow=c(2,2))
  # number vectors for each category
  if(!full.run){
    # w/out outliers
    AF <- c(1:15)
    CF <- c(16:31)
    AT <- c(32:41)
    CT <- c(42:54)
    
    # w/ outliers
    # AF <- c(1:16)
    # CF <- c(17:32)
    # AT <- c(33:45)
    # CT <- c(46:58)
    
    plot(range(dat.pca$x[,1]),range(dat.pca$x[,2]),type="n",xlab='P1',ylab='P2',main='P2 vs. P1')
    points(dat.pca$x[,1][AF], dat.pca$x[,2][AF],col='red', pch=17)
    points(dat.pca$x[,1][AT], dat.pca$x[,2][AT],col='green', pch=17)
    points(dat.pca$x[,1][CF], dat.pca$x[,2][CF],col='blue', pch=19)
    points(dat.pca$x[,1][CT], dat.pca$x[,2][CT],col='purple', pch=19)
    plot(range(dat.pca$x[,3]),range(dat.pca$x[,2]),type="n",xlab='P3',ylab='P2',main='P2 vs. P3')
    points(dat.pca$x[,3][AF], dat.pca$x[,2][AF],col='red', pch=17)
    points(dat.pca$x[,3][AT], dat.pca$x[,2][AT],col='green', pch=17)
    points(dat.pca$x[,3][CF], dat.pca$x[,2][CF],col='blue', pch=19)
    points(dat.pca$x[,3][CT], dat.pca$x[,2][CT],col='purple', pch=19)
    plot(range(dat.pca$x[,1]),range(dat.pca$x[,3]),type="n",xlab='P1',ylab='P3',main='P3 vs. P1')
    points(dat.pca$x[,1][AF], dat.pca$x[,3][AF],col='red', pch=17)
    points(dat.pca$x[,1][AT], dat.pca$x[,3][AT],col='green', pch=17)
    points(dat.pca$x[,1][CF], dat.pca$x[,3][CF],col='blue', pch=19)
    points(dat.pca$x[,1][CT], dat.pca$x[,3][CT],col='purple', pch=19)
    plot(c(0:1), c(0:1), type='n', xaxt = 'n', yaxt = 'n', ann = F, bty = 'n')
    legend(0,1, pch = c(17,17,19,19), col = c('red', 'green', 'blue', 'purple'), cex =1.5,
           legend = c('Frontal Cortex (Autistic)', 'Temporal Cortex (Autistic)', 
                      'Frontal Cortex (Control)','Temporal Cortex (Control)'))
  } else {
    # w/out outliers
    AC <- c(1:8)
    CC <- c(9:18)
    AF <- c(19:33)
    CF <- c(34:49)
    AT <- c(50:59)
    CT <- c(60:72)
    
    # w/ outliers
    # AC <- c(1:10)
    # CC <- c(11:21)
    # AF <- c(22:37)
    # CF <- c(38:53)
    # AT <- c(54:66)
    # CT <- c(67:79)
    plot(range(dat.pca$x[,1]),range(dat.pca$x[,2]),type="n",xlab='P1',ylab='P2',main='P2 vs. P1')
    points(dat.pca$x[,1][AF], dat.pca$x[,2][AF],col='red', pch=17)
    points(dat.pca$x[,1][AT], dat.pca$x[,2][AT],col='green', pch=17)
    points(dat.pca$x[,1][AC], dat.pca$x[,2][AC],col='orange', pch=17)
    points(dat.pca$x[,1][CF], dat.pca$x[,2][CF],col='blue', pch=19)
    points(dat.pca$x[,1][CT], dat.pca$x[,2][CT],col='yellow', pch=19)
    points(dat.pca$x[,1][CC], dat.pca$x[,2][CC],col='purple', pch=19)
    plot(range(dat.pca$x[,3]),range(dat.pca$x[,2]),type="n",xlab='P3',ylab='P2',main='P2 vs. P3')
    points(dat.pca$x[,3][AF], dat.pca$x[,2][AF],col='red', pch=17)
    points(dat.pca$x[,3][AT], dat.pca$x[,2][AT],col='green', pch=17)
    points(dat.pca$x[,3][AC], dat.pca$x[,2][AC],col='orange', pch=17)
    points(dat.pca$x[,3][CF], dat.pca$x[,2][CF],col='blue', pch=19)
    points(dat.pca$x[,3][CT], dat.pca$x[,2][CT],col='yellow', pch=19)
    points(dat.pca$x[,3][CC], dat.pca$x[,2][CC],col='purple', pch=19)
    plot(range(dat.pca$x[,1]),range(dat.pca$x[,3]),type="n",xlab='P1',ylab='P3',main='P3 vs. P1')
    points(dat.pca$x[,1][AF], dat.pca$x[,3][AF],col='red', pch=17)
    points(dat.pca$x[,1][AT], dat.pca$x[,3][AT],col='green', pch=17)
    points(dat.pca$x[,1][AC], dat.pca$x[,3][AC],col='orange', pch=17)
    points(dat.pca$x[,1][CF], dat.pca$x[,3][CF],col='blue', pch=19)
    points(dat.pca$x[,1][CT], dat.pca$x[,3][CT],col='yellow', pch=19)
    points(dat.pca$x[,1][CC], dat.pca$x[,3][CC],col='purple', pch=19)
    plot(c(0:1), c(0:1), type='n', xaxt = 'n', yaxt = 'n', ann = F, bty = 'n')
    legend(0,1, pch = c(17,17,17,19,19,19), cex = 1.5,
           legend = c('Frontal Cortex (Autistic)', 'Temporal Cortex (Autistic)', 
                      'Cerebellum (Autistic)', 'Frontal Cortex (Control)',
                      'Temporal Cortex (Control)','Cerebellum (Control)'),
           col = c('red', 'green', 'orange', 'blue', 'yellow', 'purple'))
  }
  par(mfrow=c(1,1))
}

# LDA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run.lda <- function(dat) {
  if (!full.run) {
    # w/out outliers
    # AF <- c(1:15)
    # CF <- c(16:31)
    # AT <- c(32:41)
    # CT <- c(42:54)
    
    # w/ outliers
    # AF <- c(1:16)
    # CF <- c(17:32)
    # AT <- c(33:45)
    # CT <- c(46:58)
    cl <- c(rep('AF', 15), rep('CF', 16), rep('AT', 10), rep('CT', 13))
    dat.cl <- data.frame(cl, t(dat))
    training.set <- rbind(dat.cl[1:10,], dat.cl[16:26,], dat.cl[32:37,], dat.cl[42:54,])
    test.set <- rbind(dat.cl[11:15,], dat.cl[27:31,], dat.cl[38:41,], dat.cl[50:54,])
    test.cl <- test.set[,1]
    test.set <- test.set[,-1]
    dat.lda <- lda(cl ~ ., training.set)
    dat.pred <- predict(dat.lda, test.set)
    dat.pred$orig.cl <- test.cl
    print(table(dat.pred$class, test.cl))
    return(dat.pred)
  } else {
    # w/out outliers
    # AC <- c(1:8)
    # CC <- c(9:18)
    # AF <- c(19:33)
    # CF <- c(34:49)
    # AT <- c(50:59)
    # CT <- c(60:72)
    
    # w/ outliers
    # AC <- c(1:10)
    # CC <- c(11:21)
    # AF <- c(22:37)
    # CF <- c(38:53)
    # AT <- c(54:66)
    # CT <- c(67:79)
    cl <- c(rep('AC', 8), rep('CC', 10), rep('AF', 15), rep('CF', 16), rep('AT', 10), rep('CT', 13))
    dat.cl <- data.frame(cl, t(dat))
    training.set <- rbind(dat.cl[1:5,], dat.cl[9:15,], dat.cl[19:28,], dat.cl[34:44,], dat.cl[50:56,], dat.cl[60:68,])
    test.set <- rbind(dat.cl[6:8,], dat.cl[16:18,], dat.cl[29:33,], dat.cl[45:49,], dat.cl[57:59,], dat.cl[69:72,])
    test.cl <- test.set[,1]
    test.set <- test.set[,-1]
    dat.lda <- lda(cl ~ ., training.set)
    dat.pred <- predict(dat.lda, test.set)
    dat.pred$orig.cl <- test.cl
    print(table(dat.pred$class, test.cl))
    return(dat.pred)
  }
}

plot.lda <- function(dat.lda){
  lds <- dat.lda$x
  pred.class <- dat.lda$class
  orig.class <- dat.lda$orig.cl
  par(mfrow=c(2,2))
  if (! full.run) {
    plot(lds[,1], lds[,2], col = as.numeric(pred.class), bg = as.numeric(orig.class), pch = 21, xlab = 'LD1', ylab = 'LD2', main = 'LD1 v LD2', cex = 1.5, lwd = 2.5)
    plot(lds[,1], lds[,3], col = as.numeric(pred.class), bg = as.numeric(orig.class), pch = 21, xlab = 'LD1', ylab = 'LD3', main = 'LD1 v LD3', cex = 1.5, lwd = 2.5)
    plot(lds[,3], lds[,2], col = as.numeric(pred.class), bg = as.numeric(orig.class), pch = 21, xlab = 'LD3', ylab = 'LD2', main = 'LD3 v LD2', cex = 1.5, lwd = 2.5)
    plot(c(0:1), c(0:1), type='n', xaxt = 'n', yaxt = 'n', ann = F, bty = 'n')
    legend(0,1, legend = c('Frontal Cortex (Autistic)', 'Frontal Cortex (Control)', 'Temporal Cortex (Autistic)', 'Temporal Cortex (Control)'),
           pch = 21, cex = 1.5, pt.bg = unique(as.numeric(orig.class)), col = unique(as.numeric(orig.class)))
    legend(0, 0.5, legend = c('Actual Class', 'Predicted Class'), pch = 21, pt.bg = c(1,0), col = c(0,1), cex = 1.5)
  } else{
    plot(lds[,1], lds[,2], col = as.numeric(pred.class), bg = as.numeric(orig.class), pch = 21, xlab = 'LD1', ylab = 'LD2', main = 'LD1 v LD2', cex = 1.5, lwd = 2.5)
    plot(lds[,1], lds[,3], col = as.numeric(pred.class), bg = as.numeric(orig.class), pch = 21, xlab = 'LD1', ylab = 'LD3', main = 'LD1 v LD3', cex = 1.5, lwd = 2.5)
    plot(lds[,3], lds[,2], col = as.numeric(pred.class), bg = as.numeric(orig.class), pch = 21, xlab = 'LD3', ylab = 'LD2', main = 'LD3 v LD2', cex = 1.5, lwd = 2.5)
    plot(c(0:1), c(0:1), type='n', xaxt = 'n', yaxt = 'n', ann = F, bty = 'n')
    legend(0,1, legend = c('Cerebellum (Autistic)', 'Cerebellum (Control)','Frontal Cortex (Autistic)',  
                           'Frontal Cortex (Control)', 'Temporal Cortex (Autistic)', 'Temporal Cortex (Control)'), 
           pt.bg = unique(as.numeric(orig.class)), col = unique(as.numeric(orig.class)), pch = 21, cex = 1.5)
    legend(0, 0.3, legend = c('Actual Class', 'Predicted Class'), pch = 21, pt.bg = c(1,0), col = c(0,1), cex = 1.5)
  }
  par(mfrow=c(1,1))
}

# Fold change ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
calc.fold.changes <- function(dat){
  if (!full.run) {
    # w/out outliers
    AF <- c(1:15)
    CF <- c(16:31)
    AT <- c(32:41)
    CT <- c(42:54)
    
    # w/ outliers
    # AF <- c(1:16)
    # CF <- c(17:32)
    # AT <- c(33:45)
    # CT <- c(46:58)
    # find gene mean expression by group
    AF.m <- apply(dat[,AF], 1, mean, na.rm = T)
    CF.m <- apply(dat[,CF], 1, mean, na.rm = T)
    AT.m <- apply(dat[,AT], 1, mean, na.rm = T)
    CT.m <- apply(dat[,CT], 1, mean, na.rm = T)
    # by control/autistic 
    CFT.m <- c(CF.m, CT.m)
    AFT.m <- c(AF.m, AT.m)
    
    # frontal cortex change
    F.fold <- (CF.m - AF.m)
    # temporal cortex change
    T.fold <- (CT.m - AT.m)
    # control/autistic change
    A.fold <- (CFT.m - AFT.m)
    
    # sort by absolute value
    F.fold <- F.fold[order(abs(F.fold))]^2
    T.fold <- T.fold[order(abs(T.fold))]^2
    A.fold <- A.fold[order(abs(A.fold))]^2
    
    # top 5 of each group
    top.5.F <- rev(tail(F.fold, n = 5))
    top.5.T <- rev(tail(T.fold, n = 5))
    top.5.A <- rev(tail(A.fold, n = 5))
    print(top.5.F)
    print(top.5.T)
    print(top.5.A)
    # FT.intersect <- intersect(labels(top.5.F), labels(top.5.T))
    # unique genes
    change.genes <- unique(c(labels(top.5.F), labels(top.5.T), labels(top.5.A)))
    
    # table of significant genes
    in.F <- change.genes %in% labels(top.5.F)
    in.T <- change.genes %in% labels(top.5.T)
    in.A <- change.genes %in% labels(top.5.A)
    change.table <- cbind(change.genes, in.F, in.T, in.A)
    print(change.table)
    return(change.table)
  } else {
    # w/out outliers
    AC <- c(1:8)
    CC <- c(9:18)
    AF <- c(19:33)
    CF <- c(34:49)
    AT <- c(50:59)
    CT <- c(60:72)
    
    # w/ outliers
    # AC <- c(1:10)
    # CC <- c(11:21)
    # AF <- c(22:37)
    # CF <- c(38:53)
    # AT <- c(54:66)
    # CT <- c(67:79)
    # find gene mean expression by group
    AC.m <- apply(dat[,AC], 1, mean, na.rm = T)
    CC.m <- apply(dat[,CC], 1, mean, na.rm = T)
    AF.m <- apply(dat[,AF], 1, mean, na.rm = T)
    CF.m <- apply(dat[,CF], 1, mean, na.rm = T)
    AT.m <- apply(dat[,AT], 1, mean, na.rm = T)
    CT.m <- apply(dat[,CT], 1, mean, na.rm = T)
    # by control/autistic 
    CFTC.m <- c(CC.m, CF.m, CT.m)
    AFTC.m <- c(AC.m, AF.m, AT.m)
    
    # cerebullum change
    C.fold <- (CC.m - AC.m)
    # frontal cortex change
    F.fold <- (CF.m - AF.m)
    # temporal cortex change
    T.fold <- (CT.m - AT.m)
    # control/autistic change
    A.fold <- (CFTC.m - AFTC.m)
    
    # sort by absolute value
    C.fold <- C.fold[order(abs(C.fold))]^2
    F.fold <- F.fold[order(abs(F.fold))]^2
    T.fold <- T.fold[order(abs(T.fold))]^2
    A.fold <- A.fold[order(abs(A.fold))]^2
    
    # top 5 of each group
    top.5.C <- rev(tail(C.fold, n = 5))
    top.5.F <- rev(tail(F.fold, n = 5))
    top.5.T <- rev(tail(T.fold, n = 5))
    top.5.A <- rev(tail(A.fold, n = 5))
    print(top.5.C)
    print(top.5.F)
    print(top.5.T)
    print(top.5.A)
    # FT.intersect <- intersect(labels(top.5.F), labels(top.5.T))
    # unique genes
    change.genes <- unique(c(labels(top.5.C), labels(top.5.F), labels(top.5.T), labels(top.5.A)))
    
    # table of significant genes
    in.C <- change.genes %in% labels(top.5.C)
    in.F <- change.genes %in% labels(top.5.F)
    in.T <- change.genes %in% labels(top.5.T)
    in.A <- change.genes %in% labels(top.5.A)
    change.table <- cbind(change.genes, in.C, in.F, in.T, in.A)
    print(change.table)
    return(change.table)
  }
  
}


#---------------------------------------------------------------------
#*********************************************************************
#*********************************************************************
#---------------------------------------------------------------------

#---------------------------------------------------
# SETUP
#---------------------------------------------------
library.setup()
download.gse28521()

# *********************************FLAG*********
# FLAG        **********************************
full.run <- T #***************FLAG******FLAG****
# FLAG        **********************************
# *********************************FLAG*********

dat <- gen.dat()

#---------------------------------------------------
# OUTLIERS
#---------------------------------------------------
outlier.tests(dat)
outlier.samps <- c('GSM706394_A_C','GSM706410_C_C','GSM706391_A_C',
                   'GSM706455_A_T','GSM706425_A_F','GSM706447_A_T','GSM706448_A_T')
out.cols <- which(colnames(dat) %in% outlier.samps)
dat <- dat[,-out.cols]


#---------------------------------------------------
# CLEANING
#---------------------------------------------------
# filtering low average expression
avgs <- rowMeans(dat)
hist(avgs, main = 'Histogram of Mean Gene Expression Values (GSE28521)', xlab = 'Mean Expression')
# small number of genes are expressed at < 9
# names of low expressed genes
low.genes <- labels(avgs[avgs < 9])

# filtering low cv
vars <- apply(dat, 1, cv)
hist(vars, main = 'Histogram of Coeffecients of Variation (GSE28521)', xlab = 'CV')
# there is not a lot of variation in the dataset in general. Most cv ~ 3 
invar.genes <- labels(vars[vars < 2])

# row number of low genes
low.dat <- which(rownames(dat) %in% low.genes)
# row number of invariant
invar.dat <- which(rownames(dat) %in% invar.genes)
# removing low genes from set
dat <- dat[-low.dat,]
# removing invariant genes from set
dat <- dat[-invar.dat,]
# revisualize data
outlier.tests(dat)


#---------------------------------------------------
# ANOVA
#---------------------------------------------------
aov.run <- run.anova(dat)


#---------------------------------------------------
# SIGNIFICANT VALUES
#---------------------------------------------------
sig.ps <- aov.run[aov.run <= 0.05]
hist(sig.ps, xlab = 'p-values', main = 'ANOVA of Significant (p <= 0.05) Temporal and Frontal Expression (GSE28521)')

siger.ps <- aov.run[aov.run <= 0.01]
hist(siger.ps, xlab = 'p-values', main = 'ANOVA of Significant (p <= 0.01) Temporal and Frontal Expression (GSE28521)')

sigist.ps <- aov.run[aov.run <= 0.003]
hist(sigist.ps, xlab = 'p-values', main = 'ANOVA of Significant (p <= 0.003) Temporal and Frontal Expression (GSE28521)')

# data filtered to only most significant genes
dat.sigs <- dat[which(rownames(dat) %in% labels(siger.ps)),]
# PARTIAL DATA w/ outliers removed
# length(sig.ps)
# [1] 1099
# length(siger.ps)
# [1] 387
# length(sigist.ps)
# [1] 193
# dim(dat.sigs)
# [1] 193  54

# FULL DATA w/out outliers
# length(sig.ps)
# [1] 8910
# length(siger.ps)
# [1] 8601
# length(sigist.ps)
# [1] 8414
# dim(dat.sigs)
# [1] 8414   72

# significant gene dendrogram
sig.cor <- cor(dat.sigs, use = 'pairwise.complete.obs')
sd <- dist(scale(sig.cor), method = 'euclidean')
shc <- hclust(sd, method = 'ward.D2')
plot(shc, main = 'Significant (p < 0.003) GSE28521 Genes')


#---------------------------------------------------
# PCA
#---------------------------------------------------
dat.pca <- pca.with.scree(dat.sigs)
# PARTIAL DATA w/ outliers removed
# summary(dat.pca)
# Importance of components:
#                         PC1    PC2     PC3     PC4     PC5 
# Standard deviation     2.9393 1.4543 1.12438 1.09677 0.84184 
# Proportion of Variance 0.4541 0.1112 0.06644 0.06322 0.03725 
# Cumulative Proportion  0.4541 0.5652 0.63165 0.69487 0.73212

# FULL DATA w/out outliers
# summary(dat.pca)
# Importance of components:
#                            PC1      PC2     PC3     PC4     PC5
# Standard deviation     86.1532 14.29606 8.67807 7.34298 6.70966
# Proportion of Variance  0.9055  0.02493 0.00919 0.00658 0.00549
# Cumulative Proportion   0.9055  0.93049 0.93967 0.94625 0.95175

pca.plots(dat.pca)


#---------------------------------------------------
# LDA
#---------------------------------------------------
dat.lda <- run.lda(dat)
# PARTIAL DATA w/ outliers removed
# test.cl
#     AF AT CF CT
# AF  4  2  0  0
# AT  0  0  1  0
# CF  1  1  1  0
# CT  0  1  3  5

# FULL DATA w/out outliers
# test.cl
#     AC AF AT CC CF CT
# AC  3  0  0  1  0  0
# AF  0  3  2  0  1  2
# AT  0  2  0  0  2  0
# CC  0  0  0  2  0  0
# CF  0  0  1  0  1  1
# CT  0  0  0  0  1  1

plot.lda(dat.lda)


#---------------------------------------------------
# FOLD CHANGES
#---------------------------------------------------
change.table <- calc.fold.changes(dat)

# PARTIAL DATA w/out outliers
# > top.5.F
# ILMN_2062112 ILMN_2153373 ILMN_1783333 ILMN_2388975 ILMN_1729212 
#    1.3339403    1.1249580    1.0945617    0.9995852    0.9974261 
# > top.5.T
# ILMN_1729212 ILMN_1795419 ILMN_1784749 ILMN_2153373 ILMN_2062112 
#    1.3825370    1.3036476    1.2139222    1.1632023    0.8306985
# > top.5.A
# ILMN_1729212 ILMN_2062112 ILMN_1795419 ILMN_1784749 ILMN_2153373 
#     1.382537     1.333940     1.303648     1.213922     1.163202 

#       change.genes   in.F    in.T    in.A   
# [1,] "ILMN_2062112" "TRUE"  "TRUE"  "TRUE"  *
# [2,] "ILMN_2153373" "TRUE"  "TRUE"  "TRUE"  *
# [3,] "ILMN_1783333" "TRUE"  "FALSE" "FALSE"
# [4,] "ILMN_2388975" "TRUE"  "FALSE" "FALSE"
# [5,] "ILMN_1729212" "TRUE"  "TRUE"  "TRUE"  *
# [6,] "ILMN_1795419" "FALSE" "TRUE"  "TRUE" 
# [7,] "ILMN_1784749" "FALSE" "TRUE"  "TRUE" 

# FULL DATA w/out outliers
# > top.5.C
# ILMN_1665207 ILMN_1676822 ILMN_1689515 ILMN_1671149 ILMN_2361603 
#    1.7942087    1.3097996    1.0723089    0.9952567    0.9060924 
# > top.5.F
# ILMN_1787680 ILMN_2062112 ILMN_2153373 ILMN_1783333 ILMN_1659781 
#     2.117442     1.333940     1.124958     1.094562     1.045396 
# > top.5.T
# ILMN_1787680 ILMN_1729212 ILMN_1795419 ILMN_1784749 ILMN_2153373 
#     2.690193     1.382537     1.303648     1.213922     1.163202 
# > top.5.A
# ILMN_1787680 ILMN_1787680 ILMN_1665207 ILMN_1729212 ILMN_2062112 
#     2.690193     2.117442     1.794209     1.382537     1.333940 

#       change.genes   in.C    in.F    in.T    in.A   
# [1,] "ILMN_1665207" "TRUE"  "FALSE" "FALSE" "TRUE"  2
# [2,] "ILMN_1676822" "TRUE"  "FALSE" "FALSE" "FALSE" 1
# [3,] "ILMN_1689515" "TRUE"  "FALSE" "FALSE" "FALSE" 1
# [4,] "ILMN_1671149" "TRUE"  "FALSE" "FALSE" "FALSE" 1
# [5,] "ILMN_2361603" "TRUE"  "FALSE" "FALSE" "FALSE" 1
# [6,] "ILMN_1787680" "FALSE" "TRUE"  "TRUE"  "TRUE"  3
# [7,] "ILMN_2062112" "FALSE" "TRUE"  "FALSE" "TRUE"  2
# [8,] "ILMN_2153373" "FALSE" "TRUE"  "TRUE"  "FALSE" 2
# [9,] "ILMN_1783333" "FALSE" "TRUE"  "FALSE" "FALSE" 2
# [10,]"ILMN_1659781" "FALSE" "TRUE"  "FALSE" "FALSE" 2
# [11,]"ILMN_1729212" "FALSE" "FALSE" "TRUE"  "TRUE"  2
# [12,]"ILMN_1795419" "FALSE" "FALSE" "TRUE"  "FALSE" 1
# [13,]"ILMN_1784749" "FALSE" "FALSE" "TRUE"  "FALSE" 1

#---------------------------------------------------
# GO
#---------------------------------------------------
# performed in scope of calc.fold.changes()
topgenes <- c(C.fold, F.fold, T.fold)
topgenes <- sort(topgenes)
topgenes <- rev(topgenes)
topgenes <- head(topgenes, n = 100)
GO.genes <- head(unique(labels(topgenes)), n = 10)

GO.data <- gpl6883@dataTable@table[gpl6883@dataTable@table$ID %in% GO.genes,]

keeps <- c('ID', 'Entrez_Gene_ID', 'GI', 'Accession','Symbol', 'Synonyms', 'Definition', 'Ontology_Component', 'Ontology_Process', 'Ontology_Function')
GO.clean <- GO.data[keeps]
Entrez.nos <- GO.clean$Entrez_Gene_ID
# > Entrez.nos
# [1]  56942  55854  84417 404216 594855   2621   2918    987  27013  55829