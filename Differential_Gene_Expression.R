library(limma)
library(edgeR)
library(corrplot)
library(pheatmap)   
library(gplots)
library(factoextra)
library(cluster)
library(ggpubr)
library(randomForest)
library(caret)
cortex_data <- read.table(file = ".../GSE76517_fox_cortex_read_counts_by_gene.txt", header=TRUE, row.names = 1)
forebrain_data <- read.table(file = ".../GSE76517_fox_forebrain_read_counts_by_gene.txt", header=TRUE, row.names = 1)
#---------------------------------------------------------------------------------------------------------#
###################
##CORTEX ANALYSIS##
###################
#----------------------------------------------------------------GeneExpression---------------------------#
edger.data.group <- factor(c(rep("Tame", 12), rep("Aggressive", 12)))
d <- DGEList(counts = cortex_data, group = edger.data.group)
d <- calcNormFactors(d)  
plotMDS(d, col=as.numeric(d$samples$group), labels=(c(rep("Tame",12), rep("Aggressive",12)))) ##displays pairwise similarity of  each  sample  in  two  automatically  determined  dimensions
##Compute count per million values by using the cpm  function.
##We filter out genes that do not achieve at least a
##minimum abundance (e.g., 1 read per million) in at least the 
##smallest group of replicates.
cps <- cpm(d)  
k <- rowSums(cps>1)>=2  
d <- d[k,]  
cps <- cpm(d)
design <- model.matrix(~0+group, data=d$samples)  
colnames(design) <- levels(d$samples$group) 
##   Calculate the dispersion estimates, relative to the design matrix:
d <- estimateGLMCommonDisp(d, design)  
d <- estimateGLMTrendedDisp(d, design)  
d <- estimateGLMTagwiseDisp(d, design)

##LIKELIHOOD RATIO TEST method
##Fit the GLM according to the design matrix:     
f <- glmFit(d, design) 
##  Use  plotBCV  function to plot the biological coefficient of variation
plotBCV(d)  
contrasts <- makeContrasts(TamevsAggressive = Tame-Aggressive, levels=design)
# Conduct a likelihood ratio test for each contrast and filter the top differentially expressed features. 
results <- glmLRT(f, contrast = contrasts)
tt <- topTags(results, n=nrow(d))  
a <- tt$table[tt$table$FDR < 0.05 & tt$table$logFC > 0,]  ##upregulated
sum.de <- summary(decideTestsDGE(results, adjust.method = "BH", p.value=0.05))
sig.p <- topTags(results, n=sum.de[1]+sum.de[3], adjust.method = "BH", sort.by = "p.value")
sig.p.d <- topTags(results, n=sum.de[1], adjust.method = "BH", sort.by = "p.value")
all.genes<- sig.p$table #ALL D.E. GENES
downreg <- sig.p.d$table ##downregulated
write.table(downreg, ".../d.118.sig.genes.cortex.txt", sep="\t",quote = FALSE)
##Upregulated genes correlation heatmap:
gene.index <- rownames(a)
corr.cortex.data <- cortex_data[gene.index,]
corr.cortex.data <- t(corr.cortex.data)
corr.mat <- cor(corr.cortex.data, method="spearman")
if (nrow(corr.mat) > 100) stop("Too many rows for heatmap")
fontsize_row = 10 - nrow(corr.mat) / 15
pheatmap(corr.mat, col=greenred(256), main="My Heatmap", cluster_cols=F, 
         fontsize_row=fontsize_row, border_color=NA)
##Downregulated genes correlation heatmap:
gene.index <- rownames(downreg)
corr.cortex.data <- cortex_data[gene.index,]
corr.cortex.data <- t(corr.cortex.data)
corr.mat <- cor(corr.cortex.data, method="spearman")
if (nrow(corr.mat) > 100) stop("Too many rows for heatmap")
fontsize_row = 10 - nrow(corr.mat) / 15
pheatmap(corr.mat, col=greenred(256), main="My Heatmap", cluster_cols=F, 
         fontsize_row=fontsize_row, border_color=NA)
write.table(a, ".../118.sig.genes.cortex.txt", sep="\t",quote = FALSE)
##Smearplot results
summary(de <- decideTestsDGE(results))
detags <- rownames(d)[as.logical(de)]
plotSmear(results, de.tags=detags)
abline(h=c(-1, 1), col="blue")


##LIMMA method
v <- voom(d, design, plot=TRUE) ##mean variance trend
contrasts <- makeContrasts(TamevsAggressive = Tame-Aggressive, levels=design)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contrasts)
efit <- eBayes(vfit)
plotSA(efit)
volcanoplot(efit)
top_genes <- topTable(efit, number = 49, adjust = "fdr")# Summary of results (number of differentially expressed genes)
top_genes <- subset(top_genes, adj.P.Val<0.05)
# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
hist(top_genes$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")
result <- decideTests(efit)
ct <- 1        # choose contrast of interest
volcanoplot(efit, coef=ct, main=colnames(efit)[ct], pch=20,
            highlight=length(which( result[,ct]!=0)), names=rep('+', nrow(efit)))
write.table(top_genes, ".../dlimmamixedger.sig.genes.cortex.txt", sep="\t",quote = FALSE)
plotMD(efit, column=1, status=result[,1], main=colnames(efit)[1], xlim=c(-1,13))
##up/down regulated
up <- which(result[,1] == 1)
upGenes <- efit[up, ]
upTable <- topTable(upGenes, number = 50, adjust.method= "BH",
                    sort.by="p")
write.table(upTable, ".../u.dlimmamixedger.sig.genes.cortex.txt", sep="\t",quote = FALSE)
down <- which(result[,1] == -1)
downGenes <- efit[down, ]
downTable <- topTable(downGenes, number = 50, adjust.method= "BH",
                    sort.by="p")
write.table(downTable, ".../d.dlimmamixedger.sig.genes.cortex.txt", sep="\t",quote = FALSE)
##Upregulated genes correlation heatmap:
gene.index <- rownames(upGenes)
corr.cortex.data <- cortex_data[gene.index,]
corr.cortex.data <- t(corr.cortex.data)
corr.mat <- cor(corr.cortex.data, method="spearman")
if (nrow(corr.mat) > 100) stop("Too many rows for heatmap")
fontsize_row = 10 - nrow(corr.mat) / 15
pheatmap(corr.mat, col=greenred(256), main="My Heatmap", cluster_cols=F, 
         fontsize_row=fontsize_row, border_color=NA)
##Downregulated genes correlation heatmap:
gene.index <- rownames(downGenes)
corr.cortex.data <- cortex_data[gene.index,]
corr.cortex.data <- t(corr.cortex.data)
corr.mat <- cor(corr.cortex.data, method="spearman")
if (nrow(corr.mat) > 100) stop("Too many rows for heatmap")
fontsize_row = 10 - nrow(corr.mat) / 15
pheatmap(corr.mat, col=greenred(256), main="My Heatmap", cluster_cols=F, 
         fontsize_row=fontsize_row, border_color=NA)


##EDGER WITH FILTERING method
edger.data.group <- factor(c(rep("Tame", 12), rep("Aggressive", 12)))
dge <- DGEList(counts = cortex_data, group = edger.data.group)
dge <- calcNormFactors(dge)
##Compute count per million values by using the  cpm  function.
##We recommend filtering out genes that do not achieve at least a
##minimum abundance (e.g., 1 read per million) in at least the 
##smallest group of replicates.
keep <- rowSums(cpm(dge)>1) >= 2
dge <- dge[keep,]
dim(dge)
dge$samples$lib.size <- colSums(dge$counts)
dge <- estimateCommonDisp(dge)
de.com <- exactTest(dge)
summary(decideTestsDGE(de.com, adjust.method = "BH", p.value=0.05))
sum.de <- summary(decideTestsDGE(de.com, adjust.method = "BH", p.value=0.05))
sum.de
sig.p <- topTags(de.com, n=sum.de[1]+sum.de[3], adjust.method = "BH", sort.by = "p.value")
sig.counts <- cortex_data[row.names(sig.p),]
sig.pre <- dge$pseudo.counts[row.names(sig.p),]
sig.genes <- cbind(sig.p$table,sig.pre)
write.table(sig.genes, ".../sig.genes.cortex.txt", sep="\t",quote = FALSE)
##Upregulated
sig.p <- topTags(de.com, n=sum.de[1], adjust.method = "BH", sort.by = "p.value")
row.names(sig.p)
sig.counts <- cortex_data[row.names(sig.p),]
sig.pre <- dge$pseudo.counts[row.names(sig.p),]
sig.genes.down <- cbind(sig.p$table,sig.pre)
write.table(sig.genes.down, ".../down.sig.genes.cortex.txt", sep="\t",quote = FALSE)
##Downregulated
sig.p <- topTags(de.com, n=sum.de[3], adjust.method = "BH", sort.by = "p.value")
row.names(sig.p)
sig.counts <- cortex_data[row.names(sig.p),]
sig.pre <- dge$pseudo.counts[row.names(sig.p),]
sig.genes.up <- cbind(sig.p$table,sig.pre)
write.table(sig.genes.up, ".../up.sig.genes.cortex.txt", sep="\t",quote = FALSE)

#--------------------------------------------------------------------------PCA---------------------------#
Group <- c(rep("Tame",12),rep("Aggressive",12))
cortex_datat <- t(cortex_data)
df_f <- cortex_datat[,apply(cortex_data, 2, var, na.rm=TRUE) != 0]
set.seed(123)
res.km <- kmeans(df_f, 2)
res.pca <- prcomp(df_f,  scale = TRUE)
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(res.km$cluster)
# Add Species groups from the original data sett
ind.coord$Group <- Group
# Data inspection
head(ind.coord)
eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)
ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "Group", size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  stat_mean(aes(color = cluster), size = 4)

#------------------------------------------------------------------RandomForest---------------------------#
##GLM genes
set.seed(123)
index <- row.names(all.genes)
rf.data.cortex <- cortex_data[index,]
cortex.rf <- t(rf.data.cortex)
Type <- c(rep("Tame",12), rep("Aggressive",12))
cortex.rf <- as.data.frame(cortex.rf)
cortex.rf <- cbind(cortex.rf,Type)
cortex.rf$Type <- factor(cortex.rf$Type)
Type <- factor(cortex.rf$Type)
str(cortex.rf)
##MODEL
rf_100 <- randomForest(Type~., data=cortex.rf, ntree=500, importance=TRUE) ##train with 100% of the data
rf_100
# Importance of each predictor.
i_scores <- varImp(rf_100, conditional=TRUE)
rf_imp <- (importance(rf_100,type = 2)) 
varImpPlot(rf_100) 

##voom genes
set.seed(123)
index <- row.names(top_genes)
rf.data.cortex <- cortex_data[index,]
cortex.rf <- t(rf.data.cortex)
Type <- c(rep("Tame",12), rep("Aggressive",12))
cortex.rf <- as.data.frame(cortex.rf)
cortex.rf <- cbind(cortex.rf,Type)
cortex.rf$Type <- factor(cortex.rf$Type)
Type <- factor(cortex.rf$Type)
str(cortex.rf)
##MODEL
rf_100 <- randomForest(Type~., data=cortex.rf, ntree=500, importance=TRUE) ##train with 100% of the data
rf_100
# Importance of each predictor.
i_scores <- varImp(rf_100, conditional=TRUE)
rf_imp <- (importance(rf_100,type = 2)) 
varImpPlot(rf_100) 

#------------------------------------------------------------------------PCA(2)---------------------------#
Group <- c(rep("Tame",12),rep("Aggressive",12))
index <- row.names(all.genes)
rf.data.cortex <- cortex_data[index,]
cortex_datat <- t(rf.data.cortex)
df_f <- cortex_datat[,apply(cortex_datat, 2, var, na.rm=TRUE) != 0]
set.seed(123)
res.km <- kmeans(df_f, 2)
res.pca <- prcomp(df_f,  scale = TRUE)
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(res.km$cluster)
# Add Species groups from the original data sett
ind.coord$Group <- Group
# Data inspection
head(ind.coord)
eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)
ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "Group", size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  stat_mean(aes(color = cluster), size = 4)

#---------------------------------------------------------------------------------------------------------#
######################
##FOREBRAIN ANALYSIS##
######################
#----------------------------------------------------------------GeneExpression---------------------------#
edger.data.group <- factor(c(rep("Tame", 11), rep("Aggressive", 11)))
d <- DGEList(counts = forebrain_data, group = edger.data.group)
d <- calcNormFactors(d)  
plotMDS(d, col=as.numeric(d$samples$group), labels=(c(rep("Tame",11), rep("Aggressive",11)))) ##displays pairwise similarity of  each  sample  in  two  automatically  determined  dimensions
##Compute count per million values by using the cpm  function.
##We filter out genes that do not achieve at least a
##minimum abundance (e.g., 1 read per million) in at least the 
##smallest group of replicates.
cps <- cpm(d)  
k <- rowSums(cps>1)>=2  
d <- d[k,]  
cps <- cpm(d)
design <- model.matrix(~0+group, data=d$samples)  
colnames(design) <- levels(d$samples$group) 
##   Calculate the dispersion estimates, relative to the design matrix:
d <- estimateGLMCommonDisp(d, design)  
d <- estimateGLMTrendedDisp(d, design)  
d <- estimateGLMTagwiseDisp(d, design)

##LIKELIHOOD RATIO TEST method
##Fit the GLM according to the design matrix:     
f <- glmFit(d, design) 
##  Use  plotBCV  function to plot the biological coefficient of variation
plotBCV(d)  
contrasts <- makeContrasts(TamevsAggressive = Tame-Aggressive, levels=design)
# Conduct a likelihood ratio test for each contrast and filter the top differentially expressed features. 
results <- glmLRT(f, contrast = contrasts)
tt <- topTags(results, n=nrow(d))  
a <- tt$table[tt$table$FDR < 0.05 & tt$table$logFC > 0,]  
write.table(a, ".../u.sig.genes.forebrain.LR.txt", sep="\t",quote = FALSE)
sum.de <- summary(decideTestsDGE(results, adjust.method = "BH", p.value=0.05))
sig.p <- topTags(results, n=sum.de[1]+sum.de[3], adjust.method = "BH", sort.by = "p.value")
sig.p.d <- topTags(results, n=sum.de[1], adjust.method = "BH", sort.by = "p.value")
all.genes<- sig.p$table
downreg <- sig.p.d$table
write.table(downreg, ".../d.sig.genes.FB.LR.txt", sep="\t",quote = FALSE)
##Smearplot all genes
summary(de <- decideTestsDGE(results))
detags <- rownames(d)[as.logical(de)]
plotSmear(results, de.tags=detags)
abline(h=c(-1, 1), col="blue")
##Upregulated genes correlation heatmap:
gene.index <- rownames(a)
corr.cortex.data <- forebrain_data[gene.index,]
corr.cortex.data <- t(corr.cortex.data)
corr.mat <- cor(corr.cortex.data, method="spearman")
if (nrow(corr.mat) > 100) stop("Too many rows for heatmap")
fontsize_row = 10 - nrow(corr.mat) / 15
pheatmap(corr.mat, col=greenred(256), main="My Heatmap", cluster_cols=F, 
         fontsize_row=fontsize_row, border_color=NA)
##Downregulated genes correlation heatmap:
gene.index <- rownames(downreg)
corr.cortex.data <- forebrain_data[gene.index,]
corr.cortex.data <- t(corr.cortex.data)
corr.mat <- cor(corr.cortex.data, method="spearman")
if (nrow(corr.mat) > 100) stop("Too many rows for heatmap")
fontsize_row = 10 - nrow(corr.mat) / 15
pheatmap(corr.mat, col=greenred(256), main="My Heatmap", cluster_cols=F, 
         fontsize_row=fontsize_row, border_color=NA)

##LIMMA method
v <- voom(d, design, plot=TRUE) ##mean variance trend
contrasts <- makeContrasts(TamevsAggressive = Tame-Aggressive, levels=design)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contrasts)
efit <- eBayes(vfit)
plotSA(efit)
volcanoplot(efit)
top_genes <- topTable(efit, number = 50, adjust = "fdr")# Summary of results (number of differentially expressed genes)
# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
hist(top_genes$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")
result <- decideTests(efit)
ct <- 1        # choose contrast of interest
volcanoplot(efit, coef=ct, main=colnames(efit)[ct], pch=20,
            highlight=length(which( result[,ct]!=0)), names=rep('+', nrow(efit)))
write.table(top_genes, ".../limmamixedger.sig.genes.forebrain.txt", sep="\t",quote = FALSE)
plotMD(efit, column=1, status=result[,1], main=colnames(efit)[1], xlim=c(-1,13))
##up/down regulated
up <- which(result[,1] == 1)
upGenes <- efit[up, ]
upTable <- topTable(upGenes, number = 50, adjust.method= "BH",
                    sort.by="p")
write.table(upTable, ".../u.dlimmamixedger.sig.genes.cortex.txt", sep="\t",quote = FALSE)
down <- which(result[,1] == -1)
downGenes <- efit[down, ]
downTable <- topTable(downGenes, number = 50, adjust.method= "BH",
                      sort.by="p")
write.table(downTable, ".../d.dlimmamixedger.sig.genes.cortex.txt", sep="\t",quote = FALSE)

##EDGER WITH FILTERING method
edger.data.group <- factor(c(rep("Tame", 11), rep("Aggressive", 11)))
dge <- DGEList(counts = forebrain_data, group = edger.data.group)
dge <- calcNormFactors(dge)
##Compute count per million values by using the  cpm  function.
##We recommend filtering out genes that do not achieve at least a
##minimum abundance (e.g., 1 read per million) in at least the 
##smallest group of replicates.
keep <- rowSums(cpm(dge)>1) >= 2
dge <- dge[keep,]
dim(dge)
dge$samples$lib.size <- colSums(dge$counts)
dge <- estimateCommonDisp(dge)
de.com <- exactTest(dge)
summary(decideTestsDGE(de.com, adjust.method = "BH", p.value=0.05))
sum.de <- summary(decideTestsDGE(de.com, adjust.method = "BH", p.value=0.05))
sum.de
sig.p <- topTags(de.com, n=sum.de[1]+sum.de[3], adjust.method = "BH", sort.by = "p.value")
sig.counts <- cortex_data[row.names(sig.p),]
sig.pre <- dge$pseudo.counts[row.names(sig.p),]
sig.genes <- cbind(sig.p$table,sig.pre)
write.table(sig.genes, ".../sig.genes.forebrain.txt", sep="\t",quote = FALSE)

#--------------------------------------------------------------------------PCA---------------------------#
Group <- c(rep("Tame",11),rep("Aggressive",11))
forebrain_data <- t(forebrain_data)
df_f <- forebrain_data[,apply(forebrain_data, 2, var, na.rm=TRUE) != 0]
set.seed(123)
res.km <- kmeans(df_f, 2)
res.pca <- prcomp(df_f,  scale = TRUE)
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(res.km$cluster)
# Add Species groups from the original data sett
ind.coord$Group <- Group
# Data inspection
head(ind.coord)
eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)
ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "Group", size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  stat_mean(aes(color = cluster), size = 4)

#------------------------------------------------------------------RandomForest---------------------------#
##GLM genes
set.seed(123)
index <- row.names(all.genes)
rf.data.fb <- forebrain_data[index,]
fb.rf <- t(rf.data.fb)
Type <- c(rep("Tame",11), rep("Aggressive",11))
fb.rf <- as.data.frame(fb.rf)
fb.rf <- cbind(fb.rf,Type)
fb.rf$Type <- factor(fb.rf$Type)
Type <- factor(fb.rf$Type)
str(fb.rf)
##MODEL
fb.rf_100 <- randomForest(Type~., data=fb.rf, ntree=500, importance=TRUE) ##train with 100% of the data
fb.rf_100
# Importance of each predictor.
i_scores <- varImp(fb.rf_100, conditional=TRUE)
rf_imp <- (importance(fb.rf_100,type = 2)) 
varImpPlot(fb.rf_100) 

#------------------------------------------------------------------------PCA(2)---------------------------#
Group <- c(rep("Tame",11),rep("Aggressive",11))
index <- row.names(all.genes)
rf.data.fb <- forebrain_data[index,]
forebrain_datat <- t(rf.data.fb)
df_f <- forebrain_datat[,apply(forebrain_datat, 2, var, na.rm=TRUE) != 0]
set.seed(123)
res.km <- kmeans(df_f, 2)
res.pca <- prcomp(df_f,  scale = TRUE)
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(res.km$cluster)
# Add Species groups from the original data sett
ind.coord$Group <- Group
# Data inspection
head(ind.coord)
eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)
ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "Group", size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  stat_mean(aes(color = cluster), size = 4)
