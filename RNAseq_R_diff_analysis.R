# source("http://bioconductor.org/biocLite.R") 
# this will run the installation script
# biocLite('edgeR')
# biocLite('DESeq2')

library(DESeq2)
library(edgeR)
library(limma)
library(Biobase)
library(gplots)

load('/Users/wjidea/Google\ Drive/Graduate_Study/9_Stats_Bioinfo/9_Bioinfo_Study/NGS_2015/NGS_Reproducibility/NGS_Reproducibility.rdata')

par(mfrow=c(1,3))
# Plot 2 replicate dataset
eset <- bottomly.2reps
cpm.mat <- log(cpm(exprs(eset)))
mean.vec <- apply(cpm.mat, 1, mean)
sdvec <- apply(cpm.mat, 1, sd)
plot(mean.vec, sdvec, pch=".", main="2 replicates", ylab="sd", xlab="Average logCPM",ylim=c(0,2))

# Plot 5 replicate dataset
eset <- bottomly.5reps
cpm.mat <- log(cpm(exprs(eset)))
mean.vec <- apply(cpm.mat, 1, mean)
sdvec <- apply(cpm.mat, 1, sd)
plot(mean.vec, sdvec, pch=".", main="5 replicates", ylab="sd", xlab="Average logCPM",ylim=c(0,2))

# Plot 10 replicate dataset
eset <- bottomly.eset
cpm.mat <- log(cpm(exprs(eset)))
mean.vec <- apply(cpm.mat, 1, mean)
sdvec <- apply(cpm.mat, 1, sd)
plot(mean.vec, sdvec, pch=".", main="10 replicates", ylab="sd", xlab="Average logCPM",ylim=c(0,2))

par(mfrow=c(1,3))

# DE analyses

# Create DESeq2 datasets
dds <- DESeqDataSetFromMatrix(countData = exprs(bottomly.eset), 
                              colData = pData(bottomly.eset), 
                              design = ~ strain )
dds <- DESeq(dds)

dds.5rep <- DESeqDataSetFromMatrix(countData = exprs(bottomly.5reps), 
                                   colData = pData(bottomly.5reps), 
                                   design = ~ strain )
dds.5rep <- DESeq(dds.5rep)

dds.2rep <- DESeqDataSetFromMatrix(countData = exprs(bottomly.2reps), 
                                   colData = pData(bottomly.2reps), 
                                   design = ~ strain )
dds.2rep <- DESeq(dds.2rep)


# Plot dispersion estimates
plotDispEsts(dds.2rep)
plotDispEsts(dds.5rep)
plotDispEsts(dds)


# EdgeR
dge <- DGEList(counts=exprs(bottomly.eset), group=pData(bottomly.eset)$strain)
# Normalize by total count
dge <- calcNormFactors(dge)

# Create the contrast matrix
design.mat <- model.matrix(~ 0 + dge$samples$group)
colnames(design.mat) <- levels(dge$samples$group)

# Estimate dispersion parameter for GLM
dge <- estimateGLMCommonDisp(dge, design.mat)
dge <- estimateGLMTrendedDisp(dge, design.mat, method="power")
dge <- estimateGLMTagwiseDisp(dge,design.mat)

# Do it all over again for 5 replicates
dge.5reps <- DGEList(counts=exprs(bottomly.5reps), group=pData(bottomly.5reps)$strain)
dge.5reps <- calcNormFactors(dge.5reps)
design.mat <- model.matrix(~ 0 + dge.5reps$samples$group)
colnames(design.mat) <- levels(dge.5reps$samples$group)

dge.5reps <- estimateGLMCommonDisp(dge.5reps, design.mat)
dge.5reps <- estimateGLMTrendedDisp(dge.5reps, design.mat, method="power")
dge.5reps <- estimateGLMTagwiseDisp(dge.5reps,design.mat)

# Do it all over again for 2 replicates
dge.2reps <- DGEList(counts=exprs(bottomly.2reps), group=pData(bottomly.2reps)$strain)
dge.2reps <- calcNormFactors(dge.2reps)
design.mat <- model.matrix(~ 0 + dge.2reps$samples$group)
colnames(design.mat) <- levels(dge.2reps$samples$group)

dge.2reps <- estimateGLMCommonDisp(dge.2reps, design.mat)
dge.2reps <- estimateGLMTrendedDisp(dge.2reps, design.mat, method="power")
dge.2reps <- estimateGLMTagwiseDisp(dge.2reps,design.mat)

# Plot mean-variance
plotBCV(dge.2reps,main = "2reps", ylim=c(0,1.5))
plotBCV(dge.5reps, main = "5reps", ylim=c(0,1.5))
plotBCV(dge, main = "10reps", ylim=c(0,1.5))

# Limma voom

par(mfrow=c(1,1))

# Create design matrix
design <- model.matrix(~ pData(bottomly.eset)$strain)

# Apply voom transformation
nf <- calcNormFactors(bottomly.eset)
v <- voom(exprs(bottomly.eset), design, lib.size=colSums(exprs(bottomly.eset))*nf, normalize.method="quantile", plot=TRUE)

# Do same for 5 replicate dataset
design <- model.matrix(~ pData(bottomly.5reps)$strain)
nf <- calcNormFactors(bottomly.5reps)
v.5reps <- voom(exprs(bottomly.5reps), design, lib.size=colSums(exprs(bottomly.5reps))*nf, 
                normalize.method="quantile", plot=TRUE)

# Do same for 2 replicates dataset
design <- model.matrix(~ pData(bottomly.2reps)$strain)
nf <- calcNormFactors(bottomly.2reps)
v.2reps <- voom(exprs(bottomly.2reps), design, lib.size=colSums(exprs(bottomly.2reps))*nf, 
                normalize.method="quantile", plot=TRUE)

# define the threshold
p.threshold <- 0.05

## edgeR ##
# Design matrix
design.mat <- model.matrix(~ 0 + dge$samples$group)
colnames(design.mat) <- c("C57BL", "DBA")

# Model fitting
fit.edgeR <- glmFit(dge, design.mat)

# Differential expression
contrasts.edgeR <- makeContrasts(C57BL - DBA, levels=design.mat)
lrt.edgeR <- glmLRT(fit.edgeR, contrast=contrasts.edgeR)

# Access results tables
edgeR_results <- lrt.edgeR$table
sig.edgeR <- decideTestsDGE(lrt.edgeR, adjust.method="BH", p.value = p.threshold)
genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR != 0)]

## DESeq2 ##
contrast.deseq2 <- list("strainC57BL.6J", "strainDBA.2J")
deseq2_results <- results(dds, contrast=contrast.deseq2)
deseq2_results$threshold <- as.logical(deseq2_results$padj < p.threshold)
genes.deseq <- row.names(deseq2_results)[which(deseq2_results$threshold)]

## voom-limma ##
# Create design matrix
design <- model.matrix(~ pData(bottomly.eset)$strain)

# Usual limma pipeline
fit.voom <- lmFit(v, design)
fit.voom <- eBayes(fit.voom)

voom_results <- topTable(fit.voom, coef=2,  adjust="BH", number = nrow(exprs(bottomly.eset)))
voom_results$threshold <- as.logical(voom_results$adj.P.Val < p.threshold)
genes.voom <- row.names(voom_results)[which(voom_results$threshold)]

intersect(genes.edgeR,genes.voom)
# Venn diagram
venn(list(edgeR = genes.edgeR, DESeq2 = genes.deseq, voom = genes.voom))

## edgeR ##
# Design matrix
design.mat <- model.matrix(~ 0 + dge.2reps$samples$group)
colnames(design.mat) <- c("C57BL", "DBA")

# Model fitting
fit.edgeR <- glmFit(dge.2reps, design.mat)

# Differential expression
contrasts.edgeR <- makeContrasts(C57BL - DBA, levels=design.mat)
lrt.edgeR <- glmLRT(fit.edgeR, contrast=contrasts.edgeR)

# Access results tables
edgeR_results_2reps <- lrt.edgeR$table
sig.edgeR.2reps <- decideTestsDGE(lrt.edgeR, adjust.method="BH", p.value = p.threshold)
genes.edgeR.2reps <- row.names(edgeR_results_2reps)[which(sig.edgeR.2reps != 0)]

## DESeq2 ##
contrast.deseq2 <- list("strainC57BL.6J", "strainDBA.2J")
deseq2_results_2reps <- results(dds.2rep, contrast=contrast.deseq2)
deseq2_results_2reps$threshold <- as.logical(deseq2_results_2reps$padj < p.threshold)
genes.deseq.2reps <- row.names(deseq2_results_2reps)[which(deseq2_results_2reps$threshold)]

## voom-limma ##
# Create design matrix
design <- model.matrix(~ pData(bottomly.2reps)$strain)

# Usual limma pipeline
fit.voom <- lmFit(v.2reps, design)
fit.voom <- eBayes(fit.voom)

voom_results_2reps <- topTable(fit.voom, coef=2,  adjust="BH", number = nrow(exprs(bottomly.2reps)))
voom_results_2reps$threshold <- as.logical(voom_results_2reps$adj.P.Val < p.threshold)
genes.voom.2reps <- row.names(voom_results_2reps)[which(voom_results_2reps$threshold)]

# Overlapping genes
length(genes.deseq.2reps)
length(genes.edgeR.2reps)
length(genes.voom.2reps)

venn(list(edgeR = genes.edgeR.2reps, DESeq2 = genes.deseq.2reps, voom = genes.voom.2reps))
