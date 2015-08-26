# source("http://bioconductor.org/biocLite.R") 
# this will run the installation script
# biocLite('edgeR')
# biocLite('DESeq2')

library(DESeq2)
library(edgeR)
library(limma)
library(Biobase)

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

