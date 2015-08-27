# gene annotation pathway analysis
biocLite(c("pathview", "edgeR", "gage"))
library(pathview)
library(edgeR)
library(gage)
rnaseq_URL <- "https://github.com/vsbuffalo/rna-seq-example/raw/master/results/raw-counts.txt"
rnaseq_data <- read.table(rnaseq_URL, header=TRUE)
colnames(rnaseq_data) <- c("WT_1", "WT_2", "hy_1", "hy_2")


data(korg)
head(korg)

# Finding your species in KEGG
org <- "arabidopsis thaliana"
species <- unlist(sapply(1:ncol(korg), function(i) {
  agrep(org, korg[, i])
}))
korg[species, 1, drop = F]

# org2 <- "fusarium oxysporum"
# species <- unlist(sapply(1:ncol(korg), function(i) {
#   agrep(org2, korg[, i])
# }))
# korg[species, 1, drop = F]


# Creating the KEGG dataset for GAGE analysis
kegg_arab <- kegg.gsets("ath")
kegg.gs <- kegg_arab$kg.sets[kegg_arab$sigmet.idx]


# 1 remove zero counts in the dataset
rnaseq_counts <- rnaseq_data # rename the dataset 

dim(rnaseq_counts)
non_zero <- rowSums(rnaseq_counts) != 0
rnaseq_counts <- rnaseq_counts[non_zero,]
# number of genes left after removing zero counts 

dim(rnaseq_counts)

# 2. Normalize library sizes

libsizes <- colSums(rnaseq_counts)
size_factor <- libsizes / exp(mean(log(libsizes)))
norm_counts <- t(t(rnaseq_counts) / size_factor)
range(norm_counts)

# 3. Variance stabilizing transformation
norm_counts <- log2(norm_counts + 8)
range(norm_counts)

# 4. Define reference and experiment samples
ref_idx <- 1:2 #wt
samp_idx <- 3:4 #hy mutant 

# 5. Enrichment analysis
native_kegg_fc <- gage(norm_counts, 
                       gsets = kegg.gs, 
                       ref = ref_idx, 
                       samp = samp_idx, 
                       compare ="unpaired")
# 6 output 
wt_hy_sig_kegg <-sigGeneSet(native_kegg_fc, outname="wt_hy.kegg")
write.table(rbind(wt_hy_sig_kegg$greater,
                  wt_hy_sig_kegg$less), 
            file = "wt_hy_sig_kegg.txt", sep = "\t")

# Combine GAGE and EdgeR analysis

# 1. Make the study design as follows:
targets <- data.frame(sample_name=colnames(rnaseq_data),Group=rep(c("WT","hy"),each=2))

# 2. Create DGElist object.
d <- DGEList(counts=rnaseq_data, group=targets$Group) # Constructs DGEList object

# 3. Remove lowly expressed genes
#before filtering:
dim(d)
keep <- rowSums(cpm(d)> 10) >= 2
d <- d[keep,]
#after filtering:
dim(d)

# 4. Normalize the data
d$samples$lib.size <- colSums(d$counts)
d <- calcNormFactors(d)
d

# 5. multidimentioanl scaling plot
plotMDS(d, labels = targets$Group, 
        col = c("darkgreen","blue")[factor(targets$Group)])

# 6. Estimate the dispersion

d1 <- estimateCommonDisp(d, verbose=T)

d1 <- estimateTagwiseDisp(d1)

plotBCV(d1)

# 7. Differentail gene expression

de.com <- exactTest(d1, pair = c("WT", "hy"))
topTags(de.com, n = 5)

FDR <- p.adjust(de.com$table$PValue, method="BH")
de1 <- decideTestsDGE(de.com, adjust.method="BH", p.value=0.05)
summary(de1)

# 8. Making the data suitable for pathway analysis

isDE <- as.logical(de1) # covert to DE set to true/false set

DEnames <- rownames(d)[isDE] # get the DE gene names
edger.fc <- de.com$table$logFC # get the log fold change
names(edger.fc) <- rownames(de.com$table) # assign row names to fold change 
exp.fc <- edger.fc
head(exp.fc)

# Getting DE gene set and fold change

DE_foldchange <- edger.fc[DEnames]

length(DE_foldchange)

# 9. GAGE analysis

fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
head(fc.kegg.p$greater[,1:5], 4)

# Significance test

wt_hy_sig_kegg <-sigGeneSet(fc.kegg.p, outname="wt_hy.kegg")

# Visualize the data using Pathview

log_fc= norm_counts[, samp_idx]-rowMeans(norm_counts[, ref_idx])

# 1. Selecting the path ids for the upregulated set.

greater_set <- native_kegg_fc$greater[, "q.val"] < 0.1 &
  !is.na(native_kegg_fc$greater[, "q.val"])
greater_ids <- rownames(native_kegg_fc$greater)[greater_set]
head(greater_ids)
# 2. Selecting path ids for down-regulated set

less_set <- native_kegg_fc$less[, "q.val"] < 0.1 &
  !is.na(native_kegg_fc$less[,"q.val"])
less_ids <- rownames(native_kegg_fc$less)[less_set]
# 3. Combine up and down-regulated path ids.

combine_ids <- substr(c(greater_ids, less_ids), 1, 8)
head(combine_ids)
# 4. Visualization

# Here we are going to get the first three pathways and visualize using pathview.

pv.out.list <- sapply(combine_ids[1:3], function(pid) pathview(
  gene.data =  exp.fc, pathway.id = pid,
  species = "ath", out.suffix="Wt_hy",
  gene.idtype="KEGG"))
# 5. Mapping DE genes to individual pathways

# We can map genes that show statistically significant changes to individual pathways. As an example, we can take DNA replication pathway and see which of the DE genes from EdgeR is present.

pv_replication <- pathview(gene.data = DE_foldchange, 
                           gene.idtype = "KEGG", 
                           pathway.id = combine_ids[3], 
                           species = "ath", 
                           out.suffix = "DNA_replication", 
                           keys.align = "y", 
                           kegg.native = T, 
                           match.data = T, 
                           key.pos = "topright")
# See how the output look like.

head(pv_replication)



