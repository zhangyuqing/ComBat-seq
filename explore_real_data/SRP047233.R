rm(list=ls())
setwd("~/../Dropbox/Work/ComBat_Seq/real_data_example/")
sapply(c("recount", "DESeq2", "edgeR","dendextend", "ggplot2"), require, character.only=TRUE)


####  Load data
# abstract_search('batch')$project
study_SRP <- "SRP047233"
download_study(study_SRP)
load(file.path(study_SRP, 'rse_gene.Rdata'))
meta_info <- as.data.frame(colData(rse_gene))

# Obtain scaled read counts (from provided coverage counts) to a desired library size
rse_gene <- read_counts(rse_gene, round=TRUE) #scale_counts(rse_gene)  

# Remove genes with all 0 counts
rse_gene <- rse_gene[apply(assays(rse_gene)$counts, 1, function(x){!all(x==0)}), ]

# Batch indicator
# split_info <- strsplit(meta_info$characteristics," ")
# batch <- rep(NA, length(split_info))
# for(i in 1:length(split_info)){
#   batch[i] <- split_info[[i]][grep("batch:", split_info[[i]])+1]
# }
# batch <- as.numeric(gsub(')', '', batch))
batch <- sapply(meta_info$characteristics, function(s){s[grep("batch",s)]})

# Condition indicator
group <- rep(NA, nrow(meta_info))
group[grep("CHD8", meta_info$title)] <- "CHD8"
group[grep("GFP", meta_info$title)] <- "GFP"
group[grep("LacZ", meta_info$title)] <- "LacZ"
group_num <- rep(0, nrow(meta_info)); group_num[group=="CHD8"] <- 1
#group <- sapply(meta_info$characteristics, function(s){s[grep("shRNA target gene",s)]})
  
# Add batch and condition in meta info
colData(rse_gene)$batch <- as.factor(batch)
colData(rse_gene)$condition <- as.factor(group)
colData(rse_gene)$condition_num <- as.factor(group_num)


####  Visualize to see if there is batch effect
p <- plotPCA(DESeqTransform(rse_gene), intgroup="batch")
ggsave(file.path(study_SRP, 'PCA_batch.png'), p, width=11, height=5, units="in")
dev.off(); rm(p)

p2 <- plotPCA(DESeqTransform(rse_gene), intgroup="condition")
ggsave(file.path(study_SRP, 'PCA_cond.png'), p2, width=11, height=5, units="in")
dev.off(); rm(p2)

counts_norm <- scale(t(assays(rse_gene)$counts), center=TRUE, scale=TRUE)
hc <- hclust(dist(counts_norm))
dend <- as.dendrogram(hc)

png(file.path(study_SRP, 'hclust_batch.png'), width=7, height=7, units="in", res=300)
colseq <- sapply(1:nrow(counts_norm), function(i){which(levels(as.factor(batch))==batch[i])})
dend_batch <- color_branches(dend, col=colseq[order.dendrogram(dend)]+1)
plot(dend_batch) 
dev.off(); rm(dend_batch)

png(file.path(study_SRP, 'hclust_cond.png'), width=7, height=7, units="in", res=300)
dend_cond <- color_branches(dend, col=group_num[order.dendrogram(dend)]+2)
plot(dend_cond) 
dev.off(); rm(dend_cond)


####  DE analysis with edgeR
# differential genes relative to biological condition
# dds <- DESeqDataSet(rse_gene, design=~batch+condition_num)
# dds <- DESeq(dds)
# res <- results(dds, name="condition_num_1_vs_0")



# differential genes relative to batch




