rm(list=ls())
setwd("~/../Dropbox/Work/ComBat_Seq/real_data_example/")
sapply(c("recount", "DESeq2", "dendextend"), require, character.only=TRUE)


####  Load data
# abstract_search('batch')$project
study_SRP <- "SRP047233"

# Download dataset
download_study(study_SRP)
load(file.path(study_SRP, 'rse_gene.Rdata'))
meta_info <- read.table(file.path(study_SRP, 'SRP047233.tsv'), sep='\t', header=TRUE, as.is=TRUE)
identical(colnames(assay(rse_gene)), meta_info$run)

# Obtain read counts with/without scaling
rse_gene <- read_counts(rse_gene, round=TRUE)  #scale_counts(rse_gene)

# Remove genes with all 0 counts
rse_gene <- rse_gene[apply(assay(rse_gene), 1, function(x){!all(x==0)}), ]

# Batch indicator
split_info <- strsplit(meta_info$characteristics," ")
batch <- rep(NA, length(split_info))
for(i in 1:length(split_info)){
  batch[i] <- split_info[[i]][grep("batch:", split_info[[i]])+1]
}
batch <- as.numeric(gsub(')', '', batch))

# Condition indicator
group <- rep(NA, nrow(meta_info))
group[grep("CHD8", meta_info$title)] <- "CHD8"
group[grep("GFP", meta_info$title)] <- "GFP"
group[grep("LacZ", meta_info$title)] <- "LacZ"
group_num <- rep(0, nrow(meta_info)); group_num[group=="CHD8"] <- 1

# Add batch and condition in meta info
colData(rse_gene)$batch <- as.factor(batch)
colData(rse_gene)$condition <- as.factor(group)
colData(rse_gene)$condition_num <- as.factor(group_num)


####  Visualize to see if there is batch effect
plotPCA(DESeqTransform(rse_gene), intgroup="batch")
#plotPCA(DESeqTransform(rse_gene), intgroup="condition_num")

counts_norm <- scale(t(assay(rse_gene)), center=TRUE, scale=TRUE)
hc <- hclust(dist(counts_norm))
dend <- as.dendrogram(hc)
dend <- color_branches(dend, col=batch[order.dendrogram(dend)]+1)
plot(dend) 


####  DE analysis with DESeq2
# differential genes relative to biological condition




