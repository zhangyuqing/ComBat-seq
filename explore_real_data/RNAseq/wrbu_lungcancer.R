rm(list=ls())
sapply(c("dplyr", "plyr", "DESeq2", "ggplot2", "reshape2", "gridExtra", "SummarizedExperiment", "scales"), require, character.only=TRUE)
setwd("~/Documents/ComBat_seq/real_data_app/")
data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/WRBU_lungcancer"

## Load data
cts <- readRDS(file.path(data_dir, "count_matrix_BU_WR.rds"))
meta_info <- readRDS(file.path(data_dir, "BU_WR_annotation_file.rds"))
# drop 4790-024 (single sample as its own batch)
cts <- cts[, -grep("4790-024", colnames(cts))]
meta_info <- filter(meta_info, kitnumber!="4790-024")
batch <- factor(meta_info$site)
group <- factor(meta_info$smoking_status)
# remove all 0 genes
tmp <- t(apply(cts,1,floor)); colnames(tmp) <- colnames(cts); cts <- tmp; rm(tmp)
cts <- cts[rowVars(cts)>0, ]


## PCA on logCPM
col_data <- meta_info; rownames(col_data) <- meta_info$sample
seobj <- SummarizedExperiment(assays=list(counts=cts), colData=col_data)
pca_obj <- plotPCA(DESeqTransform(seobj), intgroup=c("site", "smoking_status"))
plt_batch <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=site)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
       title="Color by batch (location)") +
  theme(legend.position="bottom")
plt_cond <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=smoking_status)) + 
  geom_point() + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
       title="Color by condition") +
  theme(legend.position="bottom")
png("wrbu_smoke_PCA.png", width=9, height=5, units="in", res=300)
grid.arrange(plt_batch, plt_cond, ncol=2)
dev.off()

