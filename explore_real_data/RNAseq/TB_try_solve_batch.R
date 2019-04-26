rm(list=ls())
output_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/TB"
setwd(output_dir)
sapply(c("SummarizedExperiment", "sva", "DESeq2", "ggplot2", "reshape2", "gridExtra", "scales"), require, character.only=TRUE)



####  Load data
seobj <- readRDS(file.path(output_dir, "combined.rds"))
cts <- assays(seobj)$counts
cts_combat <- assays(seobj)$combat
# batch indicator
batch <- colData(seobj)$SequencingBatch
# condition indicator
group <- colData(seobj)$Label



#### On combat
seobj_combat <- SummarizedExperiment(assays=list(counts=cts_combat), colData=colData(seobj))
pca_obj_combat <- plotPCA(DESeqTransform(seobj_combat), intgroup=c("SequencingBatch", "Label")) 
p_batch <- ggplot(pca_obj_combat$data, aes(x=PC1, y=PC2, color=SequencingBatch)) +
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_combat$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_combat$plot_env$percentVar[2])),
       title="PCA: samples colored by Batch")
p_cond <- ggplot(pca_obj_combat$data, aes(x=PC1, y=PC2, color=Label)) +
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_combat$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_combat$plot_env$percentVar[2])),
       title="PCA: samples colored by Condition")
plt_combat <- grid.arrange(p_batch, p_cond, ncol=2)
dev.off()



#### Try separate studies & ref-combat
## separate studies
seobj_africa <- seobj[, seobj$Dataset=="Africa"]
seobj_brazil <- seobj[, seobj$Dataset=="Brazil"]
seobj_G6 <- seobj[, seobj$Dataset=="G6"]
seobj_india <- seobj[, seobj$Dataset=="India"]


## in Africa
seobj_africa_logtpm <- SummarizedExperiment(assays=list(logtpm=assays(seobj_africa)$logtpm), colData=colData(seobj_africa))
pca_obj_africa <- plotPCA(DESeqTransform(seobj_africa_logtpm), intgroup=c("SequencingBatch", "Label", "Sex")) 
ggplot(pca_obj_africa$data, aes(x=PC1, y=PC2, color=Sex)) +
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_africa$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_africa$plot_env$percentVar[2])))
#clearly separate by sex; include sex as covariate


## in Brazil
seobj_brazil_logtpm <- SummarizedExperiment(assays=list(logtpm=assays(seobj_brazil)$logtpm), colData=colData(seobj_brazil))
pca_obj_brazil <- plotPCA(DESeqTransform(seobj_brazil_logtpm), intgroup=c("SequencingBatch", "Label", "Sex")) 
ggplot(pca_obj_brazil$data, aes(x=PC1, y=PC2, color=SequencingBatch, shape=Sex)) +
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_brazil$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_brazil$plot_env$percentVar[2])))
# separate both by batch and by sex

brazil_logtpm <- assays(seobj_brazil)$logtpm
brazil_batch <- as.factor(as.character(seobj_brazil$SequencingBatch))
brazil_group <- seobj_brazil$Label
brazil_sex_cov <- seobj_brazil$Sex

nonprg_batch1 <- brazil_logtpm[, brazil_batch=="Brazil_1" & brazil_group=="Non-progressor"]
nonprg_batch2 <- brazil_logtpm[, brazil_batch=="Brazil_2" & brazil_group=="Non-progressor"]
nonprg_merged <- cbind(nonprg_batch1, nonprg_batch2)
nonprg_cov <- as.factor(c(as.character(brazil_sex_cov[brazil_batch=="Brazil_1" & brazil_group=="Non-progressor"]),
                          as.character(brazil_sex_cov[brazil_batch=="Brazil_2" & brazil_group=="Non-progressor"])))
nonprg_batch <- c(rep(1, ncol(nonprg_batch1)), rep(2, ncol(nonprg_batch2)))

keep1 <- apply(nonprg_batch1, 1, function(x){var(x)!=0})
keep2 <- apply(nonprg_batch2, 1, function(x){var(x)!=0})
keep <- keep1 & keep2
nonprg_merged <- nonprg_merged[keep, ]
brazil_logtpm <- brazil_logtpm[keep, ]

nonprg_merged_adj <- ComBat(nonprg_merged, batch=nonprg_batch, mod=model.matrix(~nonprg_cov), ref.batch=1)
brazil_logtpm[, brazil_batch=="Brazil_2" & brazil_group=="Non-progressor"] <- nonprg_merged_adj[, nonprg_batch==2]

seobj_brazil_logtpm_adj <- SummarizedExperiment(assays=list(logtpm=brazil_logtpm), colData=colData(seobj_brazil))
pca_obj_brazil_adj <- plotPCA(DESeqTransform(seobj_brazil_logtpm_adj), intgroup=c("SequencingBatch", "Label", "Sex")) 
ggplot(pca_obj_brazil_adj$data, aes(x=PC1, y=PC2, color=SequencingBatch, shape=Sex)) +
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_brazil_adj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_brazil_adj$plot_env$percentVar[2])))


## in G6
seobj_G6_logtpm <- SummarizedExperiment(assays=list(logtpm=assays(seobj_G6)$logtpm), colData=colData(seobj_G6))
pca_obj_G6 <- plotPCA(DESeqTransform(seobj_G6_logtpm), intgroup=c("SequencingBatch", "Label", "Sex")) 
ggplot(pca_obj_G6$data, aes(x=PC1, y=PC2))+#, color=Sex)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_G6$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_G6$plot_env$percentVar[2])))
# seems okay, slightly separate by sex


## in India
seobj_india_logtpm <- SummarizedExperiment(assays=list(logtpm=assays(seobj_india)$logtpm), colData=colData(seobj_india))
pca_obj_india <- plotPCA(DESeqTransform(seobj_india_logtpm), intgroup=c("SequencingBatch", "Label", "Sex")) 
ggplot(pca_obj_india$data, aes(x=PC1, y=PC2)) + #, color=Label)) +
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_india$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_india$plot_env$percentVar[2])))
# fine as-is


## Merge all 4 studies
seobj_merged <- cbind(seobj_africa_logtpm[keep, ], seobj_brazil_logtpm_adj, seobj_G6_logtpm[keep, ], seobj_india_logtpm[keep, ])

merged_logtpm_adj <- ComBat(assays(seobj_merged)$logtpm, batch=seobj_merged$Dataset, 
                            mod=model.matrix(~seobj_merged$Label))
merged_logtpm_adj <- ComBat(merged_logtpm_adj, batch=seobj_merged$Sex, 
                            mod=model.matrix(~seobj_merged$Label))

  
rds_obj_sep_ref <- SummarizedExperiment(assays=list(merged_logtpm_adj), colData=colData(seobj_merged))
pca_obj_sep_ref <- plotPCA(DESeqTransform(rds_obj_sep_ref), intgroup=c("SequencingBatch", "Label", "Sex", "Dataset")) 
p2_batch <- ggplot(pca_obj_sep_ref$data, aes(x=PC1, y=PC2, color=SequencingBatch)) +
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_sep_ref$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_sep_ref$plot_env$percentVar[2])),
       title="PCA: samples colored by Batch")
p2_cond <- ggplot(pca_obj_sep_ref$data, aes(x=PC1, y=PC2, color=Label)) +
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_sep_ref$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_sep_ref$plot_env$percentVar[2])),
       title="PCA: samples colored by Condition")
p2_studies <- ggplot(pca_obj_sep_ref$data, aes(x=PC1, y=PC2, color=Dataset)) +
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_sep_ref$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_sep_ref$plot_env$percentVar[2])),
       title="PCA: samples colored by Studies")
p2_sex <- ggplot(pca_obj_sep_ref$data, aes(x=PC1, y=PC2, color=Sex)) +
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_sep_ref$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_sep_ref$plot_env$percentVar[2])),
       title="PCA: samples colored by Sex")
plt_sep_ref <- grid.arrange(p2_batch, p2_cond, p2_studies, p2_sex, ncol=2, nrow=2)
dev.off()

### Final dataset separate by Sex. Not really a batch effect. If you want to fix it, treat sex as a "batch" indicator.

grid.arrange(plt_combat, plt_sep_ref, nrow=2)

pca_obj_test <- plotPCA(DESeqTransform(seobj_merged), intgroup=c("SequencingBatch", "Label", "Sex", "Dataset")) 

