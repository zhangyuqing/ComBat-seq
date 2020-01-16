rm(list=ls())
sapply(c("sva", "dplyr", "DESeq2", "ggplot2", "reshape2", "gridExtra", "scales", 
         "RUVSeq", "ggpubr", "BatchQC"), require, character.only=TRUE)

## Parameters (change paths when necessary)
data_dir <- "~/Desktop/ComBat-seq/real_data_application"  # path to the signature data (.rds)
source("~/Desktop/ComBat-seq/real_data_application/gfrn_helpers.R")  # path to gfrn_helpers.R
source("~/Desktop/ComBat-seq/ComBat_seq.R"); source("~/Desktop/ComBat-seq/helper_seq.R")   
# path to the combat-seq scripts (or use the sva package on github, in which case comment out the above line)

pathway_regex <- c("her2", "^egfr", "kraswt")  
set.seed(1)


## Load data
sigdata <- readRDS(file.path(data_dir, "signature_data.rds"))
cts_mat <- assay(sigdata, "counts")   # count matrix (also have tpm and fpkm in there)
rownames(cts_mat) <- paste0("gene", 1:nrow(cts_mat))
batch <- colData(sigdata)$batch
group <- colData(sigdata)$group

# Take the subset (all controls & 1 condition /batch specified by pathway_regex)
pathway_condition_ind <- grep(paste(pathway_regex, collapse="|"), group)
ctrl_ind <- which(group %in% c("gfp_for_egfr", "gfp18", "gfp30"))
subset_ind <- sort(c(ctrl_ind, pathway_condition_ind))  
cts_sub <- cts_mat[, subset_ind]
batch_sub <- batch[subset_ind]
group_sub <- group[subset_ind]
# remove genes with only 0 counts in the subset & in any batch
keep1 <- apply(cts_sub[, batch_sub==1],1,function(x){!all(x==0)})
keep2 <- apply(cts_sub[, batch_sub==2],1,function(x){!all(x==0)})
keep3 <- apply(cts_sub[, batch_sub==3],1,function(x){!all(x==0)})
cts_sub <- cts_sub[keep1 & keep2 & keep3, ]


## Use ComBatSeq to adjust data
group_sub <- factor(as.character(group_sub), levels=c("gfp_for_egfr", "gfp18", "gfp30",  gsub("^", "", pathway_regex, fixed=T)))
group_sub <- plyr::revalue(group_sub, c("gfp_for_egfr"="gfp", "gfp18"="gfp", "gfp30"="gfp"))

start_time <- Sys.time()
combatseq_sub <- ComBat_seq(counts=cts_sub, batch=batch_sub, group=group_sub, shrink=FALSE)
end_time <- Sys.time()
print(end_time - start_time)


## Use original ComBat on logCPM
combat_sub <- sva::ComBat(cpm(cts_sub, log=TRUE), batch=batch_sub, mod=model.matrix(~group_sub))


## RUVseq
group1 <- plyr::revalue(as.factor(as.character(group_sub[batch_sub==1])), c("gfp"="0", "her2"="1"))
deres1 <- edgeR_DEpipe(cts_sub[, batch_sub==1], batch=NULL, group=group1,
                       include.batch=FALSE, alpha.unadj=1, alpha.fdr=1)
group2 <- plyr::revalue(as.factor(as.character(group_sub[batch_sub==2])), c("gfp"="0", "egfr"="1"))
deres2 <- edgeR_DEpipe(cts_sub[, batch_sub==2], batch=NULL, group=group2,
                       include.batch=FALSE, alpha.unadj=1, alpha.fdr=1)
group3 <- plyr::revalue(as.factor(as.character(group_sub[batch_sub==3])), c("gfp"="0", "kraswt"="1"))
deres3 <- edgeR_DEpipe(cts_sub[, batch_sub==3], batch=NULL, group=group3,
                       include.batch=FALSE, alpha.unadj=1, alpha.fdr=1)
null_genes <- Reduce(intersect, list(which(deres1$de_res$FDR>0.95), which(deres2$de_res$FDR>0.95), which(deres3$de_res$FDR>0.95)))
group_obj <- makeGroups(group_sub)
ruvseq_sub <- RUVs(cts_sub, cIdx=null_genes, scIdx=group_obj, k=1)$normalizedCounts


## Normalize library size
cts_norm <- apply(cts_sub, 2, function(x){x/sum(x)})
cts_adj_norm <- apply(combatseq_sub, 2, function(x){x/sum(x)})
cts_adjori_norm <- apply(combat_sub, 2, function(x){x/sum(x)})
cts_ruvseq_norm <- apply(ruvseq_sub, 2, function(x){x/sum(x)})


## PCA 
col_data <- data.frame(Batch=factor(batch_sub), Group=group_sub) 
rownames(col_data) <- colnames(cts_sub)

seobj <- SummarizedExperiment(assays=cts_norm, colData=col_data)
pca_obj <- plotPCA(DESeqTransform(seobj), intgroup=c("Batch", "Group"))
plt <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
       title="Unadjusted") 

seobj_adj <- SummarizedExperiment(assays=cts_adj_norm, colData=col_data)
pca_obj_adj <- plotPCA(DESeqTransform(seobj_adj), intgroup=c("Batch", "Group"))
plt_adj <- ggplot(pca_obj_adj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj$plot_env$percentVar[2])),
       title="ComBat-Seq") 

seobj_adjori <- SummarizedExperiment(assays=cts_adjori_norm, colData=col_data)
pca_obj_adjori <- plotPCA(DESeqTransform(seobj_adjori), intgroup=c("Batch", "Group"))
plt_adjori <- ggplot(pca_obj_adjori$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[2])),
       title="Original ComBat") 

## ruvseq
seobj_ruvseq <- SummarizedExperiment(assays=cts_ruvseq_norm, colData=col_data)
pca_obj_ruvseq <- plotPCA(DESeqTransform(seobj_ruvseq), intgroup=c("Batch", "Group"))
plt_ruvseq <- ggplot(pca_obj_ruvseq$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) +
  geom_point() +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_ruvseq$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_ruvseq$plot_env$percentVar[2])),
       title="RUV-Seq")

plt_PCA_full <- ggarrange(plt, plt_ruvseq, plt_adjori, plt_adj, ncol=1, nrow=4, common.legend=TRUE, legend="right")

varexp_full <- list(
  unadjusted=batchqc_explained_variation(cpm(cts_sub, log=TRUE), condition=col_data$Group, batch=col_data$Batch)$explained_variation,
  combatseq=batchqc_explained_variation(cpm(combatseq_sub, log=TRUE), condition=col_data$Group, batch=col_data$Batch)$explained_variation,
  combat=batchqc_explained_variation(combat_sub, condition=col_data$Group, batch=col_data$Batch)$explained_variation,
  ruvseq=batchqc_explained_variation(cpm(ruvseq_sub, log=TRUE), condition=col_data$Group, batch=col_data$Batch)$explained_variation
)
varexp_full_df <- melt(varexp_full)
varexp_full_df$L1 <- factor(varexp_full_df$L1, levels=c("unadjusted", "ruvseq", "combat","combatseq"))
varexp_full_df$L1 <- plyr::revalue(varexp_full_df$L1, c("unadjusted"="Unadjusted", "combatseq"="ComBat-Seq",
                                                        "ruvseq"="RUV-Seq", "combat"="Original ComBat"))
varexp_full_df$Var2 <- plyr::revalue(varexp_full_df$Var2, c("Full (Condition+Batch)"="Condition+Batch"))

plt_varexp_full <- ggplot(varexp_full_df, aes(x=Var2, y=value)) +
  geom_boxplot() +
  facet_wrap(~L1, nrow=4, ncol=1) +
  labs(y="Explained variation") +
  theme(axis.title.x = element_blank())

ggarrange(plt_PCA_full, plt_varexp_full, ncol=2, widths=c(0.55, 0.45))

