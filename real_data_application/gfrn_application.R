rm(list=ls())
sapply(c("sva", "dplyr", "DESeq2", "ggplot2", "reshape2", "gridExtra", "scales", "dendextend", 
         "RUVSeq", "ggpubr", "BatchQC"), require, character.only=TRUE)
demo <- TRUE   # if testing code, set as TRUE; if running simulations, set as FALSE
if(demo){
  setwd("~/Documents/ComBat_seq/real_data_app/")
  script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
  data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/gfrn_signature"
  source("~/Dropbox/Work/ComBat_Seq/ComBat-Seq/real_data_app/gfrn_sig_helpers.R")
}else{
  setwd("~/yuqingz/ComBat_seq/real_data_app/")
  script_dir <- ".."; data_dir <- "."
  source("gfrn_sig_helpers.R")
}
set.seed(1)

## Load ComBat-seq functions
source(file.path(script_dir, "ComBat_seq.R"))
source(file.path(script_dir, "helper_seq.R"))

## Parameters
pathway_regex <- c("her2", "^egfr", "kraswt")  #"^egfr"  

## Load data
sigdata <- readRDS(file.path(data_dir, "signature_data.rds"))
cts_mat <- assay(sigdata, "counts")   # count matrix (also have tpm and fpkm in there)
# filter out genes with 0 counts
#cts_mat <- cts_mat[apply(cts_mat, 1, function(x){all(x!=0)}), ]
rownames(cts_mat) <- paste0("gene", 1:nrow(cts_mat))
batch <- colData(sigdata)$batch
group <- colData(sigdata)$group


## Take the subset (all controls & 1 condition specified by pathway_regex)
pathway_condition_ind <- grep(paste(pathway_regex, collapse="|"), group)
#pathway_batch <- unique(batch[pathway_condition_ind])  # which batch does this pathway belong to
ctrl_ind <- which(group %in% c("gfp_for_egfr", "gfp18", "gfp30"))
subset_ind <- sort(c(ctrl_ind, pathway_condition_ind))  
cts_sub <- cts_mat[, subset_ind]
batch_sub <- batch[subset_ind]
group_sub <- group[subset_ind]
# remove genes with only 0 counts in the subset & in any batch
#cts_sub <- cts_sub[apply(cts_sub,1,function(x){!all(x==0)}), ]
keep1 <- apply(cts_sub[, batch_sub==1],1,function(x){!all(x==0)})
keep2 <- apply(cts_sub[, batch_sub==2],1,function(x){!all(x==0)})
keep3 <- apply(cts_sub[, batch_sub==3],1,function(x){!all(x==0)})
cts_sub <- cts_sub[keep1 & keep2 & keep3, ]


## Use ComBatSeq to adjust data
# group_sub <- as.character(group_sub)
# group_sub_num <- rep(1, length(group_sub))
# group_sub_num[grep("gfp", group_sub)] <- 0
group_sub <- factor(as.character(group_sub), levels=c("gfp_for_egfr", "gfp18", "gfp30",  gsub("^", "", pathway_regex, fixed=T)))
group_sub <- plyr::revalue(group_sub, c("gfp_for_egfr"="gfp", "gfp18"="gfp", "gfp30"="gfp"))
#counts=cts_sub;batch=batch_sub;group=group_sub_num;gene.subset.n=1000;Cpp=FALSE;covar_mod=NULL;full_mod=TRUE
start_time <- Sys.time()
combatseq_sub <- ComBat_seq(counts=cts_sub, batch=batch_sub, group=group_sub, gene.subset.n=1000,
                            shrink=FALSE, shrink.disp=FALSE)
# combatseq_sub <- ComBat_seq(counts=cts_sub, batch=batch_sub, group=group_sub, 
#                             ctrl=list(full_mod=TRUE, shrink=TRUE, gene.subset.n=100))$adjust_counts
end_time <- Sys.time()
print(end_time - start_time)
#combatseq_sub <- readRDS("~/Google Drive/ComBat_seq/real_data_example/RNAseq/gfrn_signature/combatseq_test.rds")


## Use original ComBat
combat_sub <- ComBat(cpm(cts_sub, log=TRUE), batch=batch_sub, mod=model.matrix(~group_sub))


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


## SVAseq
mod1 <- model.matrix(~group_sub); mod0 <- cbind(mod1[,1])
svseq <- svaseq(cts_sub, mod=mod1, mod0=mod0, n.sv=3); cat("\n")



## Normalize library size
cts_norm <- apply(cts_sub, 2, function(x){x/sum(x)})
cts_adj_norm <- apply(combatseq_sub, 2, function(x){x/sum(x)})
cts_adjori_norm <- apply(combat_sub, 2, function(x){x/sum(x)})
cts_ruvseq_norm <- apply(ruvseq_sub, 2, function(x){x/sum(x)})
cts_svaseq_norm <- rbind(t(svseq$sv), cts_norm)


## PCA 
col_data <- data.frame(Batch=factor(batch_sub), Group=group_sub) 
rownames(col_data) <- colnames(cts_sub)

seobj <- SummarizedExperiment(assays=cts_norm, colData=col_data)
pca_obj <- plotPCA(DESeqTransform(seobj), intgroup=c("Batch", "Group"))
plt <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
       title="Unadjusted") # +
  #theme(legend.position="bottom")

seobj_adj <- SummarizedExperiment(assays=cts_adj_norm, colData=col_data)
pca_obj_adj <- plotPCA(DESeqTransform(seobj_adj), intgroup=c("Batch", "Group"))
plt_adj <- ggplot(pca_obj_adj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj$plot_env$percentVar[2])),
       title="ComBat-Seq") # +
  #theme(legend.position="bottom")

seobj_adjori <- SummarizedExperiment(assays=cts_adjori_norm, colData=col_data)
pca_obj_adjori <- plotPCA(DESeqTransform(seobj_adjori), intgroup=c("Batch", "Group"))
plt_adjori <- ggplot(pca_obj_adjori$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[2])),
       title="Original ComBat") # +
  #theme(legend.position="bottom")

## ruvseq
seobj_ruvseq <- SummarizedExperiment(assays=cts_ruvseq_norm, colData=col_data)
pca_obj_ruvseq <- plotPCA(DESeqTransform(seobj_ruvseq), intgroup=c("Batch", "Group"))
plt_ruvseq <- ggplot(pca_obj_ruvseq$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) +
  geom_point() +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_ruvseq$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_ruvseq$plot_env$percentVar[2])),
       title="RUV-Seq")

## svaseq
seobj_svaseq <- SummarizedExperiment(assays=cts_svaseq_norm, colData=col_data)
pca_obj_svaseq <- plotPCA(DESeqTransform(seobj_svaseq), intgroup=c("Batch", "Group"))
plt_svaseq <- ggplot(pca_obj_svaseq$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) +
  geom_point() +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_svaseq$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_svaseq$plot_env$percentVar[2])),
       title="SVA-Seq")

png("gfrn_PCA_shrinkOff_multipath.png", width=8, height=7, units="in", res=300)
#grid.arrange(plt, plt_adj, plt_adjori, plt_ruvseq, ncol=2, nrow=2)
ggarrange(plt, plt_adj, plt_adjori, plt_ruvseq, ncol=2, nrow=2, common.legend=TRUE, legend="right")
dev.off()



## Variation explained
# varexp_lst <- list(
#   unadjusted=batchqc_explained_variation(cts_norm, condition=col_data$Group, batch=col_data$Batch)$explained_variation,
#   combatseq=batchqc_explained_variation(cts_adj_norm, condition=col_data$Group, batch=col_data$Batch)$explained_variation,
#   combat=batchqc_explained_variation(cts_adjori_norm, condition=col_data$Group, batch=col_data$Batch)$explained_variation,
#   ruvseq=batchqc_explained_variation(cts_ruvseq_norm, condition=col_data$Group, batch=col_data$Batch)$explained_variation
# )
varexp_lst <- list(
  unadjusted=batchqc_explained_variation(cpm(cts_sub, log=TRUE), condition=col_data$Group, batch=col_data$Batch)$explained_variation,
  combatseq=batchqc_explained_variation(cpm(combatseq_sub, log=TRUE), condition=col_data$Group, batch=col_data$Batch)$explained_variation,
  combat=batchqc_explained_variation(combat_sub, condition=col_data$Group, batch=col_data$Batch)$explained_variation,
  ruvseq=batchqc_explained_variation(cpm(ruvseq_sub, log=TRUE), condition=col_data$Group, batch=col_data$Batch)$explained_variation
)
varexp_df <- melt(varexp_lst)
varexp_df$L1 <- factor(varexp_df$L1, levels=c("unadjusted", "combatseq", "combat", "ruvseq"))
varexp_df$L1 <- plyr::revalue(varexp_df$L1, c("unadjusted"="Unadjusted", "combatseq"="ComBat-Seq",
                                              "combat"="Original ComBat", "ruvseq"="RUV-Seq"))
png("gfrn_varexp_shrinkOff_multipath.png", width=7, height=7, units="in", res=300)
ggplot(varexp_df, aes(x=Var2, y=value)) +
  geom_boxplot() +
  facet_wrap(~L1, nrow=2, ncol=2) +
  labs(y="Explained variation") +
  theme(axis.title.x = element_blank())
dev.off()
