rm(list=ls())
sapply(c("dplyr", "edgeR", "DESeq2", "sva", "RUVSeq", "MASS", "ggplot2", "reshape2", "gridExtra", "scales", "dendextend", "ggpubr"), 
       require, character.only=TRUE)
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
ComBat_seq_old <- ComBat_seq; rm(ComBat_seq)
source(file.path(script_dir, "ComBat_seq_newEB_sep.R")) 
source(file.path(script_dir, "helper_seq.R"))


## Load data
sigdata <- readRDS(file.path(data_dir, "signature_data.rds"))
cts_mat <- assay(sigdata, "counts")   # count matrix (also have tpm and fpkm in there)
# filter out genes with 0 counts
#cts_mat <- cts_mat[apply(cts_mat, 1, function(x){all(x!=0)}), ]
rownames(cts_mat) <- paste0("gene", 1:nrow(cts_mat))
batch <- colData(sigdata)$batch
group <- colData(sigdata)$group


## Take the subset (all controls)
ctrl_ind <- which(group %in% c("gfp_for_egfr", "gfp18", "gfp30"))
cts_ctrl <- cts_mat[, ctrl_ind]
batch_ctrl <- batch[ctrl_ind]
#remove genes with only 0 counts in all control samples
cts_ctrl <- cts_ctrl[apply(cts_ctrl,1,function(x){!all(x==0)}), ]


## Use old ComBatSeq
combatseq_ctrl_old <- ComBat_seq_old(counts=cts_ctrl, batch=batch_ctrl, group=NULL, full_mod=FALSE, gene.subset.n=1000, 
                                     Cpp=FALSE, shrink=FALSE, shrink.disp=FALSE)

## Use ComBatSeq to adjust data
#counts=cts_ctrl;batch=batch_ctrl;group=NULL;full_mod=FALSE;covar_mod=NULL
#ctrl_obj <- ComBat_seq_control_gene(cts_ctrl, batch=batch_ctrl, group=NULL)
start_time <- Sys.time()
combatseq_ctrl <- ComBat_seq(counts=cts_ctrl, batch=batch_ctrl, group=NULL, 
                             ctrl=list(full_mod=FALSE, shrink=TRUE, gene.subset.n=100))
end_time <- Sys.time()
print(end_time - start_time)
#combatseq_ctrl <- readRDS("~/Google Drive/ComBat_seq/real_data_example/RNAseq/gfrn_signature/combatseq_ctrlonly.rds")


## Use original ComBat
combat_ctrl <- ComBat(cts_ctrl, batch=batch_ctrl, mod=NULL)


## Normalize library size
cts_ctrl_norm <- apply(cts_ctrl, 2, function(x){x/sum(x)})
cts_ctrl_adj_norm <- apply(combatseq_ctrl$adjust_counts, 2, function(x){x/sum(x)})
cts_ctrl_adj_norm_old <- apply(combatseq_ctrl_old, 2, function(x){x/sum(x)})
cts_ctrl_adjori_norm <- apply(combat_ctrl, 2, function(x){x/sum(x)})


## PCA 
col_data <- data.frame(Batch=factor(batch_ctrl)); rownames(col_data) <- colnames(cts_ctrl)

seobj <- SummarizedExperiment(assays=cts_ctrl_norm, colData=col_data)
pca_obj <- plotPCA(DESeqTransform(seobj), intgroup="Batch")
plt <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=Batch)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
       title="Color by Batch - Unadjusted") +
  theme(legend.position="bottom")

seobj_adj_newEB <- SummarizedExperiment(assays=cts_ctrl_adj_norm, colData=col_data)
pca_obj_adj_newEB <- plotPCA(DESeqTransform(seobj_adj_newEB), intgroup="Batch")
plt_adj_newEB <- ggplot(pca_obj_adj_newEB$data, aes(x=PC1, y=PC2, color=Batch)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj_newEB$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj_newEB$plot_env$percentVar[2])),
       title="Color by Batch - ComBat-Seq (new EB)") +
  theme(legend.position="bottom")

seobj_adj <- SummarizedExperiment(assays=cts_ctrl_adj_norm_old, colData=col_data)
pca_obj_adj <- plotPCA(DESeqTransform(seobj_adj), intgroup="Batch")
plt_adj <- ggplot(pca_obj_adj$data, aes(x=PC1, y=PC2, color=Batch)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj$plot_env$percentVar[2])),
       title="Color by Batch - ComBat-Seq (shrink off)") +
  theme(legend.position="bottom")

seobj_adjori <- SummarizedExperiment(assays=cts_ctrl_adjori_norm, colData=col_data)
pca_obj_adjori <- plotPCA(DESeqTransform(seobj_adjori), intgroup="Batch")
plt_adjori <- ggplot(pca_obj_adjori$data, aes(x=PC1, y=PC2, color=Batch)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[2])),
       title="Color by Batch - Original ComBat") +
  theme(legend.position="bottom")

png("gfrn_newEB_PCA_GFPs.png", width=10, height=10, units="in", res=300)
grid.arrange(plt, plt_adjori, plt_adj, plt_adj_newEB, ncol=2, nrow=2)
dev.off()


## DE
source("~/Dropbox/Work/ComBat_Seq/ComBat-Seq/real_data_app/exprep_helpers.R")
DE_method <- "DESeq2"
de_measure_name <- "pvalue"
if(DE_method=="edgeR"){DE_function <- edgeR_DEpipe}else{DE_function <- DESeq2_DEpipe}
fpr_lst <- list()
set.seed(12345)

for(itr in 1:100){
  print(itr)
  rand_group <- factor(sample(c(0,1), ncol(cts_ctrl), replace=TRUE))
  de_unadj <- DE_function(cts_ctrl, batch=batch_ctrl, group=rand_group, include.batch=FALSE, 
                          alpha.unadj=1, alpha.fdr=1, covar_incl=NULL, covar=NULL)
  de_onestep <-  DE_function(cts_ctrl, batch=batch_ctrl, group=rand_group, include.batch=TRUE, 
                             alpha.unadj=1, alpha.fdr=1, covar_incl=NULL, covar=NULL)
  de_combatseq_newEB <- DE_function(combatseq_ctrl$adjust_counts, batch=batch_ctrl, group=rand_group, include.batch=FALSE, 
                                    alpha.unadj=1, alpha.fdr=1, covar_incl=NULL, covar=NULL)
  de_combatseq_old <- DE_function(combatseq_ctrl_old, batch=batch_ctrl, group=rand_group, include.batch=FALSE, 
                                  alpha.unadj=1, alpha.fdr=1, covar_incl=NULL, covar=NULL)
  
  # de_lst <- list(unadjusted=de_unadj$de_res[, de_measure_name], #length(which(de_unadj$de_res$PValue < 0.05)), 
  #                one.step=de_onestep$de_res[, de_measure_name],
  #                combatseq.old=de_combatseq_old$de_res[, de_measure_name],
  #                combatseq.newEB=de_combatseq_newEB$de_res[, de_measure_name])
  # boxplot(de_lst)
  
  fpr_lst[[itr]] <- c(unadjusted=length(which(de_unadj$de_res[, de_measure_name] < 0.05))/nrow(cts_ctrl), #length(which(de_unadj$de_res$PValue < 0.05)), 
                      one.step=length(which(de_onestep$de_res[, de_measure_name] < 0.05))/nrow(cts_ctrl),
                      combatseq.old=length(which(de_combatseq_old$de_res[, de_measure_name] < 0.05))/nrow(cts_ctrl),
                      combatseq.newEB=length(which(de_combatseq_newEB$de_res[, de_measure_name] < 0.05))/nrow(cts_ctrl))
}

fpr_res <- do.call(rbind, fpr_lst)
png(sprintf("gfrn_newEB_fpr_GFPs_%s.png", DE_method), width=5, height=5, units="in", res=300)
ggplot(melt(data.frame(fpr_res), variable.name="Method"), aes(x=Method, y=value)) +
  geom_boxplot() +
  geom_hline(yintercept=0.05, color="red") +
  ylim(0, 0.1) +
  labs(y="FPR", title=sprintf("GFRN GFP controls - %s", DE_method)) +
  theme(axis.title.x=element_blank())
dev.off()

