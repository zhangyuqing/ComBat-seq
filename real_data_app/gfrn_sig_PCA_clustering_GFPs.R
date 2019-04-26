rm(list=ls())
sapply(c("dplyr", "edgeR", "DESeq2", "sva", "RUVSeq", "MASS", "ggplot2", "reshape2", "gridExtra", "scales", "dendextend"), 
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


## Use ComBatSeq to adjust data
#counts=cts_ctrl;batch=batch_ctrl;group=NULL;full_mod=FALSE;covar_mod=NULL
start_time <- Sys.time()
combatseq_ctrl <- ComBat_seq(counts=cts_ctrl, batch=batch_ctrl, group=NULL, full_mod=FALSE, gene.subset.n=1000, Cpp=FALSE, 
                             shrink=FALSE, shrink.disp=FALSE)
end_time <- Sys.time()
print(end_time - start_time)
#combatseq_ctrl <- readRDS("~/Google Drive/ComBat_seq/real_data_example/RNAseq/gfrn_signature/combatseq_ctrlonly.rds")


## Use original ComBat
combat_ctrl <- ComBat(cts_ctrl, batch=batch_ctrl, mod=NULL)


## Normalize library size
cts_ctrl_norm <- apply(cts_ctrl, 2, function(x){x/sum(x)})
cts_ctrl_adj_norm <- apply(combatseq_ctrl, 2, function(x){x/sum(x)})
cts_ctrl_adjori_norm <- apply(combat_ctrl, 2, function(x){x/sum(x)})


## PCA of log data
col_data <- data.frame(Batch=factor(batch_ctrl)); rownames(col_data) <- colnames(cts_ctrl)

seobj <- SummarizedExperiment(assays=cts_ctrl_norm, colData=col_data)
pca_obj <- plotPCA(DESeqTransform(seobj), intgroup="Batch")
plt <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=Batch)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
       title="Color by Batch - Unadjusted") +
  theme(legend.position="bottom")
seobj_adj <- SummarizedExperiment(assays=cts_ctrl_adj_norm, colData=col_data)
pca_obj_adj <- plotPCA(DESeqTransform(seobj_adj), intgroup="Batch")
plt_adj <- ggplot(pca_obj_adj$data, aes(x=PC1, y=PC2, color=Batch)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj$plot_env$percentVar[2])),
       title="Color by Batch - ComBat-Seq") +
  theme(legend.position="bottom")
seobj_adjori <- SummarizedExperiment(assays=cts_ctrl_adjori_norm, colData=col_data)
pca_obj_adjori <- plotPCA(DESeqTransform(seobj_adjori), intgroup="Batch")
plt_adjori <- ggplot(pca_obj_adjori$data, aes(x=PC1, y=PC2, color=Batch)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[2])),
       title="Color by Batch - Original ComBat") +
  theme(legend.position="bottom")
png("gfrn_PCA_GFPs_shrinkOff.png", width=14, height=5, units="in", res=300)
grid.arrange(plt, plt_adj, plt_adjori, ncol=3)
dev.off()


## Clustering
hc <- hclust(dist(t(cts_ctrl_norm)))
dend <- as.dendrogram(hc)
dend <- color_branches(dend, col=batch_ctrl[order.dendrogram(dend)]+1)
plot(dend)  

adj_hc <- hclust(dist(t(cts_ctrl_adj_norm)))
adj_dend <- as.dendrogram(adj_hc)
adj_dend <- color_branches(adj_dend, col=batch_ctrl[order.dendrogram(adj_dend)]+1)
plot(adj_dend)  

