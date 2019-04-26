rm(list=ls())
sapply(c("sva", "dplyr", "DESeq2", "ggplot2", "reshape2", "gridExtra", "scales", "dendextend"), require, character.only=TRUE)
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
pathway_regex <- "akt"  #"^egfr"  

## Load data
sigdata <- readRDS(file.path(data_dir, "signature_data.rds"))
cts_mat <- assay(sigdata, "counts")   # count matrix (also have tpm and fpkm in there)
# filter out genes with 0 counts
#cts_mat <- cts_mat[apply(cts_mat, 1, function(x){all(x!=0)}), ]
rownames(cts_mat) <- paste0("gene", 1:nrow(cts_mat))
batch <- colData(sigdata)$batch
group <- colData(sigdata)$group


## Take the subset (all controls & 1 condition specified by pathway_regex)
pathway_condition_ind <- grep(pathway_regex, group)
pathway_batch <- unique(batch[pathway_condition_ind])  # which batch does this pathway belong to
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
group_sub <- as.character(group_sub)
group_sub_num <- rep(1, length(group_sub))
group_sub_num[grep("gfp", group_sub)] <- 0
#counts=cts_sub;batch=batch_sub;group=group_sub_num;gene.subset.n=1000;Cpp=FALSE;covar_mod=NULL;full_mod=TRUE
start_time <- Sys.time()
combatseq_sub <- ComBat_seq(counts=cts_sub, batch=batch_sub, group=group_sub_num, gene.subset.n=1000, Cpp=FALSE, 
                            shrink=FALSE, shrink.disp=FALSE)
end_time <- Sys.time()
print(end_time - start_time)
#combatseq_sub <- readRDS("~/Google Drive/ComBat_seq/real_data_example/RNAseq/gfrn_signature/combatseq_test.rds")


## Use original ComBat
combat_sub <- ComBat(cts_sub, batch=batch_sub, mod=model.matrix(~factor(group_sub_num)))


## Normalize library size
cts_norm <- apply(cts_sub, 2, function(x){x/sum(x)})
cts_adj_norm <- apply(combatseq_sub, 2, function(x){x/sum(x)})
cts_adjori_norm <- apply(combat_sub, 2, function(x){x/sum(x)})


## PCA of log data
col_data <- data.frame(Batch=factor(batch_sub), Group=factor(group_sub_num)) 
rownames(col_data) <- colnames(cts_sub)

seobj <- SummarizedExperiment(assays=cts_norm, colData=col_data)
pca_obj <- plotPCA(DESeqTransform(seobj), intgroup=c("Batch", "Group"))
plt <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
       title="Color by Batch - Unadjusted") +
  theme(legend.position="bottom")

seobj_adj <- SummarizedExperiment(assays=cts_adj_norm, colData=col_data)
pca_obj_adj <- plotPCA(DESeqTransform(seobj_adj), intgroup=c("Batch", "Group"))
plt_adj <- ggplot(pca_obj_adj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj$plot_env$percentVar[2])),
       title="Color by Batch - ComBat-Seq") +
  theme(legend.position="bottom")

seobj_adjori <- SummarizedExperiment(assays=cts_adjori_norm, colData=col_data)
pca_obj_adjori <- plotPCA(DESeqTransform(seobj_adjori), intgroup=c("Batch", "Group"))
plt_adjori <- ggplot(pca_obj_adjori$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[2])),
       title="Color by Batch - Original ComBat") +
  theme(legend.position="bottom")

png(sprintf("gfrn_PCA_shrinkOff_%s.png", gsub("^", "", pathway_regex)), 
    width=14, height=5, units="in", res=300)
grid.arrange(plt, plt_adj, plt_adjori, ncol=3)
dev.off()


## Clustering
hc <- hclust(dist(t(cts_norm)))
dend <- as.dendrogram(hc)
dend <- color_branches(dend, col=batch_sub[order.dendrogram(dend)]+1)#, groupLabels=group_sub_num)
plot(dend)#, horiz=TRUE)  

adj_hc <- hclust(dist(t(cts_adj_norm)))
adj_dend <- as.dendrogram(adj_hc)
adj_dend <- color_branches(adj_dend, col=batch_sub[order.dendrogram(adj_dend)]+1)#, groupLabels=group_sub_num)
plot(adj_dend)#, horiz=TRUE)  

