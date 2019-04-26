rm(list=ls())
sapply(c("dplyr", "plyr", "edgeR", "DESeq2", "ggplot2", "reshape2", "gridExtra", "scales", "dendextend"), 
       require, character.only=TRUE)
demo <- TRUE   # if testing code, set as TRUE; if running simulations, set as FALSE
if(demo){
  setwd("~/Documents/ComBat_seq/real_data_app/")
  script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
  data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/WRBU_lungcancer"
  source("~/Dropbox/Work/ComBat_Seq/ComBat-Seq/real_data_app/gfrn_sig_helpers.R")
}else{
  setwd("~/yuqingz/ComBat_seq/real_data_app/")
  script_dir <- ".."; data_dir <- "."
  source("gfrn_sig_helpers.R")
}
#set.seed(1)
source(file.path(script_dir, "ComBat_seq.R"))
source(file.path(script_dir, "helper_seq.R"))


## Load data
cts <- readRDS(file.path(data_dir, "count_matrix_BU_WR.rds"))
meta_info <- readRDS(file.path(data_dir, "BU_WR_annotation_file.rds"))
# drop 4790-024 (single sample as its own batch)
cts <- cts[, -grep("4790-024", colnames(cts))]
meta_info <- filter(meta_info, kitnumber!="4790-024")
batch <- factor(meta_info$site)
group <- factor(meta_info$smoking_status)
group <- as.factor(as.character(revalue(group, c("Former smoker"="0", "Current smoker"="1"))))
# remove all 0 genes
tmp <- t(apply(cts,1,floor)); colnames(tmp) <- colnames(cts); cts <- tmp; rm(tmp)
cts <- cts[rowVars(cts)>0, ]


## Run combat-seq
start_time <- Sys.time()
combatseq_cts <- ComBat_seq(counts=cts, batch=batch, group=group, gene.subset.n=1000, Cpp=FALSE,
                            shrink=FALSE, shrink.disp=FALSE)
end_time <- Sys.time()
print(end_time - start_time)


## Run original combat
combat_cts <- ComBat(cts, batch=batch, mod=model.matrix(~group))


## Normalizing library size
cts_norm <- apply(cts, 2, function(x){x/sum(x)})
cts_adj_norm <- apply(combatseq_cts, 2, function(x){x/sum(x)})
cts_adjori_norm <- apply(combat_cts, 2, function(x){x/sum(x)})


## PCA 
col_data <- meta_info; rownames(col_data) <- meta_info$sample

seobj <- SummarizedExperiment(assays=cts_norm, colData=col_data)
pca_obj <- plotPCA(DESeqTransform(seobj), intgroup=c("site", "smoking_status"))
plt_batch <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=site)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
       title="Color by batch (location) - Unadjusted") +
  theme(legend.position="bottom")
plt_cond <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=smoking_status)) + 
  geom_point() + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
       title="Color by condition - Unadjusted") +
  theme(legend.position="bottom")
rm(seobj, pca_obj)

#seobj_adj <- SummarizedExperiment(assays=list(counts=combatseq_cts), colData=col_data)
seobj_adj <- SummarizedExperiment(assays=cts_adj_norm, colData=col_data)
pca_obj_adj <- plotPCA(DESeqTransform(seobj_adj), intgroup=c("site", "smoking_status"))
plt_batch_adj <- ggplot(pca_obj_adj$data, aes(x=PC1, y=PC2, color=site)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj$plot_env$percentVar[2])),
       title="Color by batch (location) - ComBat-seq") +
  theme(legend.position="bottom")
plt_cond_adj <- ggplot(pca_obj_adj$data, aes(x=PC1, y=PC2, color=smoking_status)) + 
  geom_point() + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj$plot_env$percentVar[2])),
       title="Color by condition - ComBat-seq") +
  theme(legend.position="bottom")

seobj_adjori <- SummarizedExperiment(assays=cts_adjori_norm, colData=col_data)
pca_obj_adjori <- plotPCA(DESeqTransform(seobj_adjori), intgroup=c("site", "smoking_status"))
plt_batch_adjori <- ggplot(pca_obj_adjori$data, aes(x=PC1, y=PC2, color=site)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[2])),
       title="Color by batch (location) - Original ComBat") +
  theme(legend.position="bottom")
plt_cond_adjori <- ggplot(pca_obj_adjori$data, aes(x=PC1, y=PC2, color=smoking_status)) + 
  geom_point() + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[2])),
       title="Color by condition - Original ComBat") +
  theme(legend.position="bottom")

png("wrbu_smoke_PCA_combatseq.png", width=8, height=8, units="in", res=300)
grid.arrange(plt_batch, plt_cond, plt_batch_adj, plt_cond_adj, ncol=2, nrow=2)
dev.off()
png("wrbu_smoke_PCA_combatseq_full.png", width=8, height=13, units="in", res=300)
grid.arrange(plt_batch, plt_cond, plt_batch_adj, plt_cond_adj, plt_batch_adjori, plt_cond_adjori, ncol=2, nrow=3)
dev.off()



## Clustering
hc <- hclust(dist(t(cts_norm)))
dend <- as.dendrogram(hc)
dend <- color_branches(dend, groupLabels=batch[order.dendrogram(dend)], 
                       col=as.numeric(as.character(group))[order.dendrogram(dend)] + 2)
plot(dend)  

adj_hc <- hclust(dist(t(combatseq_cts_norm)))
adj_dend <- as.dendrogram(adj_hc)
adj_dend <- color_branches(adj_dend, groupLabels=batch[order.dendrogram(adj_dend)], 
                           col=as.numeric(as.character(group))[order.dendrogram(adj_dend)]+2)
plot(adj_dend)
