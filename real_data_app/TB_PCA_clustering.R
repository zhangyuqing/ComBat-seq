rm(list=ls())
sapply(c("edgeR", "sva", "ggplot2", "reshape2", "gridExtra", "SummarizedExperiment", "DESeq2", "scales", "BatchQC", "dendextend"), 
       require, character.only=TRUE)
demo <- TRUE  # if testing code, set as TRUE; if running simulations, set as FALSE
if(demo){
  setwd("~/Documents/ComBat_seq/real_data_app/")
  script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
  source(file.path(script_dir, "real_data_app/exprep_helpers.R"))
}else{
  setwd("~/yuqingz/ComBat_seq/real_data_app/")
  #setwd("/restricted/projectnb/combat/work/yuqingz/ComBat-seq/real_data_app")
  script_dir <- ".."
  source(file.path(script_dir, "exprep_helpers.R"))
}
source(file.path(script_dir, "ComBat_seq.R"))
source(file.path(script_dir, "helper_seq.R"))
set.seed(123)


## Load data
load("TB_ExpRep.RData")
cts <- cts[setdiff(rownames(cts),"HBA2"), ]  ##!!!!!! need to deal with this...
covar_mod <- model.matrix(~covar)


## Adjust data with ComBat-seq
start_time <- Sys.time()
cts_adj <- ComBat_seq(counts=cts, batch=batch, group=group, gene.subset.n=1000, Cpp=FALSE, covar_mod=covar_mod)
end_time <- Sys.time()
print(end_time - start_time)
#cts_adj <- readRDS("TB_combatseq.rds")


## Adjust data with Original ComBat
cts_adjori <- ComBat(cts, batch=batch, mod=model.matrix(~group+covar))


## Normalize library size
cts_norm <- apply(cts, 2, function(x){x/sum(x)})
cts_adj_norm <- apply(cts_adj, 2, function(x){x/sum(x)})
cts_adjori_norm <- apply(cts_adjori, 2, function(x){x/sum(x)})


## PCA
col_data <- data.frame(Batch=factor(batch), Group=factor(group), Cov=covar)

seobj <- SummarizedExperiment(assays=cts_norm, colData=col_data)
pca_obj <- plotPCA(DESeqTransform(seobj), intgroup=c("Batch", "Group"))
plt_batch <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=Batch)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
       title="Color by batch - Unadjusted") +
  theme(legend.position="bottom")
plt_cond <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=Group)) + 
  geom_point() + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
       title="Color by condition - Unadjusted") +
  theme(legend.position="bottom")
rm(seobj, pca_obj)

seobj_adj <- SummarizedExperiment(assays=cts_adj_norm, colData=col_data)
pca_obj_adj <- plotPCA(DESeqTransform(seobj_adj), intgroup=c("Batch", "Group"))
plt_batch_adj <- ggplot(pca_obj_adj$data, aes(x=PC1, y=PC2, color=Batch)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj$plot_env$percentVar[2])),
       title="Color by batch - ComBat-seq") +
  theme(legend.position="bottom")
plt_cond_adj <- ggplot(pca_obj_adj$data, aes(x=PC1, y=PC2, color=Group)) + 
  geom_point() + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj$plot_env$percentVar[2])),
       title="Color by condition - ComBat-seq") +
  theme(legend.position="bottom")

seobj_adjori <- SummarizedExperiment(assays=cts_adjori_norm, colData=col_data)
pca_obj_adjori <- plotPCA(DESeqTransform(seobj_adjori), intgroup=c("Batch", "Group"))
plt_batch_adjori <- ggplot(pca_obj_adjori$data, aes(x=PC1, y=PC2, color=Batch)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[2])),
       title="Color by batch - Original ComBat") +
  theme(legend.position="bottom")
plt_cond_adjori <- ggplot(pca_obj_adjori$data, aes(x=PC1, y=PC2, color=Group)) + 
  geom_point() + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[2])),
       title="Color by condition - Original ComBat") +
  theme(legend.position="bottom")

png("TB_PCA_combatseq.png", width=7, height=12, units="in", res=300)
grid.arrange(plt_batch, plt_cond, plt_batch_adj, plt_cond_adj, plt_batch_adjori, plt_cond_adjori, ncol=2, nrow=3)
dev.off()



## Explained Variance Analysis
# normalize samples by library size
var_exp_ori <- batchqc_explained_variation(cts_norm, condition=factor(group), batch=factor(batch))
var_exp_adj <- batchqc_explained_variation(cts_adj_norm, condition=factor(group), batch=factor(batch))
var_exp_adjori <- batchqc_explained_variation(cts_adjori_norm, condition=factor(group), batch=factor(batch))

plt_ori <- ggplot(melt(var_exp_ori$explained_variation), aes(x=Var2, y=value)) +
  geom_boxplot() +
  labs(y="Explained Variance", title="Unadjusted") +
  theme(axis.title.x = element_blank())
plt_adj <- ggplot(melt(var_exp_adj$explained_variation), aes(x=Var2, y=value)) +
  geom_boxplot() +
  labs(y="Explained Variance", title="ComBat-Seq") +
  theme(axis.title.x = element_blank())
plt_adjori <- ggplot(melt(var_exp_adjori$explained_variation), aes(x=Var2, y=value)) +
  geom_boxplot() +
  labs(y="Explained Variance", title="Original ComBat") +
  theme(axis.title.x = element_blank())
png("TB_explained_var.png", width=13, height=5, units="in", res=300)
grid.arrange(plt_ori, plt_adj, plt_adjori, ncol=3)
dev.off()



## Clustering
hc <- hclust(dist(t(cts_norm)))
dend <- as.dendrogram(hc)
dend <- color_branches(dend, groupLabels=batch[order.dendrogram(dend)], col=group[order.dendrogram(dend)]+2)
plot(dend)  

adj_hc <- hclust(dist(t(cts_adj_norm)))
adj_dend <- as.dendrogram(adj_hc)
adj_dend <- color_branches(adj_dend, groupLabels=batch[order.dendrogram(adj_dend)], col=group[order.dendrogram(adj_dend)]+2)
plot(adj_dend)



