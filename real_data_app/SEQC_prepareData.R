rm(list=ls())
setwd("~/Documents/ComBat_seq/real_data_app/")
sapply(c("edgeR", "DESeq2", "sva", "RUVSeq", "MASS", "ggplot2", "reshape2", "gridExtra", "seqc", "ddCt", "scales", "BatchQC"), 
       require, character.only=TRUE)
demo <- TRUE  # if testing code, set as TRUE; if running simulations, set as FALSE
if(demo){
  setwd("~/Documents/ComBat_seq/real_data_app/")
  script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
  source(file.path(script_dir, "real_data_app/libcomp_helpers.R"))
}else{
  #setwd("~/yuqingz/ComBat_seq/real_data_app/")
  setwd("/restricted/projectnb/combat/work/yuqingz/ComBat-seq/real_data_app")
  script_dir <- ".."
  source(file.path(script_dir, "libcomp_helpers.R"))
}
source(file.path(script_dir, "ComBat_seq.R")); source(file.path(script_dir, "helper_seq.R"))
#set.seed(123)


##  Load data
rnaseq_eset <- seqc.eSet("gene", "refseq")
rnaseq_eset_sub <- rnaseq_eset[, (rnaseq_eset$platform=="ILM") & (rnaseq_eset$sample %in% c("A", "B"))]
#rnaseq_eset_sub <- rnaseq_eset_sub[, rnaseq_eset_sub$center %in% c("BGI", "CNL")]
cts <- exprs(rnaseq_eset_sub)
batch <- as.character(rnaseq_eset_sub$center)
group <- as.character(rnaseq_eset_sub$sample)

# remove genes with only 0 counts in any batch
keep_lst <- lapply(levels(factor(batch)), function(batch_level){
  which(apply(cts[, batch==batch_level], 1, function(x){!all(x==0)}))
})
keep <- Reduce(intersect, keep_lst)
cts <- cts[keep, ]


##  Look for batch effect 
# - plot PCA
logcpm <- cpm(cts, log=TRUE)
pca_obj <- prcomp(t(logcpm), center=TRUE, scale.=TRUE)
plt_df <- data.frame(Batch=batch, Group=group, PC1=pca_obj$x[,1], PC2=pca_obj$x[,2])
plt_df_mlt <- melt(plt_df, id.vars=c("Batch", "Group"))
p_group <- ggplot(plt_df, aes(x=PC1, y=PC2, color=Group)) +
  geom_point() +
  labs(x=sprintf("PC1 (%s variance explained)", percent(pca_obj$sdev[1]/sum(pca_obj$sdev))),
       y=sprintf("PC2 (%s variance explained)", percent(pca_obj$sdev[2]/sum(pca_obj$sdev))),
       title="Colored by group (individuals)")
p_batch <- ggplot(plt_df, aes(x=PC1, y=PC2, color=Batch)) +
  geom_point() +
  labs(x=sprintf("PC1 (%s variance explained)", percent(pca_obj$sdev[1]/sum(pca_obj$sdev))),
       y=sprintf("PC2 (%s variance explained)", percent(pca_obj$sdev[2]/sum(pca_obj$sdev))),
       title="Colored by batch (locations)")
png("SEQC_PCA.png", width=11, height=5, units="in", res=300)
grid.arrange(p_group, p_batch, ncol=2)
dev.off()

# - variance explained by batch / group
va_res <- batchqc_explained_variation(logcpm, condition=group, batch=batch)$explained_variation
png("SEQC_explained_variation.png", width=5, height=5, units="in", res=300)
ggplot(melt(as.data.frame(va_res)), aes(x=variable, y=value)) +
  geom_boxplot() +
  labs(y="Explained Variation") +
  theme(axis.title.x = element_blank()) +
  ylim(0, max(va_res))
dev.off()

# - dispersion difference across batch
disp_lst <- lapply(levels(factor(batch)), function(batch_level){
  
})


##  DE after batch adjustment methods
dat_lst <- list(counts=cts, batch=batch, group=group, covar_incl=NULL)
DE_res <- DEpipe(dat_lst, DE_method=edgeR_DEpipe, alpha.unadj=1, alpha.fdr=1)
DE_tables <- lapply(DE_res, function(de){de$de_res})


##  Ground truth
taqman_A <- taqman[, c("A1_detection", "A2_detection", "A3_detection", "A4_detection")]
taqman_B <- taqman[, c("B1_detection", "B2_detection", "B3_detection", "B4_detection")]
de_genes_ind <- which(sapply(1:nrow(taqman), function(i){
  (all(taqman_A[i, ]=="P") & all(taqman_B[i, ]=="A")) | (all(taqman_A[i, ]=="A") & all(taqman_B[i, ]=="P"))
}))
de_genes <- as.character(taqman[de_genes_ind, "EntrezID"])
de_genes <- de_genes[!is.na(de_genes)]


##  FDR of DE genes
#boxplot(lapply(DE_tables, function(de){match(de_genes, rownames(de))}))
DE_genes <- lapply(DE_tables, function(de_tb){rownames(de_tb)[abs(de_tb$logFC) > 2]})
lapply(DE_genes, function(de_pred){length(intersect(de_pred, de_genes))/length(de_genes)})


