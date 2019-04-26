rm(list=ls())
sapply(c("edgeR", "DESeq2", "sva", "RUVSeq", "MASS", "ggplot2", "reshape2", "gridExtra"), require, character.only=TRUE)
demo <- FALSE  # if testing code, set as TRUE; if running simulations, set as FALSE
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

# Parameter
FC_cutoff <- 5

# Load data
load("TB_ExpRep.RData")
cts <- cts[setdiff(rownames(cts),c("HBA2", "HBB")), ]  ##!!!!!! need to deal with this...
#method_names <- c("raw.count", "one.step", "curr.combat", "ruvseq", "svaseq", "combatseq")

# take subset of genes to reduce computation time
cts <- cts[1:2000, ]

# Remove all-0 genes in any batch
keep1 <- apply(cts[, batch==1], 1, function(x){!all(x==0)})
keep2 <- apply(cts[, batch==2], 1, function(x){!all(x==0)})
#cts <- cts[keep1 & keep2, ]
sum(!keep1); sum(!keep2)

# Find significant difference WRT batch - library composition batch effect
# within each condition then take intersect
#cts_norm <- apply(cts, 2, function(x){x/sum(x)})  # normalize library size
libcomp_genes_list <- lapply(levels(factor(group)), function(group_level){
  batch_sub <- batch[group==group_level]
  cts_sub <- cts[, group==group_level]
  
  y <- DGEList(counts=cts_sub, group=batch_sub)
  y <- calcNormFactors(y, method="TMM")
  batch_mod <- model.matrix(~factor(batch_sub))
  y <- estimateDisp(y, batch_mod)
  fit <- glmQLFit(y, batch_mod)
  qlf <- glmQLFTest(fit, coef=2)
  de_res <- topTags(qlf, n=nrow(cts_sub))$table
  de_called <- rownames(de_res)[de_res$logFC > log(FC_cutoff)]
  return(de_called)
})

libcomp_genes <- Reduce(intersect, libcomp_genes_list)

# Find significant difference WRT condition - library composition batch effect
# within each batch then take union
DEbatch_genes_list <- lapply(levels(factor(batch)), function(batch_level){
  group_sub <- group[batch==batch_level]
  cts_sub <- cts[, batch==batch_level]
  
  y <- DGEList(counts=cts_sub, group=group_sub)
  y <- calcNormFactors(y, method="TMM")
  design <- model.matrix(~factor(group_sub))
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef=2)
  de_res <- topTags(qlf, n=nrow(cts_sub))$table
  #de_called <- rownames(de_res)[abs(de_res$logFC) > log(FC_cutoff)]
  de_called <- rownames(de_res)[de_res$FDR < 0.05]
  return(de_called)
})

libcomp_genes <- setdiff(libcomp_genes, Reduce(union, DEbatch_genes_list))


# DE after batch correction
dat_lst <- list(counts=cts, batch=batch, group=group, covar_incl=covar)
DE_res <- DEpipe(dat_lst, DE_method=edgeR_DEpipe, alpha.unadj=1, alpha.fdr=1)
DE_res <- lapply(DE_res, function(de){de$de_res})
save(DE_res, file="TB_DEpipe_results.RData")

load("TB_DEpipe_results.RData")
# Locate genes affected by library composition effect in DE results
libcomp_fdr <- lapply(DE_res, function(de){de[libcomp_genes, "FDR"]})
#boxplot(libcomp_fdr)
libcomp_fdr_mat <- melt(libcomp_fdr)
colnames(libcomp_fdr_mat) <- c("fdr", "Method")
libcomp_fdr_mat$Method <- factor(libcomp_fdr_mat$Method, levels=c("raw.count", "one.step", "ruvseq", "svaseq", "combatseq"))
png(sprintf("TB_fdr_FC%s.png", FC_cutoff), width=6, height=6, units="in", res=300)
ggplot(libcomp_fdr_mat, aes(x=Method, y=fdr, fill=Method)) +
  geom_violin() +
  #geom_boxplot(width=0.2) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_hline(yintercept=0.05, color="blue", linetype="dashed") +
  labs(y="FDR of non-DE genes affected by batch effect",
       title=paste("batch effect FC cutoff =", FC_cutoff)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

libcomp_ranks <- lapply(DE_res, function(de){match(libcomp_genes, rownames(de))})
#boxplot(libcomp_ranks)
libcomp_ranks_mat <- melt(libcomp_ranks)
colnames(libcomp_ranks_mat) <- c("ranks", "Method")
libcomp_ranks_mat$Method <- factor(libcomp_ranks_mat$Method, levels=c("raw.count", "one.step", "ruvseq", "svaseq", "combatseq"))
png(sprintf("TB_ranks_FC%s.png", FC_cutoff), width=6, height=6, units="in", res=300)
ggplot(libcomp_ranks_mat, aes(x=Method, y=ranks, fill=Method)) +
  geom_violin() +
  #geom_boxplot(width=0.2) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(y="Ranks of non-DE genes affected by batch effect",
       title=paste("batch effect FC cutoff =", FC_cutoff)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
