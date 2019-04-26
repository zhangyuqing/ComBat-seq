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
ComBat_seq_old <- ComBat_seq; rm(ComBat_seq)
source(file.path(script_dir, "ComBat_seq_newEB_logParams.R")) 
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


## Use ComBatSeq (old) to adjust data
group_sub <- as.character(group_sub)
group_sub_num <- rep(1, length(group_sub))
group_sub_num[grep("gfp", group_sub)] <- 0
#counts=cts_sub;batch=batch_sub;group=group_sub_num;gene.subset.n=1000;Cpp=FALSE;covar_mod=NULL;full_mod=TRUE
start_time <- Sys.time()
combatseq_sub_old <- ComBat_seq_old(counts=cts_sub, batch=batch_sub, group=group_sub_num, gene.subset.n=1000, Cpp=FALSE, 
                                    shrink=FALSE, shrink.disp=FALSE)
end_time <- Sys.time()
print(end_time - start_time)
#combatseq_sub <- readRDS("~/Google Drive/ComBat_seq/real_data_example/RNAseq/gfrn_signature/combatseq_test.rds")


## Use ComBatSeq - new EB to adjust data
ctrl_obj <- ComBat_seq_control_gene(cts_sub, batch=batch_sub, group=group_sub_num)
start_time <- Sys.time()
combatseq_sub <- ComBat_seq(cts_sub, batch=batch_sub, group=group_sub_num, ctrl=ctrl_obj)
end_time <- Sys.time()
print(end_time - start_time)


## Use original ComBat
combat_sub <- ComBat(cts_sub, batch=batch_sub, mod=model.matrix(~factor(group_sub_num)))


## Normalize library size
cts_norm <- apply(cts_sub, 2, function(x){x/sum(x)})
cts_adj_norm <- apply(combatseq_sub$adjust_counts, 2, function(x){x/sum(x)})
cts_adj_norm_old <- apply(combatseq_sub_old, 2, function(x){x/sum(x)})
cts_adjori_norm <- apply(combat_sub, 2, function(x){x/sum(x)})


## PCA 
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

seobj_adj <- SummarizedExperiment(assays=cts_adj_norm_old, colData=col_data)
pca_obj_adj <- plotPCA(DESeqTransform(seobj_adj), intgroup=c("Batch", "Group"))
plt_adj <- ggplot(pca_obj_adj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj$plot_env$percentVar[2])),
       title="Color by Batch - ComBat-Seq (shrink off)") +
  theme(legend.position="bottom")

seobj_adj_newEB <- SummarizedExperiment(assays=cts_adj_norm, colData=col_data)
pca_obj_adj_newEB <- plotPCA(DESeqTransform(seobj_adj_newEB), intgroup=c("Batch", "Group"))
plt_adj_newEB <- ggplot(pca_obj_adj_newEB$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj_newEB$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj_newEB$plot_env$percentVar[2])),
       title="Color by Batch - ComBat-Seq (new EB)") +
  theme(legend.position="bottom")

seobj_adjori <- SummarizedExperiment(assays=cts_adjori_norm, colData=col_data)
pca_obj_adjori <- plotPCA(DESeqTransform(seobj_adjori), intgroup=c("Batch", "Group"))
plt_adjori <- ggplot(pca_obj_adjori$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[2])),
       title="Color by Batch - Original ComBat") +
  theme(legend.position="bottom")

png(sprintf("gfrn_newEB_PCA_%s.png", gsub("^", "", pathway_regex)), 
    width=10, height=10, units="in", res=300)
grid.arrange(plt, plt_adjori, plt_adj, plt_adj_newEB, ncol=2, nrow=2)
dev.off()


##  Parameter shrinkage plot
ShrinkParamPlt <- function(cts, combatseq_res, batch_name, eb_viz_genes=NULL, outlier_gene_ind=NULL){
  is_outlier_gene <- rep(FALSE, nrow(cts)); is_outlier_gene[outlier_gene_ind] <- TRUE
  ebON_outlier_type <- rep("EB On", nrow(cts)); ebON_outlier_type[is_outlier_gene] <- "EB On - outliers"
  ebOFF_outlier_type <- rep("EB Off", nrow(cts)); ebOFF_outlier_type[is_outlier_gene] <- "EB Off - outliers"
  
  if(!is.null(eb_viz_genes)){
    params_est <- rbind(data.frame(genes=rownames(cts)[eb_viz_genes],
                                   gamma=combatseq_res$gamma_hat[eb_viz_genes, paste0('batch',batch_name)],
                                   phi=combatseq_res$phi_hat[eb_viz_genes, paste0('batch',batch_name)],
                                   EB.outlier.type=ebOFF_outlier_type[eb_viz_genes]),
                        data.frame(genes=rownames(cts)[eb_viz_genes],
                                   gamma=combatseq_res$gamma_star[eb_viz_genes, paste0('batch',batch_name)],
                                   phi=combatseq_res$phi_star[eb_viz_genes, paste0('batch',batch_name)],
                                   EB.outlier.type=ebON_outlier_type[eb_viz_genes]))
  }else{
    params_est <- rbind(data.frame(genes=rownames(cts), 
                                   gamma=combatseq_res$gamma_hat[, paste0('batch',batch_name)], 
                                   phi=combatseq_res$phi_hat[, paste0('batch',batch_name)], 
                                   EB.outlier.type=ebOFF_outlier_type),
                        data.frame(genes=rownames(cts), 
                                   gamma=combatseq_res$gamma_star[, paste0('batch',batch_name)], 
                                   phi=combatseq_res$phi_star[, paste0('batch',batch_name)], 
                                   EB.outlier.type=ebON_outlier_type))
  }
  
  plt <- ggplot(params_est) +
    geom_point(aes(x=gamma, y=phi, color=EB.outlier.type)) +
    scale_color_manual(values=c("light green", "orange", "dark green", "red")) + 
    labs(x="Mean batch effect estimates (gamma)", y="Dispersion estimates (phi)", 
         title=sprintf("Parameter estimates within batch %s", batch_name))
  if(!is.null(eb_viz_genes)){
    plt <- plt + geom_line(aes(x=gamma, y=phi, group=genes), color="grey")
  }
  return(plt)
}

# eb_plts <- lapply(sort(unique(batch)), function(b){
#   curr_out <- unique(do.call(c, ctrl_obj$outlier.genes[[b]]))
#   out_plt <- sample(curr_out,50)
#   curr_nonout <- setdiff(1:nrow(cts_sub), curr_out)
#   nonout_plt <- sample(curr_nonout,50)
#   ShrinkParamPlt(cts_sub, combatseq_sub, b, eb_viz_genes=c(out_plt, nonout_plt), outlier_gene_ind=out_plt)
# })
# png(sprintf("gfrn_ShrinkParams_newEB_%s.png", gsub("^", "", pathway_regex)), 
#     width=15, height=5, units="in", res=300)
# grid.arrange(eb_plts[[1]], eb_plts[[2]], eb_plts[[3]], ncol=3)
# dev.off()

eb_plts_sub <- lapply(sort(unique(batch)), function(b){
  curr_out <- unique(do.call(c, ctrl_obj$outlier.genes[[b]]))
  ShrinkParamPlt(cts_sub, combatseq_sub, b, eb_viz_genes=curr_out,100)
})
png(sprintf("gfrn_newEB_ShrinkParams_subgenes_%s.png", gsub("^", "", pathway_regex)), 
    width=6, height=5, units="in", res=300)
#grid.arrange(eb_plts_sub[[1]], eb_plts_sub[[2]], eb_plts_sub[[3]], ncol=3)
eb_plts_sub[[1]]
dev.off()


## DE
source("~/Dropbox/Work/ComBat_Seq/ComBat-Seq/real_data_app/exprep_helpers.R")
DE_method <- "edgeR"
de_measure_name <- "PValue"
if(DE_method=="edgeR"){DE_function <- edgeR_DEpipe}else{DE_function <- DESeq2_DEpipe}

de_unadj <- DE_function(cts_sub, batch=batch_sub, group=group_sub_num, include.batch=FALSE, 
                        alpha.unadj=1, alpha.fdr=1, covar_incl=NULL, covar=NULL)
de_onestep <-  DE_function(cts_sub, batch=batch_sub, group=group_sub_num, include.batch=TRUE, 
                           alpha.unadj=1, alpha.fdr=1, covar_incl=NULL, covar=NULL)
de_combatseq_newEB <- DE_function(combatseq_sub$adjust_counts, batch=batch_sub, group=group_sub_num, include.batch=FALSE, 
                                  alpha.unadj=1, alpha.fdr=1, covar_incl=NULL, covar=NULL)
de_combatseq_old <- DE_function(combatseq_sub_old, batch=batch_sub, group=group_sub_num, include.batch=FALSE, 
                                alpha.unadj=1, alpha.fdr=1, covar_incl=NULL, covar=NULL)
  
n_DE <- c(unadjusted=length(which(de_unadj$de_res[, de_measure_name] < 0.05))/nrow(cts_sub), #length(which(de_unadj$de_res$PValue < 0.05)), 
          one.step=length(which(de_onestep$de_res[, de_measure_name] < 0.05))/nrow(cts_sub),
          combatseq.old=length(which(de_combatseq_old$de_res[, de_measure_name] < 0.05))/nrow(cts_sub),
          combatseq.newEB=length(which(de_combatseq_newEB$de_res[, de_measure_name] < 0.05))/nrow(cts_sub))
n_DE_df <- data.frame(Method=factor(names(n_DE), levels=names(n_DE)), N.DE=n_DE)
png(sprintf("gfrn_newEB_nDE_%s_%s.png", pathway_regex, DE_method), width=5, height=5, units="in", res=300)
ggplot(n_DE_df, aes(x=Method, y=n_DE)) +
  geom_bar(stat='identity') +
  labs(y="Percentage of called DE genes (unadjusted P < 0.05)", 
       title=sprintf("GFRN GFPs + %s - %s", toupper(gsub("^", "", pathway_regex, fixed=T)), DE_method)) +
  theme(axis.title.x=element_blank())
dev.off()
