rm(list=ls())
sapply(c("edgeR", "sva", "RUVSeq", "MASS", "ggplot2", "reshape2", "gridExtra"), require, character.only=TRUE)
demo <- FALSE #TRUE  # if testing code, set as TRUE; if running simulations, set as FALSE
if(demo){
  setwd("~/Documents/ComBat_seq/real_data_app/")
  script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
  source(file.path(script_dir, "real_data_app/exprep_helpers.R"))
}else{
  setwd("~/yuqingz/ComBat_seq/real_data_app/")
  script_dir <- ".."
  source(file.path(script_dir, "exprep_helpers.R"))
}
source(file.path(script_dir, "ComBat_seq.R")); source(file.path(script_dir, "helper_seq.R"))
set.seed(123)
# load data
load("TB_ExpRep.RData")
cts <- cts[setdiff(rownames(cts),"HBA2"), ]


####  Parameters
p_eval <- 1/3  # percentage of samples splitting into the evaluation set
iterations <- 20
alpha <- 0.1


####  Run pipeline
tpr_mat_lst <- fpr_mat_lst <- tprADJ_mat_lst <- fdrADJ_mat_lst <- list()
iter <- 1
for(iter in 1:iterations){
  cat(sprintf("Simulation: %s\n\n", iter))
  
  ## Randomly split whole dataset into evaluation set & larger verification set
  splt_res <- splitData(batch=batch, group=group, p_eval=p_eval)
  eval_lst <- list(counts=cts[, splt_res$eval], batch=batch[splt_res$eval], 
                   group=group[splt_res$eval], covar_incl=covar[splt_res$eval])
  veri_lst <- list(counts=cts[, splt_res$veri], batch=batch[splt_res$veri], 
                   group=group[splt_res$veri], covar_incl=covar[splt_res$veri])
  
  ## Remove genes with only 0s in any batch in eval set or in veri set
  keep_eval1 <- apply(eval_lst$counts[, eval_lst$batch==1], 1, function(x1){!all(x1==0)})
  keep_eval2 <- apply(eval_lst$counts[, eval_lst$batch==2], 1, function(x2){!all(x2==0)})
  keep_veri1 <- apply(veri_lst$counts[, veri_lst$batch==1], 1, function(y1){!all(y1==0)})
  keep_veri2 <- apply(veri_lst$counts[, veri_lst$batch==2], 1, function(y2){!all(y2==0)})
  eval_lst$counts <- eval_lst$counts[keep_eval1 & keep_eval2 & keep_veri1 & keep_veri2, ]
  veri_lst$counts <- veri_lst$counts[keep_eval1 & keep_eval2 & keep_veri1 & keep_veri2, ]
  
  ## cut the size of gene set to save computation time of ComBat-seq
  genes_sel <- sample(1:nrow(eval_lst$counts), 5000, replace=FALSE)
  eval_lst$counts <- eval_lst$counts[genes_sel, ]
  veri_lst$counts <- veri_lst$counts[genes_sel, ]
  
  
  ## DE 
  #dat_lst=eval_lst; DE_method=edgeR_DEpipe
  DE_eval <- try(DEpipe(dat_lst=eval_lst, DE_method=edgeR_DEpipe, alpha))
  if(class(DE_eval)=="try-error"){next}
  DEcalled_eval_unadj <- lapply(DE_eval, function(de_res){de_res$unadj})
  DEcalled_eval_fdr <- lapply(DE_eval, function(de_res){de_res$fdr})
  
  DE_veri <- try(DEpipe(dat_lst=veri_lst, DE_method=edgeR_DEpipe, alpha))
  if(class(DE_veri)=="try-error"){next}
  DEcalled_veri_unadj <- lapply(DE_veri, function(de_res){de_res$unadj})
  DEcalled_veri_fdr <- lapply(DE_veri, function(de_res){de_res$fdr})
  
  
  ## Compute performances
  tpr_mat_lst[[iter]] <- fpr_mat_lst[[iter]] <- matrix(NA, nrow=length(DE_eval), ncol=length(DE_veri), 
                                                       dimnames=list(names(DE_eval), names(DE_veri)))
  tprADJ_mat_lst[[iter]] <- fdrADJ_mat_lst[[iter]] <- matrix(NA, nrow=length(DE_eval), ncol=length(DE_veri), 
                                                             dimnames=list(names(DE_eval), names(DE_veri)))
  for(i in 1:length(DE_eval)){
    for(j in 1:length(DE_veri)){
      curr_stats <- perfStats(called_vec=DEcalled_eval_unadj[[i]], ground_truth_vec=DEcalled_veri_unadj[[j]], 
                              N_genes=nrow(eval_lst$counts))
      tpr_mat_lst[[iter]][i, j] <- curr_stats["tpr"]
      fpr_mat_lst[[iter]][i, j] <- curr_stats["fpr"]
      
      curr_stats_fdr <- perfStats(called_vec=DEcalled_eval_fdr[[i]], ground_truth_vec=DEcalled_veri_fdr[[j]], 
                                  N_genes=nrow(eval_lst$counts))
      tprADJ_mat_lst[[iter]][i, j] <- curr_stats_fdr["tpr"]
      fdrADJ_mat_lst[[iter]][i, j] <- 1 - curr_stats_fdr["prec"]
    }
  }
  
  tpr_mat_lst[[iter]] <- melt(tpr_mat_lst[[iter]])
  fpr_mat_lst[[iter]] <- melt(fpr_mat_lst[[iter]])
  tprADJ_mat_lst[[iter]] <- melt(tprADJ_mat_lst[[iter]])
  fdrADJ_mat_lst[[iter]] <- melt(fdrADJ_mat_lst[[iter]])
}

save(tpr_mat_lst, fpr_mat_lst, tprADJ_mat_lst, fdrADJ_mat_lst, file="perf_stats.RData")



####  Visualize performances with boxplots
load("perf_stats.RData")

## adjusted FDR 
# remove NA iterations
tprADJ_mat_lst <- tprADJ_mat_lst[sapply(tprADJ_mat_lst, function(tpr_m){!is.null(tpr_m)})]
fdrADJ_mat_lst <- fdrADJ_mat_lst[sapply(fdrADJ_mat_lst, function(fdr_m){!is.null(fdr_m)})]
# clean up performance over iterations
tprADJ_mat <- do.call(rbind, tprADJ_mat_lst)
fdrADJ_mat <- do.call(rbind, fdrADJ_mat_lst)
colnames(tprADJ_mat) <- colnames(fdrADJ_mat) <- c("Evaluation", "Verification", "Stats")
# generate boxplots
png("perf_stats_adj.png", width=13, height=6.5, units="in", res=300)
sens_plt <- ggplot(tprADJ_mat, aes(x=Evaluation, y=Stats, color=Evaluation)) +
  geom_boxplot() +
  facet_wrap(~Verification, ncol=2) +
  labs(y="Sensitivity") +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        axis.title.x=element_blank())
fdr_plt <- ggplot(fdrADJ_mat, aes(x=Evaluation, y=Stats, color=Evaluation)) +
  geom_boxplot() +
  geom_hline(yintercept=alpha, color="blue", linetype="dashed") +
  facet_wrap(~Verification, ncol=2) +
  labs(y="1 - Precision (FDR)") +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        axis.title.x=element_blank())
grid.arrange(sens_plt, fdr_plt, ncol=2)
dev.off()

## unadjusted P values
# remove NA iterations
tpr_mat_lst <- tpr_mat_lst[sapply(tpr_mat_lst, function(tpr_m){!is.null(tpr_m)})]
fpr_mat_lst <- fpr_mat_lst[sapply(fpr_mat_lst, function(fpr_m){!is.null(fpr_m)})]
# clean up performance over iterations
tpr_mat <- do.call(rbind, tpr_mat_lst)
fpr_mat <- do.call(rbind, fpr_mat_lst)
colnames(tpr_mat) <- colnames(fpr_mat) <- c("Evaluation", "Verification", "Stats")
# generate boxplots
png("perf_stats.png", width=13, height=6.5, units="in", res=300)
tpr_plt <- ggplot(tpr_mat, aes(x=Evaluation, y=Stats, color=Evaluation)) +
  geom_boxplot() +
  facet_wrap(~Verification, ncol=2) +
  labs(y="TPR") +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        axis.title.x=element_blank())
fpr_plt <- ggplot(fpr_mat, aes(x=Evaluation, y=Stats, color=Evaluation)) +
  geom_boxplot() +
  geom_hline(yintercept=alpha, color="blue", linetype="dashed") +
  facet_wrap(~Verification, ncol=2) +
  labs(y="FPR") +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        axis.title.x=element_blank())
grid.arrange(tpr_plt, fpr_plt, ncol=2)
dev.off()
