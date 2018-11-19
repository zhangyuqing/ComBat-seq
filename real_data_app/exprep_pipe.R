rm(list=ls())
sapply(c("edgeR", "sva", "RUVSeq", "MASS", "ggplot2", "reshape2", "gridExtra"), require, character.only=TRUE)
demo <- FALSE  # if testing code, set as TRUE; if running simulations, set as FALSE
if(demo){
  setwd("~/Documents/ComBat_seq/real_data_app/")
  script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
  source(file.path(script_dir, "real_data_app/exprep_helpers.R"))
}else{
  #setwd("~/yuqingz/ComBat_seq/real_data_app/")
  setwd("/restricted/projectnb/combat/work/yuqingz/ComBat-seq/real_data_app")
  script_dir <- ".."
  source(file.path(script_dir, "exprep_helpers.R"))
}
source(file.path(script_dir, "ComBat_seq.R")); source(file.path(script_dir, "helper_seq.R"))
set.seed(123)
# load data
load("TB_ExpRep.RData")
#cts <- cts[setdiff(rownames(cts),"HBA2"), ]  ##!!!!!! need to deal with this...


####  Parameters
#command_args <- commandArgs(trailingOnly=TRUE)
p_eval <- 1/3  # percentage of samples splitting into the evaluation set
iterations <- 20
alpha_unadj <- 0.1 
alpha_fdr_seq <- seq(from=0, to=0.2, by=0.005)[-1]
  

####  Run pipeline
tpr_res <- fpr_res <- tprADJ_res <- fdrADJ_res <- list()
iter <- ii <- 1
for(ii in seq_along(alpha_fdr_seq)){
  alpha_fdr <- alpha_fdr_seq[ii]  # FDR cutoff level for evaluation set
  tpr_mat_lst <- fpr_mat_lst <- tprADJ_mat_lst <- fdrADJ_mat_lst <- list()
  
  for(iter in 1:iterations){
    print(paste("Iteration:", iter))
    
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
    genes_sel <- sample(1:nrow(eval_lst$counts), 2000, replace=FALSE)
    eval_lst$counts <- eval_lst$counts[genes_sel, ]
    veri_lst$counts <- veri_lst$counts[genes_sel, ]
    
    
    ## DE 
    #dat_lst=eval_lst; DE_method=edgeR_DEpipe; alpha.unadj=alpha_unadj; alpha.fdr=alpha_fdr
    cat("\n####  DE on Evaluation set  ####\n")
    DE_eval <- try(DEpipe(dat_lst=eval_lst, DE_method=edgeR_DEpipe, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr))
    if(class(DE_eval)=="try-error"){next}
    DEcalled_eval_unadj <- lapply(DE_eval, function(de_res){de_res$unadj})
    DEcalled_eval_fdr <- lapply(DE_eval, function(de_res){de_res$fdr})
    
    cat("\n####  DE on Verification set  ####\n")
    DE_veri <- try(DEpipe(dat_lst=veri_lst, DE_method=edgeR_DEpipe, alpha.unadj=alpha_unadj, alpha.fdr=0.1))
    if(class(DE_veri)=="try-error"){next}
    DEcalled_veri_unadj <- lapply(DE_veri, function(de_res){de_res$unadj})
    DEcalled_veri_fdr <- lapply(DE_veri, function(de_res){de_res$fdr})
    
    
    ## Compute performances
    tpr_mat_lst[[iter]] <- fpr_mat_lst[[iter]] <- matrix(NA, nrow=length(DE_eval), ncol=length(DE_veri), 
                                                         dimnames=list(names(DE_eval), names(DE_veri)))
    tprADJ_mat_lst[[iter]] <- fdrADJ_mat_lst[[iter]] <- matrix(NA, nrow=length(DE_eval), ncol=length(DE_veri), 
                                                               dimnames=list(names(DE_eval), names(DE_veri)))
    for(i in seq_along(DE_eval)){
      for(j in seq_along(DE_veri)){
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
  
  tpr_res <- c(tpr_res, list(CleanOutput(tpr_mat_lst, alpha=alpha_unadj)))
  fpr_res <- c(fpr_res, list(CleanOutput(fpr_mat_lst, alpha=alpha_unadj)))
  tprADJ_res <- c(tprADJ_res, list(CleanOutput(tprADJ_mat_lst, alpha=alpha_fdr)))
  fdrADJ_res <- c(fdrADJ_res, list(CleanOutput(fdrADJ_mat_lst, alpha=alpha_fdr)))
}

stats_out <- list(TPR=tpr_res, FPR=fpr_res, SENS=tprADJ_res, FDR=fdrADJ_res)
stats_out <- lapply(stats_out, function(s_out){do.call(rbind, s_out)})

save(stats_out, file="perf_stats.RData")
