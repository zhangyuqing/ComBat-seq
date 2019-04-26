rm(list=ls())
sapply(c("edgeR", "sva", "RUVSeq", "MASS", "ggplot2", "reshape2", "gridExtra", "ffpe"), require, character.only=TRUE)
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
source(file.path(script_dir, "ComBat_seq.R")); source(file.path(script_dir, "helper_seq.R"))
set.seed(123)
# load data
load("TB_ExpRep.RData")
cts <- cts[setdiff(rownames(cts),"HBA2"), ]  ##!!!!!! need to deal with this...
method_names <- c("raw.count", "one.step", "curr.combat", "ruvseq", "svaseq", "combatseq")

####  Parameters
#command_args <- commandArgs(trailingOnly=TRUE)
p_eval <- 1/3  # percentage of samples splitting into the evaluation set
iterations <- 300

####  Run pipeline
iter <- 1
concord_collect <- list()
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
  #dat_lst=eval_lst; DE_method=edgeR_DEpipe; alpha.unadj=1; alpha.fdr=1
  cat("\n####  DE on Evaluation set  ####\n")
  DE_eval <- try(DEpipe(dat_lst=eval_lst, DE_method=edgeR_DEpipe, alpha.unadj=1, alpha.fdr=1))
  if(class(DE_eval)=="try-error"){next}
  toptable_eval <- lapply(DE_eval, function(de){de$de_res})
  
  cat("\n####  DE on Verification set  ####\n")
  DE_veri <- try(DEpipe(dat_lst=veri_lst, DE_method=edgeR_DEpipe, alpha.unadj=1, alpha.fdr=1))
  if(class(DE_veri)=="try-error"){next}
  toptable_veri <- lapply(DE_veri, function(de){de$de_res})
  
  
  ## Rank genes with absolute logFC (or P value if logFC unavailable)
  eval_genes_ordered <- lapply(1:length(toptable_eval), function(i){
    if(names(toptable_eval)[i]=="curr.combat"){
      gordered <- rownames(eval_lst$counts)[order(toptable_eval[[i]]$PValue, decreasing=FALSE)]
    }else{
      tptb <- toptable_eval[[i]]
      gordered <- rownames(tptb)[order(abs(tptb$logFC), decreasing=TRUE)]
    }
    return(gordered)
  })
  
  veri_genes_ordered <- lapply(1:length(toptable_veri), function(i){
    if(names(toptable_veri)[i]=="curr.combat"){
      gordered <- rownames(veri_lst$counts)[order(toptable_veri[[i]]$PValue, decreasing=FALSE)]
    }else{
      tptb <- toptable_veri[[i]]
      gordered <- rownames(tptb)[order(abs(tptb$logFC), decreasing=TRUE)]
    }
    return(gordered)
  })
  
  
  ## Compute concordance
  concord_res <- list()
  for(jj in seq_along(veri_genes_ordered)){  ## jj - veri set
    concord_res[[jj]] <- list()
    for(ii in seq_along(eval_genes_ordered)){  ## ii - eval set
      concord_res[[jj]][[ii]] <- CATplot(eval_genes_ordered[[ii]], veri_genes_ordered[[jj]], make.plot=FALSE)
    }
    names(concord_res[[jj]]) <- names(toptable_eval)
  }
  names(concord_res) <- names(toptable_veri)
  
  concord_res_mlt <- melt(concord_res, id.vars=c("rank", "concordance"))
  colnames(concord_res_mlt)[3:4] <- c("Evaluation", "Verification")
  concord_res_mlt$Evaluation <- factor(concord_res_mlt$Evaluation, levels=method_names)
  concord_res_mlt$Verification <- factor(concord_res_mlt$Verification, levels=method_names)
  
  concord_collect[[iter]] <- concord_res_mlt
}

# remove NULL values (skipped iterations)
concord_collect <- concord_collect[!sapply(concord_collect, is.null)]

concord_mlt_summary <- data.frame(rank=concord_collect[[1]]$rank, Evaluation=concord_collect[[1]]$Evaluation, Verification=concord_collect[[1]]$Verification)
concord_mlt_summary$concordance <- rowMeans(do.call(cbind, lapply(concord_collect, function(cres){cres$concordance})))

save(concord_mlt_summary, concord_collect, file="CATdata.RData")

png("CATplot_TB.png", width=9, height=8, units="in", res=300)
ggplot(concord_mlt_summary, aes(x=rank, y=concordance, group=Evaluation, color=Evaluation)) +
  geom_line() +
  facet_wrap(~Verification, nrow=3) +
  labs(x="Rank", y="Concordance")
dev.off()
