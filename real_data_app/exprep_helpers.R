## Randomly split whole dataset into evaluation set & larger verification set
splitData <- function(batch, group, p_eval){
  n_samples <- length(batch)
  
  splt <- split(1:n_samples, factor(batch))
  splt <- lapply(levels(factor(batch)), function(b){split(splt[[b]], factor(group[batch==b]))})
  
  eval_ind <- veri_ind <- c()
  for(i in 1:nlevels(factor(batch))){
    for(j in 1:nlevels(factor(group))){
      curr_sep <- sample(c("eval","veri"), length(splt[[i]][[j]]), prob=c(p_eval, 1-p_eval), replace=TRUE)
      eval_ind <- c(eval_ind, splt[[i]][[j]][curr_sep=="eval"])
      veri_ind <- c(veri_ind, splt[[i]][[j]][curr_sep=="veri"])
    }
  }
  
  out <- list(eval=eval_ind, veri=veri_ind)
  return(out)
}


## DE pipelines
DEpipe <- function(dat_lst, DE_method, alpha.unadj, alpha.fdr){
  cts <- dat_lst$counts 
  batch <- dat_lst$batch
  group <- dat_lst$group
  covar_incl <- dat_lst$covar_incl
  
  # data with batch effect
  cat("## Raw count\n")
  de_called1 <- DE_method(cts, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha.unadj, alpha.fdr=alpha.fdr, covar_incl=covar_incl)  
  # one-step
  cat("## One-step edgeR\n")
  de_called2 <- DE_method(cts, batch=batch, group=group, include.batch=TRUE, alpha.unadj=alpha.unadj, alpha.fdr=alpha.fdr, covar_incl=covar_incl)  
  # current ComBat 
  cat("## Current ComBat & lm\n")
  de_called3 <- currComBat_lm_DEpipe(cts=cts, batch=batch, group=group, alpha.unadj=alpha.unadj, alpha.fdr=alpha.fdr, covar_incl=covar_incl)
  # ComBat-seq 
  cat("## ComBat-seq\n")
  covar_mod <- model.matrix(~factor(covar_incl))
  adj_counts_combatseq <- ComBat_seq(counts=cts, batch=batch, group=group, covar_mod=covar_mod)
  de_called4 <- DE_method(cts=adj_counts_combatseq, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha.unadj, alpha.fdr=alpha.fdr, covar_incl=covar_incl)
  # RUVseq 
  cat("## RUV-seq\n")
  emps <- tail(rownames(de_called1$de_res), n=10)
  uvseq <- RUVg(cts, cIdx=emps, k=1)
  de_called5 <- DE_method(cts, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha.unadj, alpha.fdr=alpha.fdr, covar_incl=covar_incl, covar=uvseq$W)  
  # SVAseq 
  cat("## SVA-seq\n")
  if(length(unique(covar_incl))>1){
    mod1 <- model.matrix(~as.factor(group) + covar_incl)
    mod0 <- model.matrix(~covar_incl)
  }else{
    mod1 <- model.matrix(~as.factor(group))
    mod0 <- cbind(mod1[,1])
  }
  svseq <- svaseq(cts, mod1, mod0, n.sv=1);cat("\n")
  de_called6 <- DE_method(cts=cts, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha.unadj, alpha.fdr=alpha.fdr, covar_incl=covar_incl, covar=svseq$sv)  
  
  out <- list(raw.count=de_called1, one.step=de_called2, curr.combat=de_called3, 
              combatseq=de_called4, ruvseq=de_called5, svaseq=de_called6)
  return(out)
}


# edgeR DE pipeline
edgeR_DEpipe <- function(cts, batch, group, include.batch, alpha.unadj, alpha.fdr, covar_incl, covar=NULL){
  y <- DGEList(counts=cts)
  y <- calcNormFactors(y, method="TMM")
  if(include.batch){
    cat("Including batch as covariate\n")
    design <- model.matrix(~ as.factor(group) + as.factor(batch))
  }else{
    cat("Default group as model matrix\n")
    design <- model.matrix(~as.factor(group))
  }
  if(!is.null(covar)){
    cat("Including surrogate variables or unwanted variation variables\n")
    design <- cbind(design, covar)
  }
  # include gender variable
  if(length(unique(covar_incl))>1){
    design <- cbind(design, Sex=model.matrix(~covar_incl)[,2])
  }
  
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef=2)
  de_res <- topTags(qlf, n=nrow(cts))$table
  
  de_called <- rownames(de_res)[de_res$PValue < alpha.unadj]
  de_called_fdr <- rownames(de_res)[de_res$FDR < alpha.fdr]
  return(list(unadj=de_called, fdr=de_called_fdr, de_res=de_res, design=design))
}


##  DESeq2 pipeline
DESeq2_DEpipe <- function(counts_mat, batch, group, include.batch, alpha.unadj, alpha.fdr, covar_incl, covar=NULL){
  if(include.batch){
    cat("Including batch as covariate\n")
    col_data <- data.frame(Batch=as.factor(batch), Group=as.factor(group))
    design_formula <- ~Batch+Group
  }else if(!is.null(covar)){
    cat("Including surrogate variables or unwanted variation variables\n")
    col_data <- cbind(model.matrix(~as.factor(group)), covar)
    colnames(col_data)[2:ncol(col_data)] <- c("Group", paste0("Covar", 1:(ncol(col_data)-2)))
    rownames(col_data) <- colnames(counts_mat)
    col_data <- as.data.frame(col_data); col_data$Group <- as.factor(col_data$Group); col_data <- col_data[, -1]
    design_formula <- as.formula(paste("~", paste(colnames(col_data), collapse="+")))
  }else{
    cat("Default group as model matrix\n")
    col_data <- data.frame(Group=as.factor(group))
    design_formula <- ~Group
  }
  
  # deal with covar_incl
  if(!is.null(covar_incl)){
    covar_mod <- model.matrix(~factor(covar_incl))
    colnames(covar_mod)[-1] <- paste0("cov",1:ncol(covar_mod))
    col_data <- data.frame(col_data, covar_mod)
    design_formula 
  }
  
  dds <- DESeqDataSetFromMatrix(countData=counts_mat, colData=col_data, design=design_formula)
  dds <- DESeq(dds)
  de_res <- results(dds, name="Group_1_vs_0")
  
  de_called <- rownames(de_res)[which(de_res$pvalue < alpha.unadj)]
  de_called_fdr <- rownames(de_res)[which(de_res$padj < alpha.fdr)]
  return(list(unadj=de_called, fdr=de_called_fdr, de_res=de_res, design=design))
}


# current ComBat + lm DE pipeline
currComBat_lm_DEpipe <- function(cts, batch, group, alpha.unadj, alpha.fdr, covar_incl){
  log_counts <- cpm(cts, log=TRUE)  # use logCPM to make data more normal
  adj_counts <- ComBat(log_counts, batch=batch, mod=model.matrix(~as.factor(group) + covar_incl))
  pval_seq <- apply(adj_counts, 1, function(x, group){
    x_norm <- scale(x, center=TRUE, scale=TRUE)
    fit3 <- lm(x_norm ~as.factor(group) + covar_incl)
    return(summary(fit3)$coefficients[2, 4])
  }, group=group)
  padj_seq <- p.adjust(pval_seq, method="fdr")
  
  de_called <- list(unadj=rownames(cts)[pval_seq < alpha.unadj], fdr=rownames(cts)[padj_seq < alpha.fdr], 
                    de_res=data.frame(PValue=pval_seq, FDR=padj_seq), design=model.matrix(~as.factor(group) + covar_incl))
  return(de_called)
}


##  Model performance
perfStats <- function(called_vec, ground_truth_vec, N_genes){
  if(length(called_vec)==0){
    tpr <- fpr <- 0
    prec <- NA
  }else{
    tp <- length(intersect(called_vec, ground_truth_vec))
    fp <- length(setdiff(called_vec, ground_truth_vec)) 
    N_DE <- length(ground_truth_vec)
    N_nonDE <- N_genes - N_DE
    
    tpr <- tp / N_DE
    fpr <- fp / N_nonDE
    prec <- tp / length(called_vec)
  }
  return(c(tpr=tpr, fpr=fpr, prec=prec))
}


##  Format output
CleanOutput <- function(stats_lst, alpha){
  # remove NA iterations
  stats_lst <- stats_lst[sapply(stats_lst, function(st){!is.null(st)})]
  # clean up performance over iterations
  stats_mat <- do.call(rbind, stats_lst)
  colnames(stats_mat) <- c("Evaluation", "Verification", "Stats")
  stats_mat$Cutoff <- alpha
  return(stats_mat)
}

