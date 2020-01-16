## library size difference
constructFCMatrix <- function(G, n_group, G_ups, G_downs, bioFC, batchFC){
  fold_changes <- matrix(NA, nrow=G, ncol=n_group) 
  fold_changes[G_ups, ] <- matrix(rep(c(1, bioFC, batchFC, bioFC*batchFC), length(G_ups)), ncol=n_group, byrow=TRUE)  # up-regulated
  fold_changes[G_downs, ] <- matrix(rep(c(bioFC, 1, bioFC*batchFC, batchFC), length(G_downs)), ncol=n_group, byrow=TRUE)  # down-regulated
  G_null <- setdiff(1:G, c(G_ups, G_downs))
  fold_changes[G_null, ] <- matrix(rep(c(1, 1, batchFC, batchFC), length(G_null)), ncol=n_group, byrow=TRUE) 
  return(fold_changes)
}

constructSizeMatrix <- function(G, size_vec){
  size_mat <- matrix(rep(size_vec, G), ncol=length(size_vec), byrow=TRUE)
  return(size_mat)
}

constructFCSampleMatrix <- function(fc_batch, batch, group){
  fc_batch_mat <- matrix(NA, nrow=nrow(fc_batch), ncol=length(batch))
  fc_batch_mat[, batch==1&group==0] <- fc_batch[,1]
  fc_batch_mat[, batch==1&group==1] <- fc_batch[,2]
  fc_batch_mat[, batch==2&group==0] <- fc_batch[,3]
  fc_batch_mat[, batch==2&group==1] <- fc_batch[,4]
  return(fc_batch_mat)
}


## library composition difference
constructFCMatrix_Comp <- function(G, FC_group, G_ups, G_downs, bioFC, batchFC){
  fold_changes <- matrix(NA, nrow=G, ncol=length(FC_group))
  
  # randomly split all genes into 2 groups - one increased in batch 2 vs 1, the other decreased
  batch_genes_split <- sample(c(1,2), G, replace=TRUE)
  batch_up_genes <- which(batch_genes_split==1)
  batch_down_genes <- which(batch_genes_split==2)
  
  # batch fold changes
  fold_changes[batch_up_genes, ] <- matrix(rep(c(1, 1, batchFC, batchFC), length(batch_up_genes)), 
                                           nrow=length(batch_up_genes), byrow=TRUE)
  fold_changes[batch_down_genes, ] <- matrix(rep(c(batchFC, batchFC, 1, 1), length(batch_down_genes)), 
                                             nrow=length(batch_down_genes), byrow=TRUE)
  
  # add biological fold changes
  fold_changes[G_ups, FC_group==1] <- fold_changes[G_ups, FC_group==1] * bioFC
  fold_changes[G_downs, FC_group==0] <- fold_changes[G_downs, FC_group==0] * bioFC  
    
  return(fold_changes)
}


estimateMLEDisp <- function(x, group){
  ff <- glm.nb(x~as.factor(group))
  return(1/ff$theta)
}


## simulate outliers
simOutliers <- function(cts_mat, fc_mat, s_mat, batch, group, rptx, outgprob, outalpha){
  genes_w_outliers <- which(sample(c(TRUE, FALSE), nrow(cts_mat), prob=c(outgprob, 1-outgprob), replace=TRUE))
  
  sample_categories <- rep(NA, ncol(cts_mat))
  sample_categories[batch==1 & group==0] <- 1
  sample_categories[batch==1 & group==1] <- 2
  sample_categories[batch==2 & group==0] <- 3
  sample_categories[batch==2 & group==1] <- 4
  
  new_cts <- cts_mat
  for(g in genes_w_outliers){
    outlier_sample <- sample(1:ncol(cts_mat), 1)
    lbd <- qnbinom(outalpha, mu=rptx[g]*fc_mat[g, sample_categories[outlier_sample]], 
                   size=s_mat[g, sample_categories[outlier_sample]], lower.tail=FALSE)
    new_cts[g, outlier_sample] <- sample(lbd:(2*lbd), 1)
  }
  return(new_cts)
}


quantDisp <- function(cts_meanadj, batch, group, DE_ind){
  counts_base_quant <- matrix(NA, nrow=nrow(cts_meanadj), ncol=ncol(cts_meanadj), dimnames=dimnames(cts_meanadj))
  counts_base_quant[, batch==1] <- cts_meanadj[, batch==1]
    
  disp_batch1 <- apply(cts_meanadj[, batch==1], 1, function(x){estimateMLEDisp(x, group=group[batch==1])})
  disp_batch2 <- apply(cts_meanadj[, batch==2], 1, function(y){estimateMLEDisp(y, group=group[batch==2])})
  
  for(i in DE_ind){
    counts_base_quant[i, batch==2&group==0] <- matchDisp(cts_meanadj[i, batch==2&group==0], 
                                                         disp=disp_batch2[i], target.disp=disp_batch1[i])
    counts_base_quant[i, batch==2&group==1] <- matchDisp(cts_meanadj[i, batch==2&group==1], 
                                                         disp=disp_batch2[i], target.disp=disp_batch1[i])
  }
  for(j in setdiff(1:nrow(cts_meanadj), DE_ind)){
    counts_base_quant[j, batch==2] <- matchDisp(cts_meanadj[j, batch==2], disp=disp_batch2[j], target.disp=disp_batch1[j])
  }
  
  return(counts_base_quant)
}


matchDisp <- function(x, disp, target.disp){
  mu <- mean(x)
  new.x <- sapply(x, function(ct){
    p <- pnbinom(ct, mu=mu, size=1/disp)
    if(abs(p-1)<1e-4){return(ct)}
    return(qnbinom(p, mu=mu, size=1/target.disp))
  })
  return(new.x)
}


edgeR_DEpipe <- function(counts_mat, batch, group, include.batch, alpha.unadj, alpha.fdr, covar=NULL){
  cat("DE tool: edgeR\n")
  y <- DGEList(counts=counts_mat)
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
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef=2)
  de_res <- topTags(qlf, n=nrow(counts_mat))$table
  
  de_called <- rownames(de_res)[de_res$PValue < alpha.unadj]
  de_called_fdr <- rownames(de_res)[de_res$FDR < alpha.fdr]
  return(list(unadj=de_called, fdr=de_called_fdr, de_res=de_res, design=design))
}


DESeq2_DEpipe <- function(counts_mat, batch, group, include.batch, alpha.unadj, alpha.fdr, covar=NULL){
  cat("DE tool: DESeq2\n")
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
  
  dds <- DESeqDataSetFromMatrix(countData=counts_mat, colData=col_data, design=design_formula)
  dds <- DESeq(dds, quiet=TRUE)
  de_res <- results(dds, name="Group_1_vs_0")
  
  if(length(which(de_res$pvalue < alpha.unadj))>0){
    de_called <- rownames(de_res)[which(de_res$pvalue < alpha.unadj)]
  }else{de_called <- character(0)}
  if(length(which(de_res$padj < alpha.fdr))>0){
    de_called_fdr <- rownames(de_res)[which(de_res$padj < alpha.fdr)]
  }else{de_called_fdr <- character(0)}
  
  return(list(unadj=de_called, fdr=de_called_fdr, de_res=de_res, 
              design_formula=design_formula, col_data=col_data))
}


currComBat_lm_DEpipe <- function(cts, batch, group, alpha.unadj, alpha.fdr){
  log_counts <- cpm(cts, log=TRUE)  # use logCPM to make data more normal
  adj_counts <- ComBat(log_counts, batch=batch, mod=model.matrix(~as.factor(group)))
  pval_seq <- apply(adj_counts, 1, function(x, group){
    x_norm <- scale(x, center=TRUE, scale=TRUE)
    fit3 <- lm(x_norm ~ as.factor(group))
    return(summary(fit3)$coefficients[2, 4])
  }, group=group)
  padj_seq <- p.adjust(pval_seq, method="fdr")
  
  de_called <- list(unadj=rownames(cts)[pval_seq < alpha.unadj], fdr=rownames(cts)[padj_seq < alpha.fdr], 
                    de_res=data.frame(PValue=pval_seq, FDR=padj_seq), design=model.matrix(~as.factor(group)))
  return(de_called)
}


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


cancelLibsizeEffect <- function(count_matrix){
  lib_sizes <- colSums(count_matrix)
  do.call(cbind, lapply(1:ncol(count_matrix), function(i){
    count_matrix[,i]/lib_sizes[i]
  }))
}
