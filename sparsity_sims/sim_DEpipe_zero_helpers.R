####   Functions to inflate zeros manually
sim_inflate_zeros <- function(counts_matrix, prop_gene_partition){
  # partition all genes into the three types
  gpart_res <- sim_random_gene_partition(nrow(counts_matrix), prop_gene_partition)
  
  p_zero_seq <- rep(0, nrow(counts_matrix))
  # do nothing for type 1; for type 2, randomly select p_zero: proportion of samples to set as zeros
  p_zero_seq[gpart_res[[2]]] <- runif(length(gpart_res[[2]]), min=0, max=1)
  # for type 3, randomly select a proportion of zeros (proportion at least 80%)
  p_zero_seq[gpart_res[[3]]] <- runif(length(gpart_res[[3]]), min=0.8, max=1)
  
  for(ii in 1:nrow(counts_matrix)){
    if(p_zero_seq[ii]!=0){
      # randomly select p_zero proportion of samples to set as zero
      zero_samples <- sample(1:ncol(counts_matrix), round(ncol(counts_matrix)*p_zero_seq[ii], 0), replace=FALSE)
      counts_matrix[ii, zero_samples] <- 0
    }
  }
  return(list(counts=counts_matrix, p_zero_seq=p_zero_seq))
}

sim_random_gene_partition <- function(n_genes, prop_gene_partition){
  gene_type_id <- sample(1:3, n_genes, replace=TRUE, prob=prop_gene_partition)
  gpart_res <- lapply(1:3, function(i){which(gene_type_id==i)})
}

filter_genes_zin <- function(cts, batch){
  keep1 <- apply(cts[, batch==1], 1, function(x){!all(x==0)})
  keep2 <- apply(cts[, batch==2], 1, function(x){!all(x==0)})
  keep <- keep1 & keep2
  return(cts[keep, ])
}

####   Functions to perform DE analysis
##  KS-test

##  Wilcox-rank-sum test

##  edgeR

##  DESeq2
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
  
  de_called <- rownames(de_res)[which(de_res$pvalue < alpha.unadj)]
  de_called_fdr <- rownames(de_res)[which(de_res$padj < alpha.fdr)]
  
  return(list(unadj=de_called, fdr=de_called_fdr, de_res=de_res, 
              design_formula=design_formula, col_data=col_data))
}


##  DESeq2 + zinbwave



##  Monocle

##  MAST


##  BPSC

##  SCDE

##  scDD


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
#G=length(fasta); FC_group=c(0,1,0,1); bioFC=bio_fold; batchFC=batch_fold
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
