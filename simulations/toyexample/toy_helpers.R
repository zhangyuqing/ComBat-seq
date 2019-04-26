NBvar <- function(mu, size){
   return(mu + mu^2/size)
}

simCounts <- function(curr_setup){
  y1_ctrl <- matrix(rnbinom(curr_setup$G_DE1 * curr_setup$N_ctrl1, mu=curr_setup$mu1, size=curr_setup$size1), ncol=curr_setup$N_ctrl1)
  y1_case <- matrix(rnbinom(curr_setup$G_DE1 * curr_setup$N_case1, mu=curr_setup$mu1_up, size=curr_setup$size1), ncol=curr_setup$N_case1)
  y1_null <- matrix(rnbinom(curr_setup$G_nonDE1 * curr_setup$N_1, mu=curr_setup$mu1, size=curr_setup$size1), ncol=curr_setup$N_1)
  
  y2_ctrl <- matrix(rnbinom(curr_setup$G_DE1 * curr_setup$N_ctrl2, mu=curr_setup$mu2, size=curr_setup$size2), ncol=curr_setup$N_ctrl2)
  y2_case <- matrix(rnbinom(curr_setup$G_DE1 * curr_setup$N_case2, mu=curr_setup$mu2_up, size=curr_setup$size2), ncol=curr_setup$N_case2)
  y2_null <- matrix(rnbinom(curr_setup$G_nonDE1 * curr_setup$N_2, mu=curr_setup$mu2, size=curr_setup$size2), ncol=curr_setup$N_2)
  
  y3_ctrl <- matrix(rnbinom(curr_setup$G_DE2 * curr_setup$N_ctrl1, mu=curr_setup$mu2, size=curr_setup$size1), ncol=curr_setup$N_ctrl1)
  y3_case <- matrix(rnbinom(curr_setup$G_DE2 * curr_setup$N_case1, mu=curr_setup$mu2_up, size=curr_setup$size1), ncol=curr_setup$N_case1)
  y3_null <- matrix(rnbinom(curr_setup$G_nonDE2 * curr_setup$N_1, mu=curr_setup$mu2, size=curr_setup$size1), ncol=curr_setup$N_1)
  
  y4_ctrl <- matrix(rnbinom(curr_setup$G_DE2 * curr_setup$N_ctrl2, mu=curr_setup$mu1, size=curr_setup$size2), ncol=curr_setup$N_ctrl2)
  y4_case <- matrix(rnbinom(curr_setup$G_DE2 * curr_setup$N_case2, mu=curr_setup$mu1_up, size=curr_setup$size2), ncol=curr_setup$N_case2)
  y4_null <- matrix(rnbinom(curr_setup$G_nonDE2 * curr_setup$N_2, mu=curr_setup$mu1, size=curr_setup$size2), ncol=curr_setup$N_2)
  
  counts <- rbind(cbind(y1_ctrl, y1_case, y2_ctrl, y2_case), cbind(y1_null, y2_null),
                  cbind(y3_ctrl, y3_case, y4_ctrl, y4_case), cbind(y3_null, y4_null))
  rownames(counts) <- paste0("gene", 1:nrow(counts))
  colnames(counts) <- paste0("sample", 1:ncol(counts))
  
  batch <- c(rep("B", ncol(y1_ctrl)+ncol(y1_case)), rep("A", ncol(y2_ctrl)+ncol(y2_case))); table(batch)
  group <- c(rep(0, ncol(y1_ctrl)), rep(1, ncol(y1_case)), rep(0, ncol(y2_ctrl)), rep(1, ncol(y2_case)))
  
  return(list(counts=counts, batch=batch, group=group))
}


#counts_mat=cts; include.batch=FALSE; alpha.unadj=alpha_unadj; alpha.fdr=alpha_fdr
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
