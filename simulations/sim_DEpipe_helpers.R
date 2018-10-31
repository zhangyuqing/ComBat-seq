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


estimateMLEDisp <- function(x, group){
  ff <- glm.nb(x~as.factor(group))
  return(1/ff$theta)
}


# disp_batch1 <- apply(cts_meanadj[,batch==1], 1, function(x){estimateMLEDisp(x, group=group[batch==1])})
# disp_batch2 <- apply(cts_meanadj[,batch==2], 1, function(y){estimateMLEDisp(y, group=group[batch==2])})
# median(disp_batch2) / median(disp_batch1)
# mean(rowVars(cts[,batch==1])); mean(rowVars(cts[,batch==2]))
# dispersion of batch 2 is still 3 times that of batch 1, wish they are the same, otherwise there still is variance batch effect
# map dispersion of batch 2 to batch 1 without changing the mean
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
#disp_batch1 <- apply(counts_base_quant[,batch==1], 1, function(x){estimateMLEDisp(x, group=group[batch==1])})
#disp_batch2 <- apply(counts_base_quant[,batch==2], 1, function(x){estimateMLEDisp(x, group=group[batch==2])})


matchDisp <- function(x, disp, target.disp){
  mu <- mean(x)
  new.x <- sapply(x, function(ct){
    p <- pnbinom(ct, mu=mu, size=1/disp)
    return(qnbinom(p, mu=mu, size=1/target.disp))
  })
  return(new.x)
}


edgeR_DEpipe <- function(counts_mat, batch, group, include.batch, alpha, covar=NULL){
  y <- DGEList(counts=counts_mat)
  y <- calcNormFactors(y, method="TMM")
  if(include.batch){
    design <- model.matrix(~ as.factor(group) + as.factor(batch))
  }else{
    design <- model.matrix(~as.factor(group))
  }
  if(!is.null(covar)){
    design <- cbind(design, covar)
  }
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef=2)
  de_res <- topTags(qlf, n=nrow(counts_mat))$table
  de_called <- rownames(de_res)[de_res$PValue < alpha]
  return(de_called)
}


DESeq2_DEpipe <- function(counts_mat, batch, group, include.batch, alpha, covar=NULL){
  if(include.batch){
    print("Including batch as covariate")
    col_data <- data.frame(Batch=as.factor(batch), Group=as.factor(group))
    design_formula <- ~Batch+Group
  }else if(!is.null(covar)){
    print("Including surrogate variables or unwanted variation variables")
    col_data <- cbind(model.matrix(~as.factor(group)), covar)
    colnames(col_data)[2:ncol(col_data)] <- c("Group", paste0("Covar", 1:(ncol(col_data)-2)))
    rownames(col_data) <- colnames(counts_mat)
    col_data <- as.data.frame(col_data); col_data$Group <- as.factor(col_data$Group); col_data <- col_data[, -1]
    design_formula <- as.formula(paste("~", paste(colnames(col_data), collapse="+")))
  }else{
    print("Default group as model matrix")
    col_data <- data.frame(Group=as.factor(group))
    design_formula <- ~Group
  }
  dds <- DESeqDataSetFromMatrix(countData=counts_mat, colData=col_data, design=design_formula)
  dds <- DESeq(dds)
  de_res <- results(dds, name="Group_1_vs_0")
  de_called <- rownames(de_res)[de_res$pvalue < alpha]
  return(de_called)
}

