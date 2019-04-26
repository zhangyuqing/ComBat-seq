DESeq2_DEpipe <- function(counts_mat, batch, group, include.batch, alpha.unadj, alpha.fdr, covar_incl=NULL, covar=NULL){
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
    colnames(covar_mod)[-1] <- paste0("cov",1:(ncol(covar_mod)-1))
    col_data <- data.frame(col_data, covar_mod[, -1])
    rownames(col_data) <- colnames(counts_mat)
    design_formula <- as.formula(paste0("~", paste(colnames(col_data), collapse="+")))
  }
  
  dds <- DESeqDataSetFromMatrix(countData=counts_mat, colData=col_data, design=design_formula)
  dds <- DESeq(dds)
  de_res <- results(dds, name="Group_1_vs_0")
  
  de_called <- rownames(de_res)[which(de_res$pvalue < alpha.unadj)]
  de_called_fdr <- rownames(de_res)[which(de_res$padj < alpha.fdr)]
  return(list(unadj=de_called, fdr=de_called_fdr, de_res=de_res, design=design_formula))
}



predLasso <- function(trn_set, tst_set, y_trn, normalize=TRUE, use.refcombat=FALSE, ...){
  if(normalize){
    ## normalize
    trn_set_norm <- t(apply(trn_set, 1, scale, center=TRUE, scale=TRUE))
    tst_set_norm <- t(apply(tst_set, 1, scale, center=TRUE, scale=TRUE))
  }else{
    trn_set_norm <- trn_set
    tst_set_norm <- tst_set
  }
  
  if(use.refcombat){
    batch_cmb <- c(rep(1, ncol(trn_set_norm)), rep(2, ncol(tst_set_norm)))
    cmb_dat <- sva::ComBat(cbind(trn_set_norm, tst_set_norm), batch=batch_cmb, mod=NULL, ref.batch=1)
    tst_set_norm <- cmb_dat[, batch_cmb==2]
  }
  
  ## train model and make predictions
  mod <- glmnet::cv.glmnet(x=t(trn_set_norm), y=as.numeric(as.character(y_trn)), family="binomial", ...)
  
  pred_trn_prob <- as.vector(predict(mod, newx=t(trn_set_norm), s="lambda.1se", type="response"))
  pred_tst_prob <- as.vector(predict(mod, newx=t(tst_set_norm), s="lambda.1se", type="response"))
  pred_trn_class <- as.vector(predict(mod, newx=t(trn_set_norm), s="lambda.1se", type="class"))
  pred_tst_class <- as.vector(predict(mod, newx=t(tst_set_norm), s="lambda.1se", type="class"))
  
  return(list(pred_trn_prob=pred_trn_prob, pred_tst_prob=pred_tst_prob,
              pred_trn_class=pred_trn_class, pred_tst_class=pred_tst_class))
}



evalPerfs <- function(preds, labels){
  rocr_pred <- prediction(preds, as.numeric(as.character(labels)))
  curr_perf <- performance(rocr_pred, "tpr", "fpr")
  plt_df <- data.frame(curr_perf@x.values, curr_perf@y.values)
  colnames(plt_df) <- c("fpr", "tpr")
  return(list(plt_df=plt_df, auc=performance(rocr_pred, "auc")@y.values[[1]]))
}

