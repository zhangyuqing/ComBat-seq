rm(list=ls())
setwd("~/Documents/ComBat_seq/real_data_app/")
sapply(c("glmnet", "SummarizedExperiment", "sva", "RUVSeq","DESeq2", "ROCR", "ggplot2", "gridExtra", "reshape2", "dplyr"), 
       require, character.only=TRUE)
script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/TB"
# load all data
rds_obj <- readRDS(file.path(data_dir, "combined.rds"))
# source combat-seq functions
source(file.path(script_dir, "ComBat_seq.R"))
source(file.path(script_dir, "helper_seq.R"))
source(file.path(script_dir, "real_data_app/TB_pred_helpers.R"))
set.seed(123)


#### Training data
# take brazil batch 1 & india as training set, africa as test set
batch_filter <- rds_obj$SequencingBatch %in% c("Brazil_1", "India")
group_filter <- rds_obj$Label %in% c("Non-progressor", "Active")
rds_obj_subset <- rds_obj[, batch_filter & group_filter]
# count matrix
cts <- assays(rds_obj_subset)$counts
# batch indicator
batch_org <- as.character(colData(rds_obj_subset)$SequencingBatch)
batch <- rep(1, length(batch_org)); batch[batch_org=="India"] <- 2
# condition indicator
group_org <- as.character(colData(rds_obj_subset)$Label)
group <- rep(0, length(group_org)); group[group_org=="Active"] <- 1
# remove genes with only zero in any of batches
keep1 <- apply(cts[, batch==1], 1, function(x){!all(x==0)})
keep2 <- apply(cts[, batch==2], 1, function(y){!all(y==0)})
cts <- cts[keep1 & keep2, ]  
# get covariate for gender
covar <- rds_obj_subset$Sex
rm(rds_obj_subset)


#### Test data
rds_obj_subset_test <- rds_obj[, (rds_obj$SequencingBatch=="Africa") & (rds_obj$Label %in% c("Non-progressor", "Active"))]
cts_test <- assays(rds_obj_subset_test)$counts
group_test_org <- as.character(colData(rds_obj_subset_test)$Label)
group_test <- rep(0, length(group_test_org)); group_test[group_test_org=="Active"] <- 1
covar_test <- rds_obj_subset_test$Sex
# take a subset to make test data balanced
cts_test <- cts_test[, 89:120]
group_test <- group_test[89:120]
covar_test <- covar_test[89:120]


#### Select highly variable genes in training data
genes_sel_names <- names(sort(apply(cts,1,var), decreasing=TRUE))[1:1001]
genes_sel_names <- setdiff(genes_sel_names, "HBA2")
cts <- cts[genes_sel_names, ]
cts_test <- cts_test[genes_sel_names, ]


#### Batch correction
## ComBat-seq
covar_mod <- model.matrix(~covar)
start_time <- Sys.time()
cts_adj <- ComBat_seq(counts=cts, batch=batch, group=group, gene.subset.n=NULL, Cpp=FALSE, covar_mod=covar_mod)
end_time <- Sys.time()
print(end_time - start_time)

## Original ComBat
cts_adjori <- ComBat(cts, batch=batch, mod=model.matrix(~group+covar))


#### Training prediction model
predLasso <- function(trn_set, tst_set, y_trn, covmod_trn=NULL, covmod_tst=NULL, ...){
  ## normalize
  trn_set_norm <- t(apply(trn_set, 1, scale, center=TRUE, scale=TRUE))
  tst_set_norm <- t(apply(tst_set, 1, scale, center=TRUE, scale=TRUE))
  # trn_set_norm <- trn_set
  # tst_set_norm <- tst_set
  # bind covariates
  if(!is.null(covmod_trn)){
    trn_set_norm <- rbind(t(covmod_trn), trn_set_norm)
    tst_set_norm <- rbind(t(covmod_tst), tst_set_norm)
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

res_unadj <- predLasso(trn_set=cts, tst_set=cts_test, y_trn=group)
res_combat <- predLasso(trn_set=cts_adjori, tst_set=cts_test, y_trn=group)
res_combatseq <- predLasso(trn_set=cts_adj, tst_set=cts_test, y_trn=group)


predres <- list(Unadjusted=res_unadj, OriginalComBat=res_combat, ComBatSeq=res_combatseq)
pred_trnscores <- lapply(predres, function(x){x$pred_trn_prob})
pred_tstscores <- lapply(predres, function(x){x$pred_tst_prob})

## Evaluate and visualize performances
perf_res <- lapply(pred_tstscores, function(preds, labels){
  rocr_pred <- prediction(preds, as.numeric(as.character(labels)))
  curr_perf <- performance(rocr_pred, "tpr", "fpr")
  plt_df <- data.frame(curr_perf@x.values, curr_perf@y.values)
  colnames(plt_df) <- c("fpr", "tpr")
  return(list(plt_df=plt_df, auc=performance(rocr_pred, "auc")@y.values[[1]]))
}, labels=group_test)

plt_lst <- list()
for(i in 1:length(perf_res)){
  plt_lst[[i]] <- ggplot(perf_res[[i]]$plt_df, aes(x=fpr, y=tpr)) +
    geom_line() +
    labs(title=sprintf("%s, AUC = %s", names(predres)[i], round(perf_res[[i]]$auc,3)),
         x="False positive rate", y="True positive rate") +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey")
}
png("TB_pred_ROCs.png", width=13, height=5, units="in", res=300)
grid.arrange(plt_lst[[1]], plt_lst[[2]], plt_lst[[3]], ncol=3)
dev.off()


#### Adding RUV-Seq & SVA0Seq 
de_called1 <- DESeq2_DEpipe(cts, batch=batch, group=group, include.batch=FALSE, alpha.unadj=1, alpha.fdr=1, covar_incl=covar)  
de_called1$de_res <- data.frame(genes=rownames(de_called1$de_res), de_called1$de_res)
de_called1$de_res <- arrange(de_called1$de_res, padj)
## RUV-seq
emps <- as.character(tail(de_called1$de_res$genes, n=10))
uvseq <- RUVg(cts, cIdx=emps, k=1)
# SVAseq 
mod1 <- model.matrix(~as.factor(group) + covar)
mod0 <- model.matrix(~covar)
svseq <- svaseq(cts, mod1, mod0, n.sv=1);cat("\n")

res_ruvseq <- predLasso(trn_set=uvseq$normalizedCounts, tst_set=cts_test, y_trn=group)
#res_ruvseq <- predLasso(trn_set=cts_adj, tst_set=cts_test, y_trn=group, cov=uvseq$W)
res_svaseq <- predLasso(trn_set=cts_adj, tst_set=cts_test, y_trn=group, cov=svseq$sv)
predres_full <- list(Unadjusted=res_unadj, OriginalComBat=res_combat, ComBatSeq=res_combatseq, RUVSeq=res_ruvseq)#, SVASeq=res_svaseq)
pred_tstscores_full <- lapply(predres_full, function(x){x$pred_tst_prob})

## evaluate and visualize performances
perf_res_full <- lapply(pred_tstscores_full, function(preds, labels){
  rocr_pred <- prediction(preds, as.numeric(as.character(labels)))
  curr_perf <- performance(rocr_pred, "tpr", "fpr")
  plt_df <- data.frame(curr_perf@x.values, curr_perf@y.values)
  colnames(plt_df) <- c("fpr", "tpr")
  return(list(plt_df=plt_df, auc=performance(rocr_pred, "auc")@y.values[[1]]))
}, labels=group_test)

plt_lst_full <- list()
for(i in 1:length(perf_res_full)){
  plt_lst_full[[i]] <- ggplot(perf_res_full[[i]]$plt_df, aes(x=fpr, y=tpr)) +
    geom_line() +
    labs(title=sprintf("%s, AUC = %s", names(predres_full)[i], round(perf_res_full[[i]]$auc,3)),
         x="False positive rate", y="True positive rate") +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey")
}
png("TB_pred_ROCs_full_normdata.png", width=13, height=9, units="in", res=300)
grid.arrange(plt_lst_full[[1]], plt_lst_full[[2]], 
             plt_lst_full[[3]], plt_lst_full[[4]], ncol=2)
dev.off()
