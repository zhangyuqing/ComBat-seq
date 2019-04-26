rm(list=ls()); demo <- FALSE
if(demo){
  setwd("~/Documents/ComBat_seq/real_data_app/")
  script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
  data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/TB"
}else{
  setwd("/restricted/projectnb/combat/work/yuqingz/ComBat-seq/real_data_app/")
  script_dir <- "../"
  data_dir <- "."
}
if(!dir.exists("./TB_pred_bootstrap_shrinkOff_subset")){dir.create("./TB_pred_bootstrap_shrinkOff_subset")}
sapply(c("glmnet", "SummarizedExperiment", "sva", "RUVSeq","DESeq2", "ROCR", "ggplot2", "gridExtra", "reshape2", "dplyr"), 
       require, character.only=TRUE)

rds_obj <- readRDS(file.path(data_dir, "combined.rds"))
source(file.path(script_dir, "ComBat_seq.R"))
source(file.path(script_dir, "helper_seq.R"))
source(file.path(script_dir, "real_data_app/TB_pred_helpers.R"))
set.seed(123)


####  Parameters  ####
command_args <- commandArgs(trailingOnly=TRUE)
test_sel_ind <- 89:120  # subsetting test set (Africa) to make it balanced
n_highvar_genes <- 500  # number of highly variable genes to use in feature reduction
fdr_cutoff <- 0.95  # control genes FDR cutoff in RUVSeq
B <- 1000  # number of bootstrap samples to run
norm_data <- as.logical(command_args[1]) #FALSE  # whether to z-score scaling features before training
use_ref_combat <- as.logical(command_args[2])  #FALSE  # whether to use ref combat to adjust 


####  Process original data into training & test  ####
## training
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

## test
rds_obj_subset_test <- rds_obj[, (rds_obj$SequencingBatch=="Africa") & (rds_obj$Label %in% c("Non-progressor", "Active"))]
cts_test <- assays(rds_obj_subset_test)$counts
group_test_org <- as.character(colData(rds_obj_subset_test)$Label)
group_test <- rep(0, length(group_test_org)); group_test[group_test_org=="Active"] <- 1
covar_test <- rds_obj_subset_test$Sex
# take a subset to make test data balanced
cts_test <- cts_test[, test_sel_ind]
group_test <- group_test[test_sel_ind]
covar_test <- covar_test[test_sel_ind]

## feature reduction - select highly variable genes in training data
genes_sel_names <- names(sort(apply(cts,1,var), decreasing=TRUE))[1:(n_highvar_genes+1)]
genes_sel_names <- setdiff(genes_sel_names, "HBA2")
cts <- cts[genes_sel_names, ]
cts_test <- cts_test[genes_sel_names, ]


####  Batch correction  ####
# ComBat-seq
covar_mod <- model.matrix(~covar)
start_time <- Sys.time()
cts_combatseq <- ComBat_seq(counts=cts, batch=batch, group=group, covar_mod=covar_mod, 
                            shrink=FALSE, shrink.disp=FALSE)  #, gene.subset.n=NULL, Cpp=FALSE)
end_time <- Sys.time()
print(end_time - start_time)

## Original ComBat
cts_combat <- ComBat(cts, batch=batch, mod=model.matrix(~group+covar))

## RUV-seq
de_called <- DESeq2_DEpipe(cts, batch=batch, group=group, include.batch=FALSE, alpha.unadj=1, alpha.fdr=1, covar_incl=covar)  
de_called$de_res <- data.frame(genes=rownames(de_called$de_res), de_called$de_res)
emps <- as.character(de_called$de_res$genes[which(de_called$de_res$padj > fdr_cutoff)])
uvseq <- RUVg(cts, cIdx=emps, k=1)
cts_ruvseq <- uvseq$normalizedCounts

## SVAseq 
mod1 <- model.matrix(~as.factor(group) + covar)
mod0 <- model.matrix(~covar)
svseq <- svaseq(cts, mod1, mod0, n.sv=1);cat("\n")
mod_tst <- model.matrix(~covar_test)
svseq_tst <- svaseq(cts_test, mod_tst, n.sv=1);cat("\n")
cts_svaseq <- rbind(t(as.matrix(svseq$sv)), cts)
cts_test_svaseq <- rbind(t(as.matrix(svseq_tst$sv)), cts_test)


####  Training prediction models  ####
res_unadj <- predLasso(trn_set=cts, tst_set=cts_test, y_trn=group, normalize=norm_data, use.refcombat=use_ref_combat)
res_combat <- predLasso(trn_set=cts_combat, tst_set=cts_test, y_trn=group, normalize=norm_data, use.refcombat=use_ref_combat)
res_combatseq <- predLasso(trn_set=cts_combatseq, tst_set=cts_test, y_trn=group, normalize=norm_data, use.refcombat=use_ref_combat)
res_ruvseq <- predLasso(trn_set=cts_ruvseq, tst_set=cts_test, y_trn=group, normalize=norm_data, use.refcombat=use_ref_combat)
#res_ruvseq <- predLasso(trn_set=cts_adj, tst_set=cts_test, y_trn=group, cov=uvseq$W, normalize=norm_data, use.refcombat=use_ref_combat)
res_svaseq <- predLasso(trn_set=cts_svaseq, tst_set=cts_test_svaseq, y_trn=group, normalize=norm_data, use.refcombat=use_ref_combat)

# collect prediction scores
predres <- list(Unadjusted=res_unadj, OriginalComBat=res_combat, ComBatSeq=res_combatseq, 
                RUVSeq=res_ruvseq, SVASeq=res_svaseq)
pred_trnscores <- lapply(predres, function(x){x$pred_trn_prob})
pred_tstscores <- lapply(predres, function(x){x$pred_tst_prob})

# evaluate performances
perf_res <- lapply(pred_tstscores, evalPerfs, labels=group_test)

# visualize performances
plt_lst <- list()
for(i in 1:length(perf_res)){
  plt_lst[[i]] <- ggplot(perf_res[[i]]$plt_df, aes(x=fpr, y=tpr)) +
    geom_line() +
    labs(title=sprintf("%s, AUC = %s", names(predres)[i], round(perf_res[[i]]$auc,3)),
         x="False positive rate", y="True positive rate") +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey")
}
png(sprintf("./TB_pred_bootstrap_shrinkOff_subset/TB_pred_ROCs%s%s.png", 
            ifelse(norm_data, "_normdata", ""), ifelse(use_ref_combat, "_refcombat", "")), 
    width=13, height=10, units="in", res=300)
grid.arrange(plt_lst[[1]], plt_lst[[2]], plt_lst[[3]], 
             plt_lst[[4]], plt_lst[[5]], nrow=2, ncol=3)
dev.off()

# save performance
saveRDS(perf_res, sprintf("./TB_pred_bootstrap_shrinkOff_subset/TB_pred_perfs%s%s.rds", 
                          ifelse(norm_data, "_normdata", ""), ifelse(use_ref_combat, "_refcombat", "")))


####  Bootstrap  ####
perf_boot_auc <- list()
for(iter in 1:B){
  print(paste("Bootstrap:", iter))
  
  # resample with replacement
  unadj_boot_ind <- sample(1:ncol(cts), ncol(cts), replace=TRUE)
  combat_boot_ind <- sample(1:ncol(cts_combat), ncol(cts_combat), replace=TRUE)
  combatseq_boot_ind <- sample(1:ncol(cts_combatseq), ncol(cts_combatseq), replace=TRUE)
  ruvseq_boot_ind <- sample(1:ncol(cts_ruvseq), ncol(cts_ruvseq), replace=TRUE)
  svaseq_boot_ind <- sample(1:ncol(cts_svaseq), ncol(cts_svaseq), replace=TRUE)
  
  boot_ind_lst <- list(Unadjusted=unadj_boot_ind, OriginalComBat=combat_boot_ind, ComBatSeq=combatseq_boot_ind, 
                            RUVSeq=ruvseq_boot_ind, SVASeq=svaseq_boot_ind)
  # train
  boot_unadj <- predLasso(trn_set=cts[, boot_ind_lst$Unadjusted], y_trn=group[boot_ind_lst$Unadjusted], 
                          tst_set=cts_test, normalize=norm_data)
  boot_combat <- predLasso(trn_set=cts_combat[, boot_ind_lst$OriginalComBat], y_trn=group[boot_ind_lst$OriginalComBat], 
                           tst_set=cts_test, normalize=norm_data)
  boot_combatseq <- predLasso(trn_set=cts_combatseq[, boot_ind_lst$ComBatSeq], y_trn=group[boot_ind_lst$ComBatSeq], 
                              tst_set=cts_test, normalize=norm_data)
  boot_ruvseq <- predLasso(trn_set=cts_ruvseq[, boot_ind_lst$RUVSeq], y_trn=group[boot_ind_lst$RUVSeq], 
                           tst_set=cts_test, normalize=norm_data)
  boot_svaseq <- predLasso(trn_set=cts_svaseq[, boot_ind_lst$SVASeq], y_trn=group[boot_ind_lst$SVASeq], 
                           tst_set=cts_test_svaseq, normalize=norm_data)
  
  # collect prediction scores
  predres_boot <- list(Unadjusted=boot_unadj, OriginalComBat=boot_combat, ComBatSeq=boot_combatseq, 
                       RUVSeq=boot_ruvseq, SVASeq=boot_svaseq)
  pred_trnscores_boot <- lapply(predres_boot, function(x){x$pred_trn_prob})
  pred_tstscores_boot <- lapply(predres_boot, function(x){x$pred_tst_prob})
  
  # evaluate performances
  perf_boot <- lapply(pred_tstscores_boot, evalPerfs, labels=group_test)
  
  # output & save bootstrap performance
  perf_boot_auc[[iter]] <- sapply(perf_boot, function(x){x$auc})
  
  first_file <- !file.exists(sprintf("./TB_pred_bootstrap_shrinkOff_subset/TB_pred_perfs_bootstrap%s%s.csv", 
                                     ifelse(norm_data, "_normdata", ""), ifelse(use_ref_combat, "_refcombat", "")))
  write.table(data.frame(t(perf_boot_auc[[iter]])), 
              sprintf("./TB_pred_bootstrap_shrinkOff_subset/TB_pred_perfs_bootstrap%s%s.csv", 
                      ifelse(norm_data, "_normdata", ""), ifelse(use_ref_combat, "_refcombat", "")),
              append=!first_file, col.names=first_file, row.names=FALSE, sep=",")
}

# visualize the results
perf_boot_auc_df <- do.call(rbind, perf_boot_auc)
perf_auc <- sapply(perf_res, function(x){x$auc})

perf_boot_auc_mlt <- melt(data.frame(perf_boot_auc_df), variable.name="Method")
perf_boot_auc_mlt$Ori <- NA
for(i in 1:length(perf_auc)){perf_boot_auc_mlt$Ori[perf_boot_auc_mlt$Method==names(perf_auc)[i]] <- perf_auc[i]}
  
png(sprintf("./TB_pred_bootstrap_shrinkOff_subset/TB_pred_bootstrap%s%s.png", 
            ifelse(norm_data, "_normdata", ""), ifelse(use_ref_combat, "_refcombat", "")), 
    width=6, height=4, units="in", res=300)
ggplot(perf_boot_auc_mlt, aes(y=value)) +
  geom_boxplot(width=0.02) +
  geom_hline(aes(yintercept=Ori, group=Method), colour='red') + 
  #scale_fill_manual(values=c("white","#999999", "#E69F00", "#56B4E9")) +
  facet_wrap(~Method, nrow=1) + #, ncol=3) +
  labs(y="Bootstrap AUC") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()

