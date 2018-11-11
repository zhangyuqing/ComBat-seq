rm(list=ls())
demo <- TRUE  # if testing code, set as TRUE; if running simulations, set as FALSE
if(demo){
  setwd("~/Documents/ComBat_seq/DE_analysis_tmp/")
  script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
  source(file.path(script_dir, "simulations/sim_DEpipe_helpers.R"))
}else{
  setwd("~/yuqingz/ComBat_seq/DE_disp_confounding_base/")
  script_dir <- ".."
  source(file.path(script_dir, "sim_DEpipe_helpers.R"))
}
sapply(c("polyester", "Biostrings", "limma", "edgeR", "DESeq2", "sva", "RUVSeq", "MASS"), require, character.only=TRUE)
source(file.path(script_dir, "ComBat_seq.R")); source(file.path(script_dir, "helper_seq.R"))
set.seed(123)


####  Parameters
command_args <- commandArgs(trailingOnly=TRUE)
disp_fold_level <- as.numeric(command_args[1])  # dispersion of batch 2 is how many times that of batch 1, 1-10
confounding_level <- as.numeric(command_args[2])  # level of confounding, 0-0.5
N_total_sample <- as.numeric(command_args[3])  # 20 / 60
#disp_fold_level <- 3; confounding_level <- 0.3; N_total_sample <- 20
factor_exam <- ifelse(demo, "demo", "DCbase")  #command_args[1]  
bio_fold <- 2  #as.numeric(command_args[2])  
batch_fold <- 1.5  #as.numeric(command_args[3])  
size_1 <- 1/0.15  #as.numeric(command_args[4])  # 1/dispersion in batch 1 
size_2 <- 1/(0.15*disp_fold_level)  #as.numeric(command_args[5])   # 1/dispersion in batch 2 
balanced <- FALSE  #as.logical(command_args[7]) 
coverage <- 5 #as.numeric(command_args[8])  

iterations <- 100 #5  #number of simulations to run
alpha <- 0.05
# exp_name <- paste0("sim", factor_exam, "_bio", bio_fold, "_batch", batch_fold, "_sizes", size_1, '_', size_2,
#                    "_N", N_total_sample, ifelse(balanced, "_B", "_U"), "_depth", coverage)
exp_name <- paste0("sim", factor_exam, "_N", N_total_sample, "_dispFC", disp_fold_level, "_percent", confounding_level)
exp_name <- gsub('.', '', exp_name, fixed=TRUE)

# FASTA annotation
read_length <- 100
fasta_file <- system.file('extdata', 'chr22.fa', package='polyester')
fasta <- readDNAStringSet(fasta_file)
# reads per transcript = transcriptlength/readlength * coverage
readspertx <- round(coverage * width(fasta) / read_length)

# study design
if(balanced & N_total_sample==10){
  N_samples <- c(2, 3, 2, 3)
}else if(balanced){
  N_samples <- rep(N_total_sample/4, 4)
}else{
  N_samples <- N_total_sample*(c(confounding_level, 1-confounding_level, 1-confounding_level, confounding_level)/2) 
}
if(sum(N_samples)!=N_total_sample){stop("ERROR in Nsamples!")}
#batch & biological vectors
batch <- c(rep(1, sum(N_samples[1:2])), rep(2, sum(N_samples[3:4])))
group <- c(rep(0, N_samples[1]), rep(1, N_samples[2]), rep(0, N_samples[3]), rep(1, N_samples[4]))

# DE
de_ground_truth_ind <- sample(1:length(fasta), 100, replace=FALSE)
G_ups <- de_ground_truth_ind[1:50]; G_downs <- de_ground_truth_ind[51:100]
G_nulls <- setdiff(1:length(fasta), c(G_ups, G_downs))

#for baseline datasets without batch effect
fold_changes_base <- constructFCMatrix(G=length(fasta), n_group=4, G_ups=G_ups, G_downs=G_downs, bioFC=bio_fold, batchFC=1)
size_mat_base <- constructSizeMatrix(G=length(fasta), size_vec=c(size_1, size_1, size_1, size_1))
#for data with batch effect
fold_changes <- constructFCMatrix(G=length(fasta), n_group=4, G_ups=G_ups, G_downs=G_downs, bioFC=bio_fold, batchFC=batch_fold)
size_mat <- constructSizeMatrix(G=length(fasta), size_vec=c(size_1, size_1, size_2, size_2))
  


####  Run pipeline
predDEgenes_edgeR_unadj <- predDEgenes_edgeR_fdr <- predDE_edgeR_resdfs <- list()
for(iter in 1:iterations){
  cat(paste("\nSimulation", iter, '\n'))
  if(dir.exists(exp_name)){unlink(exp_name, recursive=TRUE)}
  
  ####  Simulate datasets
  ## simulate data with batch effect
  simulate_experiment(fasta_file, reads_per_transcript=readspertx, size=size_mat, 
                      num_reps=N_samples, fold_changes=fold_changes, outdir=exp_name) 
  #remove fasta files to save space
  f_rm <- file.remove(file.path(exp_name, dir(exp_name)[grep(".fasta", dir(exp_name))]))
  if(!all(f_rm)){warning("Something went wrong when deleting fasta files.")}
  #load count matrix and remove fastq files to save space
  load(file.path(exp_name, "sim_counts_matrix.rda"))
  cts <- counts_matrix; rm(counts_matrix)
  rownames(cts) <- paste0("gene", 1:nrow(cts))
  
  #DE ground truth vectors
  de_ground_truth <- rownames(cts)[de_ground_truth_ind]
  true_ups <- rownames(cts)[G_ups]
  true_downs <- rownames(cts)[G_downs]
  true_nulls <- rownames(cts)[G_nulls]
  
  ## simulate baseline (no batch effect) data - independently using theoretical values for parameters
  simulate_experiment(fasta_file, reads_per_transcript=readspertx, size=size_mat_base, 
                      num_reps=N_samples, fold_changes=fold_changes_base, outdir=exp_name) 
  #remove fasta files to save space
  f_rm <- file.remove(file.path(exp_name, dir(exp_name)[grep(".fasta", dir(exp_name))]))
  if(!all(f_rm)){warning("Something went wrong when deleting fasta files.")}
  #load count matrix and remove fastq files to save space
  load(file.path(exp_name, "sim_counts_matrix.rda"))
  counts_base_indi <- counts_matrix; rm(counts_matrix)
  rownames(counts_base_indi) <- paste0("gene", 1:nrow(counts_base_indi))
  
  ## simulate baseline (no batch effect) data - match quantiles
  # fc_batch <- constructFCMatrix(G=length(fasta), n_group=4, G_ups=G_ups, G_downs=G_downs, bioFC=1, batchFC=batch_fold)
  # fc_batch_mat <- constructFCSampleMatrix(fc_batch, batch=batch, group=group)
  # cts_meanadj <- round(cts/fc_batch_mat)
  # counts_base_quant <- try(quantDisp(cts_meanadj, batch=batch, group=group, DE_ind=de_ground_truth_ind))
  
  rm(list=ls())
  load("test.RData")
  
  
  ####  DE analysis 
  ## On baseline dataset without batch effect - independent baseline 
  # edgeR
  de_called01 <- edgeR_DEpipe(counts_mat=counts_base_indi, batch=batch, group=group, include.batch=FALSE, alpha=alpha)  
  perf_stats01 <- perfStats(called_vec=de_called01$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  perf_stats01_fdr <- perfStats(called_vec=de_called01$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
  # DESeq2
  de_called01_deseq <- DESeq2_DEpipe(counts_mat=counts_base_indi, batch=batch, group=group, include.batch=FALSE, alpha=alpha)
  perf_stats01_deseq <- perfStats(called_vec=de_called01_deseq$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  perf_stats01_deseq_fdr <- perfStats(called_vec=de_called01_deseq$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
  
  ## On baseline dataset without batch effect - quantile baseline 
  # if(class(counts_base_quant)=="matrix" & !any(is.infinite(counts_base_quant))){
  #   # edgeR
  #   de_called02 <- edgeR_DEpipe(counts_mat=counts_base_quant, batch=batch, group=group, include.batch=FALSE, alpha=alpha)  
  #   perf_stats02 <- perfStats(called_vec=de_called02$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  #   perf_stats02_fdr <- perfStats(called_vec=de_called02$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  #   
  #   # DESeq2
  #   de_called02_deseq <- DESeq2_DEpipe(counts_mat=counts_base_quant, batch=batch, group=group, include.batch=FALSE, alpha=alpha)
  #   perf_stats02_deseq <- perfStats(called_vec=de_called02_deseq$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  #   perf_stats02_deseq_fdr <- perfStats(called_vec=de_called02_deseq$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  # }else{perf_stats02 <- perf_stats02_fdr <- perf_stats02_deseq <- perf_stats02_deseq_fdr <- rep(NA, 3)}
  
  
  ## On counts with batch effect 
  # edgeR
  de_called1 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha=alpha)  
  perf_stats1 <- perfStats(called_vec=de_called1$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  perf_stats1_fdr <- perfStats(called_vec=de_called1$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
  # DESeq2
  de_called1_deseq <- DESeq2_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha=alpha)
  perf_stats1_deseq <- perfStats(called_vec=de_called1_deseq$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  perf_stats1_deseq_fdr <- perfStats(called_vec=de_called1_deseq$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
  
  ## One-step - include batch as covariate
  # edgeR
  de_called2 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=TRUE, alpha=alpha)  
  perf_stats2 <- perfStats(called_vec=de_called2$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  perf_stats2_fdr <- perfStats(called_vec=de_called2$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
  # DESeq2
  de_called2_deseq <- DESeq2_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=TRUE, alpha=alpha)
  perf_stats2_deseq <- perfStats(called_vec=de_called2_deseq$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  perf_stats2_deseq_fdr <- perfStats(called_vec=de_called2_deseq$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
  
  ## Current ComBat + linear model for DE
  log_counts <- cpm(cts, log=TRUE)  # use logCPM to make data more normal
  adj_counts <- ComBat(log_counts, batch=batch, mod=model.matrix(~as.factor(group)))
  pval_seq <- apply(adj_counts, 1, function(x, group){
    x_norm <- scale(x, center=TRUE, scale=TRUE)
    fit3 <- lm(x_norm ~ as.factor(group))
    return(summary(fit3)$coefficients[2, 4])
  }, group=group)
  padj_seq <- p.adjust(pval_seq, method="fdr")
    
  de_called3 <- rownames(cts)[pval_seq < alpha]
  perf_stats3 <- perfStats(called_vec=de_called3, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  de_called3_fdr <- rownames(cts)[padj_seq < alpha]
  perf_stats3_fdr <- perfStats(called_vec=de_called3_fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
 
  de_called3 <- list(unadj=rownames(cts)[pval_seq < alpha], fdr=rownames(cts)[padj_seq < alpha], 
                     de_res=data.frame(PValue=pval_seq, FDR=padj_seq), design=model.matrix(~as.factor(group)))
  de_called3_deseq <- list(unadj=rownames(cts)[pval_seq < alpha], fdr=rownames(cts)[padj_seq < alpha], 
                           de_res=data.frame(pvalue=pval_seq, padj=padj_seq), design=model.matrix(~as.factor(group)))
  
  
  ## On adjusted count - ComBat-seq 
  adj_counts_combatseq <- ComBat_seq(counts=cts, batch=batch, group=group)
  # edgeR
  de_called5 <- edgeR_DEpipe(counts_mat=adj_counts_combatseq, batch=batch, group=group, include.batch=FALSE, alpha=alpha)  
  perf_stats5 <- perfStats(called_vec=de_called5$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  perf_stats5_fdr <- perfStats(called_vec=de_called5$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
  # DESeq2
  de_called5_deseq <- DESeq2_DEpipe(counts_mat=adj_counts_combatseq, batch=batch, group=group, include.batch=FALSE, alpha=alpha)
  perf_stats5_deseq <- perfStats(called_vec=de_called5_deseq$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  perf_stats5_deseq_fdr <- perfStats(called_vec=de_called5_deseq$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
  
  ## Compare with RUVseq 
  uvseq <- RUVg(cts, cIdx=sample(G_nulls,10,replace=F), k=1)
  # edgeR
  de_called6 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha=alpha, covar=uvseq$W)  
  perf_stats6 <- perfStats(called_vec=de_called6$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  perf_stats6_fdr <- perfStats(called_vec=de_called6$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
  # DESeq2
  de_called6_deseq <- DESeq2_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha=alpha, covar=uvseq$W)
  perf_stats6_deseq <- perfStats(called_vec=de_called6_deseq$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  perf_stats6_deseq_fdr <- perfStats(called_vec=de_called6_deseq$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
  
  ## Compare with SVAseq 
  mod1 <- model.matrix(~as.factor(group)); mod0 <- cbind(mod1[,1])
  svseq <- svaseq(cts, mod1, mod0, n.sv=1)
  # edgeR
  de_called7 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha=alpha, covar=svseq$sv)  
  perf_stats7 <- perfStats(called_vec=de_called7$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  perf_stats7_fdr <- perfStats(called_vec=de_called7$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
  # DESeq2
  de_called7_deseq <- DESeq2_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha=alpha, covar=svseq$sv)
  perf_stats7_deseq <- perfStats(called_vec=de_called7_deseq$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  perf_stats7_deseq_fdr <- perfStats(called_vec=de_called7_deseq$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
  
  
  ####  Collect and output results
  ## Unadjusted P values
  DE_res <- matrix(c(perf_stats01, perf_stats01_deseq, # baseline - indipendent simulation
                     #perf_stats02, perf_stats02_deseq, # baseline - quantile fix dispersion
                     perf_stats1, perf_stats1_deseq,  # data with batch effect
                     perf_stats2, perf_stats2_deseq,  # one-step approaches
                     perf_stats3,   # current combat
                     perf_stats5, perf_stats5_deseq, # ComBat-seq
                     perf_stats6, perf_stats6_deseq,  # RUVseq
                     perf_stats7, perf_stats7_deseq),  # SVAseq
                   nrow=3, byrow=FALSE)
  rownames(DE_res) <- names(perf_stats01)
  ## FDR adjusted Q values
  DE_res_fdr <- matrix(c(perf_stats01_fdr, perf_stats01_deseq_fdr, # baseline - indipendent simulation
                         #perf_stats02_fdr, perf_stats02_deseq_fdr, # baseline - quantile fix dispersion
                         perf_stats1_fdr, perf_stats1_deseq_fdr,  # data with batch effect
                         perf_stats2_fdr, perf_stats2_deseq_fdr,  # one-step approaches
                         perf_stats3_fdr,   # current combat
                         perf_stats5_fdr, perf_stats5_deseq_fdr, # ComBat-seq
                         perf_stats6_fdr, perf_stats6_deseq_fdr,  # RUVseq
                         perf_stats7_fdr, perf_stats7_deseq_fdr),  # SVAseq
                       nrow=3, byrow=FALSE)
  rownames(DE_res_fdr) <- names(perf_stats01_fdr)
  colnames(DE_res) <- colnames(DE_res_fdr) <- c("BaseIndi.edgeR", "BaseIndi.DESeq2", #"BaseQuant.edgeR", "BaseQuant.DESeq2",
                        "Batch.edgeR", "Batch.DESeq2", "OneStep.edgeR", "OneStep.DESeq2", "ComBat.lm", 
                        "ComBatseq.edgeR", "ComBatseq.DESeq2", 
                        "RUVseq.edgeR", "RUVseq.DESeq2", "SVAseq.edgeR", "SVAseq.DESeq2")
  DE_res <- as.data.frame(DE_res); DE_res_fdr <- as.data.frame(DE_res_fdr)
  
  
  first.file <- !file.exists(sprintf('fpr_%s.csv', exp_name))
  # for unadjusted p values, write out TPR & FPR
  write.table(DE_res["fpr", ], sprintf('fpr_%s.csv', exp_name),
              append=!first.file, col.names=first.file, row.names=FALSE, sep=",") # type 1 error rate (false positive rate)
  write.table(DE_res["tpr", ], sprintf('tpr_%s.csv', exp_name),
              append=!first.file, col.names=first.file, row.names=FALSE, sep=",") # power (true positive rate)
  # for FDR adjusted values, write out TPR & Precision
  write.table(DE_res_fdr["tpr", ], sprintf('tprADJ_%s.csv', exp_name),
              append=!first.file, col.names=first.file, row.names=FALSE, sep=",") # power (true positive rate)
  write.table(DE_res_fdr["prec", ], sprintf('precADJ_%s.csv', exp_name),
              append=!first.file, col.names=first.file, row.names=FALSE, sep=",") # precision (1-FDR:false discovery rate)
  
  
  ## Organize and cache DE results
  predDEgenes_edgeR_unadj[[iter]] <- list(base=de_called01$unadj, batch=de_called1$unadj, one.step=de_called2$unadj, 
                                          combat.lm=de_called3$unadj, combatseq=de_called5$unadj, 
                                          ruvseq=de_called6$unadj, svaseq=de_called7$unadj)
  predDEgenes_edgeR_fdr[[iter]] <- list(base=de_called01$fdr, batch=de_called1$fdr, one.step=de_called2$fdr, 
                                        combat.lm=de_called3$fdr, combatseq=de_called5$fdr, 
                                        ruvseq=de_called6$fdr, svaseq=de_called7$fdr)
  predDE_edgeR_resdfs[[iter]] <- list(base=de_called01$de_res, batch=de_called1$de_res, one.step=de_called2$de_res, 
                                      combat.lm=de_called3$de_res, combatseq=de_called5$de_res, 
                                      ruvseq=de_called6$de_res, svaseq=de_called7$de_res)
}

save(predDEgenes_edgeR_unadj, predDEgenes_edgeR_fdr, predDE_edgeR_resdfs, file="DE_results.RData")
