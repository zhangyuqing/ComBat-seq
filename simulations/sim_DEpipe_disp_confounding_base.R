rm(list=ls())
# wdir <- "~/Documents/ComBat_seq/DE_analysis_tmp/"
# script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
wdir <- "~/yuqingz/ComBat_seq/DE_disp_confounding_base/"
script_dir <- ".."
setwd(wdir)
sapply(c("polyester", "Biostrings", "limma", "edgeR", "DESeq2", "sva", "RUVSeq", "MASS"), require, character.only=TRUE)
source(file.path(script_dir, "ComBat_seq.R")); source(file.path(script_dir, "helper_seq.R"))
source(file.path(script_dir, "sim_DEpipe_helpers.R"))
set.seed(123)


####  Parameters
command_args <- commandArgs(trailingOnly=TRUE)
disp_fold_level <- as.numeric(command_args[1])  # dispersion of batch 2 is how many times that of batch 1, 1-10
confounding_level <- as.numeric(command_args[2])  # level of confounding, 0-0.5
#disp_fold_level <- 3; confounding_level <- 0.4
factor_exam <- "DCbase"  #command_args[1]  
bio_fold <- 2  #as.numeric(command_args[2])  
batch_fold <- 1.5  #as.numeric(command_args[3])  
size_1 <- 1/0.15  #as.numeric(command_args[4])  # 1/dispersion in batch 1 
size_2 <- 1/(0.15*disp_fold_level)  #as.numeric(command_args[5])   # 1/dispersion in batch 2 
N_total_sample <- 20  #as.numeric(command_args[6])  
balanced <- FALSE  #as.logical(command_args[7]) 
coverage <- 5 #as.numeric(command_args[8])  

iterations <- 100 #5  #number of simulations to run
alpha <- 0.05
# exp_name <- paste0("sim", factor_exam, "_bio", bio_fold, "_batch", batch_fold, "_sizes", size_1, '_', size_2,
#                    "_N", N_total_sample, ifelse(balanced, "_B", "_U"), "_depth", coverage)
exp_name <- paste0("sim", factor_exam, "_dispFC", disp_fold_level, "_percent", confounding_level)
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

#for baseline datasets without batch effect
fold_changes_base <- constructFCMatrix(G=length(fasta), n_group=4, G_ups=1:50, G_downs=51:100, bioFC=bio_fold, batchFC=1)
size_mat_base <- constructSizeMatrix(G=length(fasta), size_vec=c(size_1, size_1, size_1, size_1))
#for data with batch effect
fold_changes <- constructFCMatrix(G=length(fasta), n_group=4, G_ups=1:50, G_downs=51:100, bioFC=bio_fold, batchFC=batch_fold)
size_mat <- constructSizeMatrix(G=length(fasta), size_vec=c(size_1, size_1, size_2, size_2))
  
# DE
de_ground_truth_ind <- 1:100
N_DE <- length(de_ground_truth_ind)
N_nonDE <- length(fasta) - N_DE


####  Run pipeline
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
  fc_batch <- constructFCMatrix(G=length(fasta), n_group=4, G_ups=1:50, G_downs=51:100, bioFC=1, batchFC=batch_fold)
  fc_batch_mat <- constructFCSampleMatrix(fc_batch, batch=batch, group=group)
  cts_meanadj <- round(cts/fc_batch_mat)
  counts_base_quant <- try(quantDisp(cts_meanadj, batch=batch, group=group, DE_ind=de_ground_truth_ind))
  
  
  ####  DE analysis 
  ## On baseline dataset without batch effect - independent baseline 
  # edgeR
  de_called01 <- edgeR_DEpipe(counts_mat=counts_base_indi, batch=batch, group=group, include.batch=FALSE, alpha=alpha)  
  tpr01 <- length(intersect(de_called01, de_ground_truth)) / N_DE
  fpr01 <- length(setdiff(de_called01, de_ground_truth)) / N_nonDE
  # DESeq2
  de_called01_deseq <- DESeq2_DEpipe(counts_mat=counts_base_indi, batch=batch, group=group, include.batch=FALSE, alpha=alpha)
  tpr01_deseq <- length(intersect(de_called01_deseq, de_ground_truth)) / N_DE
  fpr01_deseq <- length(setdiff(de_called01_deseq, de_ground_truth)) / N_nonDE
  
  ## On baseline dataset without batch effect - quantile baseline 
  if(class(counts_base_quant)!="try-error"){
    # edgeR
    de_called02 <- edgeR_DEpipe(counts_mat=counts_base_quant, batch=batch, group=group, include.batch=FALSE, alpha=alpha)  
    tpr02 <- length(intersect(de_called02, de_ground_truth)) / N_DE
    fpr02 <- length(setdiff(de_called02, de_ground_truth)) / N_nonDE
    # DESeq2
    de_called02_deseq <- DESeq2_DEpipe(counts_mat=counts_base_quant, batch=batch, group=group, include.batch=FALSE, alpha=alpha)
    tpr02_deseq <- length(intersect(de_called02_deseq, de_ground_truth)) / N_DE
    fpr02_deseq <- length(setdiff(de_called02_deseq, de_ground_truth)) / N_nonDE
  }else{tpr02 <- fpr02 <- tpr02_deseq <- fpr02_deseq <- NA}
  
  ## On counts with batch effect 
  # edgeR
  de_called1 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha=alpha)  
  tpr1 <- length(intersect(de_called1, de_ground_truth)) / N_DE
  fpr1 <- length(setdiff(de_called1, de_ground_truth)) / N_nonDE
  # DESeq2
  de_called1_deseq <- DESeq2_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha=alpha)
  tpr1_deseq <- length(intersect(de_called1_deseq, de_ground_truth)) / N_DE
  fpr1_deseq <- length(setdiff(de_called1_deseq, de_ground_truth)) / N_nonDE
  
  ## One-step - include batch as covariate
  # edgeR
  de_called2 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=TRUE, alpha=alpha)  
  tpr2 <- length(intersect(de_called2, de_ground_truth)) / N_DE
  fpr2 <- length(setdiff(de_called2, de_ground_truth)) / N_nonDE
  # DESeq2
  de_called2_deseq <- DESeq2_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=TRUE, alpha=alpha)
  tpr2_deseq <- length(intersect(de_called2_deseq, de_ground_truth)) / N_DE
  fpr2_deseq <- length(setdiff(de_called2_deseq, de_ground_truth)) / N_nonDE
  
  ## Current ComBat + linear model for DE
  log_counts <- cpm(cts, log=TRUE)  # use logCPM to make data more normal
  adj_counts <- ComBat(log_counts, batch=batch, mod=model.matrix(~as.factor(group)))
  pval_seq <- apply(adj_counts, 1, function(x, group){
    x_norm <- scale(x, center=TRUE, scale=TRUE)
    fit3 <- lm(x_norm ~ as.factor(group))
    return(summary(fit3)$coefficients[2, 4])
  }, group=group)
  de_called3 <- rownames(cts)[pval_seq < alpha]
  tpr3 <- length(intersect(de_called3, de_ground_truth)) / N_DE
  fpr3 <- length(setdiff(de_called3, de_ground_truth)) / N_nonDE

  ## On adjusted count - ComBat-seq 
  adj_counts_combatseq <- ComBat_seq(counts=cts, batch=batch, group=group)
  # edgeR
  de_called5 <- edgeR_DEpipe(counts_mat=adj_counts_combatseq, batch=batch, group=group, include.batch=FALSE, alpha=alpha)  
  tpr5 <- length(intersect(de_called5, de_ground_truth)) / N_DE
  fpr5 <- length(setdiff(de_called5, de_ground_truth)) / N_nonDE
  # DESeq2
  de_called5_deseq <- DESeq2_DEpipe(counts_mat=adj_counts_combatseq, batch=batch, group=group, include.batch=FALSE, alpha=alpha)
  tpr5_deseq <- length(intersect(de_called5_deseq, de_ground_truth)) / N_DE
  fpr5_deseq <- length(setdiff(de_called5_deseq, de_ground_truth)) / N_nonDE
  
  ## Compare with RUVseq 
  uvseq <- RUVg(cts, cIdx=tail(1:nrow(cts),n=10), k=1)
  # edgeR
  de_called6 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha=alpha, covar=uvseq$W)  
  tpr6 <- length(intersect(de_called6, de_ground_truth)) / N_DE
  fpr6 <- length(setdiff(de_called6, de_ground_truth)) / N_nonDE
  # DESeq2
  de_called6_deseq <- DESeq2_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha=alpha, covar=uvseq$W)
  tpr6_deseq <- length(intersect(de_called6_deseq, de_ground_truth)) / N_DE
  fpr6_deseq <- length(setdiff(de_called6_deseq, de_ground_truth)) / N_nonDE
  
  ## Compare with SVAseq 
  mod1 <- model.matrix(~as.factor(group)); mod0 <- cbind(mod1[,1])
  svseq <- svaseq(cts, mod1, mod0, n.sv=1)
  # edgeR
  de_called7 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha=alpha, covar=svseq$sv)  
  tpr7 <- length(intersect(de_called7, de_ground_truth)) / N_DE
  fpr7 <- length(setdiff(de_called7, de_ground_truth)) / N_nonDE
  # DESeq2
  de_called7_deseq <- DESeq2_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha=alpha, covar=svseq$sv)
  tpr7_deseq <- length(intersect(de_called7_deseq, de_ground_truth)) / N_DE
  fpr7_deseq <- length(setdiff(de_called7_deseq, de_ground_truth)) / N_nonDE
  
  
  ####  Collect and output results
  DE_res <- matrix(c(fpr01, tpr01, fpr01_deseq, tpr01_deseq,  # baseline - indipendent simulation
                     fpr02, tpr02, fpr02_deseq, tpr02_deseq,  # baseline - quantile fix dispersion
                     fpr1, tpr1, fpr1_deseq, tpr1_deseq, # data with batch effect
                     fpr2, tpr2, fpr2_deseq, tpr2_deseq, # one-step approaches
                     fpr3, tpr3,     # current combat
                     fpr5, tpr5, fpr5_deseq, tpr5_deseq, # ComBat-seq
                     fpr6, tpr6, fpr6_deseq, tpr6_deseq, # RUVseq
                     fpr7, tpr7, fpr7_deseq, tpr7_deseq), # SVAseq
                   nrow=2)
  colnames(DE_res) <- c("BaseIndi.edgeR", "BaseIndi.DESeq2", "BaseQuant.edgeR", "BaseQuant.DESeq2",
                        "Batch.edgeR", "Batch.DESeq2", "OneStep.edgeR", "OneStep.DESeq2", "ComBat.lm", 
                        "ComBatseq.edgeR", "ComBatseq.DESeq2", 
                        "RUVseq.edgeR", "RUVseq.DESeq2", "SVAseq.edgeR", "SVAseq.DESeq2")
  rownames(DE_res) <- c("fpr", "tpr")
  DE_res <- as.data.frame(DE_res)
  
  first.file <- !file.exists(sprintf('fpr_%s.csv', exp_name))
  write.table(DE_res["fpr", ], sprintf('fpr_%s.csv', exp_name), 
              append=!first.file, col.names=first.file, row.names=FALSE, sep=",") # type 1 error rate (false positive rate)
  write.table(DE_res["tpr", ], sprintf('tpr_%s.csv', exp_name), 
              append=!first.file, col.names=first.file, row.names=FALSE, sep=",") # power (true positive rate)
}
