rm(list=ls())
#setwd("~/Google Drive/ComBat_seq/DE_analysis/")
setwd("~/yuqingz/ComBat_seq/DE_analysis/")
sapply(c("ggplot2", "reshape2", "gridExtra", "dendextend", "edgeR", "DESeq2", "polyester", "Biostrings"), 
       require, character.only=TRUE)
#script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
script_dir <- ".."
source(file.path(script_dir, "ComBat_seq.R")); source(file.path(script_dir, "helper_seq.R"))
set.seed(123)


####  Parameters
command_args <- commandArgs(trailingOnly=TRUE)
factor_exam <- command_args[1]  # "FC", "Disp", "Design", "Nsample", "BaseCts"  #the factor in examination in this sim
bio_fold <- as.numeric(command_args[2])  #2  #fold change for biological condition
batch_fold <- as.numeric(command_args[3])  #3  #fold change for batch effect - assuming that batch behaves as fold change on counts
size_1 <- as.numeric(command_args[4])   # 1/dispersion in batch 1 ("Bernard")
size_2 <- as.numeric(command_args[5])   # 1/dispersion in batch 2 ("Arnold")
N_total_sample <- as.numeric(command_args[6])  #20  #total number of samples in the study
balanced <- as.logical(command_args[7])  #TRUE  #logical, if TRUE balanced design
coverage <- as.numeric(command_args[8])  #20  
#factor_exam="FC"; bio_fold=1.5; batch_fold=2; size_1=10; size_2=1000; N_total_sample=20; balanced=FALSE; coverage=1
  
iterations <- 100  #number of simulations to run
alpha <- 0.05
exp_name <- paste0("sim", factor_exam, "_bio", bio_fold, "_batch", batch_fold, "_sizes", size_1, '_', size_2,
                   "_N", N_total_sample, ifelse(balanced, "_B", "_U"), "_depth", coverage)
exp_name <- gsub('.', '', exp_name, fixed=TRUE)

# FASTA annotation
read_length <- 100
fasta_file <- system.file('extdata', 'chr22.fa', package='polyester')
fasta <- readDNAStringSet(fasta_file)
# subset the FASTA file to first N_genes transcripts
#writeXStringSet(small_fasta, 'chr22_small.fa')
# reads per transcript = transcriptlength/readlength * coverage
readspertx <- round(coverage * width(fasta) / read_length)

# study design
if(balanced){
  N_samples <- rep(N_total_sample/4, 4)
}else{
  N_samples <- N_total_sample * c(0.1, 0.4, 0.4, 0.1)  #  20%-80% unbalanced
}

fold_changes <- matrix(NA, nrow=length(fasta), ncol=length(N_samples))  
# fold_changes[1:50, ] <- matrix(rep(c(1,bio_fold,batch_fold,bio_fold*batch_fold),50), ncol=length(N_samples), byrow=TRUE)
# fold_changes[51:450, ] <- matrix(rep(c(1,1,batch_fold,batch_fold),400), ncol=length(N_samples), byrow=TRUE)
# fold_changes[451:500, ] <- matrix(rep(c(batch_fold,bio_fold*batch_fold,1,bio_fold),50), ncol=length(N_samples), byrow=TRUE)
# fold_changes[501:length(fasta), ] <- matrix(rep(c(batch_fold,batch_fold,1,1),length(fasta)-500), ncol=length(N_samples), byrow=TRUE)
fold_changes[1:50, ] <- matrix(rep(c(1, bio_fold, batch_fold, bio_fold*batch_fold), 50),
                                ncol=length(N_samples), byrow=TRUE)  # up-regulated
fold_changes[51:100, ] <- matrix(rep(c(bio_fold, 1, bio_fold*batch_fold, batch_fold), 50),
                               ncol=length(N_samples), byrow=TRUE)  # down-regulated
fold_changes[101:length(fasta), ] <- matrix(rep(c(1, 1, batch_fold, batch_fold), length(fasta)-100),
                                            ncol=length(N_samples), byrow=TRUE)   # "null" genes

size_mat <- matrix(rep(c(size_1, size_1, size_2, size_2), length(fasta)),
                   ncol=length(N_samples), byrow=TRUE)

#de_ground_truth_ind <- c(1:50, 451:500)  # true differentially expressed genes
de_ground_truth_ind <- 1:100
N_DE <- length(de_ground_truth_ind)
N_nonDE <- length(fasta) - N_DE


####  Run pipeline
for(iter in 1:iterations){
  cat(paste("Simulation", iter, '\n'))
  
  ####  Simulate datasets
  if(dir.exists(exp_name)){unlink(exp_name, recursive=TRUE)}
  simulate_experiment(fasta_file, reads_per_transcript=readspertx, size=size_mat, 
                      num_reps=N_samples, fold_changes=fold_changes, outdir=exp_name) 
  
  # load count matrix and remove fastq files to save space
  load(file.path(exp_name, "sim_counts_matrix.rda"))
  de_ground_truth <- rownames(counts_matrix)[de_ground_truth_ind]
  
  # batch and biological vectors
  batch <- c(rep(1, sum(N_samples[1:2])), rep(2, sum(N_samples[3:4])))
  group <- c(rep(0, N_samples[1]), rep(1, N_samples[2]), rep(0, N_samples[3]), rep(1, N_samples[4]))
  
  
  ####  DE analysis (with edgeR)
  # On original counts with batch effect
  y1 <- DGEList(counts=counts_matrix, group=as.factor(group))
  y1 <- calcNormFactors(y1, method="TMM")
  design <- model.matrix(~as.factor(group))
  y1 <- estimateDisp(y1, design)
  fit1 <- glmQLFit(y1, design)
  qlf1 <- glmQLFTest(fit1, coef=2)
  de_res1 <- topTags(qlf1, n=nrow(counts_matrix))$table
  de_called1 <- rownames(de_res1)[de_res1$PValue < alpha]
  
  tpr1 <- length(intersect(de_called1, de_ground_truth)) / N_DE
  fpr1 <- length(setdiff(de_called1, de_ground_truth)) / N_nonDE

    
  # On original count - include batch as covariate
  y2 <- DGEList(counts=counts_matrix, group=as.factor(group))
  y2 <- calcNormFactors(y2, method="TMM")
  design2 <- model.matrix(~ as.factor(group) + as.factor(batch))
  y2 <- estimateDisp(y2, design2)
  fit2 <- glmQLFit(y2, design2)
  qlf2 <- glmQLFTest(fit2, coef=2)
  de_res2 <- topTags(qlf2, n=nrow(counts_matrix))$table
  de_called2 <- rownames(de_res2)[de_res2$PValue < alpha]
  
  tpr2 <- length(intersect(de_called2, de_ground_truth)) / N_DE
  fpr2 <- length(setdiff(de_called2, de_ground_truth)) / N_nonDE
    
    
  # On adjusted count - ComBat-seq
  adj_counts_combatseq <- ComBat_seq(counts=counts_matrix, batch=batch, group=group)
  
  y3 <- DGEList(counts=adj_counts_combatseq, group=as.factor(group))
  y3 <- calcNormFactors(y3, method="TMM")
  design3 <- model.matrix(~as.factor(group))
  y3 <- estimateDisp(y3, design3)
  fit3 <- glmQLFit(y3, design3)
  qlf3 <- glmQLFTest(fit3, coef=2)
  de_res3 <- topTags(qlf3, n=nrow(counts_matrix))$table
  de_called3 <- rownames(de_res3)[de_res3$PValue < alpha]
  
  tpr3 <- length(intersect(de_called3, de_ground_truth)) / N_DE
  fpr3 <- length(setdiff(de_called3, de_ground_truth)) / N_nonDE
    
  
  ####  Collect and output results
  DE_res <- matrix(c(fpr1, tpr1, fpr2, tpr2, fpr3, tpr3), nrow=2)
  colnames(DE_res) <- c("Raw counts", "One-step", "ComBat-seq")
  rownames(DE_res) <- c("fpr", "tpr")
  DE_res <- as.data.frame(DE_res)
  
  first.file <- !file.exists(sprintf('fpr_%s.csv', exp_name))
  # type 1 error rate (false positive rate)
  write.table(DE_res["fpr", ], sprintf('fpr_%s.csv', exp_name), 
              append=!first.file, col.names=first.file, row.names=FALSE, sep=",")
  # power (true positive rate)
  write.table(DE_res["tpr", ], sprintf('tpr_%s.csv', exp_name), 
              append=!first.file, col.names=first.file, row.names=FALSE, sep=",")
  
  rm(counts_matrix)
}

if(dir.exists(exp_name)){unlink(exp_name, recursive=TRUE)}
