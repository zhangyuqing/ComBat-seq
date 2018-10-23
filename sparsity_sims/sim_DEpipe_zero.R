rm(list=ls())
new_wd <- "~/Google Drive/ComBat_seq/DE_analysis_tmp/"  #"~/yuqingz/ComBat_seq/DE_analysis_ZIN/"
script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"  #".."
setwd(new_wd)
sapply(c("polyester", "Biostrings", "limma", "edgeR", "DESeq2", "sva"), require, character.only=TRUE)
source(file.path(script_dir, "NbZnCombat.R")); source(file.path(script_dir, "helper_seq.R"))   # load NbZnCombat
source(file.path(script_dir, "sparsity_sims/sim_DEpipe_zero_helpers.R"))
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
#factor_exam="FC"; bio_fold=2; batch_fold=2; size_1=50; size_2=10; N_total_sample=20; balanced=TRUE; coverage=5

prop_gene_partition <- as.numeric(c("0.6", "0.25", "0.15"))  #as.numeric(command_args[9:11])
if(sum(prop_gene_partition)!=1){stop("Wrong gene partition probability input: must sum to 1!")}
# this means that: 
# 60% genes - no zeros, 
# 25% genes - random zero fraction (ranging 0%-100%)
# 15% genes - "zin" genes, zero fraction > 80%
zero.fracs.cutoff <- 0.5   #as.numeric(command_args[12])

iterations <- 100 #5  #number of simulations to run
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
if(balanced & N_total_sample==10){
  N_samples <- c(2, 3, 2, 3)
}else if(balanced){
  N_samples <- rep(N_total_sample/4, 4)
}else{
  N_samples <- N_total_sample * c(0.1, 0.4, 0.4, 0.1)  #  20%-80% unbalanced
}

fold_changes <- matrix(NA, nrow=length(fasta), ncol=length(N_samples))  
fold_changes[1:50, ] <- matrix(rep(c(1, bio_fold, batch_fold, bio_fold*batch_fold), 50),
                               ncol=length(N_samples), byrow=TRUE)  # up-regulated
fold_changes[51:100, ] <- matrix(rep(c(bio_fold, 1, bio_fold*batch_fold, batch_fold), 50),
                                 ncol=length(N_samples), byrow=TRUE)  # down-regulated
fold_changes[101:length(fasta), ] <- matrix(rep(c(1, 1, batch_fold, batch_fold), length(fasta)-100),
                                            ncol=length(N_samples), byrow=TRUE)   # "null" genes
de_ground_truth_ind <- 1:100
#N_DE <- length(de_ground_truth_ind)
#N_nonDE <- length(fasta) - N_DE

size_mat <- matrix(rep(c(size_1, size_1, size_2, size_2), length(fasta)),
                   ncol=length(N_samples), byrow=TRUE)



####  Run pipeline
p_zeros_list <- list()
for(iter in 1:iterations){
  cat(paste("\nSimulation", iter, '\n'))

  ####  Simulate datasets
  # use polyester to simulate dataset
  if(dir.exists(exp_name)){unlink(exp_name, recursive=TRUE)}
  simulate_experiment(fasta_file, reads_per_transcript=readspertx, size=size_mat,
                      num_reps=N_samples, fold_changes=fold_changes, outdir=exp_name)

  # load count matrix 
  load(file.path(exp_name, "sim_counts_matrix.rda"))
  rownames(counts_matrix) <- paste0("gene", 1:nrow(counts_matrix))
  de_ground_truth <- rownames(counts_matrix)[de_ground_truth_ind]
  
  # batch and biological vectors
  batch <- c(rep(1, sum(N_samples[1:2])), rep(2, sum(N_samples[3:4])))
  group <- c(rep(0, N_samples[1]), rep(1, N_samples[2]), rep(0, N_samples[3]), rep(1, N_samples[4]))
  
  # inflate zeros in simulated dataset
  zin_res <- sim_inflate_zeros(counts_matrix, prop_gene_partition)
  counts_matrix <- zin_res$counts
  p_zeros_list[[iter]] <- zin_res$p_zero_seq
  
  # remove genes that are flat to all-zero due to inflation
  new_allzero_genes <- rownames(counts_matrix)[apply(counts_matrix, 1, function(x){all(x==0)})]
  if(!all(counts_matrix[new_allzero_genes, ]==0)){stop("Error in looking for new all zero genes!")}
  
  counts_matrix <- counts_matrix[setdiff(rownames(counts_matrix), new_allzero_genes), ]
  de_ground_truth <- setdiff(de_ground_truth, new_allzero_genes)
  N_DE <- length(de_ground_truth)
  N_nonDE <- nrow(counts_matrix) - N_DE
  
  
    
  ####  DE analysis (with edgeR)
  # On original counts with batch effect
  y1 <- DGEList(counts=counts_matrix)
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
  y2 <- DGEList(counts=counts_matrix)
  y2 <- calcNormFactors(y2, method="TMM")
  design2 <- model.matrix(~ as.factor(group) + as.factor(batch))
  y2 <- estimateDisp(y2, design2)
  fit2 <- glmQLFit(y2, design2)
  qlf2 <- glmQLFTest(fit2, coef=2)
  de_res2 <- topTags(qlf2, n=nrow(counts_matrix))$table
  de_called2 <- rownames(de_res2)[de_res2$PValue < alpha]
  
  tpr2 <- length(intersect(de_called2, de_ground_truth)) / N_DE
  fpr2 <- length(setdiff(de_called2, de_ground_truth)) / N_nonDE
    
  
  # On adjusted count - current ComBat + linear model for DE
  log_counts <- cpm(counts_matrix, log=TRUE)  # use logCPM to make data more normal
  adj_counts <- ComBat(log_counts, batch=batch, mod=model.matrix(~as.factor(group)))
  pval_seq <- apply(adj_counts, 1, function(x, group){
    x_norm <- scale(x, center=TRUE, scale=TRUE)
    fit3 <- lm(x_norm ~ as.factor(group))
    return(summary(fit3)$coefficients[2, 4])
  }, group=group)
  de_called3 <- rownames(counts_matrix)[pval_seq < alpha]
  
  tpr3 <- length(intersect(de_called3, de_ground_truth)) / N_DE
  fpr3 <- length(setdiff(de_called3, de_ground_truth)) / N_nonDE
  
  
  # On adjusted count - current ComBat + voom
  # design4 <- model.matrix(~as.factor(group))
  # v <- voom(adj_counts, design=design4)
  # fit4 <- lmFit(v, design4)
  # fit4 <- eBayes(fit4)
  # de_res4 <- topTable(fit4, coef=2, number=nrow(counts_matrix))
  # de_called4 <- rownames(de_res4)[de_res4$P.Value < alpha]
  # 
  # tpr4 <- length(intersect(de_called4, de_ground_truth)) / N_DE
  # fpr4 <- length(setdiff(de_called4, de_ground_truth)) / N_nonDE
  
  
  # On adjusted count - NbZnCombat + edgeR
  #cts=counts_matrix; zin.opt=TRUE; full_mod=TRUE
  adj_counts_combatseq <- NbZnCombat(cts=counts_matrix, batch=batch, group=group, 
                                     zin.opt=TRUE, zero.fracs.cutoff=zero.fracs.cutoff)
  
  y5 <- DGEList(counts=adj_counts_combatseq)
  y5 <- calcNormFactors(y5, method="TMM")
  design5 <- model.matrix(~as.factor(group))
  y5 <- estimateDisp(y5, design5)
  fit5 <- glmQLFit(y5, design5)
  qlf5 <- glmQLFTest(fit5, coef=2)
  de_res5 <- topTags(qlf5, n=nrow(counts_matrix))$table
  de_called5 <- rownames(de_res5)[de_res5$PValue < alpha]
  
  tpr5 <- length(intersect(de_called5, de_ground_truth)) / N_DE
  fpr5 <- length(setdiff(de_called5, de_ground_truth)) / N_nonDE
  
  
  # On adjusted count - NbZnCombat + voom
  # design6 <- model.matrix(~as.factor(group))
  # v2 <- voom(adj_counts_combatseq, design=design6)
  # fit6 <- lmFit(v2, design6)
  # fit6 <- eBayes(fit6)
  # de_res6 <- topTable(fit6, coef=2, number=nrow(counts_matrix))
  # de_called6 <- rownames(de_res6)[de_res6$P.Value < alpha]
  # 
  # tpr6 <- length(intersect(de_called6, de_ground_truth)) / N_DE
  # fpr6 <- length(setdiff(de_called6, de_ground_truth)) / N_nonDE
    
  
  ####  Collect and output results
  # DE_res <- matrix(c(fpr1, tpr1, fpr2, tpr2, fpr3, tpr3, fpr4, tpr4, fpr5, tpr5, fpr6, tpr6), nrow=2)
  # colnames(DE_res) <- c("RawCounts.edgeR", "OneStep.edgeR", "ComBat.lm",
  #                      "ComBat.voom", "NbZnCombat.edgeR", "NbZnCombat.voom")
  DE_res <- matrix(c(fpr1, tpr1, fpr2, tpr2, fpr3, tpr3, fpr5, tpr5), nrow=2)
  colnames(DE_res) <- c("RawCounts", "OneStep", "ComBat.lm", "NbZnCombat")
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

p_zero_mat <- do.call(cbind, p_zeros_list)
save(p_zero_mat, file=sprintf("cache_%s.RData", exp_name))

if(dir.exists(exp_name)){unlink(exp_name, recursive=TRUE)}
