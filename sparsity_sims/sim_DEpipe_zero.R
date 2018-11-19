rm(list=ls())
demo <- FALSE  # if testing code, set as TRUE; if running simulations, set as FALSE
if(demo){
  setwd("~/Documents/ComBat_seq/sparsity_sims/")
  script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
  source(file.path(script_dir, "sparsity_sims/sim_DEpipe_zero_helpers.R"))
}else{
  setwd("/restricted/projectnb/combat/work/yuqingz/ComBat-seq/sparsity_sims/")
  script_dir <- ".."
  source(file.path(script_dir, "sim_DEpipe_zero_helpers.R"))
}
sapply(c("polyester", "Biostrings", "limma", "edgeR", "DESeq2", "sva", "MASS"), require, character.only=TRUE)
source(file.path(script_dir, "NbZnCombat.R")); source(file.path(script_dir, "helper_seq.R"))
set.seed(123)


####  Parameters
command_args <- commandArgs(trailingOnly=TRUE)
disp_fold_level <- as.numeric(command_args[1])  # dispersion of batch 2 is how many times that of batch 1, 1-10
confounding_level <- as.numeric(command_args[2])  # level of confounding, 0-0.5
#disp_fold_level <- 4; confounding_level <- 0.3
N_total_sample <- 20
coverage <- 1 
factor_exam <- ifelse(demo, "demo", "Sparse")  #command_args[1]  
bio_fold <- 2  #as.numeric(command_args[2])  
batch_fold <- 1.5  #as.numeric(command_args[3])  
size_1 <- 1/0.15  #as.numeric(command_args[4])  # 1/dispersion in batch 1 
size_2 <- 1/(0.15*disp_fold_level)  #as.numeric(command_args[5])   # 1/dispersion in batch 2 
balanced <- FALSE  #as.logical(command_args[7]) 
iterations <- 20 #number of simulations to run
alpha_unadj <- 0.1
alpha_fdr_seq <- seq(from=0, to=0.2, by=0.01)[-1]
exp_name <- paste0("sim", factor_exam, "_N", N_total_sample, "_dispFC", disp_fold_level,
                   "_cnfnd", confounding_level)  #, "_depth", coverage)
exp_name <- gsub('.', '', exp_name, fixed=TRUE)

# for zero-inflation
prop_gene_partition <- as.numeric(c("0.3", "0.3", "0.4"))  #as.numeric(command_args[9:11])
if(sum(prop_gene_partition)!=1){stop("Wrong gene partition probability input: must sum to 1!")}
# this means that: 60% genes - no zeros; 25% genes - random zero fraction (ranging 0%-100%), 
# 15% genes - "zin" genes, zero fraction > 80%
zero.fracs.cutoff <- 0.6   #as.numeric(command_args[12])

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
#ground truth vectors
gene_names <- paste0("gene", 1:length(fasta))
de_ground_truth <- gene_names[de_ground_truth_ind]
true_ups <- gene_names[G_ups]
true_downs <- gene_names[G_downs]
true_nulls <- gene_names[G_nulls]

#for baseline datasets without batch effect
fold_changes_base <- constructFCMatrix_Comp(G=length(fasta), FC_group=c(0,1,0,1), G_ups=G_ups, G_downs=G_downs, bioFC=bio_fold, batchFC=1)
fold_changes_base <- fold_changes_base * sqrt(batch_fold)
size_mat_base <- constructSizeMatrix(G=length(fasta), size_vec=rep(1/mean(c(1/size_1, 1/size_2)), 4))
#for data with batch effect
fold_changes <- constructFCMatrix_Comp(G=length(fasta), FC_group=c(0,1,0,1), G_ups=G_ups, G_downs=G_downs, bioFC=bio_fold, batchFC=batch_fold)
size_mat <- constructSizeMatrix(G=length(fasta), size_vec=c(size_1, size_1, size_2, size_2))



####  Run pipeline
#iter=1; ii=1
for(ii in seq_along(alpha_fdr_seq)){
  alpha_fdr <- alpha_fdr_seq[ii]
  
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
    #load count matrix 
    load(file.path(exp_name, "sim_counts_matrix.rda"))
    cts <- counts_matrix; rm(counts_matrix)
    rownames(cts) <- gene_names
    f_rm <- file.remove(file.path(exp_name, "sim_counts_matrix.rda"))
    
    ## simulate baseline (no batch effect) data - independently using theoretical values for parameters
    simulate_experiment(fasta_file, reads_per_transcript=readspertx, size=size_mat_base, 
                        num_reps=N_samples, fold_changes=fold_changes_base, outdir=exp_name) 
    #remove fasta files to save space
    f_rm <- file.remove(file.path(exp_name, dir(exp_name)[grep(".fasta", dir(exp_name))]))
    if(!all(f_rm)){warning("Something went wrong when deleting fasta files.")}
    #load count matrix 
    load(file.path(exp_name, "sim_counts_matrix.rda"))
    counts_base_indi <- counts_matrix; rm(counts_matrix)
    rownames(counts_base_indi) <- gene_names
    f_rm <- file.remove(file.path(exp_name, "sim_counts_matrix.rda"))
    
    
    ####  Inflate zeros in both datasets
    zin_res <- sim_inflate_zeros(cts, prop_gene_partition); cts <- zin_res$counts
    zin_res_base <- sim_inflate_zeros(counts_base_indi, prop_gene_partition); counts_base_indi <- zin_res_base$counts
    # remove genes with only zeros across all samples
    cts <- cts[apply(cts, 1, function(x){!all(x==0)}), ]
    counts_base_indi <- counts_base_indi[apply(counts_base_indi, 1, function(x){!all(x==0)}), ]
    
    
    ####  DE analysis 
    # On baseline dataset without batch effect - independent baseline 
    de_called01_deseq <- try(DESeq2_DEpipe(counts_base_indi, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr))
    # On counts with batch effect 
    de_called1_deseq <- try(DESeq2_DEpipe(cts, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr))
    # One-step - include batch as covariate
    de_called2_deseq <- try(DESeq2_DEpipe(cts, batch=batch, group=group, include.batch=TRUE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr))
    # Current ComBat + linear model for DE
    #remove genes with only zero values in any batch, due to inflation
    cts_filtered <- filter_genes_zin(cts, batch)
    de_called3 <- try(currComBat_lm_DEpipe(cts_filtered, batch=batch, group=group, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr))
    # On adjusted count - ComBat-seq 
    adj_counts_combatseq <- try(NbZnCombat(cts, batch=batch, group=group, zero.fracs.cutoff=zero.fracs.cutoff))
    if(class(adj_counts_combatseq)=="try-error" | any(is.na(adj_counts_combatseq))){next}
    de_called4_deseq <- try(DESeq2_DEpipe(adj_counts_combatseq, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr))
    
    
    ####  Collect and output results
    DE_objs <- list(BaseIndi.DESeq2=de_called01_deseq,
                    Batch.DESeq2=de_called1_deseq,
                    OneStep.DESeq2=de_called2_deseq,
                    ComBat.lm=de_called3, 
                    ComBatseq.DESeq2=de_called4_deseq)
    if(any(sapply(DE_objs, function(de.obj){class(de.obj)=="try-error"}))){next}
    
    DEgenes_unadj <- lapply(DE_objs, function(de_obj){de_obj$unadj})
    DEgenes_fdr <- lapply(DE_objs, function(de_obj){de_obj$fdr})
    
    cts_lst <- list(counts_base_indi, cts, cts, cts_filtered, adj_counts_combatseq)
    de_ground_truth_lst <- lapply(cts_lst, function(curr_cts){intersect(rownames(curr_cts), de_ground_truth)})
    n_genes_vec <- sapply(cts_lst, nrow)
    rm(cts_lst)
    
    DE_res <- lapply(1:length(DEgenes_unadj), function(i){
      perfStats(called_vec=DEgenes_unadj[[i]], ground_truth_vec=de_ground_truth_lst[[i]], N_genes=n_genes_vec[i])
    })
    names(DE_res) <- names(DEgenes_unadj)
    DE_res <- data.frame(PValue.cutoff=alpha_unadj, do.call(cbind, DE_res))
    
    DE_res_fdr <- lapply(1:length(DEgenes_fdr), function(j){
      perfStats(called_vec=DEgenes_fdr[[j]], ground_truth_vec=de_ground_truth_lst[[j]], N_genes=n_genes_vec[j])
    })
    names(DE_res_fdr) <- names(DEgenes_fdr)
    DE_res_fdr <- data.frame(FDR.cutoff=alpha_fdr, do.call(cbind, DE_res_fdr))
    
    
    ## Write out DE performance
    first.file <- !file.exists(sprintf('fpr_%s.csv', exp_name))
    # for unadjusted p values, write out TPR & FPR
    write.table(DE_res["fpr", ], sprintf('fpr_%s.csv', exp_name),
                append=!first.file, col.names=first.file, row.names=FALSE, sep=",") # type 1 error rate (false positive rate)
    write.table(DE_res["tpr", ], sprintf('tpr_%s.csv', exp_name),
                append=!first.file, col.names=first.file, row.names=FALSE, sep=",") # power (true positive rate)
    # for FDR adjusted values, write out TPR & Precision
    write.table(DE_res_fdr["tpr", ], sprintf('tprADJ_%s.csv', exp_name),
                append=!first.file, col.names=first.file, row.names=FALSE, sep=",") # sensitivity (true positive rate)
    write.table(DE_res_fdr["prec", ], sprintf('precADJ_%s.csv', exp_name),
                append=!first.file, col.names=first.file, row.names=FALSE, sep=",") # precision (1-FDR:false discovery rate)
  }
}
