rm(list=ls())
demo <- FALSE  # if testing code, set as TRUE; if running simulations, set as FALSE
if(demo){  # local
  setwd("~/Documents/ComBat_seq/DE_analysis_tmp/")
  script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
  source(file.path(script_dir, "simulations/sim_DEpipe_helpers.R"))
}else{  # running on server
  #setwd("~/yuqingz/ComBat_seq/DE_compFDR/")
  setwd("/restricted/projectnb/combat/work/yuqingz/ComBat-seq/DE_libcomp_EBcomplete/")
  script_dir <- ".."
  source(file.path(script_dir, "sim_DEpipe_helpers.R"))
}
sapply(c("polyester", "Biostrings", "limma", "edgeR", "DESeq2", "sva", "RUVSeq", "MASS"), require, character.only=TRUE)
source(file.path(script_dir, "ComBat_seq.R")); source(file.path(script_dir, "helper_seq.R"))
set.seed(123)


####  Parameters
command_args <- commandArgs(trailingOnly=TRUE)
disp_fold_level <- as.numeric(command_args[1])  # dispersion of batch 2 is how many times that of batch 1, 1-5
confounding_level <- as.numeric(command_args[2])  # level of confounding, 0-0.5
N_total_sample <- as.numeric(command_args[3])  # 20 / 60
coverage <- as.numeric(command_args[4])  # 1 / 5 / 10
#disp_fold_level <- 3; confounding_level <- 0.4; N_total_sample <- 20; coverage <- 5; 
factor_exam <- ifelse(demo, "LibcompDemo", "Libcomp")  #command_args[1]  
bio_fold <- 2  #as.numeric(command_args[2])  
batch_fold <- 1.5  #as.numeric(command_args[3])  
size_1 <- 1/0.15  #as.numeric(command_args[4])  # 1/dispersion in batch 1 
size_2 <- 1/(0.15*disp_fold_level)  #as.numeric(command_args[5])   # 1/dispersion in batch 2 
balanced <- FALSE  #as.logical(command_args[7]) 
iterations <- 100 #number of simulations to run
alpha_unadj <- 0.05
alpha_fdr_seq <- seq(from=0, to=0.2, by=0.01)[-1]
exp_name <- paste0("sim", factor_exam, "_N", N_total_sample, "_dispFC", disp_fold_level, 
                   "_cnfnd", confounding_level, "_depth", coverage)
exp_name <- gsub('.', '', exp_name, fixed=TRUE)

# FASTA annotation & baseline read counts
read_length <- 100
fasta_file <- system.file('extdata', 'chr22.fa', package='polyester')
fasta <- readDNAStringSet(fasta_file)
# reads per transcript = transcriptlength/readlength * coverage
readspertx <- round(coverage * width(fasta) / read_length)
gene_names <- paste0("gene", 1:length(fasta))

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



####  Run pipeline
#iter=1; ii=1
collect_DE <- collect_ground_truth <- collect_fold_change <- list()
for(ii in seq_along(alpha_fdr_seq)){
  #collect_DE_objs[[ii]] <- list()
  alpha_fdr <- alpha_fdr_seq[ii]
  
  for(iter in 1:iterations){
    cat(paste("\nFDR cutoff", alpha_fdr, "Simulation", iter, '\n'))
    if(dir.exists(exp_name)){unlink(exp_name, recursive=TRUE)}
    
    ####  Simulate datasets
    ## true DE genes
    de_ground_truth_ind <- sample(1:length(fasta), 100, replace=FALSE)
    G_ups <- de_ground_truth_ind[1:50]; G_downs <- de_ground_truth_ind[51:100]
    G_nulls <- setdiff(1:length(fasta), c(G_ups, G_downs))
    #ground truth vectors
    de_ground_truth <- gene_names[de_ground_truth_ind]
    true_ups <- gene_names[G_ups]
    true_downs <- gene_names[G_downs]
    true_nulls <- gene_names[G_nulls]
    
    ## Fold change matrix and size matrix
    #for baseline datasets without batch effect
    fold_changes_base <- constructFCMatrix_Comp(G=length(fasta), FC_group=c(0,1,0,1), G_ups=G_ups, G_downs=G_downs, bioFC=bio_fold, batchFC=1)
    #fold_changes_base <- fold_changes_base * sqrt(batch_fold)
    size_mat_base <- constructSizeMatrix(G=length(fasta), size_vec=c(size_1, size_1, size_1, size_1))
    #size_mat_base <- constructSizeMatrix(G=length(fasta), size_vec=rep(1/mean(c(1/size_1, 1/size_2)), 4))
    #for data with batch effect
    fold_changes <- constructFCMatrix_Comp(G=length(fasta), FC_group=c(0,1,0,1), G_ups=G_ups, G_downs=G_downs, bioFC=bio_fold, batchFC=batch_fold)
    size_mat <- constructSizeMatrix(G=length(fasta), size_vec=c(size_1, size_1, size_2, size_2))
    
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
    # - reference batch 1
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
    
    
    ####  DE analysis 
    # On baseline dataset without batch effect - independent baseline 
    de_called01 <- edgeR_DEpipe(counts_mat=counts_base_indi, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr)  
    de_called01_deseq <- DESeq2_DEpipe(counts_mat=counts_base_indi, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr)
    # On counts with batch effect 
    de_called1 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr)  
    de_called1_deseq <- DESeq2_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr)
    # One-step - include batch as covariate
    de_called2 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=TRUE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr)  
    de_called2_deseq <- DESeq2_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=TRUE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr)
    # Current ComBat + linear model for DE
    de_called3 <- currComBat_lm_DEpipe(cts=cts, batch=batch, group=group, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr)
    # On adjusted count - ComBat-seq (with EB shrinkage)
    adj_counts_combatseq_ebON <- ComBat_seq(counts=cts, batch=batch, group=group, shrink=TRUE, shrink.disp=TRUE)
    de_called5_ebON <- edgeR_DEpipe(counts_mat=adj_counts_combatseq_ebON, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr)  
    de_called5_deseq_ebON <- DESeq2_DEpipe(counts_mat=adj_counts_combatseq_ebON, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr)
    # On adjusted count - ComBat-seq (shrinkage OFF)
    adj_counts_combatseq_ebOFF <- ComBat_seq(counts=cts, batch=batch, group=group, shrink=FALSE, shrink.disp=FALSE)
    de_called5_ebOFF <- edgeR_DEpipe(counts_mat=adj_counts_combatseq_ebOFF, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr)  
    de_called5_deseq_ebOFF <- DESeq2_DEpipe(counts_mat=adj_counts_combatseq_ebOFF, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr)
    # Compare with RUVseq 
    uvseq <- RUVg(cts, cIdx=sample(G_nulls, 10, replace=F), k=1)
    de_called6 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, covar=uvseq$W, include.batch=FALSE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr)  
    de_called6_deseq <- DESeq2_DEpipe(counts_mat=cts, batch=batch, group=group, covar=uvseq$W, include.batch=FALSE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr)
    # Compare with SVAseq 
    mod1 <- model.matrix(~as.factor(group)); mod0 <- cbind(mod1[,1])
    svseq <- svaseq(cts, mod=mod1, mod0=mod0, n.sv=1); cat("\n")
    de_called7 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, covar=svseq$sv, include.batch=FALSE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr)  
    de_called7_deseq <- DESeq2_DEpipe(counts_mat=cts, batch=batch, group=group, covar=svseq$sv, include.batch=FALSE, alpha.unadj=alpha_unadj, alpha.fdr=alpha_fdr)
    
    
    ####  Collect and output results
    DE_objs <- list(Base.edgeR=de_called01, Base.DESeq2=de_called01_deseq,
                    Batch.edgeR=de_called1, Batch.DESeq2=de_called1_deseq,
                    OneStep.edgeR=de_called2, OneStep.DESeq2=de_called2_deseq,
                    ComBat.lm=de_called3, 
                    ComBatseq.ebON.edgeR=de_called5_ebON, ComBatseq.ebON.DESeq2=de_called5_deseq_ebON,
                    ComBatseq.ebOFF.edgeR=de_called5_ebOFF, ComBatseq.ebOFF.DESeq2=de_called5_deseq_ebOFF,
                    RUVseq.edgeR=de_called6, RUVseq.DESeq2=de_called6_deseq,
                    SVAseq.edgeR=de_called7, SVAseq.DESeq2=de_called7_deseq)
    DEtables <- lapply(DE_objs, function(de_obj){de_obj$de_res})
    
    DEgenes_unadj <- lapply(DE_objs, function(de_obj){de_obj$unadj})
    DEgenes_fdr <- lapply(DE_objs, function(de_obj){de_obj$fdr})

    DE_res <- lapply(DEgenes_unadj, perfStats, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
    DE_res <- data.frame(PValue.cutoff=alpha_unadj, do.call(cbind, DE_res))
    DE_res_fdr <- lapply(DEgenes_fdr, perfStats, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
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

    ## Cache DE results (randomly chosen) for sanity checks
    # if(ii==3 & iter==17){
    #   collect_DE_objs[[ii]][[iter]] <- DE_objs
    # }
    
    # collect_DE[[iter]] <- DEtables
    # collect_ground_truth[[iter]] <- de_ground_truth
    # rownames(fold_changes) <- gene_names
    # collect_fold_change[[iter]] <- fold_changes
    
    rm(de_ground_truth, de_ground_truth_ind, true_downs, true_nulls, true_ups, G_downs, G_nulls, G_ups,
       fold_changes, size_mat, fold_changes_base, size_mat_base, DEtables, DE_objs)
  }
}

#save(collect_DE, collect_ground_truth, collect_fold_change, file=sprintf("DEout_%s.RData", exp_name))
