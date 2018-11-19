rm(list=ls())
demo <- TRUE  # if testing code, set as TRUE; if running simulations, set as FALSE
if(demo){
  setwd("~/Documents/ComBat_seq/sparsity_sims/")
  script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
  data_dir <- "~/Google Drive/ComBat_seq/real_data_example/SingleCell/Tung2017Batch_ScientificReports"
  source(file.path(script_dir, "sparsity_sims/sim_DEpipe_zero_helpers.R"))
}else{
  setwd("~/yuqingz/ComBat_seq/DE_analysis_ZIN/")
  script_dir <- ".."
  data_dir <- "."
  source(file.path(script_dir, "sim_DEpipe_zero_helpers.R"))
}
sapply(c("splatter", "DESeq2"), require, character.only=TRUE)
#"zinbwave", "Seurat", "MAST", "monocle", "edgeR", "BPSC", "scde"
source(file.path(script_dir, "NbZnComBat.R")); source(file.path(script_dir, "helper_seq.R"))
set.seed(123)


####  Load real data
# cts <- read.table(file.path(data_dir, "molecules.txt"), sep = "\t")  # use molecule counts rather than read counts
# anno <- read.table(file.path(data_dir, "annotation.txt"), sep = "\t", header = TRUE)
# ckeep <- anno[, "individual"] == "NA19101" | anno[, "individual"] == "NA19239"
# cts <- as.matrix(cts[, ckeep])
# group <- as.factor(as.character(anno[ckeep, "individual"]))
# batch <- as.factor(as.character(anno[ckeep, "batch"]))
# #batch <- anno[ckeep, "replicate"]
# # remove genes that aren't expressed in at least 6 cells
# gkeep <- rowSums(cts > 0) > 0
# cts <- cts[gkeep, ]
# DE_genes <- read.table(file.path(data_dir, "TPs.txt"), as.is=TRUE)[, 1]
# DE_genes <- intersect(DE_genes, rownames(cts))
# nonDE_genes <- read.table(file.path(data_dir, "TNs.txt"), as.is=TRUE)[, 1]
# nonDE_genes <- intersect(nonDE_genes, rownames(cts))
# params_real <- splatEstimate(cts, params=newSplatParams(group.prob=rep(1/2, 2)))


####  Set parameters
## For combat-seq
zero.fracs.cutoff <- 0.4

## For DE
alpha.unadj <- 0.1
alpha.fdr <- 0.1

## For splatter
data("sc_example_counts")
params <- splatEstimate(sc_example_counts)
# base
params_base <- setParams(params, update=list(mean.rate=0.5))
# batch 1
params1 <- setParams(params, update=list(mean.rate=0.5))
# batch 2
params2 <- setParams(params, update=list(mean.rate=0.5))

# For pipeline
iterations <- 100


####  Run pipeline
for(iter in 1:iterations){
  cat(paste("\nSimulation", iter, '\n'))
  
  ####  Simulate datasets
  # baseline
  sim_base <- splatSimulate(params_base, method="groups", verbose=FALSE)
  # with batch effect
  sim_batch1 <- splatSimulate(params1, method="groups", verbose=FALSE)
  sim_batch2 <- splatSimulate(params2, method="groups", verbose=FALSE)
  batch 
  
  
  ####  Adjust the data with NbZnComBat
  cts_adj <- NbZnCombat(cts, batch, group, full_mod=TRUE, zin.opt=TRUE, zero.fracs.cutoff=zero.fracs.cutoff)
  
  
  ####  DE analysis using DESeq2 
  ##  Baseline
  de0 <- DESeq2_DEpipe(cts_base, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha.unadj, alpha.fdr=alpha.fdr)
  ##  With batch effect
  de1 <- DESeq2_DEpipe(cts, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha.unadj, alpha.fdr=alpha.fdr)
  ##  One-step
  de2 <- DESeq2_DEpipe(cts, batch=batch, group=group, include.batch=TRUE, alpha.unadj=alpha.unadj, alpha.fdr=alpha.fdr)
  ##  ComBat-seq
  de3 <- DESeq2_DEpipe(cts_adj, batch=batch, group=group, include.batch=FALSE, alpha.unadj=alpha.unadj, alpha.fdr=alpha.fdr)
 
   
  ####  Compute performance
  de_lst <- list(Baseline=de0, Batch=de1, One.Step=de2, ComBatseq=de3)
  de_called_unadj <- lapply(de_lst, function(x){x$unadj})
  de_called_fdr <- lapply(de_lst, function(x){x$fdr})
  
  DE_res <- lapply(DEgenes_unadj, perfStats, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  DE_res_fdr <- lapply(DEgenes_fdr, perfStats, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
  
  ####  Write out DE performance
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
