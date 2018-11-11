rm(list=ls())
setwd("~/Documents/ComBat_seq/DE_analysis_tmp/")
script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
sapply(c("polyester", "Biostrings", "limma", "edgeR", "DESeq2", "sva", "RUVSeq", "MASS"), require, character.only=TRUE)
source(file.path(script_dir, "simulations/sim_DEpipe_helpers.R"))
source(file.path(script_dir, "ComBat_seq.R")); source(file.path(script_dir, "helper_seq.R"))
set.seed(123)

outBAL <- readRDS("test_balanced.rds")
outUN <- readRDS("test_unbalanced.rds")

outBAL$adj_cts <- ComBat_seq(outBAL$cts, batch=outBAL$batch, group=outBAL$group)
outUN$adj_cts <- ComBat_seq(outUN$cts, batch=outUN$batch, group=outUN$group)

mean(cancelLibsizeEffect(outBAL$adj_cts)[outBAL$true_nulls, outBAL$batch==1])
mean(cancelLibsizeEffect(outBAL$adj_cts)[outBAL$true_nulls, outBAL$batch==2])
mean(cancelLibsizeEffect(outUN$adj_cts)[outUN$true_nulls, outUN$batch==1])
mean(cancelLibsizeEffect(outUN$adj_cts)[outUN$true_nulls, outUN$batch==2])

mean(cancelLibsizeEffect(outBAL$cts)[, outBAL$batch==1])
mean(cancelLibsizeEffect(outBAL$cts)[, outBAL$batch==2])
mean(cancelLibsizeEffect(outUN$cts)[, outUN$batch==1])
mean(cancelLibsizeEffect(outUN$cts)[, outUN$batch==2])

yBAL <- DGEList(counts=outBAL$cts, group=outBAL$group)
yBAL <- calcNormFactors(yBAL, method="TMM")
yBAL$samples

yUN <- DGEList(counts=outUN$cts, group=outUN$group)
yUN <- calcNormFactors(yUN, method="TMM")
yUN$samples



####  DE analysis 
# On baseline dataset without batch effect - independent baseline 
de_called01 <- edgeR_DEpipe(counts_mat=counts_base_indi, batch=batch, group=group, include.batch=FALSE, alpha=alpha)  
perf_stats01 <- perfStats(called_vec=de_called01$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
perf_stats01_fdr <- perfStats(called_vec=de_called01$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
# On counts with batch effect 
de_called1 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha=alpha)  
perf_stats1 <- perfStats(called_vec=de_called1$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
perf_stats1_fdr <- perfStats(called_vec=de_called1$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
 
# One-step - include batch as covariate
de_called2 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=TRUE, alpha=alpha)  
perf_stats2 <- perfStats(called_vec=de_called2$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
perf_stats2_fdr <- perfStats(called_vec=de_called2$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
# Current ComBat + linear model for DE
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
 
# On adjusted count - ComBat-seq 
adj_counts_combatseq <- ComBat_seq(counts=cts, batch=batch, group=group)

de_called5 <- edgeR_DEpipe(counts_mat=adj_counts_combatseq, batch=batch, group=group, include.batch=FALSE, alpha=alpha)  
perf_stats5 <- perfStats(called_vec=de_called5$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
perf_stats5_fdr <- perfStats(called_vec=de_called5$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))

# Compare with RUVseq 
uvseq <- RUVg(cts, cIdx=sample(G_nulls,10,replace=F), k=1)
de_called6 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha=alpha, covar=uvseq$W)  
perf_stats6 <- perfStats(called_vec=de_called6$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
perf_stats6_fdr <- perfStats(called_vec=de_called6$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  
# Compare with SVAseq 
mod1 <- model.matrix(~as.factor(group)); mod0 <- cbind(mod1[,1])
svseq <- svaseq(cts, mod1, mod0, n.sv=1)
de_called7 <- edgeR_DEpipe(counts_mat=cts, batch=batch, group=group, include.batch=FALSE, alpha=alpha, covar=svseq$sv)  
perf_stats7 <- perfStats(called_vec=de_called7$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
perf_stats7_fdr <- perfStats(called_vec=de_called7$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(cts))
  

# Unadjusted P values
DE_res <- matrix(c(perf_stats01,  # baseline - indipendent simulation
                   perf_stats1,   # data with batch effect
                   perf_stats2,   # one-step approaches
                   perf_stats3,   # current combat
                   perf_stats5,   # ComBat-seq
                   perf_stats6,   # RUVseq
                   perf_stats7),  # SVAseq
                 nrow=3, byrow=FALSE)
rownames(DE_res) <- names(perf_stats01)
# FDR adjusted Q values
DE_res_fdr <- matrix(c(perf_stats01_fdr, # baseline - indipendent simulation
                       perf_stats1_fdr,  # data with batch effect
                       perf_stats2_fdr,  # one-step approaches
                       perf_stats3_fdr,  # current combat
                       perf_stats5_fdr,  # ComBat-seq
                       perf_stats6_fdr, # RUVseq
                       perf_stats7_fdr),  # SVAseq
                     nrow=3, byrow=FALSE)
rownames(DE_res_fdr) <- names(perf_stats01_fdr)
colnames(DE_res) <- colnames(DE_res_fdr) <- c("BaseIndi.edgeR", "Batch.edgeR", "OneStep.edgeR", "ComBat.lm", 
                                              "ComBatseq.edgeR", "RUVseq.edgeR", "SVAseq.edgeR")
DE_res <- as.data.frame(DE_res); DE_res_fdr <- as.data.frame(DE_res_fdr)
  

