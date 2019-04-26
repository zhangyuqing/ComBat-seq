rm(list=ls())
sapply(c("dplyr", "edgeR", "DESeq2", "sva", "RUVSeq", "MASS", "ggplot2", "reshape2", "gridExtra"), require, character.only=TRUE)
demo <- FALSE   # if testing code, set as TRUE; if running simulations, set as FALSE
if(demo){
  setwd("~/Documents/ComBat_seq/real_data_app/")
  script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
  data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/gfrn_signature"
  source("~/Dropbox/Work/ComBat_Seq/ComBat-Seq/real_data_app/gfrn_sig_helpers.R")
}else{
  setwd("~/yuqingz/ComBat_seq/real_data_app/")
  script_dir <- ".."; data_dir <- "."
  source("gfrn_sig_helpers.R")
}
set.seed(1)

## Load ComBat-seq functions
source(file.path(script_dir, "ComBat_seq.R"))
source(file.path(script_dir, "helper_seq.R"))
#sourceCpp(file.path(script_dir, "mcint.cpp"))

## Parameters
DEmethod_name <- "DESeq2"
pathway_regex <- "^egfr"  
alpha_fdr_seq <- seq(0, 1, 0.005)[-1]

## Load data
sigdata <- readRDS(file.path(data_dir, "signature_data.rds"))
cts_mat <- assay(sigdata, "counts")   # count matrix (also have tpm and fpkm in there)
# filter out genes with 0 counts
#cts_mat <- cts_mat[apply(cts_mat, 1, function(x){all(x!=0)}), ]
rownames(cts_mat) <- paste0("gene", 1:nrow(cts_mat))
batch <- colData(sigdata)$batch
group <- colData(sigdata)$group

## Take the subset (all controls & 1 condition specified by pathway_regex)
pathway_condition_ind <- grep(pathway_regex, group)
pathway_batch <- unique(batch[pathway_condition_ind])  # which batch does this pathway belong to
ctrl_ind <- which(group %in% c("gfp_for_egfr", "gfp18", "gfp30"))
subset_ind <- sort(c(ctrl_ind, pathway_condition_ind))  
cts_sub <- cts_mat[, subset_ind]
batch_sub <- batch[subset_ind]
group_sub <- group[subset_ind]
# remove genes with only 0 counts in the subset
cts_sub <- cts_sub[apply(cts_sub,1,function(x){!all(x==0)}), ]


## DE within batch
if(DEmethod_name=="edgeR"){DE_method <- edgeR_DEpipe}else if(DEmethod_name=="DESeq2"){DE_method <- DESeq2_DEpipe} 

cts_sub_inbatch <- cts_sub[, batch_sub==pathway_batch]
group_sub_inbatch <- as.factor(as.character(group_sub[batch_sub==pathway_batch]))
group_sub_inbatch_num <- rep(1, length(group_sub_inbatch))
group_sub_inbatch_num[grep("gfp", group_sub_inbatch)] <- 0
# remove genes with only 0 counts in this batch
cts_sub_inbatch <- cts_sub_inbatch[apply(cts_sub_inbatch,1,function(x){!all(x==0)}), ]
  
de_res_inbatch <- DE_method(cts_sub_inbatch, batch=rep(1, ncol(cts_sub_inbatch)), group=factor(group_sub_inbatch_num), 
                            include.batch=FALSE, alpha.unadj=1, alpha.fdr=1)


## DE merging controls
#apply ComBatSeq
group_sub <- as.character(group_sub)
group_sub_num <- rep(1, length(group_sub))
group_sub_num[grep("gfp", group_sub)] <- 0

#counts=cts_sub;batch=batch_sub;group=group_sub_num;gene.subset.n=1000;Cpp=FALSE;covar_mod=NULL;full_mod=TRUE
start_time <- Sys.time()
combatseq_sub <- ComBat_seq(counts=cts_sub, batch=batch_sub, group=group_sub_num, gene.subset.n=1000, Cpp=FALSE)
end_time <- Sys.time()
print(end_time - start_time)
#combatseq_sub <- readRDS("~/Google Drive/ComBat_seq/real_data_example/RNAseq/gfrn_signature/combatseq_test.rds")
data_lst <- list(Raw.Counts=cts_sub, CombatSeq=combatseq_sub)

de_res <- DEpipe(data_lst, batch=batch_sub, DE_method=DE_method, group=group_sub_num)


## Power when controlling FDR
power_res <- list(); ii <- 1
for(alpha_fdr in alpha_fdr_seq){
  de_genes <- c(list(baseline=rownames(de_res_inbatch$de_res)[de_res_inbatch$de_res[, ncol(de_res_inbatch$de_res)] < alpha_fdr]),
                lapply(de_res, function(de){rownames(de$de_res)[de$de_res[, ncol(de$de_res)] < alpha_fdr]}))
  de_genes <- lapply(de_genes, function(x){x[!is.na(x)]})
  power_res[[ii]] <- c(FDR_cutoff=alpha_fdr, sapply(de_genes, length))
  ii <- ii + 1
}
power_res <- as.data.frame(do.call(rbind, power_res))
saveRDS(power_res, file="gfrn_nDE.rds")


## Visualize results
power_res_mlt <- melt(power_res, id.vars="FDR_cutoff", variable.name="Method")
power_res_mlt$Method <- factor(power_res_mlt$Method, levels=c("baseline","raw.count","one.step","ruvseq","svaseq","combatseq"))
png("gfrn_sig_power.png", width=7, height=7, units="in", res=300)
ggplot(power_res_mlt, aes(x=FDR_cutoff, y=value, group=Method, color=Method)) +
  geom_line() +
  labs(x="FDR threshold", y="Number of detected DE genes")
dev.off()
