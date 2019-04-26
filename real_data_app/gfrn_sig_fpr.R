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

## Parameters
DEmethod_name <- "DESeq2"
iterations <- 3
alpha_fpr_seq <- seq(0, 0.2, 0.05)[-1]

## Load data
sigdata <- readRDS(file.path(data_dir, "signature_data.rds"))
cts_mat <- assay(sigdata, "counts")   # count matrix (also have tpm and fpkm in there)
# filter out genes with 0 counts
#cts_mat <- cts_mat[apply(cts_mat, 1, function(x){all(x!=0)}), ]
rownames(cts_mat) <- paste0("gene", 1:nrow(cts_mat))
batch <- colData(sigdata)$batch
group <- colData(sigdata)$group

## Take controls out
cts_ctrl <- cts_mat[, group %in% c("gfp_for_egfr", "gfp18", "gfp30")]
batch_ctrl <- batch[group %in% c("gfp_for_egfr", "gfp18", "gfp30")]
#remove genes with only 0 counts in all control samples
cts_ctrl <- cts_ctrl[apply(cts_ctrl,1,function(x){!all(x==0)}), ]

## Use ComBatSeq to adjust data
#counts=cts_ctrl;batch=batch_ctrl;group=NULL;full_mod=FALSE;covar_mod=NULL
start_time <- Sys.time()
combatseq_ctrl <- ComBat_seq(counts=cts_ctrl, batch=batch_ctrl, group=NULL, full_mod=FALSE, gene.subset.n=1000, Cpp=FALSE)
end_time <- Sys.time()
print(end_time - start_time)
#combatseq_ctrl <- readRDS("~/Google Drive/ComBat_seq/real_data_example/RNAseq/gfrn_signature/combatseq_ctrlonly.rds")
data_lst <- list(Raw.Counts=cts_ctrl, CombatSeq=combatseq_ctrl)

## DE
if(DEmethod_name=="edgeR"){DE_method <- edgeR_DEpipe}else if(DEmethod_name=="DESeq2"){DE_method <- DESeq2_DEpipe} 
fpr_df <- list(); ii <- 1
for(iter in 1:iterations){
  for(alpha_unadj in alpha_fpr_seq){
    group_rand <- sample(c(0,1), length(batch_ctrl), replace=TRUE)  # randomly split groups for the control samples
    de_res <- DEpipe(data_lst, batch=batch_ctrl, group=group_rand, DE_method=DE_method, alpha.unadj=alpha_unadj)
    de_genes <- lapply(de_res, function(de){de$unadj})
    fpr_df[[ii]] <- c(PValue_cutoff=alpha_unadj, sapply(de_genes, length)/nrow(cts_ctrl))
    ii <- ii + 1
  }
}
fpr_df <- as.data.frame(do.call(rbind, fpr_df))

## Visualize results
fpr_df_mlt <- melt(fpr_df, id.vars="PValue_cutoff", variable.name="Method")
res_df <- fpr_df_mlt %>%
  group_by(PValue_cutoff, Method) %>%
  summarise(Avg=median(value)) 
res_df$Method <- factor(res_df$Method, levels=c("raw.count","one.step","ruvseq","svaseq","combatseq"))
png("gfrn_sig_fpr.png", width=7, height=6, units="in", res=300)
ggplot(res_df, aes(x=PValue_cutoff, y=Avg, group=Method, color=Method)) +
  geom_line() +
  #geom_abline(slope=1, intercept=0, col="black") +
  #ylim(0, 0.1) + xlim(0, 0.1) +
  labs(x="P Value threshold", y=sprintf("observed FPR"))   #average over %s random split)", iterations))
dev.off()



