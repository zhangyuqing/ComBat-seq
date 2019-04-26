rm(list=ls())
sapply(c("dplyr", "plyr", "edgeR", "DESeq2", "sva", "RUVSeq", "MASS", "ggplot2", "reshape2", "gridExtra"), require, character.only=TRUE)
demo <- FALSE   # if testing code, set as TRUE; if running simulations, set as FALSE
if(demo){
  setwd("~/Documents/ComBat_seq/real_data_app/")
  script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
  data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/WRBU_lungcancer"
  source("~/Dropbox/Work/ComBat_Seq/ComBat-Seq/real_data_app/gfrn_sig_helpers.R")
}else{
  setwd("~/yuqingz/ComBat_seq/real_data_app/")
  script_dir <- ".."; data_dir <- "."
  source("gfrn_sig_helpers.R")
}
set.seed(1)
source(file.path(script_dir, "ComBat_seq.R"))
source(file.path(script_dir, "helper_seq.R"))

## Parameters
DEmethod_name <- "DESeq2"
if(DEmethod_name=="edgeR"){DE_method <- edgeR_DEpipe}else if(DEmethod_name=="DESeq2"){DE_method <- DESeq2_DEpipe} 
alpha_fpr_seq <- seq(0, 1, 0.005)[-1]
iterations <- 100

## Load data
cts <- readRDS(file.path(data_dir, "count_matrix_BU_WR.rds"))
meta_info <- readRDS(file.path(data_dir, "BU_WR_annotation_file.rds"))
# drop 4790-024 (single sample as its own batch)
cts <- cts[, -grep("4790-024", colnames(cts))]
meta_info <- filter(meta_info, kitnumber!="4790-024")
batch <- factor(meta_info$site)
group <- factor(meta_info$smoking_status)
# remove all 0 genes
tmp <- t(apply(cts,1,floor)); colnames(tmp) <- colnames(cts); cts <- tmp; rm(tmp)
cts <- cts[rowVars(cts)>0, ]

## Take Former smoker subset
cts_ctrl <- cts[, group == "Former smoker"]
batch_ctrl <- batch[group == "Former smoker"]
#remove genes with only 0 counts in all control samples
cts_ctrl <- cts_ctrl[apply(cts_ctrl,1,function(x){!all(x==0)}), ]

## Use ComBatSeq to adjust data
#counts=cts_ctrl;batch=batch_ctrl;group=NULL;full_mod=FALSE;covar_mod=NULL
start_time <- Sys.time()
combatseq_ctrl <- ComBat_seq(counts=cts_ctrl, batch=batch_ctrl, group=NULL, full_mod=FALSE, gene.subset.n=1000, Cpp=FALSE)
end_time <- Sys.time()
print(end_time - start_time)
#combatseq_ctrl <- readRDS("~/Google Drive/ComBat_seq/real_data_example/RNAseq/WRBU_lungcancer/combatseq_ctrlonly_former.rds")
data_lst <- list(Raw.Counts=cts_ctrl, CombatSeq=combatseq_ctrl)

## DE
if(DEmethod_name=="edgeR"){DE_method <- edgeR_DEpipe}else if(DEmethod_name=="DESeq2"){DE_method <- DESeq2_DEpipe} 
fpr_df <- list(); ii <- 1
for(iter in 1:iterations){
  for(alpha_unadj in alpha_fpr_seq){
    print(sprintf("iter %s, alpha: %s", iter, alpha_unadj))
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
fpr_df_mlt$PValue_cutoff <- factor(fpr_df_mlt$PValue_cutoff)
res_df <- fpr_df_mlt %>%
  dplyr::group_by(PValue_cutoff, Method) %>%
  dplyr::summarise(Avg=mean(value)) 
res_df$Method <- factor(res_df$Method, levels=c("raw.count","one.step","ruvseq","svaseq","combatseq"))
png("wrbu_smoke_fpr_enlarged.png", width=7, height=6, units="in", res=300)
ggplot(res_df, aes(x=as.numeric(as.character(PValue_cutoff)), y=Avg, group=Method, color=Method)) +
  geom_line() +
  xlim(0, 0.3) + ylim(0, 0.25)+
  labs(x="P Value threshold", y=sprintf("observed FPR"))   #average over %s random split)", iterations))
dev.off()

