rm(list=ls())
sapply(c("dplyr", "plyr", "edgeR", "DESeq2", "sva", "RUVSeq", "MASS", "ggplot2", "reshape2", "gridExtra"), require, character.only=TRUE)
demo <- TRUE   # if testing code, set as TRUE; if running simulations, set as FALSE
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
alpha_fdr_seq <- seq(0, 1, 0.005)[-1]
covar_incl_name <- "kitnumber"

if(DEmethod_name=="edgeR"){DE_method <- edgeR_DEpipe}else if(DEmethod_name=="DESeq2"){DE_method <- DESeq2_DEpipe} 


## Load data
cts <- readRDS(file.path(data_dir, "count_matrix_BU_WR.rds"))
meta_info <- readRDS(file.path(data_dir, "BU_WR_annotation_file.rds"))
# drop 4790-024 (single sample as its own batch)
cts <- cts[, -grep("4790-024", colnames(cts))]
meta_info <- filter(meta_info, kitnumber!="4790-024")
batch <- factor(meta_info$site)
group <- factor(meta_info$smoking_status)
if(!is.null(covar_incl_name)){covar_incl <- meta_info$kitnumber}
# remove all 0 genes
tmp <- t(apply(cts,1,floor)); colnames(tmp) <- colnames(cts); cts <- tmp; rm(tmp)
cts <- cts[rowVars(cts)>0, ]



## DE within each batch
cts_WR <- cts[, batch=="WR"]
group_WR <- group[batch=="WR"]
group_WR <- as.factor(as.character(revalue(group_WR, c("Former smoker"="0", "Current smoker"="1"))))
de_WR <- DE_method(cts_WR, batch=rep(1, ncol(cts_WR)), group=group_WR, include.batch=FALSE, alpha.unadj=1, alpha.fdr=1)

cts_BU <- cts[, batch=="BU"]
group_BU <- group[batch=="BU"]
group_BU <- as.factor(as.character(revalue(group_BU, c("Former smoker"="0", "Current smoker"="1"))))
de_BU <- DE_method(cts_BU, batch=rep(1, ncol(cts_BU)), group=group_BU, include.batch=FALSE, alpha.unadj=1, alpha.fdr=1)


## Double sample size
cts_WR_dup <- cbind(cts_WR, cts_WR)
group_WR_dup <- factor(c(as.character(group_WR), as.character(group_WR)))
#counts_mat=cts_WR_dup;batch=rep(1, ncol(cts_WR_dup));group=group_WR_dup;include.batch=FALSE;alpha.unadj=1;alpha.fdr=1;covar_incl=covar_incl;covar=NULL
de_WR_dup <- DE_method(cts_WR_dup, batch=rep(1, ncol(cts_WR_dup)), group=group_WR_dup, 
                       include.batch=FALSE, alpha.unadj=1, alpha.fdr=1)

cts_BU_dup <- cbind(cts_BU, cts_BU)
group_BU_dup <- factor(c(as.character(group_BU), as.character(group_BU)))
de_BU_dup <- DE_method(cts_BU_dup, batch=rep(1, ncol(cts_BU_dup)), group=group_BU_dup, 
                       include.batch=FALSE, alpha.unadj=1, alpha.fdr=1)



## DE in combined data
#counts=cts;batch=batch;group=group;gene.subset.n=1000;Cpp=FALSE;covar_mod=NULL;full_mod=TRUE
group <- as.factor(as.character(revalue(group, c("Former smoker"="0", "Current smoker"="1"))))

# run combat-seq
start_time <- Sys.time()
combatseq_cts <- ComBat_seq(counts=cts, batch=batch, group=group, gene.subset.n=1000, Cpp=FALSE)
end_time <- Sys.time()
print(end_time - start_time)   # 20MIN
# if(demo){
#   combatseq_cts <- readRDS("~/Google Drive/ComBat_seq/real_data_example/RNAseq/WRBU_lungcancer/combatseq_cts.rds")
# }else{
#   combatseq_cts <- readRDS("combatseq_cts.rds")
# }

data_lst <- list(Raw.Counts=cts, CombatSeq=combatseq_cts)
de_res <- DEpipe(data_lst, batch=batch, DE_method=DE_method, group=group)


## Power when controlling FDR
power_res <- list(); ii <- 1
for(alpha_fdr in alpha_fdr_seq){
  de_genes <- c(list(baseline_WR=rownames(de_WR$de_res)[de_WR$de_res[, ncol(de_WR$de_res)] < alpha_fdr]),
                list(baseline_BU=rownames(de_BU$de_res)[de_BU$de_res[, ncol(de_BU$de_res)] < alpha_fdr]),
                list(base_WR_dup=rownames(de_WR_dup$de_res)[de_WR_dup$de_res[, ncol(de_WR_dup$de_res)] < alpha_fdr]),
                list(base_BU_dup=rownames(de_BU_dup$de_res)[de_BU_dup$de_res[, ncol(de_BU_dup$de_res)] < alpha_fdr]),
                lapply(de_res, function(de){rownames(de$de_res)[de$de_res[, ncol(de$de_res)] < alpha_fdr]}))
  de_genes <- lapply(de_genes, function(x){x[!is.na(x)]})
  power_res[[ii]] <- c(FDR_cutoff=alpha_fdr, sapply(de_genes, length))
  ii <- ii + 1
}
power_res <- as.data.frame(do.call(rbind, power_res))
saveRDS(power_res, file="wrbu_smoke_nDE_dup.rds")


## Visualize results
power_res_mlt <- melt(power_res, id.vars="FDR_cutoff", variable.name="Method")
power_res_mlt$Method <- factor(power_res_mlt$Method, 
                               levels=c("baseline_WR", "baseline_BU", "base_WR_dup", "base_BU_dup",
                                        "raw.count","one.step","ruvseq","svaseq","combatseq"))
png("wrbu_smoke_power_dup.png", width=7, height=7, units="in", res=300)
ggplot(power_res_mlt, aes(x=FDR_cutoff, y=value, group=Method, color=Method)) +
  geom_line() +
  labs(x="FDR threshold", y="Number of detected DE genes")
dev.off()

png("wrbu_smoke_power_enlarged_dup.png", width=7, height=7, units="in", res=300)
ggplot(power_res_mlt, aes(x=FDR_cutoff, y=value, group=Method, color=Method)) +
  geom_line() +
  xlim(0,0.2)+ylim(0, 1200)+
  labs(x="FDR threshold", y="Number of detected DE genes")
dev.off()
