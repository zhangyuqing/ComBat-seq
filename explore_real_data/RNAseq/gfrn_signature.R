rm(list=ls())
setwd("/restricted/projectnb/combat/work/yuqingz/ComBat-seq/")
sapply(c("SummarizedExperiment"), require, character.only=TRUE)

####  Load data
data_dir <- "/restricted/projectnb/pathsig/signatures/fastq/u01_combined"
sigdata <- readRDS(file.path(data_dir, "signature_data.rds"))
#count matrix (also have tpm and fpkm in there)
cts_mat <- assay(sigdata, "counts")
#batch annotation
batch <- colData(sigdata)$batch
#signature annotation
group <- colData(sigdata)$group

