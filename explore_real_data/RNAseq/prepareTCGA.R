rm(list=ls())
setwd("~/Documents/ComBat_seq/TCGA")
sapply(c("curatedTCGAData", "MultiAssayExperiment"), require, character.only=TRUE)
curatedTCGAData(diseaseCode="*", assays="RNASeqGene", dry.run=TRUE)

# download data
tcga_dat_obj <- curatedTCGAData("BRCA", "RNASeqGene", FALSE)
cts <- assays(tcga_dat_obj)[[1]]

# tmp=read.delim2("~/Downloads/gdac.broadinstitute.org_GBMLGG.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/GBMLGG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt")
# tmp=read.delim2("~/Downloads/gdac.broadinstitute.org_BRCA.Merge_Clinical.Level_1.2016012800.0.0/BRCA.clin.merged.txt", header=T)
