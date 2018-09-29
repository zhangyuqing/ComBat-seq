rm(list=ls())
setwd("")
sapply(c("ggplot2", "reshape2", "gridExtra", "dendextend", "edgeR", "DESeq2", "polyester"), require, character.only=TRUE)
source("../ComBat_seq.R"); source("../helper_seq.R")
set.seed(123)


####  Parameters


####  Simulate datasets


####  Apply ComBat-seq


####  DE

