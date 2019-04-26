rm(list=ls())
setwd("~/yuqingz/ComBat_Seq/")
sapply(c("rbenchmark", "Rcpp", "ggplot2", "reshape2", "gridExtra"), require, character.only=TRUE)
source("ComBat_seq.R")
source("helper_seq.R")
Rcpp::sourceCpp("mcint.cpp")
set.seed(123)

n_gene <- 10000
n_sample <- 100

## Simulate count matrix and batch factor
y1_ctrl <- matrix(rnbinom(100*(n_sample/2), mu=10, size=20), ncol=n_sample/2); #mean(c(y1_ctrl)); var(c(y1_ctrl))
y1_case <- matrix(rnbinom(100*(n_sample/2), mu=20, size=20), ncol=n_sample/2); #mean(c(y1_case)); var(c(y1_case))
y1_null <- matrix(rnbinom((n_gene-100)*n_sample, mu=10, size=20), ncol=n_sample)

y2_ctrl <- matrix(rnbinom(100*(n_sample/2), mu=20, size=100), ncol=n_sample/2); #mean(c(y2_ctrl)); var(c(y2_ctrl))
y2_case <- matrix(rnbinom(100*(n_sample/2), mu=40, size=100), ncol=n_sample/2); #mean(c(y2_case)); var(c(y2_case))
y2_null <- matrix(rnbinom((n_gene-100)*n_sample, mu=20, size=100), ncol=n_sample)

counts <- rbind(cbind(y1_ctrl, y1_case, y2_ctrl, y2_case),
                cbind(y1_null, y2_null))
rownames(counts) <- paste0("gene", 1:nrow(counts))
colnames(counts) <- paste0("sample", 1:ncol(counts))

batch <- c(rep("B", ncol(y1_ctrl)+ncol(y1_case)), 
           rep("A", ncol(y2_ctrl)+ncol(y2_case))); table(batch)
group <- c(rep(0, ncol(y1_ctrl)), rep(1, ncol(y1_case)),
           rep(0, ncol(y2_ctrl)), rep(1, ncol(y2_case)))


cmp_res <- benchmark(all_r = ComBat_seq(counts=counts, batch=batch, group=group, Cpp=FALSE),
                     all_c = ComBat_seq(counts=counts, batch=batch, group=group, Cpp=TRUE),
                     sub_r = ComBat_seq(counts=counts, batch=batch, group=group, gene.subset.n=1000, Cpp=FALSE),
                     sub_c = ComBat_seq(counts=counts, batch=batch, group=group, gene.subset.n=1000, Cpp=TRUE), replications=1)
