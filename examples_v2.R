####  Same library size; inspired by example in Leek 2010 Nature Reviews Genetics paper 
rm(list=ls())
setwd("C:/Users/zhang/Dropbox/Work/ComBat_Seq/ComBat-Seq/")
sapply(c("ggplot2", "reshape2"), require, character.only=TRUE)
set.seed(123)

## Simulate count matrix and batch factor
y1 <- matrix(rnbinom(1000*10,mu=10,size=10),ncol=10); mean(c(y1)); var(c(y1))
y2 <- matrix(rnbinom(1000*10,mu=100,size=5),ncol=10); mean(c(y2)); var(c(y2))
y3 <- matrix(rnbinom(1000*10,mu=100,size=10),ncol=10); mean(c(y3)); var(c(y3))
y4 <- matrix(rnbinom(1000*10,mu=10,size=5), ncol=10); mean(c(y4)); var(c(y4))
counts <- rbind(cbind(y1, y2), cbind(y3, y4))
batch <- c(rep("B", ncol(y1)), rep("A", ncol(y2))); table(batch)
group <- rep(0, ncol(y1)+ncol(y2))
full_mod <- TRUE
source("helper_seq.R")

# visualize simulated count matrix
lib_sizes <- colSums(counts)


