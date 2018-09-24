rm(list=ls())
setwd("/restricted/projectnb/combat/work/yuqingz/ComBat-seq/")
sapply(c("ggplot2", "reshape2", "gridExtra"), require, character.only=TRUE)
source("ComBat_seq.R"); source("helper_seq.R")
set.seed(123)

n_gene_seq <- seq(1000, 40000, 1000)
run_time_seq <- rep(NA, length(n_gene_seq))

for(i in 1:length(n_gene_seq)){
  ## Simulate count matrix and batch factor
  y1_ctrl <- matrix(rnbinom(100*5, mu=10, size=20), ncol=5); #mean(c(y1_ctrl)); var(c(y1_ctrl))
  y1_case <- matrix(rnbinom(100*5, mu=20, size=20), ncol=5); #mean(c(y1_case)); var(c(y1_case))
  y1_null <- matrix(rnbinom((n_gene_seq[i]-100)*10, mu=10, size=20), ncol=10)
  
  y2_ctrl <- matrix(rnbinom(100*5, mu=20, size=100), ncol=5); #mean(c(y2_ctrl)); var(c(y2_ctrl))
  y2_case <- matrix(rnbinom(100*5, mu=40, size=100), ncol=5); #mean(c(y2_case)); var(c(y2_case))
  y2_null <- matrix(rnbinom((n_gene_seq[i]-100)*10, mu=20, size=100), ncol=10)
  
  counts <- rbind(cbind(y1_ctrl, y1_case, y2_ctrl, y2_case),
                  cbind(y1_null, y2_null))
  rownames(counts) <- paste0("gene", 1:nrow(counts))
  colnames(counts) <- paste0("sample", 1:ncol(counts))
  
  batch <- c(rep("B", ncol(y1_ctrl)+ncol(y1_case)), 
             rep("A", ncol(y2_ctrl)+ncol(y2_case))); table(batch)
  group <- c(rep(0, ncol(y1_ctrl)), rep(1, ncol(y1_case)),
             rep(0, ncol(y2_ctrl)), rep(1, ncol(y2_case)))
  
  
  ## Execute ComBat-seq and record running time
  start_time <- Sys.time()
  adj_counts <- ComBat_seq(counts=counts, batch=batch, group=group)
  end_time <- Sys.time()
  run_time_seq[i] <- end_time - start_time
}

run_time_df <- data.frame(N_genes=n_gene_seq, run_time=run_time_seq)
png("running_time.png", width=6, height=6, units="in", res=300)
ggplot(data=run_time_df, aes(x=N_genes, y=run_time, group=1)) +
  geom_line()+
  geom_point() +
  labs(title="Running time of ComBat-seq", 
       x="Number of genes", y="Running time")
dev.off()

