rm(list=ls())
sapply(c("SummarizedExperiment", "edgeR", "ggplot2", "reshape2", "gridExtra"), require, character.only=TRUE)
source("../../../ComBat_seq_returnParams.R")
source("../../../helper_seq.R")

### Load data
load("~/Google Drive/ComBat_seq/real_data_example/RNAseq/CHD8/rse_gene.Rdata")
cts <- assays(rse_gene)$counts
batch <- colData(rse_gene)$batch
group <- colData(rse_gene)$condition_bi


### Apply ComBat-seq
combatseq_res <- ComBat_seq(cts, batch=batch, group=group, shrink=FALSE)
params_CHD8 <- combatseq_res[4:7]
#saveRDS(params_CHD8, file="rdata/params_CHD8.rds")


### Visualize
disp_mean_plts <- lapply(levels(factor(batch)), function(b){
  curr_gamma <- params_CHD8$gamma_hat[, paste0("batch", b)]
  curr_phi <- params_CHD8$phi_hat[, paste0("batch", b)]
  curr_df <- data.frame(gamma=curr_gamma, phi=curr_phi)
  plt <- ggplot(curr_df, aes(x=gamma, y=phi)) +
    geom_point() +
    labs(x="Mean estimates (gamma)", y="Dispersion estimates (phi)",
         title=sprintf("CHD8 - batch %s, Spearman corr = %s", 
                       b, round(cor(curr_gamma, curr_phi, method="spearman"), 3)))
  return(plt)
})

png("figures/CHD8.png", width=10, height=5, units="in", res=300)
grid.arrange(disp_mean_plts[[1]], disp_mean_plts[[2]], ncol=2)
dev.off()
