rm(list=ls())
sapply(c("SummarizedExperiment", "edgeR", "ggplot2", "reshape2", "gridExtra"), require, character.only=TRUE)
source("../../../ComBat_seq_returnParams.R")
source("../../../helper_seq.R")

### Load data
data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/TB"
rds_obj <- readRDS(file.path(data_dir, "combined.rds"))
cts <- assays(rds_obj)$counts
batch <- colData(rds_obj)$SequencingBatch
group <- colData(rds_obj)$Label
covar <- colData(rds_obj)$Sex

### Apply ComBat-seq
combatseq_res <- ComBat_seq(cts, batch=batch, group=group, covar_mod=model.matrix(~covar), shrink=FALSE)
params_TB <- combatseq_res[4:7]
#saveRDS(params_TB, file="rdata/params_TB.rds")

### Visualize
disp_mean_plts <- lapply(levels(factor(batch)), function(b){
  curr_gamma <- params_TB$gamma_hat[, paste0("batch", b)]
  curr_phi <- params_TB$phi_hat[, paste0("batch", b)]
  curr_df <- data.frame(gamma=curr_gamma, phi=curr_phi)
  plt <- ggplot(curr_df, aes(x=gamma, y=phi)) +
    geom_point() +
    labs(x="Mean estimates (gamma)", y="Dispersion estimates (phi)",
         title=sprintf("TB - batch %s, Spearman corr = %s", 
                       b, round(cor(curr_gamma, curr_phi, method="spearman"), 3)))
  return(plt)
})

png("figures/TB.png", width=15, height=10, units="in", res=300)
grid.arrange(disp_mean_plts[[1]], disp_mean_plts[[2]], disp_mean_plts[[3]],
             disp_mean_plts[[4]], disp_mean_plts[[5]], nrow=2, ncol=3)
dev.off()
