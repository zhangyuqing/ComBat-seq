rm(list=ls())
sapply(c("SummarizedExperiment", "edgeR", "ggplot2", "reshape2", "gridExtra"), require, character.only=TRUE)
source("../../../ComBat_seq_returnParams.R")
source("../../../helper_seq.R")

### Load data
data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/gfrn_signature"
sigdata <- readRDS(file.path(data_dir, "signature_data.rds"))
cts_mat <- assay(sigdata, "counts")
rownames(cts_mat) <- paste0("gene", 1:nrow(cts_mat))
batch <- colData(sigdata)$batch
group <- colData(sigdata)$group

group_num <- rep(0, ncol(cts_mat))
cond_names <- levels(group)[c(1:3,7:nlevels(group))]
for(i in 1:length(cond_names)){group_num[grep(paste0("^",cond_names[i]), group)] <- i}
colData(sigdata)$condition_bi <- as.factor(group_num)
colData(sigdata)$batch <- as.factor(batch)

### Apply ComBat-seq
combatseq_res <- ComBat_seq(cts_mat, batch=batch, group=group_num, shrink=TRUE, shrink.disp=TRUE, gene.subset.n=1000)
params_gfrn <- combatseq_res[4:7]
#saveRDS(params_gfrn, file="rdata/params_gfrn.rds")

### Visualize
disp_mean_plts <- lapply(levels(factor(batch)), function(b){
  curr_gamma <- params_gfrn$gamma_hat[, paste0("batch", b)]
  curr_phi <- params_gfrn$phi_hat[, paste0("batch", b)]
  curr_df <- data.frame(gamma=curr_gamma, phi=curr_phi)
  plt <- ggplot(curr_df, aes(x=gamma, y=phi)) +
    geom_point() +
    labs(x="Mean estimates (gamma)", y="Dispersion estimates (phi)",
         title=sprintf("GFRN - batch %s", b))
  return(plt)
})

png("figures/gfrn.png", width=13, height=5, units="in", res=300)
grid.arrange(disp_mean_plts[[1]], disp_mean_plts[[2]], disp_mean_plts[[3]], ncol=3)
dev.off()
