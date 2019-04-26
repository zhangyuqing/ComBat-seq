rm(list=ls())
sapply(c("dplyr", "plyr", "SummarizedExperiment", "ggplot2", "reshape2", "gridExtra", "ggpubr"), require, character.only=TRUE)
data_name <- "CHD8"

if(data_name == "wrbu_lungcancer"){
  data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/WRBU_lungcancer"
  cts <- readRDS(file.path(data_dir, "count_matrix_BU_WR.rds"))
  meta_info <- readRDS(file.path(data_dir, "BU_WR_annotation_file.rds"))
  # drop 4790-024 (single sample as its own batch)
  cts <- cts[, -grep("4790-024", colnames(cts))]
  meta_info <- filter(meta_info, kitnumber!="4790-024")
  batch <- factor(meta_info$site)
  group <- factor(meta_info$smoking_status)
  params <- readRDS("rdata/params_wrbu_lungcancer.rds")
  plt.w <- 8; plt.h <- 8
}else if(data_name == "TB"){
  data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/TB"
  rds_obj <- readRDS(file.path(data_dir, "combined.rds"))
  cts <- assays(rds_obj)$counts
  batch <- colData(rds_obj)$SequencingBatch
  group <- colData(rds_obj)$Label
  params <- readRDS("rdata/params_TB.rds")
  plt.w <- 8; plt.h <- 8
}else if(data_name == "CHD8"){
  load("~/Google Drive/ComBat_seq/real_data_example/RNAseq/CHD8/rse_gene.Rdata")
  cts <- assays(rse_gene)$counts
  batch <- colData(rse_gene)$batch
  group <- colData(rse_gene)$condition_bi
  params <- readRDS("rdata/params_CHD8.rds")
  plt.w <- 8; plt.h <- 8
}

batch_levels <- levels(factor(batch))
if(data_name=="TB"){batch_levels <- batch_levels[1:2]}

violin_plts <- lapply(batch_levels, function(b, iqr_gamma, iqr_phi){
  curr_gamma <- params$gamma_hat[, paste0("batch", b)]
  curr_phi <- params$phi_hat[, paste0("batch", b)]
  
  gamma_upper <- quantile(curr_gamma, 0.75) + iqr_gamma * IQR(curr_gamma)
  gamma_lower <- quantile(curr_gamma, 0.25) - iqr_gamma * IQR(curr_gamma)
  phi_upper <- quantile(curr_phi, 0.75) + iqr_phi * IQR(curr_phi)
  
  ## violin plots
  violin_gamma <- ggplot(data.frame(gamma=curr_gamma), aes(x=0, y=gamma)) +
    geom_violin() +
    geom_hline(yintercept=gamma_upper, color="red") +
    geom_hline(yintercept=gamma_lower, color="red") +
    #geom_jitter() +
    coord_flip() +
    labs(x="Mean estimates (gamma)") +
    theme(axis.title.x=element_blank())
  violin_phi <- ggplot(data.frame(phi=curr_phi), aes(x=0, y=phi)) +
    geom_violin() +
    geom_hline(yintercept=phi_upper, color="red") +
    #geom_jitter() +
    coord_flip() +
    labs(x="Dispersion estimates (phi)") +
    theme(axis.title.x=element_blank())
  return(annotate_figure(ggarrange(violin_gamma, violin_phi, ncol=2),
                         top = text_grob(sprintf("Violin plots of parameters, batch %s", b))))
}, iqr_gamma=1.5, iqr_phi=1.5)

png(sprintf("./figures/distributions_%s.png", data_name), width=plt.w, height=plt.h, units="in", res=300)
do.call(grid.arrange, c(violin_plts, ncol=1))
dev.off()


# plts <- lapply(levels(factor(batch)), function(b){
#   curr_gamma <- params_wrbu_lungcancer$gamma_hat[, paste0("batch", b)]
#   curr_phi <- params_wrbu_lungcancer$phi_hat[, paste0("batch", b)]
#   
#   ## histograms
#   hist_gamma <- ggplot(data.frame(gamma=curr_gamma), aes(x=gamma)) +
#     geom_histogram(bins=round(length(curr_gamma)/500)) +
#     labs(x="Mean estimates (gamma)")
#   hist_phi <- ggplot(data.frame(phi=curr_phi), aes(x=phi)) +
#     geom_histogram(bins=round(length(curr_phi)/500)) +
#     labs(x="Dispersion estimates (phi)")
#   hist_plts <- annotate_figure(ggarrange(hist_gamma, hist_phi, ncol=2), 
#                                top = text_grob(sprintf("Histograms of parameters, lung cancer batch %s", b)))
#     
#   ## boxplots
#   box_gamma <- ggplot(data.frame(gamma=curr_gamma), aes(x=0, y=gamma)) +
#     geom_boxplot() +
#     coord_flip() +
#     labs(x="Mean estimates (gamma)") +
#     theme(axis.title.x=element_blank())
#   box_phi <- ggplot(data.frame(phi=curr_phi), aes(x=0, y=phi)) +
#     geom_boxplot() +
#     coord_flip() +
#     labs(x="Dispersion estimates (phi)") +
#     theme(axis.title.x=element_blank())
#   box_plts <- annotate_figure(ggarrange(box_gamma, box_phi, ncol=2), 
#                               top = text_grob(sprintf("Boxplots of parameters, lung cancer batch %s", b)))
#   
#   ## violin plots
#   violin_gamma <- ggplot(data.frame(gamma=curr_gamma), aes(x=0, y=gamma)) +
#     geom_violin() +
#     coord_flip() +
#     labs(x="Mean estimates (gamma)") +
#     theme(axis.title.x=element_blank())
#   violin_phi <- ggplot(data.frame(phi=curr_phi), aes(x=0, y=phi)) +
#     geom_violin() +
#     coord_flip() +
#     labs(x="Dispersion estimates (phi)") +
#     theme(axis.title.x=element_blank())
#   violin_plts <- annotate_figure(ggarrange(violin_gamma, violin_phi, ncol=2), 
#                                  top = text_grob(sprintf("Violin plots of parameters, lung cancer batch %s", b)))
#   
#   ## combine plots
#   merged_plts <- grid.arrange(hist_plts, box_plts, nrow=2)
#   return(list(hist=hist_plts, box=box_plts, violin=violin_plts, merged=merged_plts))
# }) 

### plot counts:
###### CHD8: ENSG00000135547.8

# bins_distrib <- hist(curr_gamma, breaks=round(length(curr_gamma)/500), plot=FALSE)
# bins_distrib$counts


