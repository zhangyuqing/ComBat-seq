rm(list=ls())
sapply(c("ggplot2", "reshape2", "gridExtra", "dendextend", "edgeR", "DESeq2", "scales"), require, character.only=TRUE)
source("../../ComBat_seq.R");  source("../../helper_seq.R"); source("toy_helpers.R")
set.seed(123)


#### Experiment set up
exp_name <- "check_disp_shrinkage_N5_p10000_sub1000_seed123"
if(!dir.exists(exp_name)){
  dir.create(exp_name)
  dir.create(file.path(exp_name, "figures"))
}

setup_df <- expand.grid(index=1,
                        G_DE1=50, G_DE2=50, G_nonDE1=4950, G_nonDE2=4950,
                        N_ctrl1=5, N_case1=5, N_ctrl2=5, N_case2=5,
                        #N_ctrl1=20, N_case1=20, N_ctrl2=20, N_case2=20,
                        mu1=10, 
                        mu2=c(5, 10, 13, 15, 20), 
                        bio_fold=2,
                        size1=20,
                        size2=c(10, 20, 40, 100))
setup_df$index <- 1:nrow(setup_df)
setup_df$N_1 <- setup_df$N_ctrl1 + setup_df$N_case1  
setup_df$N_2 <- setup_df$N_ctrl2 + setup_df$N_case2
setup_df$mu1_up <- setup_df$mu1 * setup_df$bio_fold 
setup_df$mu2_up <- setup_df$mu2 * setup_df$bio_fold
setup_df$var1_ctrl <- NBvar(mu=setup_df$mu1, size=setup_df$size1)
setup_df$var1_case <- NBvar(mu=setup_df$mu1_up, size=setup_df$size1)
setup_df$var2_ctrl <- NBvar(mu=setup_df$mu2, size=setup_df$size2)
setup_df$var2_case <- NBvar(mu=setup_df$mu2_up, size=setup_df$size2)
write.csv(setup_df, file=file.path(exp_name, "exp_setup.csv"))

# or with completed simultations, read in set up csv file and specify the row to study
# setup_df <- read.csv(file.path(exp_name, "exp_setup.csv"))
# ii=1
de1_ind <- 1:50
nonde1_ind <- 51:5000
de2_ind <- 5001:5050
nonde2_ind <- 5051:10000
  
stats_df_lst <- list()
for(ii in 1:nrow(setup_df)){
  curr_setup <- setup_df[ii, ]
  curr_prefix <- paste(paste0("mean", gsub('.', '', round(curr_setup$mu2 / curr_setup$mu1, 2), fixed=T)),
                       paste0("disp", gsub('.', '', round(curr_setup$size1 / curr_setup$size2, 2), fixed=T)), sep="_")
  
  #### Simulate count matrix and batch factor
  sim_res <- simCounts(curr_setup)
  counts <- sim_res$counts; batch <- sim_res$batch; group <- sim_res$group
  de_ind <- c(rep(1, curr_setup$G_DE1), rep(0, curr_setup$G_nonDE1), rep(1, curr_setup$G_DE2), rep(0, curr_setup$G_nonDE2))
  de_ground_truth <- rownames(counts)[which(de_ind==1)]
  
  
  #### Use ComBat-seq
  #covar_mod=NULL; full_mod=TRUE; gene.subset.n=5000; shrink.disp=FALSE; Cpp=FALSE
  start_time <- Sys.time()
  adj_counts <- ComBat_seq(counts=counts, batch=batch, group=group, shrink.disp=FALSE, gene.subset.n=1000)  #, Cpp=TRUE)
  end_time <- Sys.time()
  cat("\nRunning time of ComBat-seq:\n")
  print(end_time - start_time)
  
  
  ####  Differential Expression
  # de_called <- DESeq2_DEpipe(counts_mat=counts, batch=batch, group=group, include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.05)  
  # perf_stats <- perfStats(called_vec=de_called$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(counts))
  # perf_stats_fdr <- perfStats(called_vec=de_called$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(counts))
  # 
  # de_called_onestep <- DESeq2_DEpipe(counts_mat=counts, batch=batch, group=group, include.batch=TRUE, alpha.unadj=0.05, alpha.fdr=0.05)  
  # onestep_perf_stats <- perfStats(called_vec=de_called_onestep$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(counts))
  # onestep_perf_stats_fdr <- perfStats(called_vec=de_called_onestep$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(counts))
  # 
  # de_called_adj <- DESeq2_DEpipe(counts_mat=adj_counts, batch=batch, group=group, include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.05)  
  # adj_perf_stats <- perfStats(called_vec=de_called_adj$unadj, ground_truth_vec=de_ground_truth, N_genes=nrow(counts))
  # adj_perf_stats_fdr <- perfStats(called_vec=de_called_adj$fdr, ground_truth_vec=de_ground_truth, N_genes=nrow(counts))
  # 
  # stats_df <- data.frame(Unadj=perf_stats[1:2], OneStep=onestep_perf_stats[1:2], CombatSeq=adj_perf_stats[1:2])
  # stats_fdr_df <- data.frame(Unadj=perf_stats_fdr[c(1,3)], OneStep=onestep_perf_stats_fdr[c(1,3)], CombatSeq=adj_perf_stats_fdr[c(1,3)])
  # stats_fdr_df[2, ] <- 1 - stats_fdr_df[2, ]; rownames(stats_fdr_df)[2] <- "fdr"
  # # cat("Performance using un-adjusted P values\n")
  # # print(round(stats_df,2))
  # cat("Performance using FDR corrected p values\n")
  # print(round(stats_fdr_df,2))
  # stats_df_lst[[ii]] <- stats_fdr_df
  
  
  #### Clustering
  # hc <- hclust(dist(t(counts)))
  # dend <- as.dendrogram(hc)
  # dend <- color_branches(dend, groupLabels=batch[order.dendrogram(dend)], col=group[order.dendrogram(dend)]+2)
  # dendplt_ori <- plot(dend)  
  # 
  # adj_hc <- hclust(dist(t(adj_counts)))
  # adj_dend <- as.dendrogram(adj_hc)
  # adj_dend <- color_branches(adj_dend, groupLabels=batch[order.dendrogram(adj_dend)], col=group[order.dendrogram(adj_dend)]+2)
  # dendplt_adj <- plot(adj_dend) 
  # 
  # png(paste0(exp_name, "/figures/Dend_", curr_prefix, ".png"), width=9, height=5, units="in", res=300)
  # grid.arrange(dendplt_ori, dendplt_adj, ncol=2)
  # dev.off()
  
  
  #### PCA
  col_data <- data.frame(Batch=factor(batch), Group=factor(group)); rownames(col_data) <- colnames(counts)
  seobj <- SummarizedExperiment(assays=log2(counts+1), colData=col_data)
  pca_obj <- plotPCA(DESeqTransform(seobj), intgroup=c("Batch", "Group"))
  plt <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
    geom_point() + 
    labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
         y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
         title="Color by Batch - Unadjusted") +
    theme(legend.position="bottom")
  seobj_adj <- SummarizedExperiment(assays=log2(adj_counts+1), colData=col_data)
  pca_obj_adj <- plotPCA(DESeqTransform(seobj_adj), intgroup=c("Batch", "Group"))
  plt_adj <- ggplot(pca_obj_adj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
    geom_point() + 
    labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj$plot_env$percentVar[1])),
         y=sprintf("PC2: %s Variance", percent(pca_obj_adj$plot_env$percentVar[2])),
         title="Color by Batch - ComBat-Seq") +
    theme(legend.position="bottom")
  png(paste0(exp_name, "/figures/PCA_", curr_prefix, ".png"), width=9, height=5, units="in", res=300)
  grid.arrange(plt, plt_adj, ncol=2)
  dev.off()
  
  
  #### What's happening to actually counts?
  de_1_ori <- colMeans(counts[de1_ind, ])
  nonde_1_ori <- colMeans(counts[nonde1_ind, ])
  de_2_ori <- colMeans(counts[de2_ind, ])
  nonde_2_ori <- colMeans(counts[nonde2_ind, ])
  avg_gene_df_ori <- data.frame(de1=de_1_ori, nonde1=nonde_1_ori, de2=de_2_ori, nonde2=nonde_2_ori)
  avg_gene_df_ori$Samples <- factor(rownames(avg_gene_df_ori), levels=paste0("sample",1:sum(curr_setup$N_1, curr_setup$N_2)))
  avg_gene_df_ori_mlt <- melt(avg_gene_df_ori, id.vars="Samples", variable.name="Gene.Type")
  gene_plt_ori <- ggplot(avg_gene_df_ori_mlt, aes(x=Samples, y=value, group=Gene.Type, color=Gene.Type)) +
    geom_line() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(y="Average counts", title="Unadjusted")
  
  de_1_adj <- colMeans(adj_counts[de1_ind, ])
  nonde_1_adj <- colMeans(adj_counts[nonde1_ind, ])
  de_2_adj <- colMeans(adj_counts[de2_ind, ])
  nonde_2_adj <- colMeans(adj_counts[nonde2_ind, ])
  avg_gene_df_adj <- data.frame(de1=de_1_adj, nonde1=nonde_1_adj, de2=de_2_adj, nonde2=nonde_2_adj)
  avg_gene_df_adj$Samples <- factor(rownames(avg_gene_df_adj), levels=paste0("sample",1:sum(curr_setup$N_1, curr_setup$N_2)))
  avg_gene_df_adj_mlt <- melt(avg_gene_df_adj, id.vars="Samples", variable.name="Gene.Type")
  gene_plt_adj <- ggplot(avg_gene_df_adj_mlt, aes(x=Samples, y=value, group=Gene.Type, color=Gene.Type)) +
    geom_line() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(y="Average counts", title="ComBat-Seq")
  
  png(paste0(exp_name, "/figures/avg_cts_", curr_prefix, ".png"), width=7, height=8, units="in", res=300)
  grid.arrange(gene_plt_ori, gene_plt_adj, nrow=2)
  dev.off()
}
