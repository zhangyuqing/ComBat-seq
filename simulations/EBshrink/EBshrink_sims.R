rm(list=ls()); run.on.SCC <- FALSE
sapply(c("ggplot2", "reshape2", "gridExtra", "dendextend", "edgeR", "DESeq2", "scales"), require, character.only=TRUE)
if(run.on.SCC){
  setwd("")
  source("../../ComBat_seq_returnParams.R")
  source("../../helper_seq.R")
  source("EBshrink_helpers.R")
}else{
  source("../../ComBat_seq_returnParams.R")
  source("../../helper_seq.R")
  source("EBshrink_helpers.R")
}
set.seed(123)


#### Experiment set up
exp_name <- "outlier_probs_N20" #sprintf("outlier_prob%s", gsub(".", "", p, fixed=TRUE))
if(!dir.exists(exp_name)){
  dir.create(exp_name)
  dir.create(file.path(exp_name, "figures"))
  dir.create(file.path(exp_name, "rdata"))
}

eb_viz_genes <- 1:100  

setup_df <- expand.grid(index=1,
                        G_DE1=50, G_DE2=50, G_nonDE1=950, G_nonDE2=950,
                        #N_ctrl1=5, N_case1=5, N_ctrl2=5, N_case2=5,
                        N_ctrl1=20, N_case1=20, N_ctrl2=20, N_case2=20,
                        mu1=10, 
                        mu2=15, #c(5, 10, 13, 15, 20), 
                        bio_fold=2,
                        size1=20,
                        size2=40, #c(10, 20, 40, 100))
                        outgprob=seq(0.1, 0.5, 0.1),
                        outalpha=10^c(-2,-4,-8),
                        simout=c(TRUE, FALSE))  # fraction of genes containing outliers
setup_df$index <- 1:nrow(setup_df)
setup_df$N_1 <- setup_df$N_ctrl1 + setup_df$N_case1  
setup_df$N_2 <- setup_df$N_ctrl2 + setup_df$N_case2
setup_df$mu1_up <- setup_df$mu1 * setup_df$bio_fold 
setup_df$mu2_up <- setup_df$mu2 * setup_df$bio_fold
setup_df$var1_ctrl <- NBvar(mu=setup_df$mu1, size=setup_df$size1)
setup_df$var1_case <- NBvar(mu=setup_df$mu1_up, size=setup_df$size1)
setup_df$var2_ctrl <- NBvar(mu=setup_df$mu2, size=setup_df$size2)
setup_df$var2_case <- NBvar(mu=setup_df$mu2_up, size=setup_df$size2)
write.csv(setup_df, file=file.path(exp_name, "exp_setup.csv"), row.names=FALSE)
# or with completed simultations, read in set up csv file and specify the row to study
# setup_df <- read.csv(file.path(exp_name, "exp_setup.csv"))

de1_ind <- 1:50
nonde1_ind <- 51:1000
de2_ind <- 1001:1050
nonde2_ind <- 1051:2000
de_lst <- list(de1_ind=de1_ind, nonde1_ind=nonde1_ind, de2_ind=de2_ind, nonde2_ind=nonde2_ind)


stats_df_lst <- list()
# ii=27
for(ii in 1:nrow(setup_df)){
  curr_setup <- setup_df[ii, ]
  curr_prefix <- paste(paste0("mean", gsub('.', '', round(curr_setup$mu2 / curr_setup$mu1, 2), fixed=T)),
                       paste0("disp", gsub('.', '', round(curr_setup$size1 / curr_setup$size2, 2), fixed=T)), 
                       paste0("outgprob", gsub('.', '', curr_setup$outgprob, fixed=T)),
                       paste0("outalpha", curr_setup$outalpha), sep="_")
  if(as.logical(curr_setup$simout)){curr_prefix <- paste0(curr_prefix, "_simout")}
  
  
  #### Simulate count matrix with outliers and batch factor
  sim_res <- simCounts(curr_setup, sim.outliers=as.logical(curr_setup$simout))
  counts <- sim_res$counts; batch <- sim_res$batch; group <- sim_res$group; outlier_gene_ind <- sim_res$outlier_gene_ind
  de_ind <- c(rep(1, curr_setup$G_DE1), rep(0, curr_setup$G_nonDE1), rep(1, curr_setup$G_DE2), rep(0, curr_setup$G_nonDE2))
  de_ground_truth <- rownames(counts)[which(de_ind==1)]
  
  
  #### Use ComBat-seq
  start_time <- Sys.time()
  combatseq_res_ebON <- ComBat_seq(counts=counts, batch=batch, group=group, shrink=TRUE, shrink.disp=TRUE, gene.subset.n=1000)  
  end_time <- Sys.time()
  cat("\nRunning time of ComBat-seq (EB on):\n")
  print(end_time - start_time)
  
  start_time <- Sys.time()
  combatseq_res_ebOFF <- ComBat_seq(counts=counts, batch=batch, group=group, shrink=FALSE, shrink.disp=FALSE) 
  end_time <- Sys.time()
  cat("\nRunning time of ComBat-seq (EB off):\n")
  print(end_time - start_time)
  
  
  ####  Differential Expression
  de_res <- DEpipe(counts_mat=counts, batch=batch, group=group, ground_truth_vec=de_ground_truth)
  de_onestep <- DEpipe(counts_mat=counts, batch=batch, group=group, ground_truth_vec=de_ground_truth, include.batch=TRUE)
  de_adj_ebON <- DEpipe(counts_mat=combatseq_res_ebON$adjust_counts, batch=batch, group=group, ground_truth_vec=de_ground_truth)
  de_adj_ebOFF <- DEpipe(counts_mat=combatseq_res_ebOFF$adjust_counts, batch=batch, group=group, ground_truth_vec=de_ground_truth)

  stats_df <- data.frame(Unadj=de_res$perf_stats[1:2], OneStep=de_onestep$perf_stats[1:2],
                         CombatSeq.ebON=de_adj_ebON$perf_stats[1:2], CombatSeq.ebFF=de_adj_ebOFF$perf_stats[1:2])
  stats_fdr_df <- data.frame(Unadj=de_res$perf_stats_fdr[c(1,3)], OneStep=de_onestep$perf_stats_fdr[c(1,3)],
                             CombatSeq.ebON=de_adj_ebON$perf_stats_fdr[c(1,3)], CombatSeq.ebFF=de_adj_ebOFF$perf_stats_fdr[c(1,3)])
  stats_fdr_df[2, ] <- 1 - stats_fdr_df[2, ]; rownames(stats_fdr_df)[2] <- "fdr"
  # cat("Performance using un-adjusted P values\n")
  # print(round(stats_df,2))
  # cat("Performance using FDR corrected p values\n")
  # print(round(stats_fdr_df,2))
  stats_df_lst[[ii]] <- stats_fdr_df
  
  
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
  
  
  #### Scatter plot of dispersion vs mean estimates
  eb_plts <- lapply(sort(unique(batch)), function(b){ShrinkParamPlt(counts, combatseq_res_ebON, b, eb_viz_genes=NULL, outlier_gene_ind)})
  png(paste0(exp_name, "/figures/ShrinkParams_", curr_prefix, ".png"),
      width=14, height=6, units="in", res=300)
  grid.arrange(eb_plts[[1]], eb_plts[[2]], ncol=2)
  dev.off()
  
  eb_plts_sub <- lapply(sort(unique(batch)), function(b){ShrinkParamPlt(counts, combatseq_res_ebON, b, eb_viz_genes, outlier_gene_ind)})
  png(paste0(exp_name, "/figures/ShrinkParams_subgenes_", curr_prefix, ".png"), 
      width=14, height=6, units="in", res=300)
  grid.arrange(eb_plts_sub[[1]], eb_plts_sub[[2]], ncol=2)
  dev.off()
  
  
  #### PCA
  col_data <- data.frame(Batch=factor(batch), Group=factor(group)) 
  rownames(col_data) <- colnames(counts)
  plt <- PCAPlt(counts, col_data, "Unadjusted")
  plt_adj_ebON <- PCAPlt(combatseq_res_ebON$adjust_counts, col_data, "ComBat-Seq (EB on)")
  plt_adj_ebOFF <- PCAPlt(combatseq_res_ebOFF$adjust_counts, col_data, "ComBat-Seq (EB off)")
  png(paste0(exp_name, "/figures/PCA_", curr_prefix, ".png"), width=14, height=5, units="in", res=300)
  grid.arrange(plt, plt_adj_ebON, plt_adj_ebOFF, ncol=3)
  dev.off()
  
  
  #### What's happening to actually counts?
  gene_plt_ori <- AvgExpPlt(counts, de_lst, curr_setup, "Unadjusted")
  gene_plt_adj_ebON <- AvgExpPlt(combatseq_res_ebON$adjust_counts, de_lst, curr_setup, "ComBat-Seq (EB on)")
  gene_plt_adj_ebOFF <- AvgExpPlt(combatseq_res_ebOFF$adjust_counts, de_lst, curr_setup, "ComBat-Seq (EB off)")
  png(paste0(exp_name, "/figures/AvgCts_", curr_prefix, ".png"), width=7, height=12, units="in", res=300)
  grid.arrange(gene_plt_ori, gene_plt_adj_ebON, gene_plt_adj_ebOFF, nrow=3)
  dev.off()
  
  
  #### Save plot objects
  pca_plts=list(unadjusted=plt, combatseq.ebON=plt_adj_ebON, combatseq.ebOFF=plt_adj_ebOFF)
  gene_plts=list(unadjusted=gene_plt_ori, combatseq.ebON=gene_plt_adj_ebON, combatseq.ebOFF=gene_plt_adj_ebOFF)
  save(eb_plts, eb_plts_sub, pca_plts, gene_plts, file=sprintf("./%s/rdata/pltObjs_%s.RData", exp_name, curr_prefix))
}

saveRDS(stats_df_lst, sprintf("./%s/rdata/de_stats.rds", exp_name))
