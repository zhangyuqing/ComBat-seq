NBvar <- function(mu, size){
  return(mu + mu^2/size)
}


simCounts <- function(curr_setup, sim.outliers=FALSE){
  y1_ctrl <- matrix(rnbinom(curr_setup$G_DE1 * curr_setup$N_ctrl1, mu=curr_setup$mu1, size=curr_setup$size1), ncol=curr_setup$N_ctrl1)
  y1_case <- matrix(rnbinom(curr_setup$G_DE1 * curr_setup$N_case1, mu=curr_setup$mu1_up, size=curr_setup$size1), ncol=curr_setup$N_case1)
  y1_null <- matrix(rnbinom(curr_setup$G_nonDE1 * curr_setup$N_1, mu=curr_setup$mu1, size=curr_setup$size1), ncol=curr_setup$N_1)
  
  y2_ctrl <- matrix(rnbinom(curr_setup$G_DE1 * curr_setup$N_ctrl2, mu=curr_setup$mu2, size=curr_setup$size2), ncol=curr_setup$N_ctrl2)
  y2_case <- matrix(rnbinom(curr_setup$G_DE1 * curr_setup$N_case2, mu=curr_setup$mu2_up, size=curr_setup$size2), ncol=curr_setup$N_case2)
  y2_null <- matrix(rnbinom(curr_setup$G_nonDE1 * curr_setup$N_2, mu=curr_setup$mu2, size=curr_setup$size2), ncol=curr_setup$N_2)
  
  y3_ctrl <- matrix(rnbinom(curr_setup$G_DE2 * curr_setup$N_ctrl1, mu=curr_setup$mu2, size=curr_setup$size1), ncol=curr_setup$N_ctrl1)
  y3_case <- matrix(rnbinom(curr_setup$G_DE2 * curr_setup$N_case1, mu=curr_setup$mu2_up, size=curr_setup$size1), ncol=curr_setup$N_case1)
  y3_null <- matrix(rnbinom(curr_setup$G_nonDE2 * curr_setup$N_1, mu=curr_setup$mu2, size=curr_setup$size1), ncol=curr_setup$N_1)
  
  y4_ctrl <- matrix(rnbinom(curr_setup$G_DE2 * curr_setup$N_ctrl2, mu=curr_setup$mu1, size=curr_setup$size2), ncol=curr_setup$N_ctrl2)
  y4_case <- matrix(rnbinom(curr_setup$G_DE2 * curr_setup$N_case2, mu=curr_setup$mu1_up, size=curr_setup$size2), ncol=curr_setup$N_case2)
  y4_null <- matrix(rnbinom(curr_setup$G_nonDE2 * curr_setup$N_2, mu=curr_setup$mu1, size=curr_setup$size2), ncol=curr_setup$N_2)
  
  ## simulate outliers
  if(sim.outliers){
    nG <- curr_setup$G_DE1 + curr_setup$G_nonDE1 + curr_setup$G_DE2 + curr_setup$G_nonDE2
    is_outlier_gene <- sample(c(TRUE, FALSE), nG, prob=c(curr_setup$outgprob, 1-curr_setup$outgprob), replace=TRUE)
    outlier_gene_ind <- which(is_outlier_gene)
    
    gene_groups <- c(rep("DE1", curr_setup$G_DE1), rep("nonDE1", curr_setup$G_nonDE1),
                     rep("DE2", curr_setup$G_DE2), rep("nonDE2", curr_setup$G_nonDE2))
    
    y1_ctrl <- simOutliers(y1_ctrl, intersect(outlier_gene_ind, which(gene_groups=="DE1")), 
                           curr_setup$outalpha, c(mu=curr_setup$mu1, size=curr_setup$size1))
    y1_case <- simOutliers(y1_case, intersect(outlier_gene_ind, which(gene_groups=="DE1")), 
                           curr_setup$outalpha, c(mu=curr_setup$mu1_up, size=curr_setup$size1))
    y1_null <- simOutliers(y1_null, intersect(outlier_gene_ind, which(gene_groups=="nonDE1"))-nrow(y1_ctrl), 
                           curr_setup$outalpha, c(mu=curr_setup$mu1, size=curr_setup$size1))
    
    y2_ctrl <- simOutliers(y2_ctrl, intersect(outlier_gene_ind, which(gene_groups=="DE1")), 
                           curr_setup$outalpha, c(mu=curr_setup$mu2, size=curr_setup$size2))
    y2_case <- simOutliers(y2_case, intersect(outlier_gene_ind, which(gene_groups=="DE1")), 
                           curr_setup$outalpha, c(mu=curr_setup$mu2_up, size=curr_setup$size2))
    y2_null <- simOutliers(y2_null, intersect(outlier_gene_ind, which(gene_groups=="nonDE1"))-nrow(y2_ctrl), 
                           curr_setup$outalpha, c(mu=curr_setup$mu2, size=curr_setup$size2))
    
    y3_ctrl <- simOutliers(y3_ctrl, intersect(outlier_gene_ind, which(gene_groups=="DE2"))-nrow(y1_ctrl)-nrow(y1_null), 
                           curr_setup$outalpha, c(mu=curr_setup$mu2, size=curr_setup$size1))
    y3_case <- simOutliers(y3_case, intersect(outlier_gene_ind, which(gene_groups=="DE2"))-nrow(y1_ctrl)-nrow(y1_null), 
                           curr_setup$outalpha, c(mu=curr_setup$mu2_up, size=curr_setup$size1))
    y3_null <- simOutliers(y3_null, intersect(outlier_gene_ind, which(gene_groups=="nonDE2"))-nrow(y1_ctrl)-nrow(y1_null)-nrow(y3_ctrl), 
                           curr_setup$outalpha, c(mu=curr_setup$mu2, size=curr_setup$size1))
    
    y4_ctrl <- simOutliers(y4_ctrl, intersect(outlier_gene_ind, which(gene_groups=="DE2"))-nrow(y2_ctrl)-nrow(y2_null), 
                           curr_setup$outalpha, c(mu=curr_setup$mu1, size=curr_setup$size2))
    y4_case <- simOutliers(y4_case, intersect(outlier_gene_ind, which(gene_groups=="DE2"))-nrow(y2_ctrl)-nrow(y2_null), 
                           curr_setup$outalpha, c(mu=curr_setup$mu1_up, size=curr_setup$size2))
    y4_null <- simOutliers(y4_null, intersect(outlier_gene_ind, which(gene_groups=="nonDE2"))-nrow(y2_ctrl)-nrow(y2_null)-nrow(y4_ctrl), 
                           curr_setup$outalpha, c(mu=curr_setup$mu1, size=curr_setup$size2))
  }else{
    outlier_gene_ind <- NULL
  }
  
  counts <- rbind(cbind(y1_ctrl, y1_case, y2_ctrl, y2_case), cbind(y1_null, y2_null),
                  cbind(y3_ctrl, y3_case, y4_ctrl, y4_case), cbind(y3_null, y4_null))
  rownames(counts) <- paste0("gene", 1:nrow(counts))
  colnames(counts) <- paste0("sample", 1:ncol(counts))
  
  batch <- c(rep("B", ncol(y1_ctrl)+ncol(y1_case)), rep("A", ncol(y2_ctrl)+ncol(y2_case))); table(batch)
  group <- c(rep(0, ncol(y1_ctrl)), rep(1, ncol(y1_case)), rep(0, ncol(y2_ctrl)), rep(1, ncol(y2_case)))
  
  return(list(counts=counts, batch=batch, group=group, outlier_gene_ind=outlier_gene_ind))
}


#cts_block=y4_null; outlier_genes=intersect(outlier_gene_ind, which(gene_groups=="nonDE2"))-nrow(y2_ctrl)-nrow(y2_null)-nrow(y4_ctrl); 
#outalpha=curr_setup$outalpha; params=c(mu=curr_setup$mu1, size=curr_setup$size2)
simOutliers <- function(cts_block, outlier_genes, outalpha, params){
  for(g in outlier_genes){
    # select a sample to be the outlier
    outlier_sample <- sample(1:ncol(cts_block), 1)
    cts_block[g, outlier_sample] <- qnbinom(outalpha, mu=params["mu"], size=params["size"], lower.tail=FALSE)
  }
  return(cts_block)
}


DEpipe <- function(counts_mat, batch, group, ground_truth_vec, include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.05, DE_method=DESeq2_DEpipe, covar=NULL){
  de_called <- DE_method(counts_mat, batch, group, include.batch=include.batch, 
                         alpha.unadj=alpha.unadj, alpha.fdr=alpha.fdr, covar=covar)
  perf_stats <- perfStats(called_vec=de_called$unadj, ground_truth_vec=ground_truth_vec, N_genes=nrow(counts_mat))
  perf_stats_fdr <- perfStats(called_vec=de_called$fdr, ground_truth_vec=ground_truth_vec, N_genes=nrow(counts_mat))
  return(list(perf_stats=perf_stats, perf_stats_fdr=perf_stats_fdr))
}


edgeR_DEpipe <- function(counts_mat, batch, group, include.batch, alpha.unadj, alpha.fdr, covar=NULL){
  cat("DE tool: edgeR\n")
  y <- DGEList(counts=counts_mat)
  y <- calcNormFactors(y, method="TMM")
  if(include.batch){
    cat("Including batch as covariate\n")
    design <- model.matrix(~ as.factor(group) + as.factor(batch))
  }else{
    cat("Default group as model matrix\n")
    design <- model.matrix(~as.factor(group))
  }
  if(!is.null(covar)){
    cat("Including surrogate variables or unwanted variation variables\n")
    design <- cbind(design, covar)
  }
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef=2)
  de_res <- topTags(qlf, n=nrow(counts_mat))$table
  
  de_called <- rownames(de_res)[de_res$PValue < alpha.unadj]
  de_called_fdr <- rownames(de_res)[de_res$FDR < alpha.fdr]
  return(list(unadj=de_called, fdr=de_called_fdr, de_res=de_res, design=design))
}


DESeq2_DEpipe <- function(counts_mat, batch, group, include.batch, alpha.unadj, alpha.fdr, covar=NULL){
  cat("DE tool: DESeq2\n")
  if(include.batch){
    cat("Including batch as covariate\n")
    col_data <- data.frame(Batch=as.factor(batch), Group=as.factor(group))
    design_formula <- ~Batch+Group
  }else if(!is.null(covar)){
    cat("Including surrogate variables or unwanted variation variables\n")
    col_data <- cbind(model.matrix(~as.factor(group)), covar)
    colnames(col_data)[2:ncol(col_data)] <- c("Group", paste0("Covar", 1:(ncol(col_data)-2)))
    rownames(col_data) <- colnames(counts_mat)
    col_data <- as.data.frame(col_data); col_data$Group <- as.factor(col_data$Group); col_data <- col_data[, -1]
    design_formula <- as.formula(paste("~", paste(colnames(col_data), collapse="+")))
  }else{
    cat("Default group as model matrix\n")
    col_data <- data.frame(Group=as.factor(group))
    design_formula <- ~Group
  }
  
  dds <- DESeqDataSetFromMatrix(countData=counts_mat, colData=col_data, design=design_formula)
  dds <- DESeq(dds, quiet=TRUE)
  de_res <- results(dds, name="Group_1_vs_0")
  
  if(length(which(de_res$pvalue < alpha.unadj))>0){
    de_called <- rownames(de_res)[which(de_res$pvalue < alpha.unadj)]
  }else{de_called <- character(0)}
  if(length(which(de_res$padj < alpha.fdr))>0){
    de_called_fdr <- rownames(de_res)[which(de_res$padj < alpha.fdr)]
  }else{de_called_fdr <- character(0)}
  
  return(list(unadj=de_called, fdr=de_called_fdr, de_res=de_res, 
              design_formula=design_formula, col_data=col_data))
}


currComBat_lm_DEpipe <- function(cts, batch, group, alpha.unadj, alpha.fdr){
  log_counts <- cpm(cts, log=TRUE)  # use logCPM to make data more normal
  adj_counts <- ComBat(log_counts, batch=batch, mod=model.matrix(~as.factor(group)))
  pval_seq <- apply(adj_counts, 1, function(x, group){
    x_norm <- scale(x, center=TRUE, scale=TRUE)
    fit3 <- lm(x_norm ~ as.factor(group))
    return(summary(fit3)$coefficients[2, 4])
  }, group=group)
  padj_seq <- p.adjust(pval_seq, method="fdr")
  
  de_called <- list(unadj=rownames(cts)[pval_seq < alpha.unadj], fdr=rownames(cts)[padj_seq < alpha.fdr], 
                    de_res=data.frame(PValue=pval_seq, FDR=padj_seq), design=model.matrix(~as.factor(group)))
  return(de_called)
}


perfStats <- function(called_vec, ground_truth_vec, N_genes){
  if(length(called_vec)==0){
    tpr <- fpr <- 0
    prec <- NA
  }else{
    tp <- length(intersect(called_vec, ground_truth_vec))
    fp <- length(setdiff(called_vec, ground_truth_vec)) 
    N_DE <- length(ground_truth_vec)
    N_nonDE <- N_genes - N_DE
    tpr <- tp / N_DE
    fpr <- fp / N_nonDE
    prec <- tp / length(called_vec)
  }
  return(c(tpr=tpr, fpr=fpr, prec=prec))
}


cancelLibsizeEffect <- function(count_matrix){
  lib_sizes <- colSums(count_matrix)
  do.call(cbind, lapply(1:ncol(count_matrix), function(i){
    count_matrix[,i]/lib_sizes[i]
  }))
}


PCAPlt <- function(cts, col.data, plt.title){
  seobj <- SummarizedExperiment(assays=cts, colData=col.data)
  pca_obj <- plotPCA(DESeqTransform(seobj), intgroup=c("Batch", "Group"))
  plt <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
    geom_point() + 
    labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
         y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
         title=sprintf("Color by Batch - %s", plt.title)) +
    theme(legend.position="bottom")
  return(plt)
}


AvgExpPlt <- function(cts, de.lst, curr.setup, plt.title){
  de_1 <- colMeans(cts[de.lst$de1_ind, ])
  nonde_1 <- colMeans(cts[de.lst$nonde1_ind, ])
  de_2 <- colMeans(cts[de.lst$de2_ind, ])
  nonde_2 <- colMeans(cts[de.lst$nonde2_ind, ])
  avg_gene_df <- data.frame(de1=de_1, nonde1=nonde_1, de2=de_2, nonde2=nonde_2)
  avg_gene_df$Samples <- factor(rownames(avg_gene_df), levels=paste0("sample",1:sum(curr.setup$N_1, curr.setup$N_2)))
  avg_gene_df_mlt <- melt(avg_gene_df, id.vars="Samples", variable.name="Gene.Type")
  gene_plt <- ggplot(avg_gene_df_mlt, aes(x=Samples, y=value, group=Gene.Type, color=Gene.Type)) +
    geom_line() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(y="Average counts", title=plt.title)
  return(gene_plt)
}


ShrinkParamPlt <- function(cts, combatseq_res, batch_name, eb_viz_genes=NULL, outlier_gene_ind){
  is_outlier_gene <- rep(FALSE, nrow(cts)); is_outlier_gene[outlier_gene_ind] <- TRUE
  ebON_outlier_type <- rep("EB On", nrow(cts)); ebON_outlier_type[is_outlier_gene] <- "EB On - outliers"
  ebOFF_outlier_type <- rep("EB Off", nrow(cts)); ebOFF_outlier_type[is_outlier_gene] <- "EB Off - outliers"
  
  if(!is.null(eb_viz_genes)){
    # params_est <- rbind(data.frame(genes=rownames(cts)[eb_viz_genes], 
    #                                gamma=combatseq_res$gamma_hat[eb_viz_genes, paste0('batch',batch_name)], 
    #                                phi=combatseq_res$phi_hat[eb_viz_genes, paste0('batch',batch_name)], 
    #                                EB.type="EB Off",
    #                                is.outlier=is_outlier_gene[eb_viz_genes]),
    #                     data.frame(genes=rownames(cts)[eb_viz_genes], 
    #                                gamma=combatseq_res$gamma_star[eb_viz_genes, paste0('batch',batch_name)], 
    #                                phi=combatseq_res$phi_star[eb_viz_genes, paste0('batch',batch_name)], 
    #                                EB.type="EB On",
    #                                is.outlier=is_outlier_gene[eb_viz_genes]))
    params_est <- rbind(data.frame(genes=rownames(cts)[eb_viz_genes],
                                   gamma=combatseq_res$gamma_hat[eb_viz_genes, paste0('batch',batch_name)],
                                   phi=combatseq_res$phi_hat[eb_viz_genes, paste0('batch',batch_name)],
                                   EB.outlier.type=ebOFF_outlier_type[eb_viz_genes]),
                        data.frame(genes=rownames(cts)[eb_viz_genes],
                                   gamma=combatseq_res$gamma_star[eb_viz_genes, paste0('batch',batch_name)],
                                   phi=combatseq_res$phi_star[eb_viz_genes, paste0('batch',batch_name)],
                                   EB.outlier.type=ebON_outlier_type[eb_viz_genes]))
  }else{
    params_est <- rbind(data.frame(genes=rownames(cts), 
                                   gamma=combatseq_res$gamma_hat[, paste0('batch',batch_name)], 
                                   phi=combatseq_res$phi_hat[, paste0('batch',batch_name)], 
                                   EB.outlier.type=ebOFF_outlier_type),
                        data.frame(genes=rownames(cts), 
                                   gamma=combatseq_res$gamma_star[, paste0('batch',batch_name)], 
                                   phi=combatseq_res$phi_star[, paste0('batch',batch_name)], 
                                   EB.outlier.type=ebON_outlier_type))
  }
  
  plt <- ggplot(params_est) +
    geom_point(aes(x=gamma, y=phi, color=EB.outlier.type)) +
    scale_color_manual(values=c("light green", "orange", "dark green", "red")) + 
    labs(x="Mean batch effect estimates (gamma)", y="Dispersion estimates (phi)", 
         title=sprintf("Parameter estimates within batch %s", batch_name))
  if(!is.null(eb_viz_genes)){
    plt <- plt + geom_line(aes(x=gamma, y=phi, group=genes), color="grey")
  }
  return(plt)
}


ShrinkParamPlt_newEB <- function(cts, combatseq_res, batch_name, eb_viz_genes=NULL, outlier_gene_ind){
  is_outlier_gene <- rep(FALSE, nrow(cts)); is_outlier_gene[outlier_gene_ind] <- TRUE
  ebON_outlier_type <- rep("EB On", nrow(cts)); ebON_outlier_type[is_outlier_gene] <- "EB On - outliers"
  ebOFF_outlier_type <- rep("EB Off", nrow(cts)); ebOFF_outlier_type[is_outlier_gene] <- "EB Off - outliers"
  
  # outlier_detected <- union(combatseq_res$gamma_outliers, combatseq_res$phi_outliers)
  outlier_detected <- union(combatseq_res$gamma_outliers[[paste0("batch",batch_name)]], combatseq_res$phi_outliers[[paste0("batch",batch_name)]])
  outlier_detected_ind <- rownames(cts) %in% outlier_detected
  
  if(!is.null(eb_viz_genes)){
    params_est <- rbind(data.frame(genes=rownames(cts)[eb_viz_genes],
                                   gamma=combatseq_res$gamma_hat[eb_viz_genes, paste0('batch',batch_name)],
                                   phi=combatseq_res$phi_hat[eb_viz_genes, paste0('batch',batch_name)],
                                   EB.outlier.type=ebOFF_outlier_type[eb_viz_genes],
                                   Detected=outlier_detected_ind[eb_viz_genes]),
                        data.frame(genes=rownames(cts)[eb_viz_genes],
                                   gamma=combatseq_res$gamma_star[eb_viz_genes, paste0('batch',batch_name)],
                                   phi=combatseq_res$phi_star[eb_viz_genes, paste0('batch',batch_name)],
                                   EB.outlier.type=ebON_outlier_type[eb_viz_genes],
                                   Detected=outlier_detected_ind[eb_viz_genes]))
  }else{
    params_est <- rbind(data.frame(genes=rownames(cts), 
                                   gamma=combatseq_res$gamma_hat[, paste0('batch',batch_name)], 
                                   phi=combatseq_res$phi_hat[, paste0('batch',batch_name)], 
                                   EB.outlier.type=ebOFF_outlier_type,
                                   Detected=outlier_detected_ind),
                        data.frame(genes=rownames(cts), 
                                   gamma=combatseq_res$gamma_star[, paste0('batch',batch_name)], 
                                   phi=combatseq_res$phi_star[, paste0('batch',batch_name)], 
                                   EB.outlier.type=ebON_outlier_type,
                                   Detected=outlier_detected_ind))
  }
  
  plt <- ggplot(params_est) +
    geom_point(aes(x=gamma, y=phi, color=EB.outlier.type, shape=Detected)) +
    scale_color_manual(values=c("light green", "orange", "dark green", "red")) + 
    labs(x="Mean batch effect estimates (gamma)", y="Dispersion estimates (phi)", 
         title=sprintf("Parameter estimates within batch %s", batch_name)) #+
    #theme(legend.position="bottom")
  if(!is.null(eb_viz_genes)){
    plt <- plt + geom_line(aes(x=gamma, y=phi, group=genes), color="grey")
  }
  return(plt)
}
