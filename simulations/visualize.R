rm(list=ls())
sapply(c("ggplot2", "gridExtra", "reshape2", "ggpubr"), require, character.only=TRUE)
results_dir <- "~/Documents/ComBat_seq/DE_fullpipe_complete/"  # path to the CSV result files generated from sim_DEpipe.R
source("~/Desktop/ComBat-seq/simulations/visualize_helpers.R")  # path to visualize_helpers.R

N_samples <- 20
mean_level_vec <- c(1, 1.5, 2, 3)   # level of mean batch effect to show
disp_level_vec <- 1:4   # level of dispersion batch effect to show
confounding_level_vec <- 0.5  # level of confounding
coverage_level <- 5  # sequencing coverage
sel_method_list <- list(edgeR=c("Base.edgeR", "Batch.edgeR", "OneStep.edgeR", "ComBat.lm", 
                                "RUVseq.edgeR", "SVAseq.edgeR", "ComBatseq.ebOFF.edgeR"),
                        DESeq2=c("Base.DESeq2", "Batch.DESeq2", "OneStep.DESeq2", "ComBat.lm", 
                                 "RUVseq.DESeq2", "SVAseq.DESeq2", "ComBatseq.ebOFF.DESeq2"))
sim_outliers <- FALSE

#i=1;j=1
fpr_lst_edgeR <- fpr_lst_DESeq2 <- tpr_lst_edgeR <- tpr_lst_DESeq2 <- list()  
for(i in seq_along(mean_level_vec)){
  fpr_lst_edgeR[[i]] <- fpr_lst_DESeq2[[i]] <- tpr_lst_edgeR[[i]] <- tpr_lst_DESeq2[[i]] <- list()
  for(j in seq_along(disp_level_vec)){
    ## Read in result files
    exp_name <- paste0("sim", ifelse(sim_outliers, "_simout", ""), "_N", N_samples, 
                       "_meanFC", gsub(".", "", mean_level_vec[i], fixed=T), "_dispFC", disp_level_vec[j],
                       "_cnfnd", gsub(".","", confounding_level_vec[1], fixed=T), "_depth", coverage_level)
    tpr_res <- read.csv(file.path(results_dir, sprintf('tpr_%s.csv', exp_name)))
    fpr_res <- read.csv(file.path(results_dir, sprintf('fpr_%s.csv', exp_name)))
    
    # split DE algorithm used (edgeR & DESeq2)
    tpr_res_edgeR <- tpr_res[, sel_method_list$edgeR]
    fpr_res_edgeR <- fpr_res[, sel_method_list$edgeR]
    tpr_res_DESeq2 <- tpr_res[, sel_method_list$DESeq2]
    fpr_res_DESeq2 <- fpr_res[, sel_method_list$DESeq2]
    
    fpr_lst_edgeR[[i]][[j]] <- as.list(colMeans(fpr_res_edgeR, na.rm=TRUE))
    tpr_lst_edgeR[[i]][[j]] <- as.list(colMeans(tpr_res_edgeR, na.rm=TRUE))
    fpr_lst_DESeq2[[i]][[j]] <- as.list(colMeans(fpr_res_DESeq2, na.rm=TRUE))
    tpr_lst_DESeq2[[i]][[j]] <- as.list(colMeans(tpr_res_DESeq2, na.rm=TRUE))
  }
  names(fpr_lst_edgeR[[i]]) <- names(tpr_lst_edgeR[[i]]) <- paste0("Disp", disp_level_vec)
  names(fpr_lst_DESeq2[[i]]) <- names(tpr_lst_DESeq2[[i]]) <- paste0("Disp", disp_level_vec)
}
names(fpr_lst_edgeR) <- names(tpr_lst_edgeR) <- paste0("Mean", mean_level_vec)
names(fpr_lst_DESeq2) <- names(tpr_lst_DESeq2) <- paste0("Mean", mean_level_vec)

tpr_edgeR_merged <- CleanOut(tpr_lst_edgeR, "TPR", ".edgeR")
fpr_edgeR_merged <- CleanOut(fpr_lst_edgeR, "FPR", ".edgeR")
tpr_DESeq2_merged <- CleanOut(tpr_lst_DESeq2, "TPR", ".DESeq2")
fpr_DESeq2_merged <- CleanOut(fpr_lst_DESeq2, "FPR", ".DESeq2")

if(identical(fpr_edgeR_merged[,-1], tpr_edgeR_merged[,-1])){
  perfstats_edgeR <- data.frame(TPR=tpr_edgeR_merged$TPR, FPR=fpr_edgeR_merged$FPR, tpr_edgeR_merged[,-1])
}else{stop("ERROR in merging data frames! - edgeR")}
#identical(perfstats_edgeR[,-1], fpr_edgeR_merged)
if(identical(fpr_DESeq2_merged[,-1], tpr_DESeq2_merged[,-1])){
  perfstats_DESeq2 <- data.frame(TPR=tpr_DESeq2_merged$TPR, FPR=fpr_DESeq2_merged$FPR, tpr_DESeq2_merged[,-1])
}else{stop("ERROR in merging data frames! - DESeq2")}
#identical(perfstats_DESeq2[,-1], fpr_DESeq2_merged)

perfstats_edgeR$Method <- plyr::revalue(perfstats_edgeR$Method, 
                                        c("Base"="Data without batch effect", 
                                          "Batch"="With batch effect, no adjustment",
                                          "OneStep"="Including batch as covariate in DE model",
                                          "ComBat.lm"="Adjusted by original ComBat on logCPM",
                                          "RUVseq"="Adjusted by RUVSeq",
                                          "SVAseq"="Adjusted by SVASeq",
                                          "ComBatseq.ebOFF"="Adjusted by ComBat-Seq"))
perfstats_edgeR$Mean <- factor(perfstats_edgeR$Mean, levels=c("Mean1", "Mean1.5", "Mean2", "Mean3"))
perfstats_edgeR$Mean <- plyr::revalue(perfstats_edgeR$Mean, 
                                      c("Mean1"="No mean\nbatch effect",
                                        "Mean1.5"="Mean batch\nfold change 1.5",
                                        "Mean2"="Mean batch\nfold change 2",
                                        "Mean3"="Mean batch\nfold change 3"))
perfstats_edgeR$Disp <- factor(perfstats_edgeR$Disp, levels=c("Disp1", "Disp2", "Disp3", "Disp4"))
perfstats_edgeR$Disp <- plyr::revalue(perfstats_edgeR$Disp, 
                                      c("Disp1"="No dispersion\nbatch effect",
                                        "Disp2"="Dispersion batch\nfold change 2",
                                        "Disp3"="Dispersion batch\nfold change 3",
                                        "Disp4"="Dispersion batch\nfold change 4"))

ggplot(perfstats_edgeR, aes(x=FPR, y=TPR, group=Method, color=Method, shape=Method)) +
  geom_point(size=2, stroke=1.1) +
  facet_grid(Disp~Mean) +
  scale_shape_manual(values=1:nlevels(perfstats_edgeR$Method)) +
  labs(x="False positive rate (FPR)", y="True positive rate (TPR)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.direction="vertical") 
