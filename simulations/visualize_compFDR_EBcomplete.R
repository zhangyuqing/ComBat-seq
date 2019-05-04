########  FDR corrected P values  ########
rm(list=ls())
setwd("~/Documents/ComBat_seq/")
sapply(c("plyr", "ggplot2", "gridExtra", "reshape2", "ggpubr"), require, character.only=TRUE)
results_dir <- "./DE_libcomp_EBcomplete"
source("~/Dropbox/Work/ComBat_Seq/ComBat-Seq/simulations/visualize_helpers.R")

N_samples <- 20
disp_level_vec <- 1:5
confounding_level_vec <- seq(from=0.5, to=0.2, by=-0.1)
coverage_level <- 5  #c(1, 5, 10)
alpha_fdr_sel <- 0.1
sel_method_list <- list(edgeR=c("Base.edgeR", "Batch.edgeR", "OneStep.edgeR", "ComBat.lm", 
                                "RUVseq.edgeR", "SVAseq.edgeR", "ComBatseq.ebON.edgeR", "ComBatseq.ebOFF.edgeR"),
                        DESeq2=c("Base.DESeq2", "Batch.DESeq2", "OneStep.DESeq2", "ComBat.lm", 
                                 "RUVseq.DESeq2", "SVAseq.DESeq2",  "ComBatseq.ebON.DESeq2", "ComBatseq.ebOFF.DESeq2"))

#i=1;j=1
fdr_lst_edgeR <- fdr_lst_DESeq2 <- sens_lst_edgeR <- sens_lst_DESeq2 <- list()  # out layer - confounding; inner layer - disp
for(i in seq_along(confounding_level_vec)){
  fdr_lst_edgeR[[i]] <- fdr_lst_DESeq2[[i]] <- sens_lst_edgeR[[i]] <- sens_lst_DESeq2[[i]] <- list()
  for(j in seq_along(disp_level_vec)){
    ## Read in result files
    exp_name <- paste0("simLibcomp_N", N_samples, "_dispFC", disp_level_vec[j],
                       "_cnfnd", gsub(".","", confounding_level_vec[i], fixed=T), "_depth", coverage_level)
    sens_res <- read.csv(file.path(results_dir, sprintf('tprADJ_%s.csv', exp_name)))
    prec_res <- read.csv(file.path(results_dir, sprintf('precADJ_%s.csv', exp_name)))
    
    sens_res <- sens_res[sens_res$FDR.cutoff %in% seq(0.01, 0.15, 0.01), ]
    prec_res <- prec_res[prec_res$FDR.cutoff %in% seq(0.01, 0.15, 0.01), ]
      
    # split DE algorithm used (edgeR & DESeq2)
    sens_res_edgeR <- sens_res[, c("FDR.cutoff", sel_method_list$edgeR)]
    prec_res_edgeR <- prec_res[, c("FDR.cutoff", sel_method_list$edgeR)]
    sens_res_DESeq2 <- sens_res[, c("FDR.cutoff", sel_method_list$DESeq2)]
    prec_res_DESeq2 <- prec_res[, c("FDR.cutoff", sel_method_list$DESeq2)]
      
    ## Calculate average observed FDR under each cut-off value
    fdr_lst_edgeR[[i]][[j]] <- data.frame(Disp=paste0("Disp",disp_level_vec[j]), ObsNominalFDR(prec_res_edgeR)$fdr_df_mlt)
    fdr_lst_DESeq2[[i]][[j]] <- data.frame(Disp=paste0("Disp",disp_level_vec[j]), ObsNominalFDR(prec_res_DESeq2)$fdr_df_mlt)

    ## For each method, use sensitivity when the corresponding observed FDR is closeset to alpha_fdr_sel (==0.1)
    sens_lst_edgeR[[i]][[j]] <- PrecSensMapping(prec_df=prec_res_edgeR, sens_df=sens_res_edgeR, alpha_fdr_sel=alpha_fdr_sel)
    sens_lst_DESeq2[[i]][[j]] <- PrecSensMapping(prec_df=prec_res_DESeq2, sens_df=sens_res_DESeq2, alpha_fdr_sel=alpha_fdr_sel)
  }
  fdr_lst_edgeR[[i]] <- data.frame(Cnfnd=paste0("Cnfnd",confounding_level_vec[i]), do.call(rbind, fdr_lst_edgeR[[i]]))
  fdr_lst_DESeq2[[i]] <- data.frame(Cnfnd=paste0("Cnfnd",confounding_level_vec[i]), do.call(rbind, fdr_lst_DESeq2[[i]]))
  names(sens_lst_edgeR[[i]]) <- names(sens_lst_DESeq2[[i]]) <- paste0("Disp", disp_level_vec)
}
fdr_lst_edgeR <- do.call(rbind, fdr_lst_edgeR)
fdr_lst_DESeq2 <- do.call(rbind, fdr_lst_DESeq2)
colnames(fdr_lst_edgeR)[4] <- colnames(fdr_lst_DESeq2)[4] <- "Method"
names(sens_lst_edgeR) <- names(sens_lst_DESeq2) <- paste0("Cnfnd", confounding_level_vec)


##  FDR observed VS nominal
fdr_lst_edgeR$Method <- as.character(fdr_lst_edgeR$Method)
fdr_lst_edgeR$Method <- gsub(".edgeR", "", fdr_lst_edgeR$Method)
fdr_lst_edgeR$Method <- factor(fdr_lst_edgeR$Method, levels=c("Base", "Batch", "OneStep", "ComBat.lm", "RUVseq", "SVAseq", 
                                                              "ComBatseq.ebON", "ComBatseq.ebOFF"))
fdr_lst_edgeR$Method <- plyr::revalue(fdr_lst_edgeR$Method, c("Base"="No Batch Effect"))
png("EBcomplete_edgeR_FDR.png", width=12, height=6, units="in", res=300)
plt_fdr_edgeR <- ggplot(fdr_lst_edgeR, aes(x=FDR.cutoff, y=value, group=Method, color=Method)) +
  facet_grid(Cnfnd ~ Disp) +
  geom_line() +
  geom_abline(slope=1, intercept=0, color="black") +
  labs(x="evaluation set FDR cutoff", y="observed FDR", title="edgeR") +
  scale_x_continuous(limits=c(0, 0.15)) +
  scale_y_continuous(limits=c(0, 0.3)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
print(plt_fdr_edgeR) 
dev.off()

fdr_lst_DESeq2$Method <- as.character(fdr_lst_DESeq2$Method)
fdr_lst_DESeq2$Method <- gsub(".DESeq2", "", fdr_lst_DESeq2$Method)
fdr_lst_DESeq2$Method <- factor(fdr_lst_DESeq2$Method, 
                                levels=c("Base", "Batch", "OneStep", "ComBat.lm", "RUVseq", "SVAseq", 
                                         "ComBatseq.ebON", "ComBatseq.ebOFF"))
fdr_lst_DESeq2$Method <- plyr::revalue(fdr_lst_DESeq2$Method, c("Base"="No Batch Effect"))
png("EBcomplete_DESeq2_FDR.png", width=12, height=6, units="in", res=300)
plt_fdr_DESeq2 <- ggplot(fdr_lst_DESeq2, aes(x=FDR.cutoff, y=value, group=Method, color=Method)) +
  facet_grid(Cnfnd ~ Disp) +
  geom_line() +
  geom_abline(slope=1, intercept=0, color="black") +
  labs(x="evaluation set FDR cutoff", y="observed FDR", title="DESeq2") +
  scale_x_continuous(limits=c(0, 0.15)) +
  scale_y_continuous(limits=c(0, 0.3)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
print(plt_fdr_DESeq2) 
dev.off()


##  Sensitivity boxplot with observed FDR controlled
sens_lst_edgeR <- CleanSensOut(sens_lst_edgeR, ".edgeR")
sens_lst_DESeq2 <- CleanSensOut(sens_lst_DESeq2, ".DESeq2")
sens_lst_edgeR$Method <- plyr::revalue(sens_lst_edgeR$Method, c("Base"="No Batch Effect"))
sens_lst_DESeq2$Method <- plyr::revalue(sens_lst_DESeq2$Method, c("Base"="No Batch Effect"))

png("EBcomplete_edgeR_Sens.png", width=10, height=7, units="in", res=300)
plt_sens_edgeR <- ggplot(sens_lst_edgeR, aes(x=Method, y=Sensitivity, color=Method)) +
  facet_grid(Cnfnd ~ Disp) +
  geom_boxplot() +
  labs(y="Sensitivity", title="edgeR") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1)) 
print(plt_sens_edgeR)
dev.off()

png("EBcomplete_DESeq2_Sens.png", width=10, height=7, units="in", res=300)
plt_sens_DESeq2 <- ggplot(sens_lst_DESeq2, aes(x=Method, y=Sensitivity, color=Method)) +
  facet_grid(Cnfnd ~ Disp) +
  geom_boxplot() +
  labs(y="Sensitivity", title="DESeq2") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1)) 
print(plt_sens_DESeq2)
dev.off()  

  

########  Un-adjusted P values  ########
rm(list=ls())
setwd("~/Documents/ComBat_seq/")
sapply(c("ggplot2", "gridExtra", "reshape2", "ggpubr"), require, character.only=TRUE)
results_dir <- "./DE_libcomp_EBcomplete"
source("~/Dropbox/Work/ComBat_Seq/ComBat-Seq/simulations/visualize_helpers.R")

N_samples <- 20
disp_level_vec <- 1:5
confounding_level_vec <- seq(from=0.5, to=0.2, by=-0.1)
coverage_level <- 5  #c(1, 5, 10)
sel_method_list <- list(edgeR=c("Base.edgeR", "Batch.edgeR", "OneStep.edgeR", "ComBat.lm", 
                                "RUVseq.edgeR", "SVAseq.edgeR", "ComBatseq.ebON.edgeR", "ComBatseq.ebOFF.edgeR"),
                        DESeq2=c("Base.DESeq2", "Batch.DESeq2", "OneStep.DESeq2", "ComBat.lm", 
                                 "RUVseq.DESeq2", "SVAseq.DESeq2",  "ComBatseq.ebON.DESeq2", "ComBatseq.ebOFF.DESeq2"))

#i=1;j=1
fpr_lst_edgeR <- fpr_lst_DESeq2 <- tpr_lst_edgeR <- tpr_lst_DESeq2 <- list()  # out layer - confounding; inner layer - disp
for(i in seq_along(confounding_level_vec)){
  fpr_lst_edgeR[[i]] <- fpr_lst_DESeq2[[i]] <- tpr_lst_edgeR[[i]] <- tpr_lst_DESeq2[[i]] <- list()
  for(j in seq_along(disp_level_vec)){
    ## Read in result files
    exp_name <- paste0("simLibcomp_N", N_samples, "_dispFC", disp_level_vec[j],
                       "_cnfnd", gsub(".","", confounding_level_vec[i], fixed=T), "_depth", coverage_level)
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
names(fpr_lst_edgeR) <- names(tpr_lst_edgeR) <- paste0("Cnfnd", confounding_level_vec)
names(fpr_lst_DESeq2) <- names(tpr_lst_DESeq2) <- paste0("Cnfnd", confounding_level_vec)

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

perfstats_edgeR$Method <- plyr::revalue(perfstats_edgeR$Method, c("Base"="No Batch Effect"))
png("EBcomplete_edgeR_scatter.png", width=10, height=6, units="in", res=300)
p2_edgeR <- ggplot(perfstats_edgeR, aes(x=FPR, y=TPR, group=Method, color=Method, shape=Method)) +
  geom_point(size=3) +
  facet_grid(Cnfnd~Disp) +
  scale_shape_manual(values=1:nlevels(perfstats_edgeR$Method)) +
  geom_vline(xintercept=0.05, linetype="dashed") +
  labs(title="edgeR") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
print(p2_edgeR)
dev.off()

perfstats_DESeq2$Method <- plyr::revalue(perfstats_DESeq2$Method, c("Base"="No Batch Effect"))
png("EBcomplete_DESeq2_scatter.png", width=10, height=6, units="in", res=300)
p2_DESeq2 <- ggplot(perfstats_DESeq2, aes(x=FPR, y=TPR, group=Method, color=Method, shape=Method)) +
  geom_point(size=3) +
  facet_grid(Cnfnd~Disp) +
  scale_shape_manual(values=1:nlevels(perfstats_DESeq2$Method)) +
  geom_vline(xintercept=0.05, linetype="dashed") +
  labs(title="DESeq2") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
print(p2_DESeq2)
dev.off()
