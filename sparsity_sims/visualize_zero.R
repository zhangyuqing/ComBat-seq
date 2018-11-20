########  FDR corrected P values  ########
rm(list=ls())
setwd("~/Google Drive/ComBat_seq/")
sapply(c("ggplot2", "gridExtra", "reshape2", "ggpubr"), require, character.only=TRUE)
results_dir <- "./sparsity_sims"
source("~/Dropbox/Work/ComBat_Seq/ComBat-Seq/simulations/visualize_helpers.R")

N_samples <- 20
disp_level_vec <- 1:5
confounding_level_vec <- seq(from=0.5, to=0.2, by=-0.1)
alpha_fdr_sel <- 0.1

#i=1;j=1
fdr_lst <- sens_lst <- list()  # out layer - confounding; inner layer - disp
for(i in seq_along(confounding_level_vec)){
  fdr_lst[[i]] <- sens_lst[[i]] <- list()
  for(j in seq_along(disp_level_vec)){
    ## Read in result files
    exp_name <- paste0("simSparse_N", N_samples, "_dispFC", disp_level_vec[j], 
                       "_cnfnd", gsub(".","", confounding_level_vec[i], fixed=T))
    sens_res <- read.csv(file.path(results_dir, sprintf('tprADJ_%s.csv', exp_name)))
    prec_res <- read.csv(file.path(results_dir, sprintf('precADJ_%s.csv', exp_name)))
    
    ## Calculate average observed FDR under each cut-off value
    fdr_lst[[i]][[j]] <- data.frame(Disp=paste0("Disp",disp_level_vec[j]), ObsNominalFDR(prec_res)$fdr_df_mlt)
    
    ## For each method, use sensitivity when the corresponding observed FDR is closeset to alpha_fdr_sel (==0.1)
    sens_lst[[i]][[j]] <- PrecSensMapping(prec_df=prec_res, sens_df=sens_res, alpha_fdr_sel=alpha_fdr_sel)
  }
  fdr_lst[[i]] <- data.frame(Cnfnd=paste0("Cnfnd",confounding_level_vec[i]), do.call(rbind, fdr_lst[[i]]))
  names(sens_lst[[i]]) <- paste0("Disp", disp_level_vec)
}
fdr_lst <- do.call(rbind, fdr_lst)
colnames(fdr_lst)[4] <- "Method"
names(sens_lst) <- paste0("Cnfnd", confounding_level_vec)


##  FDR observed VS nominal
png("sparse_FDR.png", width=12, height=7, units="in", res=300)
plt_fdr <- ggplot(fdr_lst, aes(x=FDR.cutoff, y=value, group=Method, color=Method)) +
  facet_grid(Cnfnd ~ Disp) +
  geom_line() +
  geom_abline(slope=1, intercept=0, color="black") +
  labs(x="evaluation set FDR cutoff", y="observed FDR", title="DESeq2") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
print(plt_fdr) 
dev.off()


##  Sensitivity boxplot with observed FDR controlled
sens_lst <- CleanSensOut(sens_lst, ".DESeq2")
png("sparse_Sens.png", width=10, height=8, units="in", res=300)
plt_sens <- ggplot(sens_lst, aes(x=Method, y=Sensitivity, color=Method)) +
  facet_grid(Cnfnd ~ Disp) +
  geom_boxplot() +
  labs(y="Sensitivity", title="DESeq2") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1)) 
print(plt_sens)
dev.off()

 


########  Un-adjusted P values  ########
rm(list=ls())
setwd("~/Google Drive/ComBat_seq/")
sapply(c("ggplot2", "gridExtra", "reshape2", "ggpubr"), require, character.only=TRUE)
results_dir <- "./sparsity_sims"
source("~/Dropbox/Work/ComBat_Seq/ComBat-Seq/simulations/visualize_helpers.R")

N_samples <- 20
disp_level_vec <- 1:5
confounding_level_vec <- seq(from=0.5, to=0.2, by=-0.1)

#i=1;j=1
fpr_lst <- tpr_lst <- list()  # out layer - confounding; inner layer - disp
for(i in seq_along(confounding_level_vec)){
  fpr_lst[[i]] <- tpr_lst[[i]] <- list()
  for(j in seq_along(disp_level_vec)){
    ## Read in result files
    exp_name <- paste0("simSparse_N", N_samples, "_dispFC", disp_level_vec[j],
                       "_cnfnd", gsub(".","", confounding_level_vec[i], fixed=T))
    tpr_res <- read.csv(file.path(results_dir, sprintf('tpr_%s.csv', exp_name)))
    fpr_res <- read.csv(file.path(results_dir, sprintf('fpr_%s.csv', exp_name)))
    
    fpr_lst[[i]][[j]] <- as.list(colMeans(fpr_res, na.rm=TRUE))
    tpr_lst[[i]][[j]] <- as.list(colMeans(tpr_res, na.rm=TRUE))
  }
  names(fpr_lst[[i]]) <- names(tpr_lst[[i]]) <- paste0("Disp", disp_level_vec)
}
names(fpr_lst) <- names(tpr_lst) <- paste0("Cnfnd", confounding_level_vec)

tpr_merged <- CleanOut(tpr_lst, "TPR", ".DESeq2")
fpr_merged <- CleanOut(fpr_lst, "FPR", ".DESeq2")
if(identical(fpr_merged[,-1], tpr_merged[,-1])){
  perfstats_df <- data.frame(TPR=tpr_merged$TPR, FPR=fpr_merged$FPR, tpr_merged[,-1])
}else{stop("ERROR in merging data frames!")}
#identical(perfstats_df[,-1], fpr_merged)

png("sparse_scatter.png", width=10, height=7, units="in", res=300)
p2 <- ggplot(perfstats_df, aes(x=FPR, y=TPR, group=Method, color=Method, shape=Method)) +
  geom_point(size=3) +
  facet_grid(Cnfnd~Disp) +
  scale_shape_manual(values=1:nlevels(perfstats_df$Method)) +
  geom_vline(xintercept=0.05, linetype="dashed") +
  labs(title="DESeq2") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
print(p2)
dev.off()

