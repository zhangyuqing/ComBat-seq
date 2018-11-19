rm(list=ls())
sapply(c("ggplot2", "gridExtra", "reshape2", "ggpubr"), require, character.only=TRUE)
results_dir <- "~/Google Drive/ComBat_seq/DE_compFDR"

N_samples <- 20
disp_level_vec <- 1:5
confounding_level_vec <- seq(from=0.5, to=0.2, by=-0.1)
coverage_vec <- c(1, 5, 10)
alpha_fdr_seq <- seq(from=0, to=0.2, by=0.005)[2:4]  #[-1]
alpha_fdr_sel <- 0.1
  
sel_method_list <- list(edgeR=c("BaseIndi.edgeR", "Batch.edgeR", "OneStep.edgeR", "ComBat.lm", "RUVseq.edgeR", "SVAseq.edgeR", "ComBatseq.edgeR"),
                        DESeq2=c("BaseIndi.DESeq2", "Batch.DESeq2", "OneStep.DESeq2", "ComBat.lm", "RUVseq.DESeq2", "SVAseq.DESeq2", "ComBatseq.DESeq2"))



########  FDR corrected P values  ########
## Read in result files
fdr_all_lst <- sens_all_lst <- list()
#i=j=k=1
for(k in seq_along(coverage_vec)){
  fdr_all_lst[[k]] <- sens_all_lst[[k]] <- list()
  for(i in seq_along(confounding_level_vec)){
    fdr_all_lst[[k]][[i]] <- sens_all_lst[[k]][[i]] <- list()
    for(j in seq_along(disp_level_vec)){
      exp_name <- paste0("simCompFDR_N", N_samples, "_dispFC", disp_level_vec[j],
                         "_cnfnd", gsub(".","", confounding_level_vec[i], fixed=T), "_depth", coverage_vec[k])
      sens_res <- read.csv(file.path(results_dir, sprintf('tprADJ_%s.csv', exp_name)))
      prec_res <- read.csv(file.path(results_dir, sprintf('precADJ_%s.csv', exp_name)))
      
      sens_res_edgeR <- sens_res[, c("FDR.cutoff", sel_method_list$edgeR)]
      prec_res_edgeR <- prec_res[, c("FDR.cutoff", sel_method_list$edgeR)]
      sens_res_DESeq2 <- sens_res[, c("FDR.cutoff", sel_method_list$DESeq2)]
      prec_res_DESeq2 <- prec_res[, c("FDR.cutoff", sel_method_list$DESeq2)]
      
      ## Calculate average observed FDR under each cut-off value
      fdr_res_edgeR <- ObsNominalFDR(prec_res_edgeR, "edgeR")
      plt_fdr_edgeR <- fdr_res_edgeR$plt
      fdr_res_DESeq2 <- ObsNominalFDR(prec_res_DESeq2, "DESeq2")
      plt_fdr_DESeq2 <- fdr_res_DESeq2$plt
      
      ## For each method, use sensitivity when the corresponding observed FDR is closeset to alpha_fdr_sel (==0.1)
      sens_edgeR_ctrled <- PrecSensMapping(prec_res=prec_res_edgeR, sens_res=sens_res_edgeR)
      sens_DESeq2_ctrled <- PrecSensMapping(prec_res=prec_res_DESeq2, sens_res=sens_res_DESeq2)
    }
    names(sens_all_lst[[k]][[i]]) <- names(fdr_all_lst[[k]][[i]]) <- paste0("Disp", disp_level_vec)
  }
  names(sens_all_lst[[k]]) <- names(fdr_all_lst[[k]]) <- paste0("Cnfnd", confounding_level_vec)
}
names(sens_all_lst) <- names(fdr_all_lst) <- paste0("Depth", coverage_vec)





###  Un-adjusted P values


