rm(list=ls())
sapply(c("dplyr", "DESeq2", "ggplot2", "reshape2", "gridExtra", "alphaOutlier"), require, character.only=TRUE)
demo <- TRUE   # if testing code, set as TRUE; if running simulations, set as FALSE
if(demo){
  setwd("~/Documents/ComBat_seq/real_data_app/")
  script_dir <- "~/Dropbox/Work/ComBat_Seq/ComBat-Seq"
  data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/gfrn_signature"
  source("~/Dropbox/Work/ComBat_Seq/ComBat-Seq/real_data_app/gfrn_sig_helpers.R")
}else{
  setwd("~/yuqingz/ComBat_seq/real_data_app/")
  script_dir <- ".."; data_dir <- "."
  source("gfrn_sig_helpers.R")
}
set.seed(123)
source(file.path(script_dir, "ComBat_seq_returnParams.R"))
source(file.path(script_dir, "helper_seq.R"))

## Parameters
pathway_regex <- "akt"  
useEBpararms <- FALSE


## Load data
sigdata <- readRDS(file.path(data_dir, "signature_data.rds"))
cts_mat <- assay(sigdata, "counts")   # count matrix (also have tpm and fpkm in there)
# filter out genes with 0 counts
#cts_mat <- cts_mat[apply(cts_mat, 1, function(x){all(x!=0)}), ]
rownames(cts_mat) <- paste0("gene", 1:nrow(cts_mat))
batch <- colData(sigdata)$batch
group <- colData(sigdata)$group

# Take the subset (all controls & 1 condition specified by pathway_regex)
pathway_condition_ind <- grep(pathway_regex, group)
pathway_batch <- unique(batch[pathway_condition_ind])  # which batch does this pathway belong to
ctrl_ind <- which(group %in% c("gfp_for_egfr", "gfp18", "gfp30"))
subset_ind <- sort(c(ctrl_ind, pathway_condition_ind))  
cts_sub <- cts_mat[, subset_ind]
batch_sub <- batch[subset_ind]
group_sub <- group[subset_ind]
# remove genes with only 0 counts in the subset & in any batch
#cts_sub <- cts_sub[apply(cts_sub,1,function(x){!all(x==0)}), ]
keep1 <- apply(cts_sub[, batch_sub==1],1,function(x){!all(x==0)})
keep2 <- apply(cts_sub[, batch_sub==2],1,function(x){!all(x==0)})
keep3 <- apply(cts_sub[, batch_sub==3],1,function(x){!all(x==0)})
cts_sub <- cts_sub[keep1 & keep2 & keep3, ]


## Use ComBatSeq to adjust data
group_sub <- as.character(group_sub)
group_sub_num <- rep(1, length(group_sub))
group_sub_num[grep("gfp", group_sub)] <- 0
#counts=cts_sub;batch=batch_sub;group=group_sub_num;gene.subset.n=1000;Cpp=FALSE;covar_mod=NULL;full_mod=TRUE
start_time <- Sys.time()
combatseq_res <- ComBat_seq(counts=cts_sub, batch=batch_sub, group=group_sub_num, gene.subset.n=1000, Cpp=FALSE)#, 
                            #shrink=FALSE, shrink.disp=FALSE)
end_time <- Sys.time()
print(end_time - start_time)
#combatseq_sub <- readRDS("combatseq_gfrn_akt_params.rds")

# take the parameters
if(useEBpararms){
  par_mu <- combatseq_res$mu_star
  par_gamma <- combatseq_res$gamma_star
  par_phi <- combatseq_res$phi_star
}else{
  par_mu <- combatseq_res$mu_hat
  par_gamma <- combatseq_res$gamma_hat
  par_phi <- combatseq_res$phi_hat
}


## Outlier detection
b=1;g=0;i=1;j=1
outlier_mat <- matrix(NA, nrow=nrow(cts_sub), ncol=ncol(cts_sub))
for(b in unique(batch_sub)){
  for(g in unique(group_sub_num)){
    # within batch b & condition g
    ind_gb <- which(batch_sub==b & group_sub_num==g)
    if(length(ind_gb)>0){
      outlier_gb <- matrix(NA, nrow=nrow(cts_sub), ncol=length(ind_gb))
      cts_gb <- cts_sub[, ind_gb]
      mu_gb <- par_mu[, ind_gb]
      #gamma_gb <- par_gamma[,paste0('batch',b)]
      phi_gb <- par_phi[,paste0('batch',b)]
      
      for(i in 1:nrow(cts_gb)){
        for(j in 1:ncol(cts_gb)){
          params <- meandisp2sizeprob(mu_gb[i, j], phi_gb[i])
          if(is.na(params[2]) & cts_gb[i,j]==0){
            outlier_gb[i, j] <- FALSE
          }else{
            outlier_gb[i, j] <- aout.nbinom(cts_sub[i, j], params, alpha=0.01)$is.outlier
          }
        }
      }
      outlier_mat[, ind_gb] <- outlier_gb 
    }
  }
}
#saveRDS(outlier_mat, file="combatseq_gfrn_akt_outliers.rds")
#outlier_mat <- readRDS("combatseq_gfrn_akt_outliers.rds")

genes_to_shrink <- which(rowSums(outlier_mat) != 0)


b=1
eb_plts <- list()
for(b in unique(batch_sub)){
  params_est <- rbind(data.frame(genes=rownames(cts_sub), gamma=combatseq_res$gamma_hat[, paste0('batch',b)], 
                                 phi=combatseq_res$phi_hat[, paste0('batch',b)], EB.type="EB Off"),
                      data.frame(genes=rownames(cts_sub), gamma=combatseq_res$gamma_star[, paste0('batch',b)], 
                                 phi=combatseq_res$phi_star[, paste0('batch',b)], EB.type="EB On"))
  eb_plts[[b]] <- ggplot(params_est) +
    geom_point(aes(x=gamma, y=phi, color=EB.type)) +
    scale_color_manual(values=c("light green", "dark green")) + 
    #geom_line(aes(x=gamma, y=phi, group=genes), color="grey") +
    labs(x="Mean batch effect estimates (gamma)", y="Dispersion estimates (phi)", 
         title=sprintf("Parameter estimates within batch %s", b)) #+
    #theme_classic()
}
png(sprintf("gfrn_ShrinkParams_%s.png", gsub("^", "", pathway_regex, fixed=TRUE)), 
    width=11, height=9, units="in", res=300)
grid.arrange(eb_plts[[1]], eb_plts[[2]], eb_plts[[3]], nrow=2, ncol=2)
dev.off()


