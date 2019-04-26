rm(list=ls())
sapply(c("dplyr", "plyr", "SummarizedExperiment", "edgeR", "ggplot2", "reshape2", "gridExtra"), require, character.only=TRUE)
source("../../../ComBat_seq_returnParams.R")
source("../../../helper_seq.R")


### Load data
data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/WRBU_lungcancer"
cts <- readRDS(file.path(data_dir, "count_matrix_BU_WR.rds"))
meta_info <- readRDS(file.path(data_dir, "BU_WR_annotation_file.rds"))
# drop 4790-024 (single sample as its own batch)
cts <- cts[, -grep("4790-024", colnames(cts))]
meta_info <- filter(meta_info, kitnumber!="4790-024")
batch <- factor(meta_info$site)
group <- factor(meta_info$smoking_status)
# remove all 0 genes
tmp <- t(apply(cts,1,floor)); colnames(tmp) <- colnames(cts); cts <- tmp; rm(tmp)
cts <- cts[rowVars(cts)>0, ]


### Apply ComBat-seq
# combatseq_res <- ComBat_seq(cts, batch=batch, group=group, shrink=FALSE)
# params_wrbu_lungcancer <- combatseq_res[4:7]
# #saveRDS(params_wrbu_lungcancer, file="rdata/params_wrbu_lungcancer.rds")


### Estimate dispersion
dge_obj <- DGEList(counts=cts)
n_batch <- nlevels(batch)  # number of batches
batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch  
n_batches <- sapply(batches_ind, length)
n_sample <- sum(n_batches)

batchmod <- model.matrix(~-1+batch)  # colnames: levels(batch)
group <- as.factor(group)
mod <- model.matrix(~group)
design <- cbind(batchmod, mod)
check <- apply(design, 2, function(x) all(x == 1))
design <- as.matrix(design[,!check])
if(qr(design)$rank<ncol(design)){
  #if(ncol(design)<=(n_batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
  if(ncol(design)==(n_batch+1)){stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")}
  if(ncol(design)>(n_batch+1)){
    if((qr(design[,-c(1:n_batch)])$rank<ncol(design[,-c(1:n_batch)]))){stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
    }else{stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")}}
}

disp_common <- sapply(1:n_batch, function(i){
  if(n_batches[i]==1){
    stop("ComBat-seq doesn't support 1 sample per batch yet!")
  }else if((n_batches[i] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)){ 
    return(estimateGLMCommonDisp(cts[, batches_ind[[i]]], design=NULL, subset=nrow(cts)))
  }else{
    return(estimateGLMCommonDisp(cts[, batches_ind[[i]]], design=mod[batches_ind[[i]], ], subset=nrow(cts)))
  }
})
genewise_disp_lst <- lapply(1:n_batch, function(j){
  if(n_batches[j]==1){
    stop("ComBat-seq doesn't support 1 sample per batch yet!")
  }else if((n_batches[j] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)){
    return(rep(disp_common[j], nrow(cts)))
  }else{
    return(estimateGLMTagwiseDisp(cts[, batches_ind[[j]]], design=mod[batches_ind[[j]], ], 
                                  dispersion=disp_common[j], prior.df=0))#, span=1000))
  }
})
names(genewise_disp_lst) <- paste0('batch', levels(batch))

phi_matrix <- matrix(NA, nrow=nrow(cts), ncol=ncol(cts))
for(k in 1:n_batch){
  phi_matrix[, batches_ind[[k]]] <- vec2mat(genewise_disp_lst[[k]], n_batches[k]) 
}
glm_f <- glmFit(dge_obj, design=design, dispersion=phi_matrix, prior.count=1e-4) 
alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches/n_sample) 
new_offset <- t(vec2mat(getOffset(dge_obj), nrow(cts))) + vec2mat(alpha_g, ncol(cts))  
glm_f2 <- glmFit.default(dge_obj$counts, design=design, dispersion=phi_matrix, offset=new_offset, prior.count=1e-4) 

params_wrbu_lungcancer <- list(gamma_hat=glm_f2$coefficients[, 1:n_batch], phi_hat=do.call(cbind, genewise_disp_lst))


### Visualize
disp_mean_plts <- lapply(levels(factor(batch)), function(b){
  curr_gamma <- params_wrbu_lungcancer$gamma_hat[, paste0("batch", b)]
  curr_phi <- params_wrbu_lungcancer$phi_hat[, paste0("batch", b)]
  curr_df <- data.frame(gamma=curr_gamma, phi=curr_phi)
  plt <- ggplot(curr_df, aes(x=gamma, y=phi)) +
    geom_point() +
    labs(x="Mean estimates (gamma)", y="Dispersion estimates (phi)",
         title=sprintf("Lung cancer - batch %s, Spearman corr = %s", 
                       b, round(cor(curr_gamma, curr_phi, method="spearman"), 3)))
  return(plt)
})

png("figures/wrbu_lungcancer.png", width=10, height=5, units="in", res=300)
grid.arrange(disp_mean_plts[[1]], disp_mean_plts[[2]], ncol=2)
dev.off()
