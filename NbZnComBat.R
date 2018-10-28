#' Adjust for batch effects using an empirical Bayes framework in single-cell RNA-seq raw counts
#' 
#' NbZnComBat is an extension to the ComBat method using Negative Binomial model, with additional twists to address zero counts.
#' 
#' @param counts Raw count matrix from genomic studies (dimensions gene x sample) 
#' @param batch Batch covariate (only one batch allowed)
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' 
#' @return data A probe x sample count matrix, adjusted for batch effects.
#' 
#' @examples 
#' 
#' rm(list=ls()); load("simdata.RData"); full_mod=TRUE; normalize="none"; source("helper_seq.R")
#' 
#' @export
#' 

NbZnCombat <- function(cts, batch, group, full_mod=TRUE, zin.opt=TRUE, zero.fracs.cutoff=NULL){  #, normalize="none"){
  ########  Handle zeros
  if(any(cts==0) & (!zin.opt)){warnings("The count matrix contains zeros but zin.opt is not on. We recommend zero-inflated versions for this data.")}
  if(zin.opt){
    # Partition genes into zero-inflated ("zin") genes & non-zin genes
    zin_genes <- search_zin_genes(cts, cut.off=zero.fracs.cutoff)    #range(rowSums(cts[zin_genes, ]==0)/ncol(cts))
    non_zin_genes <- setdiff(1:nrow(cts), zin_genes)    #range(rowSums(cts[non_zin_genes, ]==0)/ncol(cts))
    
    # # Impute 0s in non-zin genes (kNN imputation)
    # library(DMwR)
    # tmp <- cts[non_zin_genes, ]
    # tmp[tmp==0] <- NA
    # cts[non_zin_genes, ] <- t(knnImputation(t(tmp), k=min(5, ncol(cts))))
  }else{
    zin_genes <- non_zin_genes <- NULL
  }
  
  
  ########  Preparation  ########  
  library(edgeR)  # require bioconductor 3.7, edgeR 3.22.1
  dge_obj <- DGEList(counts=cts, group=group)
  
  ## Prepare characteristics on batches
  batch <- as.factor(batch)
  n_batch <- nlevels(batch)  # number of batches
  batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch  
  n_batches <- sapply(batches_ind, length)
  #if(any(n_batches==1)){mean_only=TRUE; cat("Note: one batch has only one sample, setting mean.only=TRUE\n")}
  n_sample <- sum(n_batches)
  cat("Found",n_batch,'batches\n')
  
  ## Make design matrix 
  # batch
  batchmod <- model.matrix(~-1+batch)  # colnames: levels(batch)
  # covariate
  group <- as.factor(group)
  if(full_mod & nlevels(group)>1){
    cat("Using full model in ComBat-seq.\n")
    mod <- model.matrix(~group)
  }else{
    cat("Using null model in ComBat-seq.\n")
    mod <- model.matrix(~1, data=as.data.frame(t(cts)))
  }
  # combine
  design <- cbind(batchmod, mod)
  
  ## Check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  #if(!is.null(ref)){check[ref]=FALSE} ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  cat("Adjusting for",ncol(design)-ncol(batchmod),'covariate(s) or covariate level(s)\n')
  
  ## Check if the design is confounded
  if(qr(design)$rank<ncol(design)){
    #if(ncol(design)<=(n_batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n_batch+1)){stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")}
    if(ncol(design)>(n_batch+1)){
      if((qr(design[,-c(1:n_batch)])$rank<ncol(design[,-c(1:n_batch)]))){stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      }else{stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")}}
  }

  ## Check for missing values in count matrix
  NAs = any(is.na(cts))
  if(NAs){cat(c('Found',sum(is.na(cts)),'Missing Data Values\n'),sep=' ')}

  
  ########  Estimate gene-wise dispersions within each batch  ########
  cat("Estimating dispersions for each batch.\n")
  # Estimate common dispersion within each batch as an initial value
  disp_common <- sapply(1:n_batch, function(i){
    if(n_batches[i]==1){
      stop("ComBat-seq doesn't support 1 sample per batch yet!")
    }else if(n_batches[i] <= ncol(design)-ncol(batchmod)+1){
      # not enough residual degree of freedom
      curr_design <- NULL    
    }else{
      curr_design <- mod[batches_ind[[i]], ]    
    }
    # use only genes without zeros to estimate common dispersion
    genes_with_zeros_ind <- apply(cts[, batches_ind[[i]]], 1, function(x){any(x==0)})
    if(sum(!genes_with_zeros_ind)==0){stop("No gene is complete, not support yet!")}
    curr_disp <- estimateGLMCommonDisp(cts[!genes_with_zeros_ind, batches_ind[[i]]], design=curr_design, 
                                       subset=sum(!genes_with_zeros_ind))
    return(curr_disp)
  })

  ## Estimate gene-wise dispersion within each batch
  genewise_disp_lst <- lapply(1:n_batch, function(j){
    if(n_batches[j]==1){
      stop("ComBat-seq doesn't support 1 sample per batch yet!")
    }else if(n_batches[j] <= ncol(design)-ncol(batchmod)+1){
      # not enough residual degrees of freedom - use the common dispersion
      # return(estimateGLMTagwiseDisp(cts[, batches_ind[[j]]], design=NULL,
      #                               dispersion=disp_common[j], prior.df=0))
      return(rep(disp_common[j], nrow(cts)))
      #as.matrix(design[batches_ind[[j]], (n_batch+1):ncol(design)]),
    }else{
      genewise_disp_seq <- rep(NA, nrow(cts))
      genes_with_zeros_ind <- apply(cts[, batches_ind[[j]]], 1, function(x){any(x==0)})
      genewise_disp_seq[!genes_with_zeros_ind] <- estimateGLMTagwiseDisp(cts[!genes_with_zeros_ind, batches_ind[[j]]], 
                                                                         design=mod[batches_ind[[j]], ], 
                                                                         dispersion=disp_common[j], prior.df=0)
      genewise_disp_seq[genes_with_zeros_ind] <- disp_common[j]
      return(genewise_disp_seq)
    }
  })
  names(genewise_disp_lst) <- paste0('batch', levels(batch))

  # library(descend)
  # genewise_disp_lst <- lapply(1:n_batch, function(j){
  #   cts_batch <- cts[, batches_ind[[j]]]
  #   group_batch <- as.numeric(as.character(group[batches_ind[[j]]]))
  #   batch_ctrl <- DESCEND.control(max.sparse=c((ncol(cts_batch)-2)/ncol(cts_batch), 2))  
  #   #(max fraction of 0s, min number of non-zero counts)
  #   descend_res <- try(runDescend(cts_batch, Z=group_batch, Z0=group_batch, family="Poisson",
  #                                 control=batch_ctrl, show.message=FALSE, n.cores=3))
  #   batch_disp_estimates <- try(getEstimates(descend_res))
  #   if(class(batch_disp_estimates)=="try-error"){
  #     stop("Too sparse data!!")
  #   }else{
  #     # use mean and cv to calculate dispersion
  #     mean_est <- batch_disp_estimates$Mean[, 1]
  #     cv_est <- batch_disp_estimates$CV[, 1]
  #     return(compute_disp(mean_est=mean_est, cv_est=cv_est, zin.opt=zin.opt, zin_genes=zin_genes)) 
  #   }
  # })
  # names(genewise_disp_lst) <- paste0('batch', levels(batch))
  
  ## construct dispersion matrix
  phi_matrix <- matrix(NA, nrow=nrow(cts), ncol=ncol(cts))
  for(k in 1:n_batch){
    phi_matrix[, batches_ind[[k]]] <- vec2mat(genewise_disp_lst[[k]], n_batches[k]) #matrix(rep(genewise_disp_lst[[k]], n_batches[k]), ncol=n_batches[k])
  }#round(apply(phi_matrix,2,mean),2)
  
    
  ########  Estimate parameters from NB GLM  ########
  cat("Fitting Negative Binomial GLM.\n")
  glm_f <- glmFit(dge_obj, design=design, dispersion=phi_matrix, prior.count=0) #no intercept - nonEstimable; compute offset (library sizes) within function
  alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches/n_sample) #compute intercept as batch-size-weighted average from batches
  new_offset <- t(vec2mat(getOffset(dge_obj), nrow(cts))) +   # original offset - sample (library) size
    vec2mat(alpha_g, ncol(cts))  # new offset - gene background expression
  # getOffset(dge_obj) is the same as log(dge_obj$samples$lib.size)
  glm_f2 <- glmFit.default(dge_obj$counts, design=design, dispersion=phi_matrix, 
                           offset=new_offset, prior.count=0) 
  
  gamma_hat <- glm_f2$coefficients[, 1:n_batch]
  mu_hat <- glm_f2$fitted.values
  phi_hat <- do.call(cbind, genewise_disp_lst)
  #if(!identical(colnames(gamma_hat), colnames(phi_hat))){stop("gamma and phi don't match!")}
  #tmp = mu_hat - exp(glm_f2$coefficients %*% t(design) + new_offset); tmp[1:6,1:6]; mean(tmp)
  
  # handle zeros
  if(zin.opt){
    gamma_hat[zin_genes, ] <- 0
    #phi_hat[zin_genes, ] <- t(vec2mat(colMedians(phi_hat[non_zin_genes, ]), length(zin_genes)))
  }
  
  
  ########  In each batch, compute posterior estimation through Monte-Carlo integration  ########  
  cat("Estimating posterior parameters.\n")
  monte_carlo_res <- lapply(1:n_batch, function(ii){
    monte_carlo_int_NB(dat=cts[, batches_ind[[ii]]], mu=mu_hat[, batches_ind[[ii]]], 
                       gamma=gamma_hat[, ii], phi=phi_hat[, ii])
    #dat=cts[, batches_ind[[ii]]]; mu=mu_hat[, batches_ind[[ii]]]; gamma=gamma_hat[, ii]; phi=phi_hat[, ii]
  })
  names(monte_carlo_res) <- paste0('batch', levels(batch))
  
  gamma_star_mat <- lapply(monte_carlo_res, function(res){res$gamma_star})
  gamma_star_mat <- do.call(cbind, gamma_star_mat)
  phi_star_mat <- lapply(monte_carlo_res, function(res){res$phi_star})
  phi_star_mat <- do.call(cbind, phi_star_mat)
  
  
  ########  Obtain adjusted batch-free distribution  ########
  mu_star <- matrix(NA, nrow=nrow(cts), ncol=ncol(cts))
  for(jj in 1:n_batch){
    mu_star[, batches_ind[[jj]]] <- mu_hat[, batches_ind[[jj]]] / exp(vec2mat(gamma_star_mat[, jj], n_batches[jj]))
  }
  phi_star <- rowMeans(phi_star_mat)
  
  
  ########  Adjust the data  ########  
  cat("Adjusting the data.\n")
  adjust_cts <- matrix(NA, nrow=nrow(cts), ncol=ncol(cts), dimnames=dimnames(cts))
  for(kk in 1:n_batch){
    cts_sub <- cts[, batches_ind[[kk]]]
    old_mu <- mu_hat[, batches_ind[[kk]]]
    old_phi <- phi_hat[, kk]
    new_mu <- mu_star[, batches_ind[[kk]]]
    new_phi <- phi_star
    adjust_cts[, batches_ind[[kk]]] <- match_quantiles(counts_sub=cts_sub, 
                                                       old_mu=old_mu, old_phi=old_phi, 
                                                       new_mu=new_mu, new_phi=new_phi)
  }
  adjust_cts[is.na(adjust_cts)] <- 0
  
  return(adjust_cts)
}
