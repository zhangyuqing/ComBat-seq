#' Adjust for batch effects using an empirical Bayes framework in RNA-seq raw counts
#' 
#' ComBat_seq is an extension to the ComBat method using Negative Binomial model.
#' 
#' @param counts Raw count matrix from genomic studies (dimensions gene x sample) 
#' @param batch Batch covariate (only one batch allowed)
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' 
#' @return data A probe x sample count matrix, adjusted for batch effects.
#' 
#' @examples 
y1_ctrl <- matrix(rnbinom(1000*5,mu=10,size=10),ncol=5); mean(c(y1_ctrl)); var(c(y1_ctrl))
y1_case <- matrix(rnbinom(1000*5,mu=15,size=10),ncol=5); mean(c(y1_case)); var(c(y1_case))
y2_ctrl <- matrix(rnbinom(1000*3,mu=15,size=5), ncol=3); mean(c(y2_ctrl)); var(c(y2_ctrl))
y2_case <- matrix(rnbinom(1000*3,mu=20,size=5), ncol=3); mean(c(y2_case)); var(c(y2_case))
y3_ctrl <- matrix(rnbinom(1000*1,mu=5, size=20),ncol=1); mean(c(y3_ctrl)); var(c(y3_ctrl))
y3_case <- matrix(rnbinom(1000*1,mu=10,size=20),ncol=1); mean(c(y3_case)); var(c(y3_case))
counts <- cbind(y1_ctrl, y1_case, y2_ctrl, y2_case, y3_ctrl, y3_case)
batch <- c(rep("2007", ncol(y1_ctrl)+ncol(y1_case)), 
           rep("2006", ncol(y2_ctrl)+ncol(y2_case)), 
           rep("2005", ncol(y3_ctrl)+ncol(y3_case))); table(batch)
group <- c(rep(0, ncol(y1_ctrl)), rep(1, ncol(y1_case)),
           rep(0, ncol(y2_ctrl)), rep(1, ncol(y2_case)),
           rep(0, ncol(y3_ctrl)), rep(1, ncol(y3_case)))
mod <- model.matrix(~factor(group)) # Define the design matrix for the full model
#'
#' @export
#' 

ComBat_seq <- function(counts, batch, mod=NULL){
  ########  Preparation  ########  
  library(edgeR)  # require bioconductor 3.7, edgeR 3.22.1, otherwise run  # source("glmfit.R")
  
  ## Prepare characteristics on batches
  batch <- as.factor(batch)
  n_batch <- nlevels(batch)  # number of batches
  batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch  
  n_batches <- sapply(batches_ind, length)
  #if(any(n_batches==1)){mean_only=TRUE; cat("Note: one batch has only one sample, setting mean.only=TRUE\n")}
  n_sample <- sum(n_batches)
  
  ## Make batch design matrix & combine with covariates
  batchmod <- model.matrix(~-1+batch)  
  cat("Found",n_batch,'batches\n')
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
  NAs = any(is.na(counts))
  if(NAs){cat(c('Found',sum(is.na(counts)),'Missing Data Values\n'),sep=' ')}
  
  
  ########  Estimate gene-wise dispersions within each batch  ########
  ## Estimate common dispersion within each batch as an initial value
  disp_common <- sapply(1:n_batch, function(i){
    if(n_batches[i]==1){
      stop("Not supporting 1 sample per batch yet!")
    }else if(n_batches[i] <= ncol(design)-ncol(batchmod)+1){ 
      # not enough residual degree of freedom
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=NULL, subset=nrow(counts)))
      #as.matrix(design[batches_ind[[i]], (n_batch+1):ncol(design)]),
    }else{
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=mod[batches_ind[[i]], ], subset=nrow(counts)))
    }
  })
  
  ## Estimate gene-wise dispersion within each batch 
  genewise_disp_lst <- lapply(1:n_batch, function(j){
    if(n_batches[j]==1){
      stop("Not supporting 1 sample per batch yet!")
    }else if(n_batches[j] <= ncol(design)-ncol(batchmod)+1){
      # not enough residual degrees of freedom
      return(estimateGLMTagwiseDisp(counts[, batches_ind[[j]]], design=NULL, dispersion=disp_common[j], prior.df=0))
      #as.matrix(design[batches_ind[[j]], (n_batch+1):ncol(design)]),
    }else{
      return(estimateGLMTagwiseDisp(counts[, batches_ind[[j]]], design=mod[batches_ind[[j]], ], dispersion=disp_common[j], prior.df=0))
    }
  })
  
  ## construct dispersion matrix
  phi_matrix <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(k in 1:n_batch){
    phi_matrix[, batches_ind[[k]]] <- matrix(rep(genewise_disp_lst[[k]], n_batches[k]), ncol=n_batches[k])
  }
  
    
  ########  Estimate parameters from NB GLM  ########
  glmff <- glmFit(counts, design=design, dispersion=phi_matrix, offset=apply(counts,2,sum))
  
  Beta_hat <- glmff$coefficients[, (n_batch+1):ncol(design)]
  gamma_hat <- glmff$coefficients[, 1:n_batch]
  mu_hat <- glmff$fitted.values
  phi_hat <- do.call(cbind, genewise_disp_lst)
  
  
  ########  Compute posterior estimation through Monte-Carlo integration  ########  
  
  
  
  ########  Obtain adjusted batch-free distribution  ########
  
  
  
  ########  Adjust the data  ########  
  
}
