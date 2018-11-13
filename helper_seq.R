# Expand a vector into matrix (columns as the original vector)
vec2mat <- function(vec, n_times){
  return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}


# Monte Carlo integration function
monte_carlo_int_NB <- function(dat, mu, gamma, phi){
  pos_res <- lapply(1:nrow(dat), function(i){
    ph <- phi[-i]		
    m <- mu[-i,!is.na(dat[i,])]
    x <- dat[i,!is.na(dat[i,])]
    LH <- sapply(1:(nrow(dat)-1), function(j){prod(dnbinom(x, mu=m[j,], size=1/ph[j]))})
    LH[is.nan(LH)]=0
    c(gamma.star=sum(gamma[-i]*LH)/sum(LH), phi.star=sum(phi[-i]*LH)/sum(LH))
  })
  pos_res <- do.call(rbind, pos_res)
  res <- list(gamma_star=pos_res[, "gamma.star"], phi_star=pos_res[, "phi.star"])	
  return(res)
} 


# Match quantiles
match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi){
  new_counts_sub <- matrix(NA, nrow=nrow(counts_sub), ncol=ncol(counts_sub))
  for(a in 1:nrow(counts_sub)){
    for(b in 1:ncol(counts_sub)){
      if(counts_sub[a, b]==0){
        new_counts_sub[a,b] <- 0
      }else{
        tmp_p <- pnbinom(counts_sub[a, b], mu=old_mu[a, b], size=1/old_phi[a])
        if(abs(tmp_p-1)<1e-4){
          new_counts_sub[a,b] <- counts_sub[a, b]  
          # for outlier count, if p==1, will return Inf values -> use original count instead
        }else{
          new_counts_sub[a,b] <- qnbinom(tmp_p, mu=new_mu[a, b], size=1/new_phi[a])
        }
      }
    }
  }
  return(new_counts_sub)
}



#### Functions to address zeros
#cut_off=zero.fracs.cutoff
search_zin_genes <- function(cts, cut.off=NULL){
  # calculate observed zero fractions for each gene
  obs_zeros <- rowSums(cts==0)
  obs_zeros_frac <- obs_zeros / ncol(cts)
  
  # look for genes with zero fractions larger than the cut_off threshold 
  if(is.null(cut.off)){
    # cut off not specified, look for a "turning point" in observed zero fractions
    
  }else{
    # cut off specified, use the user-defined cut-off value
    zin_genes <- as.numeric(which(obs_zeros_frac > cut.off))
  }
  
  return(zin_genes)
}


## The following function takes the mean and CV estimates from DESCEND, and calculate dispersion estimates for each gene.
compute_disp <- function(mean_est, cv_est, zin.opt, zin_genes){  
  disp_est <- rep(NA, length(mean_est))
  
  mean_na_ind <- is.na(mean_est)
  cv_na_ind <- is.na(cv_est)
  est_ind <- (!mean_na_ind) & (!cv_na_ind)
  
  disp_est[est_ind] <- meanCV2disp(m=mean_est[est_ind], cv=cv_est[est_ind])
  
  return(disp_est)
}

meanCV2disp <- function(m, cv){
  disp <- cv^2 - 1/m
  return(disp)
}


## This function is an expansion of edgeR tag-wise dispersion function with additional handling of zeros
estimateGLMTagwiseDisp_ZIN <- function(curr_cts, curr_design, disp_init, zin_genes){
  library(edgeR)
  
  
  
  estimateGLMTagwiseDisp(curr_cts, design=curr_design, dispersion=disp_init, prior.df=0)
}
