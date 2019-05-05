####  Expand a vector into matrix (columns as the original vector)
vec2mat <- function(vec, n_times){
  return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}


####  Monte Carlo integration functions
monte_carlo_int_NB <- function(dat, mu, gamma, phi, gene.subset.n){
  pos_res <- lapply(1:nrow(dat), function(i){
    m <- mu[-i,!is.na(dat[i,])]
    x <- dat[i,!is.na(dat[i,])]
    gamma_sub <- gamma[-i]
    phi_sub <- phi[-i]
    
    # take a subset of genes to do integration - save time
    if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & length(gene.subset.n)==1){
      if(i==1){cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))}
      mcint_ind <- sample(1:(nrow(dat)-1), gene.subset.n, replace=FALSE)
      m <- m[mcint_ind, ]; gamma_sub <- gamma_sub[mcint_ind]; phi_sub <- phi_sub[mcint_ind]
      G_sub <- gene.subset.n
    }else{
      if(i==1){cat("Using all genes for Monte Carlo integration; the function runs very slow for large number of genes\n")}
      G_sub <- nrow(dat)-1
    }
    
    LH <- sapply(1:G_sub, function(j){prod(dnbinom(x, mu=m[j,], size=1/phi_sub[j]))})
    LH[is.nan(LH)]=0; 
    if(sum(LH)==0 | is.na(sum(LH))){
      return(c(gamma.star=as.numeric(gamma[i]), phi.star=as.numeric(phi[i])))
    }else{
      return(c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH)))
    }
  })
  pos_res <- do.call(rbind, pos_res)
  res <- list(gamma_star=pos_res[, "gamma.star"], phi_star=pos_res[, "phi.star"])	
  return(res)
} 


####  Match quantiles
match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi){
  new_counts_sub <- matrix(NA, nrow=nrow(counts_sub), ncol=ncol(counts_sub))
  for(a in 1:nrow(counts_sub)){
    for(b in 1:ncol(counts_sub)){
      if(counts_sub[a, b] <= 1){
        new_counts_sub[a,b] <- counts_sub[a, b]
      }else{
        tmp_p <- pnbinom(counts_sub[a, b]-1, mu=old_mu[a, b], size=1/old_phi[a])
        if(abs(tmp_p-1)<1e-4){
          new_counts_sub[a,b] <- counts_sub[a, b]  
          # for outlier count, if p==1, will return Inf values -> use original count instead
        }else{
          new_counts_sub[a,b] <- 1+qnbinom(tmp_p, mu=new_mu[a, b], size=1/new_phi[a])
        }
      }
    }
  }
  return(new_counts_sub)
}




####  Experiments
monte_carlo_int_NB_sep_wrapper <- function(dat, mu, gamma, phi, gene.subset.n){
  gamma_star <- rep(NA, length(gamma))
  phi_star <- rep(NA, length(phi))
  
  type1_ind <- which(gamma > 0)
  type1_res <- monte_carlo_int_NB(dat[type1_ind, ], mu[type1_ind, ], gamma[type1_ind], phi[type1_ind], gene.subset.n)
  gamma_star[type1_ind] <- type1_res[["gamma_star"]]
  phi_star[type1_ind] <- type1_res[["phi_star"]]
  
  type2_ind <- which(gamma < 0)
  type2_res <- monte_carlo_int_NB(dat[type2_ind, ], mu[type2_ind, ], gamma[type2_ind], phi[type2_ind], gene.subset.n)
  gamma_star[type2_ind] <- type2_res[["gamma_star"]]
  phi_star[type2_ind] <- type2_res[["phi_star"]]
  
  if(any(is.na(gamma_star))){
    gamma_zero_ind <- which(is.na(gamma_star))
    gamma_star[gamma_zero_ind] <- gamma[gamma_zero_ind]
    phi_star[gamma_zero_ind] <- phi[gamma_zero_ind]
  }
  return(list(gamma_star=gamma_star, phi_star=phi_star))
} 

monte_carlo_int_NB_logParams <- function(dat, mu, gamma, phi, gene.subset.n){
  gamma_log <- log(abs(gamma))
  gamma_negative <- which(gamma < 0)
  phi_log <- log(phi)
  
  pos_res <- lapply(1:nrow(dat), function(i){
    #ph <- phi[-i]		
    m <- mu[-i,!is.na(dat[i,])]
    x <- dat[i,!is.na(dat[i,])]
    gamma_sub <- gamma[-i]
    phi_sub <- phi[-i]
    gamma_sub_log <- gamma_log[-i]
    phi_sub_log <- phi_log[-i]
    
    # take a subset of genes to do integration - save time
    if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & length(gene.subset.n)==1){
      if(i==1){cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))}
      mcint_ind <- sample(1:(nrow(dat)-1), gene.subset.n, replace=FALSE)
      m <- m[mcint_ind, ]; #ph <- ph[mcint_ind]; 
      gamma_sub <- gamma_sub[mcint_ind]; phi_sub <- phi_sub[mcint_ind]
      gamma_sub_log <- gamma_sub_log[mcint_ind]; phi_sub_log <- phi_sub_log[mcint_ind]
      G_sub <- gene.subset.n
    }else{
      if(i==1){cat("Using all genes for Monte Carlo integration; the function runs very slow for large number of genes\n")}
      G_sub <- nrow(dat)-1
    }
    
    LH <- sapply(1:G_sub, function(j){prod(dnbinom(x, mu=m[j,], size=1/phi_sub[j]))})
    LH[is.nan(LH)]=0; 
    if(sum(LH)==0 | is.na(sum(LH))){
      return(c(gamma.star=as.numeric(gamma[i]), phi.star=as.numeric(phi[i])))
    }else{
      return(c(gamma.star=exp(sum(gamma_sub_log*LH)/sum(LH)), phi.star=exp(sum(phi_sub_log*LH)/sum(LH))))
    }
  })
  pos_res <- do.call(rbind, pos_res)
  
  gamma_star <- pos_res[, "gamma.star"]
  gamma_negative <- setdiff(gamma_negative, which(gamma_star == gamma))
  gamma_star[gamma_negative] <- - gamma_star[gamma_negative]
  
  res <- list(gamma_star=gamma_star, phi_star=pos_res[, "phi.star"])	
  return(res)
} 

# C++ version (deprecated)
monte_carlo_int_NB_cpp <- function(dat, mu, gamma, phi, gene.subset.n){
  pos_res <- lapply(1:nrow(dat), function(i){
    ph <- phi[-i]		
    m <- mu[-i,!is.na(dat[i,])]
    x <- dat[i,!is.na(dat[i,])]
    gamma_sub <- gamma[-i]
    phi_sub <- phi[-i]
    
    # take a subset of genes to do integration - save time
    if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & length(gene.subset.n)==1){
      if(i==1){cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))}
      mcint_ind <- sample(1:(nrow(dat)-1), gene.subset.n, replace=FALSE)
      m <- m[mcint_ind, ]; ph <- ph[mcint_ind]; gamma_sub <- gamma_sub[mcint_ind]; phi_sub <- phi_sub[mcint_ind]
      G_sub <- gene.subset.n
    }else{
      if(i==1){cat("Using all genes for Monte Carlo integration; the function runs very slow for large number of genes\n")}
      G_sub <- nrow(dat)-1
    }
    
    LH <- sapply(1:G_sub, function(j){mcintCpp(as.numeric(x), as.numeric(m[j,]), 1/ph[j], length(x))})
    LH[is.nan(LH)]=0; LH[is.infinite(LH)]=0
    if(sum(LH)==0 | is.na(sum(LH))){
      return(c(gamma.star=as.numeric(gamma[i]), phi.star=as.numeric(phi[i])))
    }else{
      return(c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH)))
    }
    #c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH))
  })
  pos_res <- do.call(rbind, pos_res)
  res <- list(gamma_star=pos_res[, "gamma.star"], phi_star=pos_res[, "phi.star"])	
  return(res)
} 

#shrink only outliers
#(use all genes)
# monte_carlo_int_NB_outliers <- function(dat, mu, gamma, phi, gene.subset.n, outlier.genes){
#   ## posterior estimates
#   pos_res <- lapply(1:nrow(dat), function(i){
#     ph <- phi[-i]
#     m <- mu[-i,!is.na(dat[i,])]
#     x <- dat[i,!is.na(dat[i,])]
#     gamma_sub <- gamma[-i]
#     phi_sub <- phi[-i]
# 
#     # take a subset of genes to do integration - save time
#     if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & length(gene.subset.n)==1){
#       if(i==1){cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))}
#       #set.seed(123)
#       mcint_ind <- sample(1:(nrow(dat)-1), gene.subset.n, replace=FALSE)
#       m <- m[mcint_ind, ]; ph <- ph[mcint_ind]; gamma_sub <- gamma_sub[mcint_ind]; phi_sub <- phi_sub[mcint_ind]
#       G_sub <- gene.subset.n
#     }else{
#       if(i==1){cat("Using all genes for Monte Carlo integration; the function runs very slow for large number of genes\n")}
#       G_sub <- nrow(dat)-1
#     }
# 
#     LH <- sapply(1:G_sub, function(j){prod(dnbinom(x, mu=m[j,], size=1/ph[j]))})
#     LH[is.nan(LH)]=0;
#     if(sum(LH)==0 | is.na(sum(LH))){
#       return(c(gamma.star=as.numeric(gamma[i]), phi.star=as.numeric(phi[i])))
#     }else{
#       return(c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH)))
#     }
#     #c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH))
#   })
#   pos_res <- do.call(rbind, pos_res)
#   #res <- list(gamma_star=pos_res[, "gamma.star"], phi_star=pos_res[, "phi.star"])
# 
#   ## detect outliers
#   # is_gamma_outliers <- gamma > quantile(gamma, 0.999) | gamma < quantile(gamma, 0.001)
#   # gamma_star <- gamma; gamma_star[is_gamma_outliers] <- pos_res[is_gamma_outliers, "gamma.star"]
#   # is_phi_outliers <- phi > quantile(phi, 0.99)
#   # phi_star <- phi; phi_star[is_phi_outliers] <- pos_res[is_phi_outliers, "phi.star"]
#   #
#   # if(!is.null(rownames(dat))){
#   #   gamma_outliers <- rownames(dat)[is_gamma_outliers]
#   #   phi_outliers <- rownames(dat)[is_phi_outliers]
#   # }else{
#   #   gamma_outliers <- which(is_gamma_outliers)
#   #   phi_outliers <- which(is_phi_outliers)
#   # }
# 
#   is_gamma_outliers <- 1:length(gamma) %in% outlier.genes$mean.genes
#   gamma_star <- gamma
#   if(sum(is_gamma_outliers) > 0){
#     gamma_star[is_gamma_outliers] <- pos_res[is_gamma_outliers, "gamma.star"]
#   }
#   is_phi_outliers <- 1:length(phi) %in% outlier.genes$disp.genes
#   phi_star <- phi
#   if(sum(is_phi_outliers) > 0){
#     phi_star[is_phi_outliers] <- pos_res[is_phi_outliers, "phi.star"]
#   }
# 
#   res <- list(gamma_star=gamma_star, phi_star=phi_star)  #, gamma_outliers=gamma_outliers, phi_outliers=phi_outliers)
#   return(res)
# }

#(shrink only outliers)
monte_carlo_int_NB_outliers <- function(dat, mu, gamma, phi, gene.subset.n, outlier.genes){
  if(!is.null(outlier.genes) & length(outlier.genes$mean.genes)>0){
    Gshrink <- unique(union(outlier.genes$disp.genes, outlier.genes$mean.genes))
    if(is.null(Gshrink)){Gshrink <- 1:nrow(dat)}
  }else{
    Gshrink <- 1:nrow(dat)
  }
  cat(sprintf("Shrinking %s genes\n", length(Gshrink)))

  ## posterior estimates
  pos_res_shrink <- lapply(Gshrink, function(i){
    ph <- phi[-i]
    m <- mu[-i,!is.na(dat[i,])]
    x <- dat[i,!is.na(dat[i,])]
    gamma_sub <- gamma[-i]
    phi_sub <- phi[-i]

    # take a subset of genes to do integration - save time
    if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & length(gene.subset.n)==1){
      if(i==Gshrink[1]){cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))}
      #set.seed(123)
      mcint_ind <- sample(1:(nrow(dat)-1), gene.subset.n, replace=FALSE)
      m <- m[mcint_ind, ]; ph <- ph[mcint_ind]; gamma_sub <- gamma_sub[mcint_ind]; phi_sub <- phi_sub[mcint_ind]
      G_sub <- gene.subset.n
    }else{
      if(i==Gshrink[1]){cat("Using all genes for Monte Carlo integration; the function runs very slow for large number of genes\n")}
      G_sub <- nrow(dat)-1
    }

    LH <- sapply(1:G_sub, function(j){prod(dnbinom(x, mu=m[j,], size=1/ph[j]))})
    LH[is.nan(LH)]=0;
    if(sum(LH)==0 | is.na(sum(LH))){
      return(c(gamma.star=as.numeric(gamma[i]), phi.star=as.numeric(phi[i])))
    }else{
      return(c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH)))
    }
    #c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH))
  })
  pos_res_shrink <- do.call(rbind, pos_res_shrink)
  pos_res <- cbind(gamma.star=gamma, phi.star=phi)
  pos_res[Gshrink, ] <- pos_res_shrink
  #res <- list(gamma_star=pos_res[, "gamma.star"], phi_star=pos_res[, "phi.star"])

  is_gamma_outliers <- 1:length(gamma) %in% outlier.genes$mean.genes
  gamma_star <- gamma
  if(sum(is_gamma_outliers) > 0){
    gamma_star[is_gamma_outliers] <- pos_res[is_gamma_outliers, "gamma.star"]
  }
  is_phi_outliers <- 1:length(phi) %in% outlier.genes$disp.genes
  phi_star <- phi
  if(sum(is_phi_outliers) > 0){
    phi_star[is_phi_outliers] <- pos_res[is_phi_outliers, "phi.star"]
  }

  res <- list(gamma_star=gamma_star, phi_star=phi_star)  #, gamma_outliers=gamma_outliers, phi_outliers=phi_outliers)
  return(res)
}

#(shrink only outliers & use genes with low counts)
# monte_carlo_int_NB_outliers <- function(dat, mu, gamma, phi, gene.subset.n, outlier.genes, eb.g.pl){
#   if(!is.null(outlier.genes)){
#     Gshrink <- unique(union(outlier.genes$disp.genes, outlier.genes$mean.genes))
#     if(is.null(Gshrink)){Gshrink <- 1:nrow(dat)}
#   }else{
#     Gshrink <- 1:nrow(dat)
#   }
#   
#   ## posterior estimates
#   pos_res_shrink <- lapply(Gshrink, function(i){
#     ph <- phi[-i]
#     m <- mu[-i,!is.na(dat[i,])]
#     x <- dat[i,!is.na(dat[i,])]
#     gamma_sub <- gamma[-i]
#     phi_sub <- phi[-i]
#     
#     # take a subset of genes to do integration - save time
#     if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & length(gene.subset.n)==1){
#       if(!is.null(eb.g.pl)){
#         gene.subset.n <- min(gene.subset.n, length(eb.g.pl))
#         if(i==Gshrink[1]){cat(sprintf("Using %s random, lowly expressed genes for Monte Carlo integration\n", gene.subset.n))}
#         #set.seed(123)
#         mcint_ind <- sample(1:length(eb.g.pl), gene.subset.n, replace=FALSE)
#         m <- m[eb.g.pl[mcint_ind], ]; ph <- ph[eb.g.pl[mcint_ind]]; 
#         gamma_sub <- gamma_sub[eb.g.pl[mcint_ind]]; phi_sub <- phi_sub[eb.g.pl[mcint_ind]]
#         G_sub <- gene.subset.n
#       }else{
#         if(i==Gshrink[1]){cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))}
#         #set.seed(123)
#         mcint_ind <- sample(1:(nrow(dat)-1), gene.subset.n, replace=FALSE)
#         m <- m[mcint_ind, ]; ph <- ph[mcint_ind]; gamma_sub <- gamma_sub[mcint_ind]; phi_sub <- phi_sub[mcint_ind]
#         G_sub <- gene.subset.n
#       }
#     }else{
#       if(i==Gshrink[1]){cat("Using all genes for Monte Carlo integration; the function runs very slow for large number of genes\n")}
#       G_sub <- nrow(dat)-1
#     }
#     
#     LH <- sapply(1:G_sub, function(j){prod(dnbinom(x, mu=m[j,], size=1/ph[j]))})
#     LH[is.nan(LH)]=0;
#     if(sum(LH)==0 | is.na(sum(LH))){
#       return(c(gamma.star=as.numeric(gamma[i]), phi.star=as.numeric(phi[i])))
#     }else{
#       return(c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH)))
#     }
#     #c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH))
#   })
#   pos_res_shrink <- do.call(rbind, pos_res_shrink)
#   pos_res <- cbind(gamma.star=gamma, phi.star=phi)
#   pos_res[Gshrink, ] <- pos_res_shrink
#   #res <- list(gamma_star=pos_res[, "gamma.star"], phi_star=pos_res[, "phi.star"])
#   
#   is_gamma_outliers <- 1:length(gamma) %in% outlier.genes$mean.genes
#   gamma_star <- gamma
#   if(sum(is_gamma_outliers) > 0){
#     gamma_star[is_gamma_outliers] <- pos_res[is_gamma_outliers, "gamma.star"]
#   }
#   is_phi_outliers <- 1:length(phi) %in% outlier.genes$disp.genes
#   phi_star <- phi
#   if(sum(is_phi_outliers) > 0){
#     phi_star[is_phi_outliers] <- pos_res[is_phi_outliers, "phi.star"]
#   }
#   
#   res <- list(gamma_star=gamma_star, phi_star=phi_star)  #, gamma_outliers=gamma_outliers, phi_outliers=phi_outliers)
#   return(res)
# }

####  Generate control variable
# IQR of parameter distributions
ComBat_seq_control_IQR <- function(cts, batch, group=NULL, covar_mod=NULL, plt_distrib=TRUE){  
  full_mod <- FALSE
  if(!is.null(group)){if(nlevels(as.factor(group))>1){full_mod <- TRUE}}
  
  ##  estimate batch parameters
  dge_obj <- DGEList(counts=cts)
  batch <- as.factor(batch)
  n_batch <- nlevels(batch)  
  batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch  
  n_batches <- sapply(batches_ind, length)
  #if(any(n_batches==1)){mean_only=TRUE; cat("Note: one batch has only one sample, setting mean.only=TRUE\n")}
  n_sample <- sum(n_batches)
  batchmod <- model.matrix(~-1+batch)  
  group <- as.factor(group)
  if(full_mod & nlevels(group)>1){
    mod <- model.matrix(~group)
  }else{
    mod <- model.matrix(~1, data=as.data.frame(t(cts)))
  }
  if(!is.null(covar_mod)){covar_mod <- covar_mod[, !apply(covar_mod, 2, function(x){all(x==1)})]}
  mod <- cbind(mod, covar_mod)
  design <- cbind(batchmod, mod)
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[,!check])
  if(qr(design)$rank<ncol(design)){
    if(ncol(design)==(n_batch+1)){stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")}
    if(ncol(design)>(n_batch+1)){
      if((qr(design[,-c(1:n_batch)])$rank<ncol(design[,-c(1:n_batch)]))){
        stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
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
                                    dispersion=disp_common[j], prior.df=0))
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
  
  gamma_hat <- glm_f2$coefficients[, 1:n_batch]
  phi_hat <- do.call(cbind, genewise_disp_lst)
  rownames(phi_hat) <- rownames(gamma_hat)
  

  ##  suggest genes to be treated as outliers
  batch_levels <- levels(factor(batch))
  outlier_genes <- lapply(batch_levels, function(b){
    curr_gamma <- gamma_hat[, paste0("batch", b)]
    curr_phi <- phi_hat[, paste0("batch", b)]
    curr_cts <- cts[, batch==b]
    
    # candidate genes by boxplot quantile + IQR
    mean.upper.cutoff <- quantile(curr_gamma, 0.75) + 15 * IQR(curr_gamma)
    mean.lower.cutoff <- quantile(curr_gamma, 0.25) - 15 * IQR(curr_gamma)
    mean.genes <- c(which(curr_gamma > mean.upper.cutoff), which(curr_gamma < mean.lower.cutoff))
    
    disp.cutoff <- quantile(curr_phi, 0.75) + 30 * IQR(curr_phi)
    disp.genes <- which(curr_phi > disp.cutoff)
    
    # remove lowly expressed genes
    nsample_cutoff <- round(ncol(curr_cts) * 0.2)  # expressed in at least 20% samples of the current batch
    cts_cutoff <- quantile(curr_cts[curr_cts!=0], 0.25) # maximum count in this gene should be at least over 25% quantile of non-zero counts
    
    if(length(mean.genes)>0){
      mean.genes.keep <- sapply(mean.genes, function(g){
        return(sum(curr_cts[g, ]!=0) > nsample_cutoff & max(curr_cts[g, ], na.rm=T) > cts_cutoff)
      })
      if(sum(mean.genes.keep)>0){mean.genes <- mean.genes[mean.genes.keep]}else{mean.genes <- NULL}
    }
    
    if(length(disp.genes)>0){
      disp.genes.keep <- sapply(disp.genes, function(g){
        return(sum(curr_cts[g, ]!=0) > nsample_cutoff & max(curr_cts[g, ], na.rm=T) > cts_cutoff)
      })
      if(sum(disp.genes.keep)>0){disp.genes <- disp.genes[disp.genes.keep]}else{disp.genes <- NULL}
    }
    
    return(list(mean.genes=mean.genes, disp.genes=disp.genes))
  })
  names(outlier_genes) <- paste0("batch", batch_levels)
  
  
  ##  adjust shrink based on outlier detection
  mean.genes.pooled <- unique(do.call(c, lapply(outlier_genes, function(glst){glst$mean.genes})))
  disp.genes.pooled <- unique(do.call(c, lapply(outlier_genes, function(glst){glst$disp.genes})))
  if(!is.null(mean.genes.pooled) | !is.null(disp.genes.pooled)){
    shrink <- TRUE;  gene.subset.n <- min(nrow(cts)-1, 1000)
  }else{
    shrink <- FALSE;  gene.subset.n <- NULL
  }
  
  
  ##  visualize batch parameter estimates 
  if(plt_distrib){
    violin_plts <- lapply(batch_levels, function(b){
      curr_gamma <- gamma_hat[, paste0("batch", b)]
      curr_phi <- phi_hat[, paste0("batch", b)]
      curr_outliers <- outlier_genes[[paste0("batch", b)]]
      
      is.gamma.outlier <- 1:length(curr_gamma) %in% curr_outliers$mean.genes
      is.phi.outlier <- 1:length(curr_phi) %in% curr_outliers$disp.genes
      
      gamma_df <- data.frame(gamma=curr_gamma, is.outlier=is.gamma.outlier)
      violin_gamma <- ggplot(gamma_df, aes(x=0, y=gamma)) +
        geom_violin() +
        coord_flip() +
        labs(x="Mean estimates (gamma)") +
        theme(axis.title.x=element_blank())
      if(!is.null(curr_outliers$mean.genes) & length(curr_outliers$mean.genes)>0){
        violin_gamma <- violin_gamma + 
          geom_point(data=subset(gamma_df, is.outlier), position=position_jitter(width = 0.02), color="red")
      }
      
      phi_df <- data.frame(phi=curr_phi, is.outlier=is.phi.outlier)
      violin_phi <- ggplot(phi_df, aes(x=0, y=phi)) +
        geom_violin() +
        coord_flip() +
        labs(x="Dispersion estimates (phi)") +
        theme(axis.title.x=element_blank())
      if(!is.null(curr_outliers$disp.genes) & length(curr_outliers$disp.genes)>0){
        violin_phi <- violin_phi + 
          geom_point(data=subset(phi_df, is.outlier), position=position_jitter(width = 0.02), color="red")
      }
      
      return(annotate_figure(ggarrange(violin_gamma, violin_phi, ncol=2),
                             top=text_grob(sprintf("Batch effect parameters with suggested outliers, batch %s", b))))
    })
    do.call(grid.arrange, c(violin_plts, ncol=1))
  }
  
  outlier_genes_output <- outlier_genes
  if(is.null(mean.genes.pooled) & is.null(disp.genes.pooled)){outlier_genes_output <- NULL}
  return(list(ctrl.obj=list(full_mod=full_mod, shrink=shrink, gene.subset.n=gene.subset.n,outlier.genes=outlier_genes_output),
              params=list(gamma=gamma_hat, phi=phi_hat)))
}

# expression based: low-high counts definition of outlier
ComBat_seq_control_gene <- function(cts, batch, group=NULL){  #, cts_percent=0.95,){  
  full_mod <- FALSE
  if(!is.null(group)){if(nlevels(as.factor(group))>1){full_mod <- TRUE}}
  
  ##  suggest genes to be treated as outliers
  cts_norm <- apply(cts, 2, function(x){x/sum(x)}) # correct for library size
  batch_levels <- levels(factor(batch))
  outlier_genes <- lapply(batch_levels, function(b){
    curr_cts_comp <- cts_norm[, batch==b]
    low_comp_cutoff <- quantile(curr_cts_comp[curr_cts_comp!=0], 0.05) # 5% quantile of non-zero values of library compositions
    percent_low_comps_in_genes <- apply(curr_cts_comp, 1, function(x){length(which(x <= low_comp_cutoff))}) / ncol(curr_cts_comp)
    genes_in_threat <- which(percent_low_comps_in_genes > 0.8)  #genes with low counts in over 80% samples 
    if(length(genes_in_threat) > 1){
      genes_in_threat_vars <- rowVars(curr_cts_comp[genes_in_threat, ])
      out_genes <- genes_in_threat[genes_in_threat_vars > quantile(genes_in_threat_vars, 0.95)]
    }else if(length(genes_in_threat) == 1){
      out_genes <- genes_in_threat
    }else{
      out_genes <- NULL
    }
    return(list(mean.genes=out_genes, disp.genes=out_genes))
    # is.outlier.gene <- apply(curr_cts_comp, 1, function(x){any(x > cts_cutoff)})
    # return(list(mean.genes=which(is.outlier.gene), disp.genes=which(is.outlier.gene)))
  })
  names(outlier_genes) <- paste0("batch", batch_levels)
  
  
  ##  adjust shrink based on outlier detection
  mean.genes.pooled <- unique(do.call(c, lapply(outlier_genes, function(glst){glst$mean.genes})))
  if(!is.null(mean.genes.pooled) & length(mean.genes.pooled) > 0){
    shrink <- TRUE;  gene.subset.n <- min(nrow(cts)-1, 1000)
  }else{
    shrink <- FALSE;  gene.subset.n <- NULL
  }
  
  outlier_genes_output <- outlier_genes
  if(is.null(mean.genes.pooled)){outlier_genes_output <- NULL}
  return(list(full_mod=full_mod, shrink=shrink, gene.subset.n=gene.subset.n, outlier.genes=outlier_genes_output))
}

####  Functions to address zeros
# #cut_off=zero.fracs.cutoff
# search_zin_genes <- function(cts, cut.off=NULL){
#   # calculate observed zero fractions for each gene
#   obs_zeros <- rowSums(cts==0)
#   obs_zeros_frac <- obs_zeros / ncol(cts)
#   
#   # look for genes with zero fractions larger than the cut_off threshold 
#   if(is.null(cut.off)){
#     # cut off not specified, look for a "turning point" in observed zero fractions
#     
#   }else{
#     # cut off specified, use the user-defined cut-off value
#     zin_genes <- as.numeric(which(obs_zeros_frac > cut.off))
#   }
#   
#   return(zin_genes)
# }
# 
# 
# ## The following function takes the mean and CV estimates from DESCEND, and calculate dispersion estimates for each gene.
# compute_disp <- function(mean_est, cv_est, zin.opt, zin_genes){  
#   disp_est <- rep(NA, length(mean_est))
#   
#   mean_na_ind <- is.na(mean_est)
#   cv_na_ind <- is.na(cv_est)
#   est_ind <- (!mean_na_ind) & (!cv_na_ind)
#   
#   disp_est[est_ind] <- meanCV2disp(m=mean_est[est_ind], cv=cv_est[est_ind])
#   
#   return(disp_est)
# }
# 
# meanCV2disp <- function(m, cv){
#   disp <- cv^2 - 1/m
#   return(disp)
# }
# 
# 
# ## This function is an expansion of edgeR tag-wise dispersion function with additional handling of zeros
# estimateGLMTagwiseDisp_ZIN <- function(curr_cts, curr_design, disp_init, zin_genes){
#   library(edgeR)
#   
#   
#   
#   estimateGLMTagwiseDisp(curr_cts, design=curr_design, dispersion=disp_init, prior.df=0)
# }
