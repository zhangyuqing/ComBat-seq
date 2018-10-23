sim_inflate_zeros <- function(counts_matrix, prop_gene_partition){
  # partition all genes into the three types
  gpart_res <- sim_random_gene_partition(nrow(counts_matrix), prop_gene_partition)
  
  p_zero_seq <- rep(0, nrow(counts_matrix))
  # do nothing for type 1; for type 2, randomly select p_zero: proportion of samples to set as zeros
  p_zero_seq[gpart_res[[2]]] <- runif(length(gpart_res[[2]]), min=0, max=1)
  # for type 3, randomly select a proportion of zeros (proportion at least 80%)
  p_zero_seq[gpart_res[[3]]] <- runif(length(gpart_res[[3]]), min=0.8, max=1)
  
  for(ii in 1:nrow(counts_matrix)){
    if(p_zero_seq[ii]!=0){
      # randomly select p_zero proportion of samples to set as zero
      zero_samples <- sample(1:ncol(counts_matrix), round(ncol(counts_matrix)*p_zero_seq[ii],0), replace=FALSE)
      counts_matrix[ii, zero_samples] <- 0
    }
  }
  return(list(counts=counts_matrix, p_zero_seq=p_zero_seq))
}


sim_random_gene_partition <- function(n_genes, prop_gene_partition){
  gene_type_id <- sample(1:3, n_genes, replace=TRUE, prob=prop_gene_partition)
  gpart_res <- lapply(1:3, function(i){which(gene_type_id==i)})
}

