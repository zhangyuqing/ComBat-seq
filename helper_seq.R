# Expand a vector into matrix (columns as the original vector)
vec2mat <- function(vec, n_times){
  return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}


# Monte Carlo integration function
monte_carlo_int_NB <- function(dat, mu, gamma, phi){
  gamma_star <- phi_star <- rep(NA, nrow(dat))
  for(i in 1:nrow(dat)){
    ph <- phi[-i]		
    m <- mu[-i,!is.na(dat[i,])]
    x <- dat[i,!is.na(dat[i,])]
    LH <- sapply(1:(nrow(dat)-1), function(j){prod(dnbinom(x, mu=m[j,], size=1/ph[j]))})
    LH[is.nan(LH)]=0
    gamma_star[i] <- sum(gamma[-i]*LH)/sum(LH)
    phi_star[i] <- sum(phi[-i]*LH)/sum(LH)
  }
  res <- list(gamma_star=gamma_star, phi_star=phi_star)	
  return(res)
} 


# Match quantiles
match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi){
  new_counts_sub <- matrix(NA, nrow=nrow(counts_sub), ncol=ncol(counts_sub))
  for(a in 1:nrow(counts_sub)){
    for(b in 1:ncol(counts_sub)){
      tmp_p <- pnbinom(counts_sub[a, b], mu=old_mu[a, b], size=1/old_phi[a])
      new_counts_sub[a,b] <- qnbinom(tmp_p, mu=new_mu[a, b], size=1/new_phi[a])
    }
  }
  return(new_counts_sub)
}
