## Example 1
y1_ctrl <- matrix(rnbinom(1000*5,mu=10,size=10),ncol=5); mean(c(y1_ctrl)); var(c(y1_ctrl))
y1_case <- matrix(rnbinom(1000*5,mu=15,size=10),ncol=5); mean(c(y1_case)); var(c(y1_case))
y2_ctrl <- matrix(rnbinom(1000*3,mu=15,size=5), ncol=3); mean(c(y2_ctrl)); var(c(y2_ctrl))
y2_case <- matrix(rnbinom(1000*3,mu=20,size=5), ncol=3); mean(c(y2_case)); var(c(y2_case))
y3_ctrl <- matrix(rnbinom(1000*1,mu=5, size=20),ncol=1); mean(c(y3_ctrl)); var(c(y3_ctrl))
y3_case <- matrix(rnbinom(1000*1,mu=10,size=20),ncol=1); mean(c(y3_case)); var(c(y3_case))
counts <- cbind(y1_ctrl, y1_case, y2_ctrl, y2_case, y3_ctrl, y3_case)
batch <- c(rep("B", ncol(y1_ctrl)+ncol(y1_case)), 
           rep("A", ncol(y2_ctrl)+ncol(y2_case)), 
           rep("C", ncol(y3_ctrl)+ncol(y3_case))); table(batch)
group <- c(rep(0, ncol(y1_ctrl)), rep(1, ncol(y1_case)),
           rep(0, ncol(y2_ctrl)), rep(1, ncol(y2_case)),
           rep(0, ncol(y3_ctrl)), rep(1, ncol(y3_case)))
full_mod <- TRUE
source("helper_seq.R")


## Example 2 - 2 batches
rm(list=ls())
y1_ctrl <- matrix(rnbinom(100*5,mu=10,size=10),ncol=5); mean(c(y1_ctrl)); var(c(y1_ctrl))
y1_case <- matrix(rnbinom(100*5,mu=20,size=10),ncol=5); mean(c(y1_case)); var(c(y1_case))
y1_null <- matrix(rnbinom(1900*10,mu=10,size=10),ncol=10)
y2_ctrl <- matrix(rnbinom(100*3,mu=100,size=5), ncol=3); mean(c(y2_ctrl)); var(c(y2_ctrl))
y2_case <- matrix(rnbinom(100*3,mu=200,size=5), ncol=3); mean(c(y2_case)); var(c(y2_case))
y2_null <- matrix(rnbinom(1900*6,mu=100,size=5), ncol=6)
counts <- rbind(cbind(y1_ctrl, y1_case, y2_ctrl, y2_case),
                cbind(y1_null, y2_null))
batch <- c(rep("B", ncol(y1_ctrl)+ncol(y1_case)), 
           rep("A", ncol(y2_ctrl)+ncol(y2_case))); table(batch)
group <- c(rep(0, ncol(y1_ctrl)), rep(1, ncol(y1_case)),
           rep(0, ncol(y2_ctrl)), rep(1, ncol(y2_case)))
full_mod <- TRUE
ref_batch <- NULL
source("helper_seq.R")


## Example 3 - a simpler example with equal sample size
rm(list=ls())
lambda_matrix <- matrix(NA, nrow=40, ncol=20)
for(i in 1:10){
  lambda_matrix[, i] <- c(rnorm(10, mean=0.01, sd=0.001), rnorm(10, mean=0.02, sd=0.001), 
                          rnorm(10, mean=0.03, sd=0.001), rnorm(10, mean=0.04, sd=0.001))
}
for(i in 11:20){
  lambda_matrix[, i] <- c(rnorm(10, mean=0.04, sd=0.001), rnorm(10, mean=0.03, sd=0.001), 
                          rnorm(10, mean=0.02, sd=0.001), rnorm(10, mean=0.01, sd=0.001))
}
for(i in 1:ncol(lambda_matrix)){
  lambda_matrix[, i] <- lambda_matrix[, i] / sum(lambda_matrix[, i])
}
lib_size_vec <- rep(1000, 20)
mu_true <- matrix(NA, nrow=40, ncol=20)
for(k in 1:20){
  mu_true[, k] <- lambda_matrix[, k] * lib_size_vec[k]
}
counts <- matrix(0, nrow=40, ncol=20)
for(g in 1:40){
  for(k in 1:20){
    size <- ifelse(k<=10, 1, 0.5)
    counts[g, k] <- rnbinom(1, mu=mu_true[g, k], size=size)
  }
}
group <- rep(1, 20)
batch <- c(rep(1,10), rep(2,10))
full_mod <- TRUE
normalize <- "TMM"
source("helper_seq.R")


## Example 4 - larger example with equal sample size
rm(list=ls())
#### fixed parameters
N_genes <- 200
N_DE <- 100
N_sample_vec <- rep(5, 4)
N_batch <- 2
#### dispersions phi_gi
phi_true <- cbind(runif(N_genes*N_batch, min=0.05, max=0.5), 
                  runif(N_genes*N_batch, min=0.5, max=1))
#### fraction of reads pi_gij (lambda_gij)
lambda_1_ctrl <- matrix(NA, nrow=N_DE, ncol=N_sample_vec[1])
lambda_1_case <- matrix(NA, nrow=N_DE, ncol=N_sample_vec[2])
lambda_2_ctrl <- matrix(NA, nrow=N_DE, ncol=N_sample_vec[3])
lambda_2_case <- matrix(NA, nrow=N_DE, ncol=N_sample_vec[4])
lambda_1_null <- matrix(NA, nrow=N_genes-N_DE, ncol=sum(N_sample_vec[1:2]))
lambda_2_null <- matrix(NA, nrow=N_genes-N_DE, ncol=sum(N_sample_vec[3:4]))
for(g in 1:N_genes){
  if(g <= N_DE){
    lambda_1_ctrl[g, ] <- rgamma(N_sample_vec[1], shape=1/phi_true[g, 1], rate=1/phi_true[g, 1])
    lambda_1_case[g, ] <- rgamma(N_sample_vec[2], shape=1/phi_true[g, 1], rate=0.5/phi_true[g, 1])
    lambda_2_ctrl[g, ] <- rgamma(N_sample_vec[3], shape=1/phi_true[g, 2], rate=0.4/phi_true[g, 2])
    lambda_2_case[g, ] <- rgamma(N_sample_vec[4], shape=1/phi_true[g, 2], rate=0.2/phi_true[g, 2])
  }else{
    lambda_1_null[g-N_DE, ] <- rgamma(sum(N_sample_vec[1:2]), shape=1/phi_true[g, 1], rate=1/phi_true[g, 1])
    lambda_2_null[g-N_DE, ] <- rgamma(sum(N_sample_vec[3:4]), shape=1/phi_true[g, 2], rate=0.4/phi_true[g, 2])
  }
}
# combine
lambda_true <- rbind(cbind(lambda_1_ctrl, lambda_1_case, lambda_2_ctrl, lambda_2_case),
                     cbind(lambda_1_null, lambda_2_null))
# scale to probability/proportion - columns have sum 1
for(k in 1:ncol(lambda_true)){
  lambda_true[, k] <- lambda_true[, k] / sum(lambda_true[, k])
}
#### library size for each sample
lib_size_vec <- rep(10000, sum(N_sample_vec))
# calculate mu_gij
mu_true <- matrix(0, nrow=N_genes, ncol=sum(N_sample_vec))
for(k in 1:sum(N_sample_vec)){
  mu_true[, k] <- lambda_true[, k] * lib_size_vec[k]
}
#### simulate Y_gij from NB(mu_gij, phi_gi)
counts <- matrix(0, nrow=N_genes, ncol=sum(N_sample_vec))
for(g in 1:N_genes){
  for(k in 1:sum(N_sample_vec)){
    i <- ifelse(k <= sum(N_sample_vec[1:2]), 1, 2)
    counts[g, k] <- rnbinom(1, mu=mu_true[g, k], size=1/phi_true[g, i])
  }
}
#### other parameters for combat-seq
group <- c(rep(0,N_sample_vec[1]), rep(1,N_sample_vec[2]), 
           rep(0,N_sample_vec[3]), rep(1,N_sample_vec[4]))
batch <- c(rep(1, sum(N_sample_vec[1:2])), rep(2, sum(N_sample_vec[3:4])))
full_mod <- TRUE
normalize <- "TMM"
source("helper_seq.R")




####  Example - simulate difference in relative abundance & library size
rm(list=ls())
set.seed(123)
## fixed parameters
N_genes <- 2000
N_DE <- 100
N_sample_vec <- c(10, 10, 10, 10)
N_samples <- sum(N_sample_vec)
lib_size_1 <- 10  # batch 1 library size
lib_size_2 <- 10  # batch 2 library size

## calculate a few characteristics
group <- c(rep(0,N_sample_vec[1]), rep(1,N_sample_vec[2]), 
           rep(0,N_sample_vec[3]), rep(1,N_sample_vec[4]))
design <- model.matrix(~as.factor(group))
batch <- c(rep("Bernard", sum(N_sample_vec[1:2])), rep("Arnold", sum(N_sample_vec[3:4])))

## simulate parameter by linear model (form of the GLM)
alpha_seq <- rnorm(N_genes, mean=log(0.5), sd=0.5)
beta_seq <- c(rnorm(N_DE, mean=log(2), sd=0.5), rep(0, N_genes-N_DE))  # Null genes don't respond to biological condition
# batch parameter matrices
gamma_batch1 <- c(rep(log(0.5),50), rep(log(1.5), N_DE-50), rep(log(0.5),1000), rep(log(1.5), N_genes-N_DE-1000))
gamma_batch2 <- c(rep(log(1.5),50), rep(log(0.5), N_DE-50), rep(log(1.5),1000), rep(log(0.5), N_genes-N_DE-1000))
gamma_true <- cbind(gamma_batch1, gamma_batch2)

phi_batch1 <- runif(N_genes, min=0.05, max=0.5)
phi_batch2 <- runif(N_genes, min=0.5, max=1)
phi_true <- cbind(phi_batch1, phi_batch2)

# library size
lib_size_vec <- c(rep(lib_size_1, sum(N_sample_vec[1:2])), rep(lib_size_2, sum(N_sample_vec[3:4]))) 

# calculate lambda matrix
bio_matrix <- cbind(alpha_seq, beta_seq) %*% t(design)  
gamma_matrix <- phi_matrix <- matrix(NA, nrow=N_genes, ncol=N_samples)
gamma_matrix[, batch=="Bernard"] <- gamma_batch1; gamma_matrix[, batch=="Arnold"] <- gamma_batch2
phi_matrix[, batch=="Bernard"] <- phi_batch1; phi_matrix[, batch=="Arnold"] <- phi_batch2  
  
lambda_matrix <- exp(bio_matrix + gamma_matrix)  
# for(k in 1:N_samples){
#   lambda_matrix[, k] <- lambda_matrix[, k] / sum(lambda_matrix[, k])
# } 
  
mu_matrix <- matrix(NA, nrow=N_genes, ncol=N_samples)
for(k in 1:N_samples){
  mu_matrix[, k] <- lambda_matrix[, k] * lib_size_vec[k]
}  
  
#### simulate Y_gij from NB(mu_gij, phi_gi)
counts <- matrix(NA, nrow=N_genes, ncol=N_samples)
for(g in 1:N_genes){
  for(k in 1:N_samples){
    counts[g, k] <- rnbinom(1, mu=mu_matrix[g, k], size=1/phi_matrix[g, k])
  }
}

phi_true <- phi_matrix
gamma_true <- gamma_matrix
alpha_true <- alpha_seq
beta_true <- beta_seq
lambda_true <- lambda_matrix
lib_size_true <- lib_size_vec

save(counts, group, batch,
     phi_true, gamma_true, alpha_true, beta_true, lambda_true, lib_size_true, 
     file="simdata.RData")

